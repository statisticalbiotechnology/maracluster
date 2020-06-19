#!/usr/bin/env nextflow
echo true
println(params)
// Since we can add labels to the directories, the output path of the files are
// different compared to work
publish_output_path = params.output_path + params.output_label

// Check if resume folder has been set
resume_directory = file("NULL")
if (params.resume_directory != "") {
  resume_directory = file(params.resume_directory)
  println("Resuming from directory: ${resume_directory}")
}

file_def = file(params.batch_file)  // batch_file

if (params.workflow == "MSconvert") {
  spectra_in = Channel.from(1)  // This prevents name collision crash
  
  Channel  // non-mzML files with proper labeling which will be converted
    .from(file_def.readLines())
    .map { it -> it.tokenize('\t') }
    .map { it -> file(it[0]) }
    .into { spectra_convert }
} else {
  // Preprocessing file_list
  Channel  // mzML files with proper labeling
    .from(file_def.readLines())
    .map { it -> it.tokenize('\t') }
    .filter { it[0].tokenize('.')[-1] == "mzML" }  // filters any input that is not .mzML
    .map { it -> file(it[0]) }
    .into { spectra_in }  // Puts the files into spectra_in  
    
  Channel  // non-mzML files with proper labeling which will be converted
    .from(file_def.readLines())
    .map { it -> it.tokenize('\t') }
    .filter{ it[0].tokenize('.')[-1] != "mzML" }  // filters any input that is .mzML
    .map { it -> file(it[0]) }
    .into { spectra_convert }
}

// Replaces the lines of non-mzML with their corresponding converted mzML counterpart
count = 0
amount_of_non_mzML = 0
all_lines = file_def.readLines()
for( line in all_lines ){
  file_path = line.tokenize('\t')[0]
  file_label = line.tokenize('\t')[-1]
  file_name = (file_path.tokenize('/')[-1]).tokenize('.')[0]  // Split '/', take last. Split '.', take first
  file_extension = file_path.tokenize('.')[-1]
  if ( params.boxcar_convert == true ) {
    all_lines[count] = params.output_path + "/work/boxcar_converted_${params.random_hash}/mzML/boxcar_converted/" + file_name + ".mzML" + '\t' + file_label
  } else if ( file_extension != "mzML" ) {
    // Note: if you are running only msconvert on mzML files, the path will be wrong. However, since msconvert+maracluster
    // is not integrated yet, I can let it slide
    all_lines[count] = params.output_path + "/work/converted_${params.random_hash}/converted/" + file_name + ".mzML" + '\t' + file_label
    amount_of_non_mzML++
  } else {
    // add as is, no change
  }
  count++
}

// Create new batch file to use in work directory. Add the "corrected" raw to mzML paths here
file_def = file("$params.output_path/work/file_list_${params.random_hash}.txt")
file_def_publish = file("$publish_output_path/file_list.txt")
file_def.text = ""  // Clear file, if it exists
file_def_publish.text = ""  // Clear file, if it exists
total_files = 0
for( line in all_lines ){
  file_def << line + '\n'  // need to add \n
  file_def_publish << line + '\n'
  total_files++
}

println("Total files = " + total_files)
//println("Files that will be converted = " + amount_of_non_mzML)

file_params = file("$publish_output_path/params.txt")
file_params << "$params" + '\n'  // need to add \n

if( params.parallel_msconvert == true ) {
  spectra_convert_channel = spectra_convert  // No collect = parallel processing, one file in each process
} else {
  spectra_convert_channel = spectra_convert.collect()  // Collect = everything is run in one process
}

process msconvert {
  /* Note: maracluster needs the file_list with paths to the files.
  Problem: We do not know the mzML file location during work and writing large files to publishDir is slow,
  so there might be a possibility that some strange issues with incomplete mzML files being inputted to maracluster.
  The solution is to have two publishdirs, one where the true files are being inserted for future use,
  while we point in the file_list to the files in the work directory with symlinks.
  The symlinks are fast to write + they point to complete files + we know the symlinks location --> fixed :)
  */
  errorStrategy 'retry'  // If wine crashes for some reason, try again once
  publishDir "$params.output_path/work/converted_${params.random_hash}", mode: 'symlink', overwrite: true, pattern: "converted/*"
  if( params.publish_msconvert == true ){
    publishDir "$publish_output_path", mode: 'copy', overwrite: true, pattern: "converted/*"
  }
  containerOptions "$params.custom_mounts"
  maxForks params.parallel_msconvert_max_forks
  input:
    file f from spectra_convert_channel
  output:
    file("converted/*") into spectra_converted
  script:
  """
  mkdir -p converted
  python -s /usr/local/bin/command_wrapper.py 'wine msconvert ${f} --verbose --filter "peakPicking true 1-" -o converted ${params.msconvert_additional_arguments} | tee -a stdout.txt'
  """
}

// Concatenate these channels, even if one channel is empty, you get the contents either way
combined_channel = spectra_in.concat(spectra_converted)
  
if ( params.boxcar_convert == true) {
  process boxcar_convert {
    publishDir "$params.output_path/work/boxcar_converted_${params.random_hash}", mode: 'symlink', overwrite: true
    if( params.publish_boxcar_convert == true ){
      publishDir "$publish_output_path", mode: 'copy', overwrite: true
    }
    containerOptions "$params.custom_mounts"
    input:
      file('mzML/*') from combined_channel.collect()
    output:
      file("mzML/boxcar_converted/*") into boxcar_channel
    script:
    """
    python -s /usr/local/bin/boxcar_converter.py mzML/ ${params.boxcar_convert_additional_arguments} 2>&1 | tee -a stdout.txt
    """
  }
} else {
  boxcar_channel = combined_channel
}

// Clone channel we created to use in multiple processes (one channel per process, no more)
boxcar_channel.flatten().into {
  combined_channel_normal
  combined_channel_parallel_1
  combined_channel_parallel_2
}

// Normal, non-parallel process
process maracluster {
  // The normal process, no parallelization
  if( params.publish_maracluster == true ){
    publishDir publish_output_path, mode: 'copy', overwrite: true,  pattern: "maracluster_output/*"
  }
  containerOptions "$params.custom_mounts"
  input:
    file 'list.txt' from file_def
    file('mzML/*') from combined_channel_normal.collect()
    file('maracluster_output_resume') from resume_directory  // optional
  output:
    file("maracluster_output/*") into maracluster_out_normal includeInputs true
  when:
    (params.workflow == "Full" || params.workflow == "MaRaCluster") && params.parallel_maracluster == false
  script:
  """
  if [ -n "${params.resume_directory}" ]; then
    mkdir -p maracluster_output
    cp -as \$(pwd)/maracluster_output_resume/* maracluster_output/
  fi
  maracluster batch --batch list.txt ${params.maracluster_additional_arguments} 2>&1 | tee -a stdout.txt
  """
}


// read files and split into precursor m/z batches
process maracluster_parallel_1 {
  if( params.publish_maracluster == true ){
    publishDir publish_output_path, mode: 'copy', overwrite: true,  pattern: "maracluster_output/*"
  }
  containerOptions "$params.custom_mounts"
  input:
    file 'list.txt' from file_def
    file('mzML/*') from combined_channel_parallel_1.collect()
    file('maracluster_output_resume') from resume_directory  // optional
  output:
    file "maracluster_output/*.dat_file_list.txt" into dat_file_queue, dat_file_queue_for_overlaps
    file "maracluster_output/*" into maracluster_out_1_to_2, maracluster_out_1_to_3, maracluster_out_1_to_4 includeInputs true
  when:
    (params.workflow == "Full" || params.workflow == "MaRaCluster") && params.parallel_maracluster == true
  script:
  """
  if [ -n "${params.resume_directory}" ]; then
    mkdir -p maracluster_output
    cp -asf \$(pwd)/maracluster_output_resume/* maracluster_output/
  fi
  maracluster index --batch list.txt ${params.maracluster_additional_arguments} --output-folder ./maracluster_output 2>&1 | tee -a stdout.txt
  """
}

dat_file_queue
  .collectFile()  // Get file, will wait for process to finish
  .map { it.text }  // Convert file to text
  .splitText()  // Split text by line
  .map { it -> it.tokenize('\n')[0].tokenize('/')[-1] }
  .into { processing_tree }

// process each of the precursor m/z batches
process maracluster_parallel_2 {
  if( params.publish_maracluster == true ){
    publishDir publish_output_path, mode: 'copy', overwrite: true,  pattern: "maracluster_output/*"
  }
  containerOptions "$params.custom_mounts"
  input:
    file 'list.txt' from file_def
    val datFile from processing_tree
    file('maracluster_output/*') from maracluster_out_1_to_2.collect()
    file('maracluster_output_resume') from resume_directory  // optional
  output:
    file "maracluster_output/*pvalue*" into maracluster_out_2_to_3, maracluster_out_2_to_4 includeInputs true
  when:
    (params.workflow == "Full" || params.workflow == "MaRaCluster") && params.parallel_maracluster == true
  script:
  """
  if [ -n "${params.resume_directory}" ]; then
    mkdir -p maracluster_output
    cp -asf \$(pwd)/maracluster_output_resume/* maracluster_output/
  fi
  maracluster pvalue --batch list.txt --specIn maracluster_output/${datFile} --scanInfoFN maracluster_output/MaRaCluster.scan_info.dat --peakCountsFN maracluster_output/MaRaCluster.peak_counts.dat --prefix $datFile --clusteringTree maracluster_output/${datFile}.pvalue_tree.tsv --pvalOut maracluster_output/${datFile}.pvalues.dat ${params.maracluster_additional_arguments} 2>&1 | tee -a stdout.txt
  """
}

dat_file_queue_for_overlaps
  .collectFile()  // Get file, will wait for process to finish
  .map { it.text }  // Convert file to text
  .splitText()  // Split text by line
  .count()
  .flatMap { it -> 0..it }
  .into { overlap_processing_tree }

// process each of the overlap batches
process maracluster_parallel_3 {
  if( params.publish_maracluster == true ){
    publishDir publish_output_path, mode: 'copy', overwrite: true,  pattern: "maracluster_output/*"
  }
  containerOptions "$params.custom_mounts"
  input:
    file 'list.txt' from file_def
    val overlapBatchIdx from overlap_processing_tree
    file('maracluster_output_resume') from resume_directory  // optional
    file("maracluster_output/*") from maracluster_out_1_to_3.collect()
    file("maracluster_output/*") from maracluster_out_2_to_3.collect()
  output:
    file "maracluster_output/overlap.*.pvalue_tree.tsv" into maracluster_out_3_to_4 includeInputs true
  when:
    (params.workflow == "Full" || params.workflow == "MaRaCluster") && params.parallel_maracluster == true
  script:
  """
  if [ -n "${params.resume_directory}" ]; then
    mkdir -p maracluster_output
    cp -asf \$(pwd)/maracluster_output_resume/* maracluster_output/
  fi
  maracluster overlap --batch list.txt --datFNfile maracluster_output/MaRaCluster.dat_file_list.txt --scanInfoFN maracluster_output/MaRaCluster.scan_info.dat --overlapBatchIdx ${overlapBatchIdx} ${params.maracluster_additional_arguments} 2>&1 | tee -a stdout.txt
  """
}

// process each of the overlap batches
process maracluster_parallel_4 {
  if( params.publish_maracluster == true ){
    publishDir publish_output_path, mode: 'copy', overwrite: true,  pattern: "maracluster_output/*"
  }
  containerOptions "$params.custom_mounts"
  input:
    file 'list.txt' from file_def
    file('maracluster_output_resume') from resume_directory  // optional
    file("maracluster_output/*") from maracluster_out_1_to_4.collect()
    file("maracluster_output/*") from maracluster_out_2_to_4.collect()
    file("maracluster_output/*") from maracluster_out_3_to_4.collect()
  output:
    file "maracluster_output/*" into maracluster_out_parallel includeInputs true
  when:
    (params.workflow == "Full" || params.workflow == "MaRaCluster") && params.parallel_maracluster == true
  script:
  """
  if [ -n "${params.resume_directory}" ]; then
    mkdir -p maracluster_output
    cp -asf \$(pwd)/maracluster_output_resume/* maracluster_output/
  fi
  maracluster batch --batch list.txt --scanInfoFN maracluster_output/MaRaCluster.scan_info.dat ${params.maracluster_additional_arguments} 2>&1 | tee -a stdout.txt
  """
}

// Concatenate maracluster_out. We don't know which the user did, so we mix them
c1 = maracluster_out_normal
c2 = maracluster_out_parallel
maracluster_out = c1.concat(c2)

workflow.onComplete {
    println("MARACLUSTER PIPELINE COMPLETED")
    if (params.email != "" && params.sendfiles == true && params.workflow == "Full") {
      sendMail{to params.email
               subject "Workflow ${workflow.runName} output files"
               body "$params"
               attach file_def_publish}
    }
}
