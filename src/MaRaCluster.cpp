/******************************************************************************  
  Copyright 2015 Matthew The <matthew.the@scilifelab.se>
  Licensed under the Apache License, Version 2.0 (the "License");
  you may not use this file except in compliance with the License.
  You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.
  
 ******************************************************************************/
 
#include "MaRaCluster.h"

namespace maracluster {

MaRaCluster::MaRaCluster() :
    mode_(NONE), call_(""), percOutFN_(""), fnPrefix_("MaRaCluster"), 
    scanDescFN_(""), peakCountFN_(""), datFNFile_(""), 
    scanInfoFN_(""), pvaluesFN_(""), clusterFileFN_(""),
    pvalVecInFileFN_(""), pvalueVectorsBaseFN_(""), overlapBatchFileFN_(""), 
    spectrumBatchFileFN_(""), spectrumInFN_(""), spectrumOutFN_(""),
    spectrumLibraryFN_(""), matrixFN_(""), resultTreeFN_(""),
    skipFilterAndSort_(false), writeAll_(false), precursorTolerance_(20),
    precursorToleranceDa_(false), dbPvalThreshold_(-5.0), 
    chargeUncertainty_(0), minConsensusClusterSize_(1u)
{
  boost::filesystem::path outputPath = boost::filesystem::current_path() / boost::filesystem::path("maracluster_output");
  std::string outputFolder_ = outputPath.string();
}

MaRaCluster::~MaRaCluster() {}

std::string MaRaCluster::greeter() {
  std::ostringstream oss;
  oss << "MaRaCluster version " << VERSION << ", ";
  oss << "Build Date " << __DATE__ << " " << __TIME__ << std::endl;
  oss << "Copyright (c) 2015-17 Matthew The. All rights reserved.\n"
      << "Written by Matthew The (matthew.the@scilifelab.se) in the\n"
      << "School of Biotechnology at the Royal Institute of Technology in Stockholm.\n";
  return oss.str();
}

std::string MaRaCluster::extendedGreeter(time_t& startTime) {
  std::ostringstream oss;
  char* host = std::getenv("HOSTNAME");
  oss << greeter();
  oss << "Issued command:" << std::endl << call_ << std::endl;
  oss << "Started " << ctime(&startTime) << std::endl;
  oss.seekp(-1, ios_base::cur);
  if (host) oss << " on " << host << std::endl;
  return oss.str();
}

bool MaRaCluster::parseOptions(int argc, char **argv) {
  std::ostringstream callStream;
  callStream << argv[0];
  for (int i = 1; i < argc; i++) {
    callStream << " " << argv[i];
  }
  callStream << std::endl;
  call_ = callStream.str();
  call_ = call_.substr(0,call_.length()-1); // trim ending carriage return
  
  std::ostringstream intro;
  intro << greeter() << "\nUsage:\n";
  intro << "  maracluster batch -b <msfile_list> [-f <output_folder>]\n";
  intro << "    where msfile_list is a flat text file with absolute paths\n";
  intro << "    to the spectrum files to be clustered, one on each line.\n";
  intro << std::endl;
  
  // init
  CommandLineParser cmd(intro.str());
  cmd.defineOption("b",
      "batch",
      "File with spectrum files to be processed in batch, one per line",
      "filename");
  cmd.defineOption("f",
      "output-folder",
      "Writable folder for output files",
      "path");
  cmd.defineOption("a",
      "prefix",
      "Output files will be prefixed as e.g. <prefix>.clusters_p10.tsv (default: 'MaRaCluster')",
      "name");
  cmd.defineOption("i",
      "specIn",
      "File readable by ProteoWizard (e.g. ms2, mzML) with spectra to be"     
      " clustered",
      "filename");
  cmd.defineOption("z",
      "lib",
      "File readable by ProteoWizard (e.g. ms2, mzML) with spectral library",
      "filename");
  cmd.defineOption("l",
      "clusterFile",
      "File with filepaths and scannrs, separated by tabs",
      "filename");
  cmd.defineOption("q",
      "pvalOut",
      "File where p-values will be written to",
      "filename");
  cmd.defineOption("r",
      "pvecOut",
      "File basename where p-values vectors will be written to",
      "filename");
  cmd.defineOption("d",
      "percOut",
      "Tab delimited percolator output file containing peptides and qvalues",
      "filename");
  cmd.defineOption("g",
      "peakCountsFN",
      "File to write/read peak counts binary file",
      "filename");
  cmd.defineOption("s",
      "scanInfoFN",
      "File to write/read scan number list binary file",
      "filename");
  cmd.defineOption("j",
      "datFNfile",
      "File with a list of binary spectrum files, one per line",
      "filename");
  cmd.defineOption("m",
      "clusteringMatrix",
      "File containing the pvalue distance matrix input used for clustering in binary format.",
      "filename");
  cmd.defineOption("u",
      "clusteringTree",
      "File containing the clustering tree result as a list of merged scannrs with corresponding p value.",
      "filename");
  cmd.defineOption("c",
      "clusterThresholds",
      "Clustering thresholds listed as a comma separated list (default: -30.0,-25.0,-20.0,-15.0,-10.0,-5.0)",
      "string");
  cmd.defineOption("e",
      "skipFilterAndSort",
      "Skips filtering and sorting of the input matrix, only use if the input is a filtered and sorted binary list p-values.",
      "",
      TRUE_IF_SET);
  cmd.defineOption("p",
      "precursorTolerance",
      "Set precursor tolerance in units of ppm or Da. If no unit is specified ppm is assumed. The units have to be \"Da\" or \"ppm\", case sensitive. (default: 20.0ppm).",
      "string");
  cmd.defineOption("C",
      "chargeUncertainty",
      "Set charge uncertainty, i.e. if set to 1, then for a spectrum with precursor ion charge C, also precursor ion charges C-1 and C+1 are considered (default: 0).",
      "int");
  cmd.defineOption("t",
      "pvalThreshold",
      "Set log(p-value) threshold for database insertion (default: -5.0).",
      "double");
  cmd.defineOption("w",
      "overlapBatch",
      "File with 2 tab separated columns as: tail_file <tab> head_file for"
      " which overlapping p-values should be calculated",
      "filename");
  cmd.defineOption("y",
      "pvalueVecFile",
      "File with pvalue vectors",
      "filename");
  cmd.defineOption("o",
      "specOut",
      "File where you want the merged spectra to be written",
      "filename");
  cmd.defineOption("M",
      "minClusterSize",
      "Set the minimum size for a cluster for producing consensus spectra (default: 1).",
      "int");
  cmd.defineOption("S",
      "splitMassChargeStates",
      "Split mass charge states in spectrum output file into separate spectrum copies with the same peak list, as some formats (e.g. mgf) and software packages (e.g. MS-GF+) do not support multiple charge states for a single peak list (default: auto-detect from output file format).",
      "",
      TRUE_IF_SET);
  cmd.defineOption("v",
      "verbatim",
      "Set the verbatim level (lowest: 0, highest: 5, default: 3).",
      "int");
      
  // finally parse and handle return codes (display help etc...)
  cmd.parseArgs(argc, argv);
  
  // there should be one argument left specifying the mode
  if (cmd.arguments.size() == 1) {
    std::string mode = cmd.arguments[0];
    if (mode == "batch") mode_ = BATCH;
    else if (mode == "pvalue") mode_ = PVALUE;
    else if (mode == "unit-test") mode_ = UNIT_TEST;
    else if (mode == "index") mode_ = INDEX;
    else if (mode == "cluster") mode_ = CLUSTER;
    else if (mode == "consensus") mode_ = CONSENSUS;
    else if (mode == "search") mode_ = SEARCH;
    else if (mode == "profile-consensus") mode_ = PROFILE_CONSENSUS;
    else if (mode == "profile-search") mode_ = PROFILE_SEARCH;
    else {
      std::cerr << "Error: Unknown mode: " << mode << std::endl;
      std::cerr << "Invoke with -h option for help" << std::endl;
      return 0;
    }
  } else {
    std::cerr << "Error: too few or too many arguments." << std::endl;
    std::cerr << "Invoke with -h option for help" << std::endl;
    return 0; // ...error
  }
  
  // now query the parsing results
  
  // file input for maracluster batch and index (also for some other methods)
  if (cmd.optionSet("b")) spectrumBatchFileFN_ = cmd.options["b"];
  if (cmd.optionSet("f")) outputFolder_ = cmd.options["f"];
  if (cmd.optionSet("a")) fnPrefix_ = cmd.options["a"];
  
  // file output for maracluster batch and index (also for some other methods)
  if (cmd.optionSet("j")) datFNFile_ = cmd.options["j"];
  if (cmd.optionSet("g")) peakCountFN_ = cmd.options["g"];
  if (cmd.optionSet("s")) scanInfoFN_ = cmd.options["s"];
  
  // file input options for maracluster pvalue
  if (cmd.optionSet("i")) spectrumInFN_ = cmd.options["i"];
  if (cmd.optionSet("l")) clusterFileFN_ = cmd.options["l"];
  if (cmd.optionSet("w")) overlapBatchFileFN_ = cmd.options["w"];
  if (cmd.optionSet("y")) pvalVecInFileFN_ = cmd.options["y"];
  
  // file output options for maracluster pvalue
  if (cmd.optionSet("q")) pvaluesFN_ = cmd.options["q"];
  if (cmd.optionSet("r")) pvalueVectorsBaseFN_ = cmd.options["r"];
  
  // file input options for maracluster cluster
  if (cmd.optionSet("m")) matrixFN_ = cmd.options["m"];
  if (cmd.optionSet("u")) resultTreeFN_ = cmd.options["u"];
  if (cmd.optionSet("e")) skipFilterAndSort_ = true;
  if (cmd.optionSet("d")) scanDescFN_ = cmd.options["d"];
  if (cmd.optionSet("c")) {
    std::istringstream ss(cmd.options["c"]);
    std::string token;
    while(std::getline(ss, token, ',')) {
      clusterThresholds_.push_back(atof(token.c_str()));
    }
    if (clusterThresholds_.empty()) {
      std::cerr << "Error: invalid input for clusterThreshold parameter: " << 
                    cmd.options["c"] << std::endl;
      return false;
    }
  } else {
    clusterThresholds_.push_back(-30.0);
    clusterThresholds_.push_back(-25.0);
    clusterThresholds_.push_back(-20.0);
    clusterThresholds_.push_back(-15.0);
    clusterThresholds_.push_back(-10.0);
    clusterThresholds_.push_back(-5.0);
  }
  std::sort(clusterThresholds_.begin(), clusterThresholds_.end());
  
  // file output option for maracluster consensus
  if (cmd.optionSet("o")) spectrumOutFN_ = cmd.options["o"];
  if (cmd.optionSet("M")) minConsensusClusterSize_ = cmd.getInt("M", 0, 100000);
  if (cmd.optionSet("S")) MSFileHandler::splitMassChargeStates_ = true;
  
  // file input option for maracluster search
  if (cmd.optionSet("z")) spectrumLibraryFN_ = cmd.options["z"];
  
  // general options
  if (cmd.optionSet("t")) dbPvalThreshold_ = cmd.getDouble("t", -1000.0, 0.0);
  if (cmd.optionSet("p")) {
    std::string precursorToleranceString = cmd.options["p"];
    size_t unitIdx = precursorToleranceString.find("Da");
    if (unitIdx != std::string::npos) {
      precursorTolerance_ = atof(precursorToleranceString.substr(0, unitIdx).c_str());
      precursorToleranceDa_ = true;
    } else {
      unitIdx = precursorToleranceString.find("ppm");
      precursorTolerance_ = atof(precursorToleranceString.substr(0, unitIdx).c_str());
    }
  }
  if (cmd.optionSet("C")) chargeUncertainty_ = cmd.getInt("C", 0, 5);
  if (cmd.optionSet("v")) Globals::VERB = cmd.getInt("v", 0, 5);

  return true;
}

int MaRaCluster::createIndex() {
  if (spectrumBatchFileFN_.empty()) {
    std::cerr << "Error: no batch file specified with -b flag" << std::endl;
    return EXIT_FAILURE;
  }
  
  if (peakCountFN_.empty())
    peakCountFN_ = outputFolder_ + "/" + fnPrefix_ + ".peak_counts.dat";
  if (scanInfoFN_.empty())
    scanInfoFN_ = outputFolder_ + "/" + fnPrefix_ + ".scan_info.dat";
  if (datFNFile_.empty())
    datFNFile_ = outputFolder_ + "/" + fnPrefix_ + ".dat_file_list.txt";
  
  SpectrumFileList fileList;
  fileList.initFromFile(spectrumBatchFileFN_);
  
  if (!Globals::fileExists(datFNFile_) || !Globals::fileExists(scanInfoFN_)) {    
    BatchSpectrumFiles spectrumFiles(outputFolder_, chargeUncertainty_);
    spectrumFiles.splitByPrecursorMz(fileList, datFNFile_, peakCountFN_, 
        scanInfoFN_, precursorTolerance_, precursorToleranceDa_);
  } else {
    std::cerr << "Read dat-files from " << datFNFile_ << 
        " and scan numbers from " << scanInfoFN_ <<
        ". Remove these files to generate new dat-files." << std::endl;
  }
  
  return EXIT_SUCCESS;
}

int MaRaCluster::doClustering(const std::vector<std::string> pvalFNs, 
    std::vector<std::string> pvalTreeFNs, SpectrumFileList& fileList) {
  // input checks
  if (scanInfoFN_.empty()) {
    std::cerr << "Error: no scannrs file specified with -s flag" << std::endl;
    return EXIT_FAILURE;
  }
  if (matrixFN_.empty() && resultTreeFN_.empty()) {
    std::cerr << "Error: no input file specified with -m or -u flag" << std::endl;
    return EXIT_FAILURE;
  }

  // create output file paths
  if (resultTreeFN_.empty()) {
    resultTreeFN_ = outputFolder_ + "/overlap.pvalue_tree.tsv";
  }
  std::string clusterBaseFN = outputFolder_ + "/" + fnPrefix_ + ".clusters_";
  
  // start clustering
  if (!Globals::fileExists(resultTreeFN_)) {
    std::cerr << "Starting p-value clustering." << std::endl;
    
    if (!Globals::fileExists(matrixFN_)) {
      bool removeUnidirection = true;
      bool tsvInput = false;
      PvalueFilterAndSort::filterAndSort(pvalFNs, matrixFN_, tsvInput, removeUnidirection);
    } else {
      std::cerr << "Using p-values from " << matrixFN_ << 
          ". Remove this file to re-sort and filter the p-values." << std::endl;
    }
    
    SparseClustering matrix;
    matrix.setMergeOffset(fileList.getMergeOffset());
    matrix.initMatrix(matrixFN_);
    matrix.setClusterPairFN(resultTreeFN_);
    matrix.doClustering(std::min(dbPvalThreshold_, clusterThresholds_.back()));
    remove(matrixFN_.c_str());
  } else {
    std::cerr << "Previous clustering results are available in " << 
        resultTreeFN_ << ". Remove this file to redo the clustering." << std::endl;
  }
  pvalTreeFNs.push_back(resultTreeFN_);
  
  // write clusters
  BatchSpectrumClusters clustering;
  clustering.printClusters(pvalTreeFNs, clusterThresholds_, fileList, scanInfoFN_, scanDescFN_, clusterBaseFN);
  
  return EXIT_SUCCESS;
}

int MaRaCluster::run() {
  time_t startTime;
  clock_t startClock;
  time(&startTime);
  startClock = clock();
  
  if (Globals::VERB > 0) {
    std::cerr << extendedGreeter(startTime);
  }
  
  boost::filesystem::path rootPath (outputFolder_);
  boost::system::error_code returnedError;
  boost::filesystem::create_directories( rootPath, returnedError );
  
  if (!boost::filesystem::exists(rootPath)) {
    std::cerr << "Error: could not create output directory at " << outputFolder_ << std::endl;
    return EXIT_FAILURE;
  }
  
  switch (mode_) {
    case BATCH: {
      // This executes the entire pipeline in one go
      /* 
      maracluster -b /media/storage/mergespec/data/Linfeng/batchcluster/all.txt \
                  -f /media/storage/mergespec/data/Linfeng/batchcluster/output/ 
      */
      if (spectrumBatchFileFN_.empty()) {
        std::cerr << "Error: no batch file specified with -b flag" << std::endl;
        return EXIT_FAILURE;
      }
      
      int error = createIndex();
      if (error != EXIT_SUCCESS) return EXIT_FAILURE;
      
      std::vector<std::string> datFNs;
      {
        BatchSpectrumFiles spectrumFiles(outputFolder_, chargeUncertainty_);
        spectrumFiles.readDatFNsFromFile(datFNFile_, datFNs);
      }
      
      std::vector<std::string> pvalFNs;
      std::vector<std::string> pvalTreeFNs;
      std::vector< std::pair<std::string, std::string> > overlapFNs(datFNs.size() - 1);
      
      {
        for (size_t i = 0; i < datFNs.size(); ++i) {
          // make sure the file exists
          if (!Globals::fileExists(datFNs[i])) {
            std::cerr << "Ignoring missing data file " << datFNs[i] << std::endl;
            continue;
          }

          std::string datFN = datFNs[i];
          std::string pvalueVectorsBaseFN = datFN + ".pvalue_vectors";
          if (i < datFNs.size() - 1) {
            overlapFNs[i].first = pvalueVectorsBaseFN + ".tail.dat";
          }
          if (i > 0) {
            overlapFNs[i-1].second = pvalueVectorsBaseFN + ".head.dat";
          }
          std::string pvaluesFN = datFN + ".pvalues.dat";
          std::string pvalueTreeFN = datFN + ".pvalue_tree.tsv";
          
          if (!Globals::fileExists(pvalueTreeFN)) {
            BatchPvalueVectors pvecs(pvaluesFN, precursorTolerance_, precursorToleranceDa_, dbPvalThreshold_);
            {
              BatchSpectra spectra;
              spectra.readBatchSpectra(datFN);
              spectra.sortSpectraByPrecMz();
              
              PeakCounts peakCounts;
              peakCounts.readFromFile(peakCountFN_);
              pvecs.calculatePvalueVectors(spectra.getSpectra(), peakCounts);
            }
            pvecs.writePvalueVectors(pvalueVectorsBaseFN, writeAll_);
            pvecs.batchCalculateAndClusterPvalues(pvalueTreeFN, scanInfoFN_);
          } else {
            std::cerr << "Using p-value tree from " << pvalueTreeFN <<
                ". Remove this file to generate a new p-value tree." << std::endl;
          }
          if (Globals::fileExists(pvaluesFN)) {
            pvalFNs.push_back(pvaluesFN);
          }
          pvalTreeFNs.push_back(pvalueTreeFN);
        }
      }
      
      if (resultTreeFN_.empty()) {
        resultTreeFN_ = outputFolder_ + "/overlap.pvalue_tree.tsv";
      }
      
      if (!Globals::fileExists(resultTreeFN_)) {
        std::string pvaluesFN = outputFolder_ + "/overlap.pvalues.dat";
        if (overlapFNs.size() > 0) {
          if (!Globals::fileExists(pvaluesFN)) {
            BatchPvalueVectors pvecs(pvaluesFN, precursorTolerance_, precursorToleranceDa_, dbPvalThreshold_);
            pvecs.processOverlapFiles(overlapFNs);
          } else {
            std::cerr << "Using p-values from " << pvaluesFN << 
                ". Remove this file to generate new p-values." << std::endl;
          }
          
          if (Globals::fileExists(pvaluesFN)) {
            pvalFNs.push_back(pvaluesFN);
          }
        }
      }
      
      if (matrixFN_.empty()) {
        matrixFN_ = outputFolder_ + "/poisoned.pvalues.dat";
      }
      
      SpectrumFileList fileList;
      fileList.initFromFile(spectrumBatchFileFN_);
      error = doClustering(pvalFNs, pvalTreeFNs, fileList); 
      
      time_t endTime;
      clock_t endClock = clock();
      time(&endTime);
      double diff_time = difftime(endTime, startTime);
      
      std::cerr << "Running MaRaCluster took: "
        << ((double)(endClock - startClock)) / (double)CLOCKS_PER_SEC
        << " cpu seconds or " << diff_time << " seconds wall time" << std::endl;
                 
      return error;
    }
    case INDEX:
    {
      // maracluster index -b /media/storage/mergespec/data/batchcluster/Linfeng/all.txt
      return createIndex();
    }
    case PVALUE:
    {
      if (pvaluesFN_.empty())
        pvaluesFN_ = outputFolder_ + "/" + fnPrefix_ + ".pvalues.tsv";
        
      if (pvalVecInFileFN_.size() > 0) { 
        // direct input of p-value vector file
        // maracluster pvalue -y /media/storage/mergespec/data/batchtest/1300.ms2.pvalue_vectors.tsv
        
        BatchPvalueVectors pvecs(pvaluesFN_, precursorTolerance_, precursorToleranceDa_, dbPvalThreshold_);
        pvecs.parsePvalueVectorFile(pvalVecInFileFN_);
        if (resultTreeFN_.size() > 0) {
          pvecs.batchCalculateAndClusterPvalues(resultTreeFN_, scanInfoFN_);
        } else {
          pvecs.batchCalculatePvalues();
        }
        PvalueFilterAndSort::convertBinaryPvalToTsv(pvaluesFN_, pvaluesFN_);
      } else if (overlapBatchFileFN_.size() > 0) { 
        // calculate p-values in overlap between two windows
        // maracluster pvalue -w data/batchcluster/overlap_files.txt
        BatchPvalueVectors pvecs(pvaluesFN_, precursorTolerance_, precursorToleranceDa_, dbPvalThreshold_);
        std::vector< std::pair<std::string, std::string> > overlapFNs;
        pvecs.parseBatchOverlapFile(overlapBatchFileFN_, overlapFNs);
        pvecs.processOverlapFiles(overlapFNs);
      } else if (clusterFileFN_.size() > 0) {
        // calculate p-values from a cluster in a scan description list
        // maracluster pvalue -l <scan_desc_file> -g <peak_counts_file>
        if (peakCountFN_.empty()) {
          std::cerr << "Error: no peak counts file specified with -g flag" << std::endl;
          return EXIT_FAILURE;
        }
        std::string spectrumOutFN = "";
        MSFileExtractor fileExtractor(spectrumOutFN);
        
        std::vector<BatchSpectrum> batchSpectra;
        fileExtractor.parseClusterFileForExtract(clusterFileFN_);
        fileExtractor.extractToBatchSpectrumList(batchSpectra);
        
        PeakCounts peakCounts;
        peakCounts.readFromFile(peakCountFN_);
        
        std::cerr << "Read " << batchSpectra.size() << " spectra." << std::endl;
        
        BatchSpectra spectra;
        spectra.setBatchSpectra(batchSpectra);
        spectra.sortSpectraByPrecMz();
        
        BatchPvalueVectors pvecs(pvaluesFN_, precursorTolerance_, precursorToleranceDa_, dbPvalThreshold_);
        pvecs.calculatePvalueVectors(spectra.getSpectra(), peakCounts);
        pvecs.batchCalculatePvalues();
        PvalueFilterAndSort::convertBinaryPvalToTsv(pvaluesFN_, pvaluesFN_);
      } else {
        // calculate p-values
        /*
        maracluster pvalue -i /media/storage/mergespec/data/batchcluster/Linfeng/600.ms2 \
                        -b /media/storage/mergespec/data/batchcluster/Linfeng/all.txt \
                        -g /media/storage/mergespec/data/batchcluster/Linfeng/peak_counts.dat \
        maracluster pvalue -i /media/storage/mergespec/data/batchcluster/Linfeng/unit_test/601.dat \
                        -b /media/storage/mergespec/data/batchcluster/Linfeng/all.txt \
                        -g /media/storage/mergespec/data/batchcluster/Linfeng/unit_test/peak_counts.dat \
                        -t -10.0
        */
        if (peakCountFN_.empty()) {
          std::cerr << "Error: no peak counts file specified with -g flag" << std::endl;
          return EXIT_FAILURE;
        }
        if (spectrumBatchFileFN_.empty() && spectrumInFN_.empty()) {
          std::cerr << "Error: no input file specified with -i or -b flag" << std::endl;
          return EXIT_FAILURE;
        }
        std::string pvalueVectorsBaseFN = outputFolder_ + "/" + fnPrefix_ + ".pvalue_vectors";
        
        BatchPvalueVectors pvecs(pvaluesFN_, precursorTolerance_, precursorToleranceDa_, dbPvalThreshold_);
        {
          BatchSpectra spectra;
          if (spectrumInFN_.substr(spectrumInFN_.size() - 3) == "dat") {
            spectra.readBatchSpectra(spectrumInFN_);
          } else {
            SpectrumFileList fileList;
            fileList.initFromFile(spectrumBatchFileFN_);
            spectra.convertToBatchSpectra(spectrumInFN_, fileList);
          }
          spectra.sortSpectraByPrecMz();
          
          PeakCounts peakCounts;
          peakCounts.readFromFile(peakCountFN_);
          pvecs.calculatePvalueVectors(spectra.getSpectra(), peakCounts);
          if (pvalueVectorsBaseFN.size() > 0) {
            pvecs.writePvalueVectors(pvalueVectorsBaseFN, writeAll_);
          }
        }
        
        //pvecs.parsePvalueVectorFile(pvalueVectorsBaseFN + ".dat");
        if (resultTreeFN_.size() > 0) {
          pvecs.batchCalculateAndClusterPvalues(resultTreeFN_, scanInfoFN_);
        } else {
          pvecs.batchCalculatePvalues();
        }
      }
      return EXIT_SUCCESS;
    }
    case CLUSTER:
    {
      std::vector<std::string> pvalFNs, pvalTreeFNs;
      pvalFNs.push_back(matrixFN_);
      
      if (spectrumBatchFileFN_.empty()) {
        std::cerr << "Error: no batch file specified with -b flag" << std::endl;
        return EXIT_FAILURE;
      }
      SpectrumFileList fileList;
      fileList.initFromFile(spectrumBatchFileFN_);
      
      return doClustering(pvalFNs, pvalTreeFNs, fileList);
    }
    case CONSENSUS:
    {
      if (spectrumOutFN_.empty())
        spectrumOutFN_ = outputFolder_ + "/" + fnPrefix_ + ".consensus.ms2";
      
      if (!MSFileHandler::validMs2OutputFN(spectrumOutFN_)) {
        return EXIT_FAILURE;
      }
      
      if (MSFileHandler::getOutputFormat(spectrumOutFN_) == "mgf") {
        MSFileHandler::splitMassChargeStates_ = true;
      }
      
      if (clusterFileFN_.empty() || !Globals::fileExists(clusterFileFN_)) {
        std::cerr << "Error: Could not find cluster input file (-l flag) " << clusterFileFN_ << std::endl;
        return EXIT_FAILURE;
      }
      
      MSFileMerger msFileMerger(spectrumOutFN_);
      
      std::cerr << "Parsing cluster file" << std::endl;
      msFileMerger.parseClusterFileForMerge(clusterFileFN_, minConsensusClusterSize_);
      std::cerr << "Finished parsing cluster file" << std::endl;
      
      std::cerr << "Merging clusters" << std::endl;
      msFileMerger.mergeSpectra();
      std::cerr << "Finished merging clusters" << std::endl;
      
      return EXIT_SUCCESS;
    }
    case SEARCH:
    {
      if (spectrumLibraryFN_.empty()) {
        std::cerr << "Error: no spectrum library file specified with -z flag" << std::endl;
        return EXIT_FAILURE;
      } else if (spectrumBatchFileFN_.empty() && spectrumInFN_.empty()) {
        std::cerr << "Error: no query spectrum file(s) specified with -i or -b flag" << std::endl;
        return EXIT_FAILURE;
      } else if (!spectrumBatchFileFN_.empty() && !spectrumInFN_.empty()) {
        std::cerr << "Error: ambiguous query spectrum file(s) input, please use only one of the -i and -b flags" << std::endl;
        return EXIT_FAILURE;
      }
      
      if (pvaluesFN_.empty())
        pvaluesFN_ = outputFolder_ + "/" + fnPrefix_ + ".pvalues.dat";
      
      if (!Globals::fileExists(pvaluesFN_)) {
        // read in the query spectra
        BatchSpectra querySpectra;
        SpectrumFileList fileList;
        if (spectrumInFN_.size() > 0) {
          querySpectra.convertToBatchSpectra(spectrumInFN_, fileList);
        } else {
          fileList.initFromFile(spectrumBatchFileFN_);
          querySpectra.convertToBatchSpectra(fileList);
        }
        querySpectra.sortSpectraByPrecMz();
        
        if (peakCountFN_.empty()) {
          peakCountFN_ = outputFolder_ + "/" + fnPrefix_ + ".peak_counts.dat";
          
          BatchSpectrumFiles spectrumFiles(outputFolder_);
          std::vector<double> precMzsAccumulated;
          spectrumFiles.getPeakCountsAndPrecursorMzs(fileList, precMzsAccumulated, peakCountFN_);
        }
        std::cerr << "Reading peak counts" << std::endl;
        PeakCounts peakCounts;
        peakCounts.readFromFile(peakCountFN_);
        //peakCounts.setSmoothingMode(1);
        std::cerr << "Finished reading peak counts" << std::endl;
        
        // read in the library spectra
        PvalueCalculator::kMinScoringPeaks = 5u;
        BatchSpectra librarySpectra;
        librarySpectra.convertToBatchSpectra(spectrumLibraryFN_, fileList);
        librarySpectra.sortSpectraByPrecMz();
        
        BatchPvalueVectors pvecs(pvaluesFN_, precursorTolerance_, precursorToleranceDa_, dbPvalThreshold_);
        pvecs.calculatePvalueVectors(librarySpectra.getSpectra(), peakCounts);
        pvecs.batchCalculatePvaluesLibrarySearch(querySpectra.getSpectra());
      } else {
        std::cerr << "Using p-values from " << pvaluesFN_ << 
            ". Remove this file to generate new p-values." << std::endl;
      }
      
      return EXIT_SUCCESS;
    }
    case PROFILE_CONSENSUS:
    {
      if (spectrumOutFN_.empty())
        spectrumOutFN_ = outputFolder_ + "/" + fnPrefix_ + ".consensus.ms2";
      MSFileMerger msFileMerger(spectrumOutFN_);
      
      std::cerr << "Parsing cluster file" << std::endl;
      msFileMerger.parseClusterFileForMerge(clusterFileFN_, minConsensusClusterSize_);
      std::cerr << "Finished parsing cluster file" << std::endl;
      
      std::cerr << "Merging clusters" << std::endl;
      msFileMerger.mergeSpectra();
      std::cerr << "Finished merging clusters" << std::endl;
      
      return EXIT_SUCCESS;
    }
    case PROFILE_SEARCH:
    {
      if (spectrumLibraryFN_.empty()) {
        std::cerr << "Error: no spectrum library file specified with -z flag" << std::endl;
        return EXIT_FAILURE;
      } else if (peakCountFN_.empty()) {
        std::cerr << "Error: no peak counts file specified with -g flag" << std::endl;
        return EXIT_FAILURE;
      } else if (spectrumBatchFileFN_.empty() && spectrumInFN_.empty()) {
        std::cerr << "Error: no query spectrum file(s) specified with -i or -b flag" << std::endl;
        return EXIT_FAILURE;
      } else if (!spectrumBatchFileFN_.empty() && !spectrumInFN_.empty()) {
        std::cerr << "Error: ambiguous query spectrum file(s) input, please use only one of the -i and -b flags" << std::endl;
        return EXIT_FAILURE;
      }
      
      if (pvaluesFN_.empty())
        pvaluesFN_ = outputFolder_ + "/" + fnPrefix_ + ".pvalues.dat";
      
      if (!Globals::fileExists(pvaluesFN_)) {
        std::cerr << "Reading peak counts" << std::endl;
        PeakCounts peakCounts;
        peakCounts.readFromFile(peakCountFN_);
        //peakCounts.setSmoothingMode(1);
        std::cerr << "Finished reading peak counts" << std::endl;
        
        // read in the query spectra
        BatchSpectra querySpectra;
        SpectrumFileList fileList;
        if (spectrumInFN_.size() > 0) {
          querySpectra.convertToBatchSpectra(spectrumInFN_, fileList);
        } else {
          fileList.initFromFile(spectrumBatchFileFN_);
          querySpectra.convertToBatchSpectra(fileList);
        }
        querySpectra.sortSpectraByPrecMz();
        
        // read in the library spectra
        PvalueCalculator::kMinScoringPeaks = 5u;
        BatchSpectra librarySpectra;
        librarySpectra.convertToBatchSpectra(spectrumLibraryFN_, fileList);
        librarySpectra.sortSpectraByPrecMz();
        
        BatchPvalueVectors pvecs(pvaluesFN_, precursorTolerance_, precursorToleranceDa_, dbPvalThreshold_);
        pvecs.calculatePvalueVectors(librarySpectra.getSpectra(), peakCounts);
        pvecs.batchCalculatePvaluesLibrarySearch(querySpectra.getSpectra());
      } else {
        std::cerr << "Using p-values from " << pvaluesFN_ << 
            ". Remove this file to generate new p-values." << std::endl;
      }
      
      return EXIT_SUCCESS;
    }
    case UNIT_TEST:
    {
      unsigned int failures = 0u;
      if (PeakCountMatrix::peakMatrixUnitTest()) {
        std::cerr << "PeakCountMatrix unit tests succeeded" << std::endl;
      } else {
        std::cerr << "PeakCountMatrix unit tests failed" << std::endl;
        ++failures;
      }
      
      if (SpectrumCountVector::specVectorUnitTest()) {
        std::cerr << "SpecCountVector unit tests succeeded" << std::endl;
      } else {
        std::cerr << "SpecCountVector unit tests failed" << std::endl;
        ++failures;
      }
      
      if (PeakCounts::peakCountsSerializationUnitTest()) {
        std::cerr << "PeakCounts serialization unit tests succeeded" << std::endl;
      } else {
        std::cerr << "PeakCounts serialization unit tests failed" << std::endl;
        ++failures;
      }
      
      if (PvalueCalculator::binaryPeakMatchUnitTest()) {
        std::cerr << "PvalueCalculator peak matching unit tests succeeded" << std::endl;
      } else {
        std::cerr << "PvalueCalculator peak matching unit tests failed" << std::endl;
        ++failures;
      }
      
      if (PvalueCalculator::pvalUnitTest()) {
        std::cerr << "PvalueCalculator unit tests succeeded" << std::endl;
      } else {
        std::cerr << "PvalueCalculator unit tests failed" << std::endl;
        ++failures;
      }
      
      if (PvalueCalculator::pvalPolyfitUnitTest()) {
        std::cerr << "PvalueCalculator polyfit unit tests succeeded" << std::endl;
      } else {
        std::cerr << "PvalueCalculator polyfit unit tests failed" << std::endl;
        ++failures;
      }
      /*
      if (PvalueFilterAndSort::unitTest()) {
        std::cerr << "PvalueFilterAndSort unit tests succeeded" << std::endl;
      } else {
        std::cerr << "PvalueFilterAndSort unit tests failed" << std::endl;
        ++failures;
      }
      */
      /*
      if (PvalueCalculator::pvalUniformUnitTest()) {
        std::cerr << "PvalueCalculator uniform distribution unit tests succeeded" << std::endl;
      } else {
        std::cerr << "PvalueCalculator uniform distribution unit tests failed" << std::endl;
        ++failures;
      }
      */
      /*
      if (BatchSpectrumFiles::limitsUnitTest()) {
        std::cerr << "BatchSpectrumFiles limits unit tests succeeded" << std::endl;
      } else {
        std::cerr << "BatchSpectrumFiles limits unit tests failed" << std::endl;
        ++failures;
      }
      */
      if (MSClusterMerge::mergeUnitTest()) {
        std::cerr << "Consensus spectra unit tests succeeded" << std::endl;
      } else {
        std::cerr << "Consensus spectra unit tests failed" << std::endl;
        ++failures;
      }
      
      if (BatchSpectrumClusters::scanDescReadUnitTest()) {
        std::cerr << "Scan desc reading unit tests succeeded" << std::endl;
      } else {
        std::cerr << "Scan desc reading unit tests failed" << std::endl;
        ++failures;
      }
      
      if (failures > 0) {
        std::cerr << std::endl << "UNIT TESTS FAILED: " << failures << std::endl;
        return EXIT_FAILURE;
      } else {
        std::cerr << std::endl << "UNIT TESTS SUCCEEDED" << std::endl;
        return EXIT_SUCCESS;
      }
    }
    default:
    {
      return EXIT_FAILURE;
    }
  }

  return EXIT_FAILURE;  
  
}

} /* namespace maracluster */
