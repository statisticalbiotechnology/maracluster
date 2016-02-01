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
 
#include <cstdio>
#include <ctime>
#include <cstring>
#include <cstdlib>
#include <climits>
#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>

#include <boost/filesystem.hpp>

#include "pwiz/data/msdata/MSDataFile.hpp"

#include "Option.h"
#include "PvalueCalculator.h"
#include "SpectrumFileList.h"

#include "PercolatorInterface.cpp"

#include "BatchGlobals.h"
#include "BatchSpectrumFiles.h"
#include "BatchSpectrumClusters.h"
#include "BatchSpectra.h"
#include "BatchPvalueVectors.h"
#include "BatchPvalues.h"
#include "MSFileExtractor.h"
#include "MSFileMerger.h"
#include "MSClusterMerge.h"
#include "PvalueFilterAndSort.h"
#include "SparseClustering.h"

enum Mode { BATCH, PVALUE, UNIT_TEST, INDEX, CLUSTER, CONSENSUS, SEARCH };
Mode mode_;

std::string percOutFN_ = "";
std::string fnPrefix_ = "MaRaCluster";
std::string scanDescFN_ = "";
std::string peakCountFN_ = "";
std::string datFNFile_ = "";
std::string scanNrsFN_ = "";
std::string pvaluesFN_ = "";
std::string clusterFileFN_ = "";
std::string pvalVecInFileFN_ = "";
std::string pvalueVectorsBaseFN_ = "";
std::string overlapBatchFileFN_ = "";
boost::filesystem::path outputPath = boost::filesystem::current_path() / boost::filesystem::path("maracluster_output");
std::string outputFolder_ = outputPath.string();
std::string spectrumBatchFileFN_ = "";
std::string spectrumInFN_ = "";
std::string spectrumOutFN_ = "";
std::string spectrumLibraryFN_ = "";

std::string matrixFN_ = "";
std::string resultTreeFN_ = "";
bool skipFilterAndSort_ = false;
std::vector<double> clusterThresholds_;

bool parseOptions(int argc, char **argv) {
  std::ostringstream intro;
  intro << "Usage:\n";
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
      "Output files will be prefixed as e.g. <prefix>.pvalue_tree.tsv (default: 'MaRaCluster')",
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
      "scanNrsFN",
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
      "precursorTolerancePpm",
      "Set precursor ppm tolerance (default: 20.0).",
      "double");
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
  if (cmd.optionSet("s")) scanNrsFN_ = cmd.options["s"];
  
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
    if (clusterThresholds_.size() == 0) {
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
  
  // file input option for maracluster search
  if (cmd.optionSet("z")) spectrumLibraryFN_ = cmd.options["z"];
  
  // general options
  if (cmd.optionSet("t")) BatchPvalueVectors::dbPvalThreshold_ = cmd.getDouble("t", -1000.0, 0.0);
  if (cmd.optionSet("p")) BatchPvalueVectors::massRangePPM_ = cmd.getDouble("p", 0.0, 1e6);
  if (cmd.optionSet("v")) BatchGlobals::VERB = cmd.getInt("v", 0, 5);

  return true;
}

int createIndex(const std::string& outputFolder, const std::string& fnPrefix, 
    const std::string& spectrumBatchFileFN, std::string& peakCountFN,
    std::string& scanNrsFN, std::string& datFNFile) {
  if (spectrumBatchFileFN.size() == 0) {
    std::cerr << "Error: no batch file specified with -b flag" << std::endl;
    return EXIT_FAILURE;
  }
  
  if (peakCountFN.size() == 0)
    peakCountFN = outputFolder + "/" + fnPrefix + ".peak_counts.dat";
  if (scanNrsFN.size() == 0)
    scanNrsFN = outputFolder + "/" + fnPrefix + ".scannrs.dat";
  if (datFNFile.size() == 0)
    datFNFile = outputFolder + "/" + fnPrefix + ".dat_file_list.txt";
  
  SpectrumFileList fileList;
  fileList.initFromFile(spectrumBatchFileFN);
  
  if (!BatchGlobals::fileExists(datFNFile)) {    
    BatchSpectrumFiles spectrumFiles(outputFolder);
    spectrumFiles.splitByPrecursorMass(fileList, datFNFile, peakCountFN, scanNrsFN);
  } else {
    std::cerr << "Read dat-files from " << datFNFile << 
        ". Remove this file to generate new dat-files." << std::endl;
  }
  
  if (!BatchGlobals::fileExists(scanNrsFN)) {    
    BatchSpectrumFiles spectrumFiles(outputFolder);
    spectrumFiles.writeScannrs(fileList, scanNrsFN);
  } else {
    std::cerr << "Read scan numbers from " << scanNrsFN << 
        ". Remove this file to generate a new scannr list." << std::endl;
  }
  
  return EXIT_SUCCESS;
}

int doClustering(const std::string& outputFolder, const std::string& fnPrefix, 
    const std::vector<std::string> pvalFNs, const std::string& scanNrsFN, 
    const std::string& scanDescFN, const std::vector<double>& clusterThresholds, 
    SpectrumFileList& fileList, bool skipFilterAndSort, 
    std::string& matrixFN, std::string& resultTreeFN) {
  // input checks
  if (scanNrsFN.size() == 0) {
    std::cerr << "Error: no scannrs file specified with -s flag" << std::endl;
    return EXIT_FAILURE;
  }
  if (matrixFN.size() == 0 && resultTreeFN.size() == 0) {
    std::cerr << "Error: no input file specified with -m or -u flag" << std::endl;
    return EXIT_FAILURE;
  }

  // create output file paths
  if (resultTreeFN.size() == 0) {
    resultTreeFN = outputFolder + "/" + fnPrefix + ".pvalue_tree.tsv";
  }
  std::string clusterBaseFN = outputFolder + "/" + fnPrefix + ".clusters_";
  
  // start clustering
  if (!BatchGlobals::fileExists(resultTreeFN)) {
    std::cerr << "Starting p-value clustering." << std::endl;
    
    if (!skipFilterAndSort) {
      matrixFN = outputFolder + "/" + fnPrefix + ".pvalue_triplets.dat";
    }
    
    if (!BatchGlobals::fileExists(matrixFN)) {
      bool removeUnidirection = true;
      bool tsvInput = false;
      PvalueFilterAndSort::filterAndSort(pvalFNs, matrixFN, tsvInput, removeUnidirection);
    } else {
      std::cerr << "Using p-values from " << matrixFN << 
          ". Remove this file to re-sort and filter the p-values." << std::endl;
    }
    
    SparseClustering matrix;
    matrix.setMergeOffset(fileList.getMergeOffset());
    matrix.initMatrix(matrixFN);
    matrix.setClusterPairFN(resultTreeFN);
    matrix.doClustering(std::min(BatchPvalueVectors::dbPvalThreshold_, clusterThresholds.back())); 
  } else {
    std::cerr << "Previous clustering results are available in " << 
        resultTreeFN << ". Remove this file to redo the clustering." << std::endl;
  }
  
  // write clusters
  BatchSpectrumClusters clustering;
  clustering.printClusters(resultTreeFN, clusterThresholds, fileList, scanNrsFN, scanDescFN, clusterBaseFN);
  
  return EXIT_SUCCESS;
}

int main(int argc, char* argv[]) {
  try {
    if (parseOptions(argc, argv)) {
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
          // maracluster -b /media/storage/mergespec/data/batchcluster/Linfeng/all.txt \
                         -f /media/storage/mergespec/data/batchcluster/Linfeng/output/
          if (spectrumBatchFileFN_.size() == 0) {
            std::cerr << "Error: no batch file specified with -b flag" << std::endl;
            return EXIT_FAILURE;
          } else {
            SpectrumFileList fileList;
            fileList.initFromFile(spectrumBatchFileFN_);
          }
          
          int error = createIndex(outputFolder_, fnPrefix_, 
                                  spectrumBatchFileFN_, peakCountFN_, 
                                  scanNrsFN_, datFNFile_);
          if (error != EXIT_SUCCESS) return EXIT_FAILURE;
          
          std::vector<std::string> datFNs;
          {
            BatchSpectrumFiles spectrumFiles(outputFolder_);
            spectrumFiles.readDatFNsFromFile(datFNFile_, datFNs);
          }
          
          std::vector<std::string> pvalFNs;
          std::vector< std::pair<std::string, std::string> > overlapFNs(datFNs.size() - 1);
          
          {
            PeakCounts peakCounts;
            peakCounts.readFromFile(peakCountFN_);
            
            for (size_t i = 0; i < datFNs.size(); ++i) {
              // make sure the file exists
              if (!BatchGlobals::fileExists(datFNs[i])) {
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
              
              if (!BatchGlobals::fileExists(pvaluesFN)) {
                BatchSpectra spectra(pvaluesFN);
                spectra.readBatchSpectra(datFN);
                spectra.calculatePvalueVectors(peakCounts);
                spectra.writePvalueVectors(pvalueVectorsBaseFN);
                spectra.calculatePvalues();
              } else {
                std::cerr << "Using p-values from " << pvaluesFN << 
                    ". Remove this file to generate new p-values." << std::endl;
              }
              pvalFNs.push_back(pvaluesFN);
            }
          }
          
          std::string pvaluesFN = outputFolder_ + "/overlaps.pvalues.dat";
          if (overlapFNs.size() > 0) {
            if (!BatchGlobals::fileExists(pvaluesFN)) {
              BatchPvalueVectors pvecs(pvaluesFN);
              pvecs.processOverlapFiles(overlapFNs);
            } else {
              std::cerr << "Using p-values from " << pvaluesFN << 
                  ". Remove this file to generate new p-values." << std::endl;
            }
            
            if (BatchGlobals::fileExists(pvaluesFN)) {
              pvalFNs.push_back(pvaluesFN);
            }
          }
          
          if (matrixFN_.size() == 0) {
            matrixFN_ = outputFolder_ + "/" + fnPrefix_ + ".pvalue_triplets.dat";
          }
          
          SpectrumFileList fileList;
          fileList.initFromFile(spectrumBatchFileFN_);
          error = doClustering(outputFolder_, fnPrefix_, pvalFNs, scanNrsFN_, scanDescFN_,
                              clusterThresholds_, fileList, skipFilterAndSort_, 
                              matrixFN_, resultTreeFN_);          
          return error;
        }
        case INDEX:
        {
          // maracluster index -b /media/storage/mergespec/data/batchcluster/Linfeng/all.txt
          return createIndex(outputFolder_, fnPrefix_, spectrumBatchFileFN_, 
                             peakCountFN_, scanNrsFN_, datFNFile_);
        }
        case PVALUE:
        {
          if (pvaluesFN_.size() == 0)
            pvaluesFN_ = outputFolder_ + "/" + fnPrefix_ + ".pvalues.tsv";
            
          if (pvalVecInFileFN_.size() > 0) { 
            // direct input of p-value vector file
            // maracluster pvalue -y /media/storage/mergespec/data/batchtest/1300.ms2.pvalue_vectors.tsv
            
            BatchPvalueVectors pvecs(pvaluesFN_);
            pvecs.parsePvalueVectorFile(pvalVecInFileFN_);
            pvecs.batchCalculatePvalues();
          } else if (overlapBatchFileFN_.size() > 0) { 
            // calculate p-values in overlap between two windows
            // maracluster pvalue -w data/batchcluster/overlap_files.txt
            BatchPvalueVectors pvecs(pvaluesFN_);
            std::vector< std::pair<std::string, std::string> > overlapFNs;
            pvecs.parseBatchOverlapFile(overlapBatchFileFN_, overlapFNs);
            pvecs.processOverlapFiles(overlapFNs);
          } else if (clusterFileFN_.size() > 0) {
            // calculate p-values from a cluster in a scan description list
            // maracluster pvalue -l <scan_desc_file> -g <peak_counts_file>
            if (peakCountFN_.size() == 0) {
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
            
            BatchSpectra spectra(pvaluesFN_);
            spectra.setBatchSpectra(batchSpectra);
            spectra.calculatePvalueVectors(peakCounts);
            spectra.calculatePvalues();
            PvalueFilterAndSort::convertBinaryPvalToTsv(pvaluesFN_, pvaluesFN_);
          } else {
            // calculate p-values
            // maracluster pvalue -i /media/storage/mergespec/data/batchcluster/Linfeng/600.ms2 \
                            -b /media/storage/mergespec/data/batchcluster/Linfeng/all.txt \
                            -g /media/storage/mergespec/data/batchcluster/Linfeng/peak_counts.dat \
            // maracluster pvalue -i /media/storage/mergespec/data/batchcluster/Linfeng/unit_test/601.dat \
                            -b /media/storage/mergespec/data/batchcluster/Linfeng/all.txt \
                            -g /media/storage/mergespec/data/batchcluster/Linfeng/unit_test/peak_counts.dat \
                            -t -10.0
            if (peakCountFN_.size() == 0) {
              std::cerr << "Error: no peak counts file specified with -g flag" << std::endl;
              return EXIT_FAILURE;
            }
            if (spectrumBatchFileFN_.size() == 0) {
              std::cerr << "Error: no batch file specified with -b flag" << std::endl;
              return EXIT_FAILURE;
            }
            
            PeakCounts peakCounts;
            peakCounts.readFromFile(peakCountFN_);
            
            SpectrumFileList fileList;
            fileList.initFromFile(spectrumBatchFileFN_);
            
            BatchSpectra spectra(pvaluesFN_);
            if (spectrumInFN_.substr(spectrumInFN_.size() - 3) == "dat") {
              spectra.readBatchSpectra(spectrumInFN_);
            } else {
              spectra.convertToBatchSpectra(spectrumInFN_, fileList);
            }
            
            spectra.calculatePvalueVectors(peakCounts);
            spectra.calculatePvalues();
          }
          return EXIT_SUCCESS;
        }
        case CLUSTER:
        {
          std::vector<std::string> pvalFNs;
          pvalFNs.push_back(matrixFN_);
          
          if (spectrumBatchFileFN_.size() == 0) {
            std::cerr << "Error: no batch file specified with -b flag" << std::endl;
            return EXIT_FAILURE;
          }
          SpectrumFileList fileList;
          fileList.initFromFile(spectrumBatchFileFN_);
          
          return doClustering(outputFolder_, fnPrefix_, pvalFNs, scanNrsFN_, 
              scanDescFN_, clusterThresholds_, fileList, skipFilterAndSort_, 
              matrixFN_, resultTreeFN_);
        }
        case CONSENSUS:
        {
          if (spectrumOutFN_.size() == 0)
            spectrumOutFN_ = outputFolder_ + "/" + fnPrefix_ + ".consensus.ms2";
          MSFileMerger msFileMerger(spectrumOutFN_);
          
          std::cerr << "Parsing cluster file" << std::endl;
          msFileMerger.parseClusterFileForMerge(clusterFileFN_);
          std::cerr << "Finished parsing cluster file" << std::endl;
          
          std::cerr << "Merging clusters" << std::endl;
          msFileMerger.mergeSpectra();
          std::cerr << "Finished merging clusters" << std::endl;
          
          return EXIT_SUCCESS;
        }
        case SEARCH:
        {
          if (spectrumLibraryFN_.size() == 0) {
            std::cerr << "Error: no spectrum library file specified with -z flag" << std::endl;
            return EXIT_FAILURE;
          } else if (peakCountFN_.size() == 0) {
            std::cerr << "Error: no peak counts file specified with -g flag" << std::endl;
            return EXIT_FAILURE;
          } else if (spectrumBatchFileFN_.size() == 0 && spectrumInFN_.size() == 0) {
            std::cerr << "Error: no query spectrum file(s) specified with -i or -b flag" << std::endl;
            return EXIT_FAILURE;
          } else if (spectrumBatchFileFN_.size() != 0 && spectrumInFN_.size() != 0) {
            std::cerr << "Error: ambiguous query spectrum file(s) input, please use only one of the -i and -b flags" << std::endl;
            return EXIT_FAILURE;
          }
          
          if (pvaluesFN_.size() == 0)
            pvaluesFN_ = outputFolder_ + "/" + fnPrefix_ + ".pvalues.dat";
          
          if (!BatchGlobals::fileExists(pvaluesFN_)) {
            std::cerr << "Reading peak counts" << std::endl;
            PeakCounts peakCounts;
            peakCounts.readFromFile(peakCountFN_);
            peakCounts.setSmoothingMode(1);
            std::cerr << "Finished reading peak counts" << std::endl;
            
            // read in the query spectra
            BatchSpectra querySpectra("");
            SpectrumFileList fileList;
            if (spectrumInFN_.size() > 0) {
              querySpectra.convertToBatchSpectra(spectrumInFN_, fileList);
            } else {
              fileList.initFromFile(spectrumBatchFileFN_);
              querySpectra.convertToBatchSpectra(fileList);
            }
            
            // read in the library spectra
            BatchSpectra librarySpectra(pvaluesFN_);
            librarySpectra.convertToBatchSpectra(spectrumLibraryFN_, fileList);
            librarySpectra.calculatePvalueVectors(peakCounts);
            
            librarySpectra.librarySearch(querySpectra);
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
          /*
          if (PvalueCalculator::pvalUnitTest()) {
            std::cerr << "PvalueCalculator unit tests succeeded" << std::endl;
          } else {
            std::cerr << "PvalueCalculator unit tests failed" << std::endl;
            ++failures;
          }
          */
          
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
    }
  }
  catch (std::exception& e)
  {
    std::cerr << "Error: " << e.what() << std::endl;
  }
  catch (...)
  {
    std::cerr << "Caught unknown exception.\n";
  }

  return EXIT_FAILURE;  
  
}
