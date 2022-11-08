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
    peakCountFN_(""), datFNFile_(""), scanInfoFN_(""), addSpecIds_(false),
    pvaluesFN_(""), clusterFileFN_(""), pvalVecInFileFN_(""), 
    pvalueVectorsBaseFN_(""), overlapBatchFileFN_(""), overlapBatchIdx_(0u), 
    spectrumBatchFileFN_(""), spectrumInFN_(""), spectrumOutFN_(""),
    spectrumLibraryFN_(""), matrixFN_(""), resultTreeFN_(""),
    skipFilterAndSort_(false), writeAll_(false), precursorTolerance_(20),
    precursorToleranceDa_(false), dbPvalThreshold_(-5.0), 
    chargeUncertainty_(0), minConsensusClusterSize_(1u)
{
  boost::filesystem::path outputPath = boost::filesystem::current_path() / boost::filesystem::path("maracluster_output");
  outputFolder_ = outputPath.string();
  
  boost::filesystem::path outputDatPath = outputPath / boost::filesystem::path("dat_files");
  datFolder_ = outputDatPath.string();
}

MaRaCluster::~MaRaCluster() {}

std::string MaRaCluster::greeter() {
  std::ostringstream oss;
  oss << "MaRaCluster version " << VERSION << ", ";
  oss << "Build Date " << __DATE__ << " " << __TIME__ << std::endl;
  oss << "Copyright (c) 2015-19 Matthew The. All rights reserved.\n"
      << "Written by Matthew The (matthewt@kth.se) in the\n"
      << "School of Biotechnology at the Royal Institute of Technology in Stockholm.\n";
  return oss.str();
}

std::string MaRaCluster::extendedGreeter(time_t& startTime) {
  std::ostringstream oss;
  char* host = std::getenv("HOSTNAME");
  oss << greeter();
  oss << "Issued command:" << std::endl << call_ << std::endl;
  oss << "Started " << ctime(&startTime) << std::endl;
  oss.seekp(-1, std::ios_base::cur);
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
      "File with spectrum files to be processed in batch, one per line. Files should be readable by ProteoWizard (e.g. ms2, mgf, mzML).",
      "filename");
  cmd.defineOption("f",
      "output-folder",
      "Writable folder for output files (default: ./maracluster_output).",
      "path");
  cmd.defineOption("p",
      "precursorTolerance",
      "Set precursor tolerance in units of ppm or Da. The units have to be \"Da\" or \"ppm\", case sensitive; if no unit is specified ppm is assumed (default: 20.0ppm).",
      "string");
  cmd.defineOption("t",
      "pvalThreshold",
      "Set log(p-value) threshold (default: -5.0).",
      "double");
  cmd.defineOption("c",
      "clusterThresholds",
      "Clustering thresholds at which to produce cluster files; listed as a comma separated list (default: -30.0,-25.0,-20.0,-15.0,-10.0,-5.0)",
      "string");
  cmd.defineOption("v",
      "verbatim",
      "Set the verbatim level (lowest: 0, highest: 5, default: 3).",
      "int");
  cmd.defineOption("D",
      "dat-folder",
      "Writable folder for converted .dat binary files. Can be used to re-use already converted spectrum files (default: ./maracluster_output/dat_files).",
      "path");
  cmd.defineOption("I",
      "addSpecIds",
      "Add column with spectrum id/title to clustering file output.",
      "",
      TRUE_IF_SET);
  cmd.defineOption("a",
      "prefix",
      "Output files will be prefixed as e.g. <prefix>.clusters_p10.tsv (default: 'MaRaCluster')",
      "name");
  cmd.defineOption("l",
      "clusterFile",
      "Input file for generating consensus spectra containing filepaths and scan numbers, separated by tabs.",
      "filename");
  cmd.defineOption("o",
      "specOut",
      "Output file for the consensus spectra. Can be in any format supported by ProteoWizard (e.g. ms2, mzML).",
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
  cmd.defineOption("i",
      "specIn",
      "Input file readable by ProteoWizard (e.g. ms2, mzML). For multiple input files use the -b/--batch option instead.",
      "filename");
  cmd.defineOption("y",
      "pvalueVecFile",
      "Input file with pvalue vectors",
      "filename");
  cmd.defineOption("r",
      "pvecOut",
      "Output file basename for p-values vectors.",
      "filename");
  cmd.defineOption("w",
      "overlapBatch",
      "File with 2 tab separated columns as: tail_file <tab> head_file for"
      " which overlapping p-values should be calculated",
      "filename");
  cmd.defineOption("W",
      "overlapBatchIdx",
      "Index of overlap to process, requires -j/--datFNfile to be specified.",
      "filename");
  cmd.defineOption("q",
      "pvalOut",
      "File where p-values will be written to.",
      "filename");
  cmd.defineOption("d",
      "percOut",
      "Tab delimited percolator output file containing peptides and qvalues. This is meant for annotation of the clusterfile.",
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
  cmd.defineOption("e",
      "skipFilterAndSort",
      "Skips filtering and sorting of the input matrix, only use if the input is a filtered and sorted binary list p-values.",
      "",
      TRUE_IF_SET);
  cmd.defineOption("C",
      "chargeUncertainty",
      "Set charge uncertainty, i.e. if set to 1, then for a spectrum with precursor ion charge C, also precursor ion charges C-1 and C+1 are considered (default: 0).",
      "int");
  cmd.defineOption("z",
      "lib",
      "File readable by ProteoWizard (e.g. ms2, mzML) with spectral library",
      "filename");
      
  // finally parse and handle return codes (display help etc...)
  cmd.parseArgs(argc, argv);
  
  // there should be one argument left specifying the mode
  if (cmd.arguments.size() == 1) {
    std::string mode = cmd.arguments[0];
    if (mode == "batch") mode_ = BATCH;
    else if (mode == "pvalue") mode_ = PVALUE;
    else if (mode == "overlap") mode_ = OVERLAP;
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
  if (cmd.optionSet("batch")) spectrumBatchFileFN_ = cmd.options["batch"];
  if (cmd.optionSet("output-folder")) {
    outputFolder_ = cmd.options["output-folder"];
    datFolder_ = outputFolder_ + "/dat_files";
  }
  if (cmd.optionSet("dat-folder")) datFolder_ = cmd.options["dat-folder"];
  if (cmd.optionSet("prefix")) fnPrefix_ = cmd.options["prefix"];
  
  // file output for maracluster batch and index (also for some other methods)
  if (cmd.optionSet("datFNfile")) datFNFile_ = cmd.options["datFNfile"];
  if (cmd.optionSet("peakCountsFN")) peakCountFN_ = cmd.options["peakCountsFN"];
  if (cmd.optionSet("scanInfoFN")) scanInfoFN_ = cmd.options["scanInfoFN"];
  if (cmd.optionSet("addSpecIds")) addSpecIds_ = true;
  
  // file input options for maracluster pvalue
  if (cmd.optionSet("specIn")) spectrumInFN_ = cmd.options["specIn"];
  if (cmd.optionSet("clusterFile")) clusterFileFN_ = cmd.options["clusterFile"];
  if (cmd.optionSet("overlapBatch")) overlapBatchFileFN_ = cmd.options["overlapBatch"];
  if (cmd.optionSet("overlapBatchIdx")) overlapBatchIdx_ = cmd.getInt("overlapBatchIdx", 0, 100000);
  if (cmd.optionSet("pvalueVecFile")) pvalVecInFileFN_ = cmd.options["pvalueVecFile"];
  
  // file output options for maracluster pvalue
  if (cmd.optionSet("pvalOut")) pvaluesFN_ = cmd.options["pvalOut"];
  if (cmd.optionSet("pvecOut")) pvalueVectorsBaseFN_ = cmd.options["pvecOut"];
  
  // file input options for maracluster cluster
  if (cmd.optionSet("clusteringMatrix")) matrixFN_ = cmd.options["clusteringMatrix"];
  if (cmd.optionSet("clusteringTree")) resultTreeFN_ = cmd.options["clusteringTree"];
  if (cmd.optionSet("skipFilterAndSort")) skipFilterAndSort_ = true;
  
  // file output option for maracluster consensus
  if (cmd.optionSet("specOut")) spectrumOutFN_ = cmd.options["specOut"];
  if (cmd.optionSet("minClusterSize")) minConsensusClusterSize_ = cmd.getInt("minClusterSize", 0, 100000);
  if (cmd.optionSet("splitMassChargeStates")) MSFileHandler::splitMassChargeStates_ = true;
  
  // file input option for maracluster search
  if (cmd.optionSet("lib")) spectrumLibraryFN_ = cmd.options["lib"];
  
  // general options
  if (cmd.optionSet("pvalThreshold")) dbPvalThreshold_ = cmd.getDouble("pvalThreshold", -1000.0, 0.0);
  if (cmd.optionSet("clusterThresholds")) {
    std::istringstream ss(cmd.options["clusterThresholds"]);
    std::string token;
    while(std::getline(ss, token, ',')) {
      clusterThresholds_.push_back(atof(token.c_str()));
    }
    if (clusterThresholds_.empty()) {
      std::cerr << "Error: invalid input for clusterThreshold parameter: " << 
                    cmd.options["clusterThresholds"] << std::endl;
      return false;
    }
  } else {
    if (dbPvalThreshold_ > -30.0) clusterThresholds_.push_back(-30.0);
    if (dbPvalThreshold_ > -25.0) clusterThresholds_.push_back(-25.0);
    if (dbPvalThreshold_ > -20.0) clusterThresholds_.push_back(-20.0);
    if (dbPvalThreshold_ > -15.0) clusterThresholds_.push_back(-15.0);
    if (dbPvalThreshold_ > -10.0) clusterThresholds_.push_back(-10.0);
    if (dbPvalThreshold_ > -5.0) clusterThresholds_.push_back(-5.0);
    if (std::find(clusterThresholds_.begin(), clusterThresholds_.end(), dbPvalThreshold_) == clusterThresholds_.end()) {
      clusterThresholds_.push_back(dbPvalThreshold_);
    }
  }
  std::sort(clusterThresholds_.begin(), clusterThresholds_.end());
  
  if (cmd.optionSet("precursorTolerance")) {
    std::string precursorToleranceString = cmd.options["precursorTolerance"];
    size_t unitIdx = precursorToleranceString.find("Da");
    if (unitIdx != std::string::npos) {
      precursorTolerance_ = atof(precursorToleranceString.substr(0, unitIdx).c_str());
      precursorToleranceDa_ = true;
    } else {
      unitIdx = precursorToleranceString.find("ppm");
      precursorTolerance_ = atof(precursorToleranceString.substr(0, unitIdx).c_str());
    }
  }
  if (cmd.optionSet("chargeUncertainty")) chargeUncertainty_ = cmd.getInt("chargeUncertainty", 0, 5);
  if (cmd.optionSet("verbatim")) Globals::VERB = cmd.getInt("verbatim", 0, 5);

  return true;
}

int MaRaCluster::createIndex() {
  if (spectrumBatchFileFN_.empty()) {
    std::cerr << "Error: no batch file specified with -b/--batch flag" << std::endl;
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
    SpectrumFiles spectrumFiles(outputFolder_, datFolder_, chargeUncertainty_, addSpecIds_);
    spectrumFiles.convertToDat(fileList);
    spectrumFiles.splitByPrecursorMz(fileList, datFNFile_, peakCountFN_, 
        scanInfoFN_, precursorTolerance_, precursorToleranceDa_);
  } else {
    std::cerr << "Read dat-files from " << datFNFile_ << 
        " and scan numbers from " << scanInfoFN_ <<
        ". Remove these files to generate new dat-files." << std::endl;
  }
  
  return EXIT_SUCCESS;
}

int MaRaCluster::processDatFiles(const std::vector<std::string>& datFNs) {
  for (size_t i = 0; i < datFNs.size(); ++i) {
    // make sure the file exists
    if (!Globals::fileExists(datFNs[i])) {
      std::cerr << "Ignoring missing data file " << datFNs[i] << std::endl;
      continue;
    }

    std::string datFN = datFNs[i];
    std::string pvalueVectorsBaseFN = datFN + ".pvalue_vectors";
    std::string pvaluesFN = datFN + ".pvalues.dat"; // only contains pvalues from overlap regions
    std::string pvalueTreeFN = datFN + ".pvalue_tree.tsv";
    
    int error = clusterSpectra(datFN, pvaluesFN, pvalueVectorsBaseFN, pvalueTreeFN);
    if (error != EXIT_SUCCESS) return error;
  }
  
  return EXIT_SUCCESS;
}

int MaRaCluster::clusterSpectra(const std::string& spectrumInFN, 
    const std::string& pvaluesFN, const std::string& pvalueVectorsBaseFN, 
    const std::string& pvalueTreeFN) {  
  if (!Globals::fileExists(pvalueTreeFN)) {
    PvalueVectors pvecs(pvaluesFN, precursorTolerance_, precursorToleranceDa_, dbPvalThreshold_);
    {
      Spectra spectra;
      if (spectrumInFN.substr(spectrumInFN.size() - 3) == "dat") {
        spectra.readBatchSpectra(spectrumInFN);
      } else {
        SpectrumFileList fileList;
        fileList.initFromFile(spectrumBatchFileFN_);
        spectra.convertToBatchSpectra(spectrumInFN, fileList);
      }
      spectra.sortSpectraByPrecMz();
      
      PeakCounts peakCounts;
      peakCounts.readFromFile(peakCountFN_);
      pvecs.calculatePvalueVectors(spectra.getSpectra(), peakCounts);
    }
    
    if (pvalueVectorsBaseFN.size() > 0) {
      pvecs.writePvalueVectors(pvalueVectorsBaseFN, writeAll_);
    }
    
    // write to temporary file, rename to real filename when finished to prevent using incomplete results
    std::string pvalueTreeTmpFN = pvalueTreeFN + ".tmp"; 
    if (pvalueTreeFN.size() > 0) {
      pvecs.batchCalculateAndClusterPvalues(pvalueTreeTmpFN, scanInfoFN_);
    } else {
      pvecs.batchCalculatePvalues();
    }
    rename(pvalueTreeTmpFN.c_str(), pvalueTreeFN.c_str());
  } else {
    std::cerr << "Using p-value tree from " << pvalueTreeFN <<
        ". Remove this file to generate a new p-value tree." << std::endl;
  }
  
  return EXIT_SUCCESS;
}

void MaRaCluster::getOverlapBatches(const std::vector<std::string>& datFNs, 
    std::vector<std::string>& pvalTreeFNs, 
    std::vector<OverlapBatch>& overlapBatches) {
  for (size_t i = 0; i < datFNs.size(); ++i) {
    // make sure the file exists
    if (!Globals::fileExists(datFNs[i])) {
      std::cerr << "Ignoring missing data file " << datFNs[i] << std::endl;
      continue;
    }

    std::string datFN = datFNs[i];
    std::string pvalueVectorsBaseFN = datFN + ".pvalue_vectors";
    std::string pvaluesFN = datFN + ".pvalues.dat"; // only contains pvalues from overlap regions
    std::string pvalueTreeFN = datFN + ".pvalue_tree.tsv";
    
    overlapBatches[i].pvalFNs.push_back(pvaluesFN + ".head.dat");
    overlapBatches[i].headPvecFile = pvalueVectorsBaseFN + ".head.dat";
    
    overlapBatches[i+1].pvalFNs.push_back(pvaluesFN + ".tail.dat");
    overlapBatches[i+1].tailPvecFile = pvalueVectorsBaseFN + ".tail.dat";
    
    pvalTreeFNs.push_back(pvalueTreeFN);
  }
}

int MaRaCluster::processOverlapBatches(
    std::vector<OverlapBatch>& overlapBatches,
    std::vector<std::string>& pvalTreeFNs, const unsigned int mergeOffset) {
  for (size_t overlapIdx = 0; overlapIdx < overlapBatches.size(); ++overlapIdx) {
    int error = processOverlapBatch(overlapBatches[overlapIdx], overlapIdx, mergeOffset);
    if (error != EXIT_SUCCESS) return error;
    
    std::string resultTreeFN = outputFolder_ + "/overlap." + boost::lexical_cast<std::string>(overlapIdx) + ".pvalue_tree.tsv";
    pvalTreeFNs.push_back(resultTreeFN);
  }
  
  return EXIT_SUCCESS;
}

int MaRaCluster::processOverlapBatch(OverlapBatch& overlapBatch, 
    const unsigned int overlapIdx, const unsigned int mergeOffset) {
  std::string pvaluesFN = outputFolder_ + "/overlap." + boost::lexical_cast<std::string>(overlapIdx) + ".pvalues.dat";
  std::string resultTreeFN = outputFolder_ + "/overlap." + boost::lexical_cast<std::string>(overlapIdx) + ".pvalue_tree.tsv";
  std::string matrixFN = outputFolder_ + "/poisoned." + boost::lexical_cast<std::string>(overlapIdx) + ".pvalues.dat";
  if (!overlapBatch.tailPvecFile.empty() && !overlapBatch.headPvecFile.empty()) {
    if (!Globals::fileExists(pvaluesFN)) {
      PvalueVectors pvecs(pvaluesFN, precursorTolerance_, precursorToleranceDa_, dbPvalThreshold_);
      
      pvecs.batchCalculatePvaluesOverlap(overlapBatch.tailPvecFile, overlapBatch.headPvecFile);
      pvecs.clearPvalueVectors();
      remove(overlapBatch.tailPvecFile.c_str());
      remove(overlapBatch.headPvecFile.c_str());
    } else {
      std::cerr << "Using p-values from " << pvaluesFN << 
          ". Remove this file to generate new p-values." << std::endl;
    }
    
    if (Globals::fileExists(pvaluesFN)) {
      overlapBatch.pvalFNs.push_back(pvaluesFN);
    }
  }
  
  if (Globals::VERB > 2) {
    std::cerr << "Clustering overlap batch " << overlapIdx << ":" << std::endl;
    BOOST_FOREACH (const std::string& s, overlapBatch.pvalFNs) {
      std::cerr << "  " << s << std::endl;
    }
  }
  
  int error = doClustering(overlapBatch.pvalFNs, resultTreeFN, matrixFN, mergeOffset);
  return error;
}

int MaRaCluster::doClustering(const std::vector<std::string> pvalFNs, 
    std::string& resultTreeFN, const std::string& matrixFN, 
    const unsigned int mergeOffset) {
  // input checks
  if (matrixFN.empty() && resultTreeFN.empty()) {
    std::cerr << "Error: no input file specified with -m/--clusteringMatrix or -u/--clusteringTree flag" << std::endl;
    return EXIT_FAILURE;
  }

  // create output file paths
  if (resultTreeFN.empty()) {
    resultTreeFN = outputFolder_ + "/overlap.pvalue_tree.tsv";
  }
  
  // start clustering
  if (!Globals::fileExists(resultTreeFN)) {
    std::cerr << "Starting p-value clustering." << std::endl;
    
    if (!Globals::fileExists(matrixFN)) {
      bool tsvInput = false;
      PvalueFilterAndSort::filterAndSort(pvalFNs, matrixFN, tsvInput);
    } else {
      std::cerr << "Using p-values from " << matrixFN << 
          " . Remove this file to re-sort and filter the p-values." << std::endl;
    }
    
    SparseClustering matrix;
    matrix.setMergeOffset(mergeOffset);
    matrix.initMatrix(matrixFN);
    matrix.setClusterPairFN(resultTreeFN);
    matrix.doClustering(std::min(dbPvalThreshold_, clusterThresholds_.back()));
    remove(matrixFN.c_str());
  } else {
    std::cerr << "Previous clustering results are available in " << 
        resultTreeFN << " . Remove this file to redo the clustering." << std::endl;
  }
  
  return EXIT_SUCCESS;
}

int MaRaCluster::mergeSpectra() {
  if (spectrumOutFN_.empty())
    spectrumOutFN_ = outputFolder_ + "/" + fnPrefix_ + ".consensus.ms2";
  
  if (!MSFileHandler::validMs2OutputFN(spectrumOutFN_)) {
    return EXIT_FAILURE;
  }
  
  if (MSFileHandler::getOutputFormat(spectrumOutFN_) == "mgf") {
    MSFileHandler::splitMassChargeStates_ = true;
  }
  
  if (clusterFileFN_.empty() || !Globals::fileExists(clusterFileFN_)) {
    std::cerr << "Error: Could not find cluster input file (-l/--clusterFile flag) " << clusterFileFN_ << std::endl;
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

// adapted from https://github.com/crux-toolkit/crux-toolkit/blob/master/src/util/crux-utils.cpp
bool MaRaCluster::parseUrl(std::string url, std::string* host, std::string* path) {
  if (!host || !path) {
    return false;
  }
  // find protocol
  size_t protocolSuffix = url.find("://");
  if (protocolSuffix != std::string::npos) {
    url = url.substr(protocolSuffix + 3);
  }
  size_t pathBegin = url.find('/');
  if (pathBegin == std::string::npos) {
    *host = url;
    *path = "/";
  } else {
    *host = url.substr(0, pathBegin);
    *path = url.substr(pathBegin);
  }
  if (host->empty()) {
    *host = *path = "";
    return false;
  }
  return true;
}

void MaRaCluster::httpRequest(const std::string& url, const std::string& data) {
  // Parse URL into host and path components
  std::string host, path;
  if (!parseUrl(url, &host, &path)) {
    if (Globals::VERB > 2) {
      std::cerr << "Warning: Failed parsing URL " << url << std::endl;
    }
    return;
  }

  using namespace boost::asio;

  // Establish TCP connection to host on port 80
  io_service service;
  ip::tcp::resolver resolver(service);
  ip::tcp::resolver::iterator endpoint = resolver.resolve(ip::tcp::resolver::query(host, "80"));
  ip::tcp::socket sock(service);
  connect(sock, endpoint);
  
  std::size_t seed = 0;
  boost::hash_combine(seed, ip::host_name());
  boost::hash_combine(seed, sock.local_endpoint().address().to_string());
  std::stringstream stream;
  stream << std::hex << seed;
  
  std::string placeholder = "CID_PLACEHOLDER";
  std::string cid = stream.str();
  
  std::string newData(data);
  
  if (Globals::VERB > 3) {
    std::cerr << "Analytics data string: " << newData << std::endl;
  }
  
  newData.replace(newData.find(placeholder), placeholder.length(), cid);
  
  // Determine method (GET if no data; otherwise POST)
  std::string method = newData.empty() ? "GET" : "POST";
  std::ostringstream lengthString;
  lengthString << newData.length();
  
  std::string contentLengthHeader = newData.empty()
    ? ""
    : "Content-Length: " + lengthString.str() + "\r\n";
  // Send the HTTP request
  std::string request =
    method + " " + path + " HTTP/1.1\r\n"
    "Host: " + host + "\r\n" +
    contentLengthHeader +
    "Connection: close\r\n"
    "\r\n" + newData;
  sock.send(buffer(request));
}

void MaRaCluster::postToAnalytics(const std::string& appName) {
  // Post data to Google Analytics
  // For more information, see: https://developers.google.com/analytics/devguides/collection/protocol/v1/devguide
  try {
    std::stringstream paramBuilder;
    paramBuilder << "v=1"                // Protocol verison
                 << "&tid=UA-165948942-3" // Tracking ID
                 << "&cid=CID_PLACEHOLDER" // Unique device ID
                 << "&t=event"           // Hit type
                 << "&ec=maracluster"     // Event category
                 << "&ea=" << appName    // Event action
                 << "&el="               // Event label
#ifdef _MSC_VER
                      "win"
#elif __APPLE__
                      "mac"
#else
                      "linux"
#endif
                   << '-' << VERSION;
    httpRequest(
      "http://www.google-analytics.com/collect",
      paramBuilder.str());
  } catch (...) {
  }
}


int MaRaCluster::run() {
  time_t startTime;
  clock_t startClock;
  time(&startTime);
  startClock = clock();
  
  if (Globals::VERB > 0) {
    std::cerr << extendedGreeter(startTime);
  }
  
  std::string appName = "maracluster";
  postToAnalytics(appName);
  
  if (!spectrumOutFN_.empty()) {
    boost::filesystem::path spectrumOutPath (spectrumOutFN_);
    boost::filesystem::path spectrumOutFolder = spectrumOutPath.parent_path();
    if (!spectrumOutFolder.empty() && !boost::filesystem::exists(spectrumOutFolder)) {
      boost::system::error_code returnedError;
      boost::filesystem::create_directories( spectrumOutFolder, returnedError );
      
      if (!boost::filesystem::exists(spectrumOutFolder)) {
        std::cerr << "Error: could not create output directory at " << outputFolder_ << std::endl;
        return EXIT_FAILURE;
      }
    }
  }
  
  boost::filesystem::path rootPath (outputFolder_);
  if (!boost::filesystem::exists(rootPath)) {
    boost::system::error_code returnedError;
    boost::filesystem::create_directories( rootPath, returnedError );
    
    if (!boost::filesystem::exists(rootPath)) {
      std::cerr << "Error: could not create output directory at " << outputFolder_ << std::endl;
      return EXIT_FAILURE;
    }
  }
  
  switch (mode_) {
    case BATCH: {
      // This executes the entire pipeline in one go
      // maracluster -b /media/storage/mergespec/data/Linfeng/batchcluster/all.txt \
                  -f /media/storage/mergespec/data/Linfeng/batchcluster/output/
      
      /*************************************************************************
        Step 1: Read in spectra and divide them in batches by precursor m/z
       ************************************************************************/
       
      int error = createIndex();
      if (error != EXIT_SUCCESS) return error;      
            
      /*************************************************************************
        Step 2: Process each of the batches and cluster non-overlapping regions
       ************************************************************************/
      
      std::vector<std::string> datFNs;
      SpectrumFiles::readDatFNsFromFile(datFNFile_, datFNs);
      
      if (datFNs.empty()) {
        std::cerr << "Error: could not find any ms2 spectra in the input files." << std::endl;
        return EXIT_FAILURE;
      }
      
      error = processDatFiles(datFNs);
      if (error != EXIT_SUCCESS) return error;
      
      /*************************************************************************
        Step 3: Process each of the overlapping regions
       ************************************************************************/
      
      SpectrumFileList fileList;
      fileList.initFromFile(spectrumBatchFileFN_);
      unsigned int mergeOffset = fileList.getMergeOffset();
      
      std::vector<std::string> pvalTreeFNs;
      std::vector<OverlapBatch> overlapBatches(datFNs.size() + 1);
      getOverlapBatches(datFNs, pvalTreeFNs, overlapBatches);
      
      error = processOverlapBatches(overlapBatches, pvalTreeFNs, mergeOffset);
      if (error != EXIT_SUCCESS) return error;
      
      /*************************************************************************
        Step 4: Combine all p-value trees into one large tree and write clusters
       ************************************************************************/
       
      SpectrumClusters clustering;
      std::string clusterBaseFN = outputFolder_ + "/" + fnPrefix_ + ".clusters_";
      clustering.printClusters(pvalTreeFNs, clusterThresholds_, fileList, 
                               scanInfoFN_, clusterBaseFN);
      
      /*************************************************************************
        Step 5: Create and write consensus spectra
       ************************************************************************/
       
      if (!spectrumOutFN_.empty()) {
        if (clusterFileFN_.empty()) {
          clusterFileFN_ = SpectrumClusters::getClusterFN(clusterBaseFN, 
                                                          dbPvalThreshold_);
        }
        if (Globals::VERB > 1) {
          std::cerr << "Creating consensus spectra using clusters generated in "
                    << clusterFileFN_ << std::endl;
        }
        error = mergeSpectra();
        if (error != EXIT_SUCCESS) return EXIT_FAILURE;
      }
      
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
      
      if (pvalueVectorsBaseFN_.empty())
        pvalueVectorsBaseFN_ = outputFolder_ + "/" + fnPrefix_ + ".pvalue_vectors";
      
      if (pvalVecInFileFN_.size() > 0) { 
        // direct input of p-value vector file
        // maracluster pvalue -y /media/storage/mergespec/data/batchtest/1300.ms2.pvalue_vectors.tsv
        
        PvalueVectors pvecs(pvaluesFN_, precursorTolerance_, precursorToleranceDa_, dbPvalThreshold_);
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
        PvalueVectors pvecs(pvaluesFN_, precursorTolerance_, precursorToleranceDa_, dbPvalThreshold_);
        std::vector< std::pair<std::string, std::string> > overlapFNs;
        pvecs.parseBatchOverlapFile(overlapBatchFileFN_, overlapFNs);
        pvecs.processOverlapFiles(overlapFNs);
      } else if (clusterFileFN_.size() > 0) {
        // calculate p-values from a cluster in a scan description list
        // maracluster pvalue -l <scan_desc_file> -g <peak_counts_file>
        if (peakCountFN_.empty()) {
          std::cerr << "Error: no peak counts file specified with -g/--peakCountsFN flag" << std::endl;
          return EXIT_FAILURE;
        }
        std::string spectrumOutFN = "";
        MSFileExtractor fileExtractor(spectrumOutFN);
        
        std::vector<Spectrum> batchSpectra;
        fileExtractor.parseClusterFileForExtract(clusterFileFN_);
        fileExtractor.extractToBatchSpectrumList(batchSpectra);
        
        PeakCounts peakCounts;
        peakCounts.readFromFile(peakCountFN_);
        
        std::cerr << "Read " << batchSpectra.size() << " spectra." << std::endl;
        
        Spectra spectra;
        spectra.setBatchSpectra(batchSpectra);
        spectra.sortSpectraByPrecMz();
        
        PvalueVectors pvecs(pvaluesFN_, precursorTolerance_, precursorToleranceDa_, dbPvalThreshold_);
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
          std::cerr << "Error: no peak counts file specified with -g/--peakCountsFN flag" << std::endl;
          return EXIT_FAILURE;
        }
        if (spectrumBatchFileFN_.empty() && spectrumInFN_.empty()) {
          std::cerr << "Error: no input file specified with -i/--specIn or -b/--batch flag" << std::endl;
          return EXIT_FAILURE;
        }
        
        clusterSpectra(spectrumInFN_, pvaluesFN_, pvalueVectorsBaseFN_, resultTreeFN_);
      }
      return EXIT_SUCCESS;
    }
    case OVERLAP:
    {
      std::vector<std::string> datFNs;
      SpectrumFiles::readDatFNsFromFile(datFNFile_, datFNs);
      
      SpectrumFileList fileList;
      fileList.initFromFile(spectrumBatchFileFN_);
      unsigned int mergeOffset = fileList.getMergeOffset();
      
      std::vector<std::string> pvalTreeFNs;
      std::vector<OverlapBatch> overlapBatches(datFNs.size() + 1);
      getOverlapBatches(datFNs, pvalTreeFNs, overlapBatches);
      
      int error = processOverlapBatch(overlapBatches[overlapBatchIdx_], overlapBatchIdx_, mergeOffset);
      return error;
    }
    case CLUSTER:
    {
      std::vector<std::string> pvalFNs, pvalTreeFNs;
      pvalFNs.push_back(matrixFN_);
      
      if (spectrumBatchFileFN_.empty()) {
        std::cerr << "Error: no batch file specified with -b/--batch flag" << std::endl;
        return EXIT_FAILURE;
      }
      if (scanInfoFN_.empty()) {
        std::cerr << "Error: no scannrs file specified with -s/--scanInfoFN flag" << std::endl;
        return EXIT_FAILURE;
      }
      
      if (!skipFilterAndSort_) {
        bool tsvInput = false;
        PvalueFilterAndSort::filterAndSort(pvalFNs, matrixFN_, tsvInput);
      }
      
      SpectrumFileList fileList;
      fileList.initFromFile(spectrumBatchFileFN_);
      int error = doClustering(pvalFNs, resultTreeFN_, matrixFN_, fileList.getMergeOffset());
      if (error != EXIT_SUCCESS) return EXIT_FAILURE;
      
      pvalTreeFNs.push_back(resultTreeFN_);
      
      SpectrumClusters clustering;
      std::string clusterBaseFN = outputFolder_ + "/" + fnPrefix_ + ".clusters_";
      clustering.printClusters(pvalTreeFNs, clusterThresholds_, fileList, scanInfoFN_, clusterBaseFN);
      
      return EXIT_SUCCESS;
    }
    case CONSENSUS:
    {
      return mergeSpectra();
    }
    case SEARCH:
    {
      if (spectrumLibraryFN_.empty()) {
        std::cerr << "Error: no spectrum library file specified with -z/--lib flag" << std::endl;
        return EXIT_FAILURE;
      } else if (spectrumBatchFileFN_.empty() && spectrumInFN_.empty()) {
        std::cerr << "Error: no query spectrum file(s) specified with -i/--specIn or -b/--batch flag" << std::endl;
        return EXIT_FAILURE;
      } else if (!spectrumBatchFileFN_.empty() && !spectrumInFN_.empty()) {
        std::cerr << "Error: ambiguous query spectrum file(s) input, please use only one of the -i/--specIn and -b/--batch flags" << std::endl;
        return EXIT_FAILURE;
      }
      
      if (pvaluesFN_.empty())
        pvaluesFN_ = outputFolder_ + "/" + fnPrefix_ + ".pvalues.dat";
      
      if (!Globals::fileExists(pvaluesFN_)) {
        // read in the query spectra
        Spectra querySpectra;
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
          
          SpectrumFiles spectrumFiles(outputFolder_, datFolder_);
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
        Spectra librarySpectra;
        librarySpectra.convertToBatchSpectra(spectrumLibraryFN_, fileList);
        librarySpectra.sortSpectraByPrecMz();
        
        PvalueVectors pvecs(pvaluesFN_, precursorTolerance_, precursorToleranceDa_, dbPvalThreshold_);
        pvecs.calculatePvalueVectors(librarySpectra.getSpectra(), peakCounts);
        if (pvalueVectorsBaseFN_.size() > 0) {
          writeAll_ = true;
          pvecs.writePvalueVectors(pvalueVectorsBaseFN_, writeAll_);
        }
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
        Spectra querySpectra;
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
        Spectra librarySpectra;
        librarySpectra.convertToBatchSpectra(spectrumLibraryFN_, fileList);
        librarySpectra.sortSpectraByPrecMz();
        
        PvalueVectors pvecs(pvaluesFN_, precursorTolerance_, precursorToleranceDa_, dbPvalThreshold_);
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
      if (SpectrumFiles::limitsUnitTest()) {
        std::cerr << "SpectrumFiles limits unit tests succeeded" << std::endl;
      } else {
        std::cerr << "SpectrumFiles limits unit tests failed" << std::endl;
        ++failures;
      }
      */
      if (MSClusterMerge::mergeUnitTest()) {
        std::cerr << "Consensus spectra unit tests succeeded" << std::endl;
      } else {
        std::cerr << "Consensus spectra unit tests failed" << std::endl;
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
