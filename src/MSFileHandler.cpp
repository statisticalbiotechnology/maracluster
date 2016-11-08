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
 
#include "MSFileHandler.h"

using pwiz::msdata::MSData;
using pwiz::msdata::MSDataFile;
using pwiz::msdata::SpectrumListSimplePtr;
using pwiz::msdata::SpectrumListSimple;
using pwiz::msdata::SpectrumListPtr;
using pwiz::msdata::SpectrumPtr;
using pwiz::msdata::Spectrum;
using pwiz::msdata::SelectedIon;
using pwiz::msdata::Software;
using pwiz::msdata::SoftwarePtr;
using pwiz::msdata::DataProcessingPtr;
using pwiz::msdata::DataProcessing;
using pwiz::msdata::ProcessingMethod;

unsigned int MSFileHandler::chargeErrorTolerance_ = 0u;
bool MSFileHandler::splitMassChargeStates_ = false;

bool MSFileHandler::validMs2OutputFN(std::string& outputFN) {
  std::string outputFormat = getOutputFormat(outputFN);
  if (outputFormat != "ms2" && outputFormat != "mzml" && outputFormat != "mzxml" && outputFormat != "mgf") {
    std::cerr << "Unknown output format:" << outputFormat << "; valid extensions are mgf, ms2, mzml and mzxml." << std::endl;
    return false;
  } else {
    return true;
  }
}

std::string MSFileHandler::getOutputFormat(const std::string& outputFN) {
  size_t found = outputFN.find_last_of(".");
  std::string outputFormat = outputFN.substr(found+1);
  std::transform(outputFormat.begin(), outputFormat.end(), outputFormat.begin(), ::tolower);
  return outputFormat;
}

void MSFileHandler::writeMSData(MSData& msd, const std::string& outputFN) {
  std::cerr << "Writing file " << outputFN << std::endl;
  
  msd.cvs = pwiz::msdata::defaultCVList();  
  if (msd.id.empty()) {
    msd.id = msd.run.id = "MaRacluster_consensus_spectra";
  }
  
  std::string outputFormat = getOutputFormat(outputFN);
  /* see http://proteowizard.sourceforge.net/dox/_m_s_data_file_8hpp.html */
  if (outputFormat == "mgf") {
    MSDataFile::write(msd, outputFN, MSDataFile::Format_MGF); 
  } else if (outputFormat == "ms2") {
    MSDataFile::write(msd, outputFN, MSDataFile::Format_MS2); 
  } else if (outputFormat == "mzml") {
    MSDataFile::write(msd, outputFN, MSDataFile::Format_mzML); 
  } else if (outputFormat == "mzxml") {
    MSDataFile::write(msd, outputFN, MSDataFile::Format_mzXML); 
  } else {
    std::cerr << "ERROR: Could not write spectra, unknown output format: " << outputFormat << "." << std::endl;
  }
  
  std::cerr << "Finished writing file " << outputFN << std::endl;
}

void MSFileHandler::msgfFixMzML(const std::string& spectrumInFN) {
  MSReaderList readerList;
  MSDataFile msd(spectrumInFN, &readerList);
  SpectrumListPtr sl = msd.run.spectrumListPtr;
  
  SpectrumListSimplePtr slNew(new SpectrumListSimple);
  
  for (unsigned int j = 0; j < sl->size(); ++j) {
    if ((j+1) % 5000 == 0) {
      std::cerr << "Fixing spectra " << j+1 << 
                   "/" << sl->size() << std::endl;
    }
    SpectrumPtr s = sl->spectrum(j, true);
    std::vector<MassChargeCandidate> mccs;
    SpectrumHandler::getMassChargeCandidates(s, mccs);
    unsigned int scannr = SpectrumHandler::getScannr(s);
    addSpectrumWithMccs(s, mccs, scannr, slNew);
  }
  
  size_t idx = 0;
  BOOST_FOREACH (SpectrumPtr& s, slNew->spectra) {
    s->index = idx++;
  }
  
  MSDataFile msdNew(spectrumInFN);
  msdNew.run.spectrumListPtr = slNew;
  writeMSData(msdNew, spectrumOutFN_); 
}

void MSFileHandler::addSpectrumWithMccs(SpectrumPtr consensusSpec, 
     std::vector<MassChargeCandidate>& consensusMccs,
     unsigned int scannr, SpectrumListSimplePtr mergedSpectra) {
  if (splitMassChargeStates_) {
    consensusSpec->precursors.at(0).selectedIons.clear();
    int i = 0;
    BOOST_FOREACH (const MassChargeCandidate& mcc, consensusMccs) {
      SpectrumPtr consensusSpecCopy(new Spectrum(*consensusSpec));
      SpectrumHandler::setScannr(consensusSpecCopy, scannr*100+(++i));
      
      std::vector<MassChargeCandidate> mccs;
      mccs.push_back(mcc);
      SpectrumHandler::setMassChargeCandidates(consensusSpecCopy, mccs);
      
    #pragma omp critical (add_to_merged_speclist)
      {
        mergedSpectra->spectra.push_back(consensusSpecCopy);
      }
    }
  } else {
    SpectrumHandler::setScannr(consensusSpec, scannr);
    SpectrumHandler::setMassChargeCandidates(consensusSpec, consensusMccs);    
  #pragma omp critical (add_to_merged_speclist)
    {
      mergedSpectra->spectra.push_back(consensusSpec);
    }
  }
}

size_t MSFileHandler::getSpectrumIdxFromScannr(SpectrumListPtr sl, 
                                               unsigned int scannr) {
  size_t result = sl->find("scan=" + boost::lexical_cast<std::string>(scannr));
  if (result >= sl->size()) { // rescue attempt for MGF files
    pwiz::msdata::IndexList indexList = sl->findSpotID(
        "scan=" + boost::lexical_cast<std::string>(scannr));
    if (indexList.size() > 0) {
      result = indexList[0];
    } else {
      result = scannr;
    }
  }
  return result;
}

void MSFileHandler::calcPeakCount(SpectrumListPtr specList, PeakCounts& mzMap,
                                  std::string resultFN, bool normalizeXCorr, 
                                  bool addPrecursorMass, bool truncatePeaks, 
                                  unsigned int chargeFilter) {
  std::cerr << "Calculating peak count\n";
  
  if (!addPrecursorMass) {
    mzMap.setRelativeToPrecMz(0.001);
  }
  
  for (unsigned int i = 0; i < specList->size(); ++i) {
    SpectrumPtr s = specList->spectrum(i, true);
    
    std::map<unsigned int, bool> chargeSeen;
    std::vector<MassChargeCandidate> mccs;
    SpectrumHandler::getMassChargeCandidates(s, mccs);
    
    BOOST_FOREACH (const MassChargeCandidate mcc, mccs) {
      unsigned int charge = (std::min)(mcc.charge, mzMap.getMaxCharge());
      if (chargeFilter != 0 && charge != chargeFilter) continue;
      if (!chargeSeen[charge]) {
        chargeSeen[charge] = true;
        
        // in the last bin we do not truncate the spectrum
        if (charge == mzMap.getMaxCharge()) charge = 100u;
        
        std::vector<MZIntensityPair> mziPairs;
        SpectrumHandler::getMZIntensityPairs(s, mziPairs);  
        
        unsigned int numScoringPeaks = PvalueCalculator::getMaxScoringPeaks(mcc.mass);
        double precMz = mcc.precMz;
        if (normalizeXCorr) {
          std::vector<BinnedMZIntensityPair> mziPairsBinned;
          BinAndRank::binAndXCorrNormalize(mziPairs, mziPairsBinned);
          mzMap.addSpectrum(mziPairsBinned, precMz, charge, numScoringPeaks);
        } else {
          mzMap.addSpectrum(mziPairs, precMz, charge, mcc.mass, numScoringPeaks);
        }
      }
    }
    
    if (i % 5000 == 0) {
      std::cerr << "Processing spectrum " << i << "/" << specList->size() << 
                   std::endl;
    }
  }
  
  std::cerr << "Finished calculating peak count" << std::endl;
}


void MSFileHandler::calcRankDotProducts(const std::string& spectrumInFN) {
  std::cerr << "Calculating similarity scores!\n";

  // manipulate your MSData object here
  MSReaderList readerList;
  MSDataFile msd(spectrumInFN, &readerList);
  SpectrumListPtr sl = msd.run.spectrumListPtr;

  if (sl.get()) {
    std::cerr << "# of spectra in original dataset: " << sl->size() << std::endl;
  }
  
  std::map<unsigned int, std::pair<double, unsigned int> > mzMap;
  for (unsigned int i = 0; i < (std::min)(static_cast<unsigned int>(sl->size()),1000u); ++i) {
    if (i % 100 == 0) std::cerr << "Processing spectrum " << i << std::endl;
    SpectrumPtr s = sl->spectrum(i, true);
    std::vector<MZIntensityPair> mziPairs;
    SpectrumHandler::getMZIntensityPairs(s, mziPairs);
    //BinRanks brIn;
    BinRankScores brIn;
    BinAndRank::binAndRank(mziPairs, brIn);
    BinAndRank::normalizeRankScores(brIn);
  
    for (unsigned int j = i+1; j < (std::min)(static_cast<unsigned int>(sl->size()), 1000u); ++j) {
      SpectrumPtr sFrom = sl->spectrum(j, true);
      std::vector<MZIntensityPair> mziPairsFrom;
      SpectrumHandler::getMZIntensityPairs(sFrom, mziPairsFrom);
      //BinRanks brFrom;
      BinRankScores brFrom;
      BinAndRank::binAndRank(mziPairsFrom, brFrom);
      double dp = BinAndRank::rankDotProduct(brIn, brFrom);
      //std::cout << i << '\t' << j << '\t' << dp << std::endl;
      std::cout << dp << std::endl;
    }
	}
  
  for (std::map<unsigned int, std::pair<double, unsigned int> >::iterator it = mzMap.begin(); it != mzMap.end(); ++it) {
  	std::cout << it->first << '\t' << it->second.second << std::endl;
  }
}
