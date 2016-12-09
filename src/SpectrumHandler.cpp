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
 
#include "SpectrumHandler.h"

unsigned int SpectrumHandler::getScannr(pwiz::msdata::SpectrumPtr s) {
  // special check for the SCANS field in mgf files
  if (s->hasCVParam(pwiz::cv::MS_peak_list_scans)) {
    std::string mgfScansField = s->cvParam(pwiz::cv::MS_peak_list_scans).valueAs<std::string>();
    return atoi(mgfScansField.c_str());
  }
  std::stringstream ss(s->id);
  while (ss.good()) {
    std::string tmp;
    ss >> tmp;
    if (tmp.substr(0,5) == "scan=")
      return atoi(tmp.substr(5).c_str());

    if (tmp.substr(0,6) == "index=")
      return atoi(tmp.substr(6).c_str());
  }
  std::cerr << "Warning: could not extract scannr. Returning index " << s->index << " (" << s->id << ")" << std::endl;
  
  // using the spectrum index should be more robust and
  // independent of the input format's indexing system, 
  // but makes it harder to search for the original spectra 
  // in the original files, e.g. for the mgf and ms2 format
  return s->index;
}

// TODO: this could be made a bit nicer by retaining other attributes and only replacing the scan attribute
void SpectrumHandler::setScannr(pwiz::msdata::SpectrumPtr s, unsigned int scanNr) {
  std::string id = "scan=" + boost::lexical_cast<std::string>(scanNr);
  s->id = id;
  s->set(pwiz::msdata::MS_spectrum_title, id); // sets the TITLE attribute in mgf files
}

void SpectrumHandler::setScannr(pwiz::msdata::SpectrumPtr s, const ScanId& scanId) {
  setScannr(s, hash_value(scanId));
}

double SpectrumHandler::interpolateIntensity(MZIntensityPair p1, MZIntensityPair p2, double mz) {
  return ( (p2.intensity - p1.intensity) / (p2.mz - p1.mz) * (mz - p1.mz) + p1.intensity );
}

double SpectrumHandler::getMaxIntensity(const std::vector<MZIntensityPair>& mziPairs) {
	return std::max_element(mziPairs.begin(), mziPairs.end(), lessIntensity)->intensity;
}

double SpectrumHandler::getSNR(std::vector<MZIntensityPair>& mziPairs) {
  std::sort( mziPairs.begin(), mziPairs.end(), SpectrumHandler::greaterIntensity );
  double accTop = 0.0, accBottom = 0.0;
  for (std::vector<MZIntensityPair>::iterator it = mziPairs.begin(); it != mziPairs.begin() + 10; ++it) {
    accTop += (*it).intensity;
  }
  for (std::vector<MZIntensityPair>::iterator it = mziPairs.end() - 10; it != mziPairs.end(); ++it) {
    accBottom += (*it).intensity;
  }
  return accTop/accBottom;
}

void SpectrumHandler::scaleIntensities(std::vector<MZIntensityPair>& mziPairs, double scaling) { 
	BOOST_FOREACH (MZIntensityPair& mziPair, mziPairs) {
		mziPair.intensity *= scaling;
	}
}

void SpectrumHandler::normalizeIntensities(std::vector<MZIntensityPair>& mziPairs) {
	scaleIntensities(mziPairs, 1.0/getMaxIntensity(mziPairs));
}

void SpectrumHandler::normalizeIntensitiesSNR(std::vector<MZIntensityPair>& mziPairs) {
  double SNR = getSNR(mziPairs);
  double maxInt = getMaxIntensity(mziPairs);
  double multiplicity = static_cast<unsigned int>(SNR/10) + 1;
  BOOST_FOREACH (MZIntensityPair& mziPair, mziPairs) {
		mziPair.intensity *= multiplicity*SNR/maxInt;
		mziPair.multiplicity = multiplicity;
	}
}

void SpectrumHandler::normalizeIntensitiesMSCluster(std::vector<MZIntensityPair>& mziPairs) {
  double totalPeakIntensity = 0.0;
  BOOST_FOREACH (const MZIntensityPair& mziPair, mziPairs) {
		totalPeakIntensity += mziPair.intensity;
	}
	
	const float normalizationValue = 1000.0 / totalPeakIntensity;
	
	BOOST_FOREACH (MZIntensityPair& mziPair, mziPairs) {
		mziPair.intensity *= normalizationValue;
	}
}

void SpectrumHandler::printIntensities(std::vector<MZIntensityPair>& mziPairs) {
  BOOST_FOREACH (MZIntensityPair& mziPair, mziPairs) {
		std::cout << mziPair.mz << "\t" << mziPair.intensity << std::endl;
	}
}

void SpectrumHandler::convertMSDataMZIntensityPairs(std::vector<pwiz::msdata::MZIntensityPair>& MSDataMziPairs, std::vector<MZIntensityPair>& mziPairs) {
  mziPairs.clear();
  BOOST_FOREACH(pwiz::msdata::MZIntensityPair MSDataMziPair, MSDataMziPairs) {
    mziPairs.push_back(MZIntensityPair(MSDataMziPair));
  }
}

void SpectrumHandler::convertMSDataMZIntensityPairs(std::vector<MZIntensityPair>& mziPairs, std::vector<pwiz::msdata::MZIntensityPair>& MSDataMziPairs) {
  MSDataMziPairs.clear();
  BOOST_FOREACH(MZIntensityPair mziPair, mziPairs) {
    MSDataMziPairs.push_back(pwiz::msdata::MZIntensityPair(mziPair.mz, mziPair.intensity));
  }
}

void SpectrumHandler::getMZIntensityPairs(pwiz::msdata::SpectrumPtr s, std::vector<MZIntensityPair>& mziPairs) {
  std::vector<pwiz::msdata::MZIntensityPair> MSDataMziPairs;
  s->getMZIntensityPairs(MSDataMziPairs);
  convertMSDataMZIntensityPairs(MSDataMziPairs, mziPairs);
}

void SpectrumHandler::setMZIntensityPairs(pwiz::msdata::SpectrumPtr s, std::vector<MZIntensityPair>& mziPairs) {
  std::vector<pwiz::msdata::MZIntensityPair> MSDataMziPairs;
  convertMSDataMZIntensityPairs(mziPairs, MSDataMziPairs);
  s->setMZIntensityPairs(MSDataMziPairs, pwiz::cv::MS_number_of_detector_counts);
}

void SpectrumHandler::getMassChargeCandidates(pwiz::msdata::SpectrumPtr s, 
    std::vector<MassChargeCandidate>& mcc, int chargeUncertainty) {
  mcc.clear();
  if (s->precursors.size() > 0) {
    std::vector<pwiz::msdata::SelectedIon>& ions = s->precursors.at(0).selectedIons;
    
    for (std::vector<pwiz::msdata::SelectedIon>::iterator it = ions.begin(); it != ions.end(); ++it) {
      int charge = it->cvParam(pwiz::cv::MS_charge_state).valueAs<int>();
      double precMz = it->cvParam(pwiz::cv::MS_selected_ion_m_z).valueAs<double>();
      if (charge == 0) {
        charge = it->cvParam(pwiz::cv::MS_possible_charge_state).valueAs<int>();
        if (charge == 0) charge = 2;
      }
      for (int chargeState = std::max(charge - chargeUncertainty, 1); 
            chargeState <= charge + chargeUncertainty; ++chargeState) {
        double mass = calcMass(precMz, chargeState);
        /*
        if (it->hasCVParam(pwiz::cv::MS_accurate_mass_OBSOLETE)) {
          mass = it->cvParam(pwiz::cv::MS_accurate_mass_OBSOLETE).valueAs<double>();
        }
        */
        mcc.push_back(MassChargeCandidate(chargeState, precMz, mass));
      }
    }
    
    std::sort(mcc.begin(), mcc.end(), MassChargeCandidate::lessChargeMass);
  }
}

void SpectrumHandler::setMassChargeCandidates(pwiz::msdata::SpectrumPtr s,
    std::vector<MassChargeCandidate>& mccs) {
  if (s->precursors.at(0).activation.empty()) {
    s->precursors.at(0).activation.set(pwiz::cv::MS_dissociation_method);
  }
  s->precursors.at(0).selectedIons.clear();
  BOOST_FOREACH (const MassChargeCandidate& mcc, mccs) {
    s->precursors.at(0).selectedIons.push_back(
        pwiz::msdata::SelectedIon(mcc.precMz, mcc.charge));
    // used by ms2 output format for EZ lines
    s->precursors.at(0).selectedIons.back().cvParams.push_back(
        pwiz::data::CVParam(pwiz::cv::MS_accurate_mass_OBSOLETE, mcc.mass, pwiz::cv::UO_mass_unit));
    s->precursors.at(0).isolationWindow.set(
        pwiz::cv::MS_isolation_window_target_m_z, mcc.precMz, pwiz::cv::MS_m_z);
  }
}

void SpectrumHandler::updateMassChargeCandidates(pwiz::msdata::SpectrumPtr s,
    int chargeUncertainty) {
  std::vector<MassChargeCandidate> mccs;
  getMassChargeCandidates(s, mccs, chargeUncertainty);
  setMassChargeCandidates(s, mccs);
}

unsigned int SpectrumHandler::getCharge(pwiz::msdata::SpectrumPtr s) {
  if (s->precursors.size() > 0) {
    std::vector<pwiz::msdata::SelectedIon>& ions = s->precursors.at(0).selectedIons;
    return ions.begin()->cvParam(pwiz::cv::MS_charge_state).valueAs<unsigned int>();
  } else {
    return 2;
  }
}

double SpectrumHandler::getPrecMz(pwiz::msdata::SpectrumPtr s) {
  if (s->precursors.size() > 0) {
    std::vector<pwiz::msdata::SelectedIon>& ions = s->precursors.at(0).selectedIons;
    return ions.begin()->cvParam(pwiz::cv::MS_selected_ion_m_z).valueAs<double>();
  } else {
    return 0.0;
  }
}

double SpectrumHandler::getRetentionTime(pwiz::msdata::SpectrumPtr s) {
  if (s->scanList.scans.size() > 0) {
    return s->scanList.scans.back().cvParam(pwiz::cv::MS_scan_start_time).valueAs<double>();
  } else {
    return 0.0;
  }
}

bool SpectrumHandler::isMs2Scan(pwiz::msdata::SpectrumPtr s) {
  return (s->cvParam(pwiz::cv::MS_ms_level).valueAs<int>() == 2);
}

/* 
 *
 * The following two functions are copies (with small modifications) of
 * proteowizard/pwiz/analysis/spectrum_processing/SpectrumList_MetadataFixer.cpp
 *
 */
void SpectrumHandler::replaceCvParam(pwiz::msdata::ParamContainer& pc, 
    pwiz::cv::CVID cvid, double value, pwiz::cv::CVID unit = pwiz::cv::CVID_Unknown) {
  std::vector<pwiz::data::CVParam>::iterator itr;
  
  itr = std::find(pc.cvParams.begin(), pc.cvParams.end(), cvid);
  if (itr == pc.cvParams.end()) {
    pc.set(cvid, value, unit);
  } else {
    itr->value = boost::lexical_cast<std::string>(value);
    itr->units = unit;
  }
}


void SpectrumHandler::fixMetaData(pwiz::msdata::SpectrumPtr s) {
  pwiz::msdata::BinaryDataArrayPtr mzArray = s->getMZArray();
  pwiz::msdata::BinaryDataArrayPtr intensityArray = s->getIntensityArray();
  if (!mzArray.get() || !intensityArray.get())
      return;

  std::vector<double>& mzs = mzArray->data;
  std::vector<double>& intensities = intensityArray->data;

  double tic = 0;
  if (!mzs.empty()) {
    double bpmz, bpi = -1;
    for (size_t i=0, end=mzs.size(); i < end; ++i) {
      tic += intensities[i];
      if (bpi < intensities[i]) {
        bpi = intensities[i];
        bpmz = mzs[i];
      }
    }

    replaceCvParam(*s, pwiz::cv::MS_base_peak_intensity, bpi, pwiz::cv::MS_number_of_detector_counts);
    replaceCvParam(*s, pwiz::cv::MS_base_peak_m_z, bpmz, pwiz::cv::MS_m_z);
    replaceCvParam(*s, pwiz::cv::MS_lowest_observed_m_z, mzs.front(), pwiz::cv::MS_m_z);
    replaceCvParam(*s, pwiz::cv::MS_highest_observed_m_z, mzs.back(), pwiz::cv::MS_m_z);
  }

  replaceCvParam(*s, pwiz::cv::MS_TIC, tic);
}
