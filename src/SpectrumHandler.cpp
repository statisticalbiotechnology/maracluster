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
  // independent of the input format's indexing system
  return s->index;
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

void SpectrumHandler::getMassChargeCandidates(pwiz::msdata::SpectrumPtr s, std::vector<MassChargeCandidate>& mcc) {
  mcc.clear();
  std::vector<pwiz::msdata::SelectedIon>& ions = s->precursors.at(0).selectedIons;
  
  for (std::vector<pwiz::msdata::SelectedIon>::iterator it = ions.begin(); it != ions.end(); ++it) {
    unsigned int charge = it->cvParam(pwiz::cv::MS_charge_state).valueAs<unsigned int>();
    double precMz = it->cvParam(pwiz::cv::MS_selected_ion_m_z).valueAs<double>();
    if (charge == 0) {
      charge = it->cvParam(pwiz::cv::MS_possible_charge_state).valueAs<unsigned int>();
      if (charge == 0) charge = 2;
    }
    double mass = calcMass(precMz, charge);
    if (it->hasCVParam(pwiz::cv::MS_accurate_mass_OBSOLETE)) {
      mass = it->cvParam(pwiz::cv::MS_accurate_mass_OBSOLETE).valueAs<double>();
    }
    mcc.push_back(MassChargeCandidate(charge, precMz, mass));
  }
  
  std::sort(mcc.begin(), mcc.end(), MassChargeCandidate::lessChargeMass);
}

unsigned int SpectrumHandler::getCharge(pwiz::msdata::SpectrumPtr s) {
  std::vector<pwiz::msdata::SelectedIon>& ions = s->precursors.at(0).selectedIons;
  return ions.begin()->cvParam(pwiz::cv::MS_charge_state).valueAs<unsigned int>();
}

double SpectrumHandler::getPrecMz(pwiz::msdata::SpectrumPtr s) {
  std::vector<pwiz::msdata::SelectedIon>& ions = s->precursors.at(0).selectedIons;
  return ions.begin()->cvParam(pwiz::cv::MS_selected_ion_m_z).valueAs<double>();
}

double SpectrumHandler::getRetentionTime(pwiz::msdata::SpectrumPtr s) {
  if (s->scanList.scans.size() > 0) {
    return s->scanList.scans.back().cvParam(pwiz::cv::MS_scan_start_time).valueAs<double>();
  } else {
    return 0.0;
  }
}

