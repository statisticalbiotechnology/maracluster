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
 
#ifndef MARACLUSTER_SPECTRUMHANDLER_H_
#define MARACLUSTER_SPECTRUMHANDLER_H_

#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <boost/foreach.hpp>
#include "pwiz/data/msdata/MSData.hpp"
#include "pwiz/data/msdata/MSDataFile.hpp"
#include "pwiz/data/common/cv.hpp"

#include "MZIntensityPair.h"
#include "MassChargeCandidate.h"
#include "ScanId.h"

#if defined(WIN32) && defined(VENDOR_SUPPORT)
#include "pwiz/data/vendor_readers/ExtendedReaderList.hpp"
#else
#include "pwiz/data/msdata/DefaultReaderList.hpp"
#endif

namespace maracluster {

#if defined(WIN32) && defined(VENDOR_SUPPORT)
typedef pwiz::msdata::ExtendedReaderList MSReaderList;
#else
typedef pwiz::msdata::DefaultReaderList MSReaderList;
#endif

const double PROTON_MASS = 1.00727646677;

class SpectrumHandler {  
 public:
  SpectrumHandler() {}
  
  inline static bool lessMZ(const MZIntensityPair& a, const MZIntensityPair& b) { return (a.mz < b.mz) || (a.mz == b.mz && a.intensity < b.intensity); }
  inline static bool greaterIntensity(const MZIntensityPair& a, const MZIntensityPair& b) { return (a.intensity > b.intensity); }
  inline static bool greaterMultiplicity(const MZIntensityPair& a, const MZIntensityPair& b) { 
    return ( (a.multiplicity > b.multiplicity) || ( (a.multiplicity == b.multiplicity) && (a.intensity > b.intensity) ) ); 
  }
  inline static bool lessIntensity(const MZIntensityPair& a, const MZIntensityPair& b) { return (a.intensity < b.intensity); }
  
  inline static bool lessScannr(pwiz::msdata::SpectrumPtr s1, pwiz::msdata::SpectrumPtr s2) { return getScannr(s1) < getScannr(s2); }
  
  inline static double calcMass(double precMz, unsigned int charge) {
    return precMz * charge - PROTON_MASS * (charge - 1);
  }
  inline static double calcPrecMz(double mass, unsigned int charge) {
    return (mass + PROTON_MASS * (charge - 1)) / charge;
  }
  
  static unsigned int getScannr(pwiz::msdata::SpectrumPtr s);
  static void setScannr(pwiz::msdata::SpectrumPtr s, unsigned int scanNr);
  static void setScannr(pwiz::msdata::SpectrumPtr s, const ScanId& scanId);
  
  static std::string getScanTitle(pwiz::msdata::SpectrumPtr s);
  
  static double interpolateIntensity(MZIntensityPair p1, MZIntensityPair p2, double mz);
  static void scaleIntensities(std::vector<MZIntensityPair>& mziPairs, double scaling);
  
  static double getMaxIntensity(const std::vector<MZIntensityPair>& mziPairs);
  static double getSNR(std::vector<MZIntensityPair>& mziPairs);
  static void normalizeIntensities(std::vector<MZIntensityPair>& mziPairs);
  static void normalizeIntensitiesSNR(std::vector<MZIntensityPair>& mziPairs);
  static void normalizeIntensitiesMSCluster(std::vector<MZIntensityPair>& mziPairs);
  static void printIntensities(std::vector<MZIntensityPair>& mziPairs);
  
  static void getMZIntensityPairs(pwiz::msdata::SpectrumPtr s, std::vector<MZIntensityPair>& mziPairs);
  static void setMZIntensityPairs(pwiz::msdata::SpectrumPtr s, std::vector<MZIntensityPair>& mziPairs);
  
  static void getMassChargeCandidates(pwiz::msdata::SpectrumPtr s, 
    std::vector<MassChargeCandidate>& mccs, int chargeUncertainty = 0);
  static void setMassChargeCandidates(pwiz::msdata::SpectrumPtr s, 
    std::vector<MassChargeCandidate>& mccs);
  static void updateMassChargeCandidates(pwiz::msdata::SpectrumPtr s,
    int chargeUncertainty = 0);
  
  static void fixMetaData(pwiz::msdata::SpectrumPtr s);
  
  static unsigned int getCharge(pwiz::msdata::SpectrumPtr s);
  static double getPrecMz(pwiz::msdata::SpectrumPtr s);
  static double getRetentionTime(pwiz::msdata::SpectrumPtr s);
  static bool isMs2Scan(pwiz::msdata::SpectrumPtr s);
 private:
  static void convertMSDataMZIntensityPairs(std::vector<pwiz::msdata::MZIntensityPair>& MSDataMziPairs, std::vector<MZIntensityPair>& mziPairs);
  static void convertMSDataMZIntensityPairs(std::vector<MZIntensityPair>& mziPairs, std::vector<pwiz::msdata::MZIntensityPair>& MSDataMziPairs);
  
  static void replaceCvParam(pwiz::msdata::ParamContainer& pc, 
    pwiz::cv::CVID cvid, double value, pwiz::cv::CVID unit);
    
};

} /* namespace maracluster */

#endif /* MARACLUSTER_SPECTRUMHANDLER_H_ */
