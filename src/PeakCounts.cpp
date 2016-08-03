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
 
#include "PeakCounts.h"

const std::map<unsigned int, unsigned int> PeakCountMatrix::emptyRow_ = std::map<unsigned int, unsigned int>();

unsigned int PeakCountMatrix::get(unsigned int row, unsigned int col) const { 
  if (peakCountMap.find(row) != peakCountMap.end()) {
    if (peakCountMap.find(row)->second.find(col) != peakCountMap.find(row)->second.end()) {
      return peakCountMap.find(row)->second.find(col)->second;
    }
  }
  return 0u;
}

const std::map<unsigned int, unsigned int>& PeakCountMatrix::getRow(unsigned int row) const {
  if (peakCountMap.find(row) != peakCountMap.end()) {
    return peakCountMap.find(row)->second;
  } else {
    return emptyRow_;
  }
}
void PeakCountMatrix::getRowIndices(std::vector<unsigned int>& rowIndices) const {
  rowIndices.clear();
  BOOST_FOREACH(PeakCountMatrixRow row, peakCountMap) {
    rowIndices.push_back(row.first);
  }
}

unsigned int PeakCountMatrix::getRowPeakCount(unsigned int row, unsigned int maxBin) const {
  unsigned int rowPeakCount = 0;
  if (peakCountMap.find(row) != peakCountMap.end()) {
    BOOST_FOREACH(ColumnValuePair cvPair, peakCountMap.find(row)->second) {
      if (cvPair.first < maxBin) {
        rowPeakCount += cvPair.second;
      } else {
        break;
      }
    }
  }
  return rowPeakCount;
}

void PeakCountMatrix::add(PeakCountMatrix& other) {
  BOOST_FOREACH(PeakCountMatrixRow row, other.getMatrix()) {
    unsigned int rowIdx = row.first;
    BOOST_FOREACH(ColumnValuePair cvPair, row.second) {
      add(rowIdx, cvPair.first, cvPair.second);
    }
  }
}

void PeakCountMatrix::subtract(PeakCountMatrix& other) {
  BOOST_FOREACH(PeakCountMatrixRow row, other.getMatrix()) {
    unsigned int rowIdx = row.first;
    BOOST_FOREACH(ColumnValuePair cvPair, row.second) {
      subtract(rowIdx, cvPair.first, cvPair.second);
    }
  }
}

unsigned int SpectrumCountVector::get(unsigned int row) const { 
  if (specCountMap.find(row) != specCountMap.end())
    return specCountMap.find(row)->second;
  else 
    return 0u;
}

void SpectrumCountVector::getRowIndices(std::vector<unsigned int>& rowIndices) const {
  rowIndices.clear();
  BOOST_FOREACH(SpecCountMapRow row, specCountMap) {
    rowIndices.push_back(row.first);
  }
}

void SpectrumCountVector::add(SpectrumCountVector& other) {
  std::vector<unsigned int> rowIndices;
  other.getRowIndices(rowIndices);
  BOOST_FOREACH(unsigned int rowIdx, rowIndices) {
    add(rowIdx, other.get(rowIdx));
  }
}

void SpectrumCountVector::subtract(SpectrumCountVector& other) {
  std::vector<unsigned int> rowIndices;
  other.getRowIndices(rowIndices);
  BOOST_FOREACH(unsigned int rowIdx, rowIndices) {
    subtract(rowIdx, other.get(rowIdx));
  }
}

void PeakCounts::setRelativeToPrecMz(double fracPeakBinWidth) {
  relativeToPrecMz = true;
  peakBinWidth = fracPeakBinWidth;
  peakBinShift = 0.0;
}

bool PeakCounts::equalParams(const PeakCounts& otherPeakCounts) const {
  return (otherPeakCounts.getMaxCharge() == maxCharge 
           && std::abs(otherPeakCounts.getPrecBinWidth() - precBinWidth) < 1e-9 
           && std::abs(otherPeakCounts.getPrecBinShift() - precBinShift) < 1e-9  
           && std::abs(otherPeakCounts.getPeakBinWidth() - peakBinWidth) < 1e-9 
           && std::abs(otherPeakCounts.getPeakBinShift() - peakBinShift) < 1e-9
           && otherPeakCounts.isRelativeToPrecMz() == relativeToPrecMz);
}

void PeakCounts::addSpectrum(std::vector<MZIntensityPair>& mziPairs, const double precMz, const unsigned int charge, const double precMass, const unsigned int numQueryPeaks) {
  unsigned int nSpectrumPeaks = 0;
  unsigned int chargeBin = getChargeBin(charge);
  unsigned int precursorBin = getPrecBin(precMz);
  specCountVectors.at(chargeBin).add(precursorBin, 1);
  
  std::vector<unsigned int> peakBins;
  BinSpectra::binBinaryTruncated(mziPairs, peakBins, numQueryPeaks, precMass);
  BOOST_FOREACH (const unsigned int bin, peakBins) {
    if (relativeToPrecMz) {
      peakCountMatrices.at(chargeBin).add(0, getRelPeakBin(bin, precMz, peakBinWidth), 1);
    } else {
      peakCountMatrices.at(chargeBin).add(precursorBin, bin, 1);
    }
    ++nSpectrumPeaks;
  }
  
  nTotalPeaks += nSpectrumPeaks;
}

void PeakCounts::addSpectrum(std::vector<BinnedMZIntensityPair>& mziPairs, const double precMz, const unsigned int charge, const unsigned int numQueryPeaks) {
  unsigned int nSpectrumPeaks = 0;
  unsigned int chargeBin = getChargeBin(charge);
  unsigned int precursorBin = getPrecBin(precMz);
  specCountVectors.at(chargeBin).add(precursorBin, 1);
  
  if (truncatePeaks) std::sort( mziPairs.begin(), mziPairs.end(), SpectrumHandler::greaterIntensity );
  BOOST_FOREACH (const BinnedMZIntensityPair& mziPairBinned, mziPairs) {
    if (mziPairBinned.intensity > intThresh) {
      if (relativeToPrecMz) {
        peakCountMatrices.at(chargeBin).add(0, getRelPeakBin(mziPairBinned.mz, precMz, peakBinWidth), 1);
      } else {
        peakCountMatrices.at(chargeBin).add(precursorBin, mziPairBinned.binIdx, 1);
      }
      ++nSpectrumPeaks;
      if (truncatePeaks && nSpectrumPeaks >= numQueryPeaks) break;
    }
  }
  
  nTotalPeaks += nSpectrumPeaks;
}

void PeakCounts::add(PeakCounts& otherPeakCounts) {
  if (equalParams(otherPeakCounts)) {
    nTotalPeaks += otherPeakCounts.getNumTotalPeaks();
    for (unsigned int charge = 1; charge <= maxCharge; ++charge) {
      getPeakCountMatrix(charge).add(otherPeakCounts.getPeakCountMatrix(charge));
      getSpecCountVector(charge).add(otherPeakCounts.getSpecCountVector(charge));
    }
  } else {
    std::cerr << "WARNING (PeakCounts::add(otherPeakCounts)): PeakCount parameters do not match, addition cancelled." << std::endl << std::endl;
    std::cerr << otherPeakCounts.getPrecBinWidth() << " " << precBinWidth << " " << (otherPeakCounts.getPrecBinWidth() - precBinWidth) << std::endl;
    std::cerr << otherPeakCounts.getPrecBinShift() << " " << precBinShift << std::endl;
    std::cerr << otherPeakCounts.getPeakBinWidth() << " " << peakBinWidth << std::endl;
    std::cerr << otherPeakCounts.getPeakBinShift() << " " << peakBinShift << std::endl;
    std::cerr << otherPeakCounts.isRelativeToPrecMz() << " " << relativeToPrecMz << " " << (otherPeakCounts.isRelativeToPrecMz() == relativeToPrecMz) << std::endl;
  }
}

void PeakCounts::subtract(PeakCounts& otherPeakCounts) {
  if (equalParams(otherPeakCounts)) {
    nTotalPeaks -= otherPeakCounts.getNumTotalPeaks();
    for (unsigned int charge = 1; charge <= maxCharge; ++charge) {
      getPeakCountMatrix(charge).subtract(otherPeakCounts.getPeakCountMatrix(charge));
      getSpecCountVector(charge).subtract(otherPeakCounts.getSpecCountVector(charge));
    }
  } else {
    std::cerr << "WARNING (PeakCounts::subtract(otherPeakCounts)): PeakCount parameters do not match, subtraction cancelled." << std::endl << std::endl;
  }
}

double PeakCounts::getDistance(const PeakCounts& otherPeakCounts) const {
  double distance = 0.0;
  if (equalParams(otherPeakCounts)) {
    for (unsigned int charge = 1; charge <= maxCharge; ++charge) {
      PeakDistribution distribution1, distribution2;
      generateRelativePeakDistribution(charge, distribution1);
      otherPeakCounts.generateRelativePeakDistribution(charge, distribution2);
      distance += distribution1.getDistance(distribution2) / maxCharge;
    }
  } else {
    std::cerr << "WARNING (PeakCounts::getDistance(otherPeakCounts)): PeakCount parameters do not match, maximum distance returned." << std::endl << std::endl;
    distance = 1.0;
  }
  return distance;
}

// TODO: generate unit test for this
void PeakCounts::generatePeakDistribution(double precMz, unsigned int charge,
    PeakDistribution& distribution, unsigned int numQueryPeaks) {
  unsigned int precMzBin = getPrecBin(precMz);
  unsigned int chargeBin = getChargeBin(charge);
  
  if (peakDistCache_[chargeBin].find(precMzBin) != peakDistCache_[chargeBin].end()) {
    distribution = peakDistCache_[chargeBin][precMzBin];
    // TODO: find a better way to do caching. Here I needed to do an extra check, 
    // since the copying of the distribution is sometimes not completely done 
    // when pvalCalc already tries to access it.
    if (distribution.getDistribution().size() != 0) return;
  }
  
  unsigned int precWindow, priorCount, windowRange;
  if (smoothingMode_ == 1) {
    precWindow = 0;
    priorCount = 1;
    windowRange = 0;
  } else {
    precWindow = 5;
    priorCount = 5;
    windowRange = 5;
  }
  unsigned int windowBinSize = 2*windowRange + 1;

  unsigned int minPrecBin = precMzBin - precWindow;
  unsigned int maxPrecBin = precMzBin + precWindow;
  int precMzInt = static_cast<int>(std::ceil(getPrecMz(precMzBin + 1)));
  unsigned int maxPeakBin = getPeakBin(precMzInt*charge) + 1;

  unsigned int totalSpecCount = priorCount*maxPeakBin / numQueryPeaks;
  std::vector<double> peakCountSum(maxPeakBin, static_cast<double>(totalSpecCount*numQueryPeaks)/maxPeakBin);  
  for (unsigned int precBin = minPrecBin; precBin <= maxPrecBin; ++precBin) {
    unsigned int specCount = specCountVectors.at(chargeBin).get(precBin);
    totalSpecCount += specCount;
    unsigned int rowPeakCount = peakCountMatrices.at(chargeBin).getRowPeakCount(precBin, maxPeakBin);
    if (rowPeakCount > 0u) {
      double multFactor = static_cast<double>(numQueryPeaks * specCount) / rowPeakCount; // correction for spectra not containing maxScoringPeaks
      //std::cerr << multFactor << std::endl;
      BOOST_FOREACH (PeakCountMatrix::ColumnValuePair colValPair, peakCountMatrices.at(chargeBin).getRow(precBin)) {
        unsigned int col = colValPair.first;
        unsigned int value = colValPair.second;
        
        if (col < maxPeakBin) {
          peakCountSum.at(col) += value * multFactor;
        }
      }
    }
  }

  //std::cerr << "Num spectra: " << totalSpecCount << std::endl;
  BOOST_FOREACH (double & peakProb, peakCountSum) {
    peakProb /= totalSpecCount;
  }

  distribution.init(maxPeakBin);
  double runningAvg = 0.0;  
  for (unsigned int bin = 0; bin <= windowRange; ++bin) {
    runningAvg += peakCountSum.at(bin)/windowBinSize;
  }
  
  for (unsigned int bin = 0; bin < maxPeakBin; ++bin) {
	  distribution.insert(bin, runningAvg);
	  if (bin + windowRange + 1 < maxPeakBin) {
    	runningAvg += peakCountSum.at(bin + windowRange + 1)/windowBinSize;
    }
    if (bin >= windowRange) {
      runningAvg -= peakCountSum.at(bin - windowRange)/windowBinSize;
    }
  }
  
  if (peakDistCache_[chargeBin].find(precMzBin) == peakDistCache_[chargeBin].end()) {
  #pragma omp critical (save_peak_dist)
    {
      peakDistCache_[chargeBin][precMzBin] = distribution;
    }
  }
}

// TODO: generate unit test for this
void PeakCounts::generateRelativePeakDistribution(unsigned int charge, PeakDistribution& distribution) const {
  double precWindow = 5.0;
  unsigned int priorCount = 5;
  unsigned int windowRange = 5;
  unsigned int windowBinSize = 2*windowRange + 1;
  
  double fracPeakBinWidth = 0.001 * charge;
  unsigned int maxPeakBin = static_cast<unsigned int>(charge / fracPeakBinWidth);
  
  unsigned int chargeBin = getChargeBin(charge);
  
  std::vector<unsigned int> peakCountSum(maxPeakBin, priorCount);
  
  BOOST_FOREACH (PeakCountMatrix::PeakCountMatrixRow matrixRow, peakCountMatrices.at(chargeBin).getMatrixConst()) {
    double precMz = getPrecMz(matrixRow.first);
    BOOST_FOREACH (PeakCountMatrix::ColumnValuePair colValPair, matrixRow.second) {
      unsigned int fracPeakBin = getRelPeakBin(colValPair.first, precMz, fracPeakBinWidth);
      unsigned int value = colValPair.second;
      
      if (fracPeakBin < maxPeakBin) {
        peakCountSum[fracPeakBin] += value;
      }
    }
  }
  
  distribution.init(maxPeakBin);
  double normSum = 0.0;
  double runningAvg = 0.0;
  
  for (unsigned int bin = 0; bin <= windowRange; ++bin) {
    runningAvg += static_cast<double>(peakCountSum.at(bin))/windowBinSize;
  }
  
  for (unsigned int bin = 0; bin < maxPeakBin; ++bin) {
  	distribution.insert(bin, runningAvg);
  	normSum += runningAvg;
  	if (bin + windowRange + 1 < maxPeakBin) {
    	runningAvg += static_cast<double>(peakCountSum.at(bin + windowRange + 1))/windowBinSize;
    }
    if (bin >= windowRange) {
      runningAvg -= static_cast<double>(peakCountSum.at(bin - windowRange))/windowBinSize;
    }
  }
  
  distribution.rescale(1.0/normSum);
}

void PeakCounts::readFromFile(const std::string& peakCountFN) {
  std::ifstream peakCountStream(peakCountFN.c_str(), std::ios_base::binary | std::ios_base::in);
  std::string serializedPeakCounts(
      (std::istreambuf_iterator<char>(peakCountStream)),
      std::istreambuf_iterator<char>());
  
  deserializePeakCounts(serializedPeakCounts, *this);
}

void PeakCounts::serialize(std::string& peakCountsSerialized) {
  std::stringstream os(std::ios_base::binary | std::ios_base::out);
  {
    boost::archive::binary_oarchive oa(os, boost::archive::no_header);
    oa << *this;
  }
  peakCountsSerialized = os.str();
}

void PeakCounts::serializePeakCounts(PeakCounts& peakCounts, std::string& peakCountsSerialized) {
  std::stringstream os(std::ios_base::binary | std::ios_base::out);
  {
    boost::archive::binary_oarchive oa(os, boost::archive::no_header);
    oa << peakCounts;
  }
  peakCountsSerialized = os.str();
}

void PeakCounts::deserialize(std::string& peakCountsSerialized) {
  std::stringstream is(peakCountsSerialized, std::ios_base::binary | std::ios_base::in);
  {
    boost::archive::binary_iarchive ia(is, boost::archive::no_header);
    ia >> *this;
  }
}

void PeakCounts::deserializePeakCounts(std::string& peakCountsSerialized, PeakCounts& peakCounts) {
  std::stringstream is(peakCountsSerialized, std::ios_base::binary | std::ios_base::in);
  {
    boost::archive::binary_iarchive ia(is, boost::archive::no_header);
    ia >> peakCounts;
  }
}

void PeakCounts::print(const std::string& resultBaseFN) {
  for (unsigned int charge = 1; charge <= getMaxCharge(); 
       ++charge) {
    std::string resultFN = resultBaseFN +  
        boost::lexical_cast<std::string>(charge) + ".txt";
    std::ofstream resultStream(resultFN.c_str());
    print(resultStream, charge);
  }
}

void PeakCounts::print(std::ostream& os, const unsigned int charge) {
  unsigned int chargeBin = getChargeBin(charge);
  if (relativeToPrecMz) {
    double windowSize = 0.005; // window for smoothing of distribution
    double scalingFactor = 1.0/(nTotalPeaks*peakBinWidth);
    int windowRange = static_cast<int>(windowSize/peakBinWidth);
    int windowBinSize = 2*windowRange + 1;
    double runningAvg = 0.0;
    for (unsigned int bin = 0; bin <= static_cast<unsigned int>(windowRange); ++bin) {
      runningAvg += peakCountMatrices.at(chargeBin).get(0, bin)/windowBinSize;
    }
    for (unsigned int bin = 0; bin < static_cast<unsigned int>(2.5/peakBinWidth); ++bin) {
    	os << bin * peakBinWidth << '\t' << runningAvg * scalingFactor << std::endl;
    	runningAvg += peakCountMatrices.at(chargeBin).get(0, bin + windowRange + 1)/windowBinSize;
      if (static_cast<int>(bin) - windowRange >= 0) {
        runningAvg -= peakCountMatrices.at(chargeBin).get(0, bin - windowRange)/windowBinSize;
      }
    }
  } else {
    // write headers with peak bin locations
    double minPeakMz = 100;
    double maxPeakMz = 2400;
    os << 0 << '\t' << 0;
    for (unsigned int bin = BinSpectra::getBin(minPeakMz,peakBinWidth); bin < BinSpectra::getBin(maxPeakMz,peakBinWidth); ++bin) {
    	os << '\t' << BinSpectra::getMZ(bin, peakBinWidth, peakBinShift);
    }
    os << std::endl;
    
    // write peak counts
    std::vector<unsigned int> rowIndices;
    peakCountMatrices.at(chargeBin).getRowIndices(rowIndices);
    unsigned int lowerBound = rowIndices[0];
    unsigned int upperBound = rowIndices[rowIndices.size()-1];
    for (unsigned int rowIdx = lowerBound; rowIdx <= upperBound; ++rowIdx) {
      os << BinSpectra::getMZ(rowIdx, precBinWidth, precBinShift);
      os << '\t' << specCountVectors.at(chargeBin).get(rowIdx);
      for (unsigned int bin = BinSpectra::getBin(minPeakMz,peakBinWidth); bin < BinSpectra::getBin(maxPeakMz,peakBinWidth); ++bin) {
      	os << '\t' << peakCountMatrices.at(chargeBin).get(rowIdx, bin);
      }
      os << std::endl;
    }
  }
}

bool PeakCounts::peakCountsSerializationUnitTest() {
  PeakCounts pk1, pk2, pk3, pk4;
  std::string pk1serialized, pk2serialized;
  
  std::vector<MZIntensityPair> mziPairs1, mziPairs2;
  
  mziPairs1.push_back(MZIntensityPair(50,1));
  mziPairs1.push_back(MZIntensityPair(100,2));
  
  mziPairs2.push_back(MZIntensityPair(50,1));
  mziPairs2.push_back(MZIntensityPair(75,3));
  
  pk1.addSpectrum(mziPairs1, 150, 2u, 300, 100u);
  pk1.addSpectrum(mziPairs2, 150, 2u, 300, 100u);
  
  pk2.addSpectrum(mziPairs1, 150, 2u, 300, 100u);
  pk2.addSpectrum(mziPairs2, 250, 2u, 500, 100u);
  
  serializePeakCounts(pk1, pk1serialized);
  serializePeakCounts(pk2, pk2serialized);
  
  /*
  std::cerr << pk1.getPeakCountMatrix().get(9,99) << std::endl;
  std::cerr << pk3.getPeakCountMatrix().get(9,99) << std::endl;
  
  std::cerr << pk1.getPeakCountMatrix().get(10,100) << std::endl;
  std::cerr << pk3.getPeakCountMatrix().get(10,100) << std::endl;

  std::cerr << pk1.getPeakCountMatrix().get(11,101) << std::endl;
  std::cerr << pk3.getPeakCountMatrix().get(11,101) << std::endl;
  */
  
  deserializePeakCounts(pk1serialized, pk3);
  deserializePeakCounts(pk2serialized, pk4);
  
  if (pk1.getSpecCountVector(2u).get(150) != 2 ||
      pk1.getSpecCountVector(2u).get(150) != pk3.getSpecCountVector(2u).get(150)) {
    std::cerr << "Deserialization returned a false specCountVector" << std::endl;
    
    std::cerr << pk1.getSpecCountVector(2u).get(150) << std::endl;
    std::cerr << pk3.getSpecCountVector(2u).get(150) << std::endl;
    return false;
  }
  
  if (pk1.getPeakCountMatrix(2u).get(150,100) != 1 || 
      pk1.getPeakCountMatrix(2u).get(150,100) != pk3.getPeakCountMatrix(2u).get(150,100)) {
    std::cerr << "Deserialization returned a false peakCountMatrix" << std::endl;
    
    std::cerr << pk1.getPeakCountMatrix(2u).get(149,99) << std::endl;
    std::cerr << pk3.getPeakCountMatrix(2u).get(149,99) << std::endl;
    
    std::cerr << pk1.getPeakCountMatrix(2u).get(150,100) << std::endl;
    std::cerr << pk3.getPeakCountMatrix(2u).get(150,100) << std::endl;

    std::cerr << pk1.getPeakCountMatrix(2u).get(151,101) << std::endl;
    std::cerr << pk3.getPeakCountMatrix(2u).get(151,101) << std::endl;
    return false;
  }
  
  pk3.add(pk2);
  
  serializePeakCounts(pk3, pk1serialized);
  deserializePeakCounts(pk1serialized, pk4);
    
  if (pk3.getSpecCountVector(2u).get(150) != 3 ||
      pk3.getSpecCountVector(2u).get(150) != pk4.getSpecCountVector(2u).get(150)) {
    std::cerr << "Deserialization followed by addition returned a false specCountVector" << std::endl;
    
    std::cerr << pk4.getSpecCountVector(2u).get(150) << std::endl;
    std::cerr << pk3.getSpecCountVector(2u).get(150) << std::endl;
    return false;
  }
  
  if (pk3.getPeakCountMatrix(2u).get(150,50) == 2 ||
      pk3.getPeakCountMatrix(2u).get(150,50) != pk4.getPeakCountMatrix(2u).get(150,50)) {
    std::cerr << "Deserialization followed by addition returned a false peakCountMatrix" << std::endl;
    
    std::cerr << pk4.getPeakCountMatrix(2u).get(149,49) << std::endl;
    std::cerr << pk3.getPeakCountMatrix(2u).get(149,49) << std::endl;
    
    std::cerr << pk4.getPeakCountMatrix(2u).get(150,50) << std::endl;
    std::cerr << pk3.getPeakCountMatrix(2u).get(150,50) << std::endl;

    std::cerr << pk4.getPeakCountMatrix(2u).get(151,51) << std::endl;
    std::cerr << pk3.getPeakCountMatrix(2u).get(151,51) << std::endl;
    return false;
  }
  
  pk3.add(pk4);
  
  serializePeakCounts(pk3, pk1serialized);
  deserializePeakCounts(pk1serialized, pk4);
  
  return true;
}

bool PeakCountMatrix::peakMatrixUnitTest() {
  PeakCountMatrix m, n, p;
  
  m.add(3,4,1);
  m.add(3,4,1);
  m.add(4,2,3);
  
  n.add(3,4,3);
  n.add(5,1,1);
  
  p.add(3,4,4);
  p.add(4,2,1);
  
  p.subtract(m);
  m.add(n);
  
  if (m.get(3,4) == 5 && n.get(3,4) == 3 && p.get(3,4) == 2 && p.get(4,2) == 0) {
    return true;
  } else {
    std::cerr << m.get(3,4) << std::endl;
    std::cerr << n.get(3,4) << std::endl;
    std::cerr << p.get(3,4) << std::endl;
    std::cerr << p.get(4,2) << std::endl;
    return false;
  }
}

bool SpectrumCountVector::specVectorUnitTest() {
  SpectrumCountVector s, t, v;
  
  s.add(5,2);
  s.add(3,1);
  
  t.add(5,1);
  t.add(5,2);
  
  v.add(5,1);
  v.add(3,4);
  
  v.subtract(s);
  s.add(t);
  
  if (s.get(5) == 5 && t.get(5) == 3 && v.get(5) == 0 && v.get(3) == 3) {
    return true; 
  } else {
    std::cerr << s.get(5) << std::endl;
    std::cerr << t.get(5) << std::endl;
    std::cerr << v.get(5) << std::endl;
    std::cerr << v.get(3) << std::endl;
    return false;
  }
}
