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
 
#ifndef PEAK_COUNTS_H
#define PEAK_COUNTS_H


#include <vector>
#include <map>

#include <iostream>

#include <utility>
#include <algorithm>
#include <numeric>
#include <cmath>

#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>

#include <boost/serialization/serialization.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/split_member.hpp>

#include "pwiz/data/msdata/MSData.hpp"

#include "MZIntensityPair.h"
#include "BinSpectra.h"
#include "PvalueCalculator.h"
#include "PeakDistribution.h"

/**
 * Creates a histogram of spectrum peaks split out by precursor
 */

class PeakCountMatrix {
    
  public:
    typedef std::map<unsigned int, std::map<unsigned int, unsigned int> > PeakCountMatrixMap;
    typedef std::pair<unsigned int, std::map<unsigned int, unsigned int> > PeakCountMatrixRow;
    typedef std::pair<unsigned int, unsigned int> ColumnValuePair;
    
    PeakCountMatrix() { }
    
    void add(PeakCountMatrix& other);
    void subtract(PeakCountMatrix& other);
    inline void add(unsigned int row, unsigned int col, unsigned int value) { peakCountMap[row][col] += value; }
    inline void subtract(unsigned int row, unsigned int col, unsigned int value) {
      peakCountMap[row][col] -= (std::min)(value, peakCountMap[row][col]);
    }
    
    // in parallel regions only use these const getters
    unsigned int get(unsigned int row, unsigned int col) const;
    const std::map<unsigned int, unsigned int>& getRow(unsigned int row) const;
    void getRowIndices(std::vector<unsigned int>& rowIndices) const;
    unsigned int getRowPeakCount(unsigned int row, unsigned int maxBin) const;
    
    PeakCountMatrixMap& getMatrix() { return peakCountMap; }
    const PeakCountMatrixMap& getMatrixConst() const { return peakCountMap; }
    
    unsigned int size() const { return peakCountMap.size(); }
    static bool peakMatrixUnitTest();
  private:
    PeakCountMatrixMap peakCountMap;
    static const std::map<unsigned int, unsigned int> emptyRow_;
    
    // splitting the serialization member function, because "ar & peakCountMap" 
    // caused a segfault (probably due to some pointer issues with a 2D map)
    friend class boost::serialization::access;
    template<class Archive>
    void save(Archive & ar, const unsigned int version) const {
      unsigned int numRows, numCols, col, value;
      std::vector<unsigned int> rowIndices;
      getRowIndices(rowIndices);
      numRows = peakCountMap.size();
      ar & numRows;
      BOOST_FOREACH(unsigned int row, rowIndices) {
        numCols = peakCountMap.find(row)->second.size();
        ar & row;
        ar & numCols;
        BOOST_FOREACH(ColumnValuePair colVal, peakCountMap.find(row)->second) {
          col = colVal.first;
          value = colVal.second;
          ar & col;
          ar & value;
        }
      }
    }
    
    template<class Archive>
    void load(Archive & ar, const unsigned int version) {
      unsigned int numRows, numCols, row, col;
      ar & numRows;
      for (unsigned int i = 0; i < numRows; ++i) {
        ar & row;
        ar & numCols;
        for (unsigned int j = 0; j < numCols; ++j) {
          ar & col;
          ar & peakCountMap[row][col];
        } 
      }
    }
    
    BOOST_SERIALIZATION_SPLIT_MEMBER()
};

class SpectrumCountVector {
  typedef std::pair<unsigned int, unsigned int> SpecCountMapRow;
  
  public:
    SpectrumCountVector() {}
    
    void add(SpectrumCountVector& other);
    void subtract(SpectrumCountVector& other);
    inline void add(unsigned int row, unsigned int value) { specCountMap[row] += value; }
    inline void subtract(unsigned int row, unsigned int value) { 
      specCountMap[row] -= (std::min)(value, specCountMap[row]);
    }
    unsigned int get(unsigned int row) const;
    void getRowIndices(std::vector<unsigned int>& rowIndices) const;
    
    unsigned int size() const { return specCountMap.size(); }
    static bool specVectorUnitTest();
  private:
    std::map<unsigned int, unsigned int> specCountMap;
    
    // splitting the serialization member function, because "ar & specCountMap" 
    // only deserialized the first entry, while ignoring the rest
    friend class boost::serialization::access;
    template<class Archive>
    void save(Archive & ar, const unsigned int version) const {
      unsigned int numRows, value;
      std::vector<unsigned int> rowIndices;
      getRowIndices(rowIndices);
      numRows = specCountMap.size();
      ar & numRows;
      BOOST_FOREACH(unsigned int row, rowIndices) {
        value = specCountMap.find(row)->second;
        ar & row;
        ar & value;
      }
    }
    
    template<class Archive>
    void load(Archive & ar, const unsigned int version) {
      unsigned int numRows, row;
      ar & numRows;
      for (unsigned int i = 0; i < numRows; ++i) {
        ar & row;
        ar & specCountMap[row];
      }
    }
    
    BOOST_SERIALIZATION_SPLIT_MEMBER()
};
 
class PeakCounts {  
  public:
    
    PeakCounts() : smoothingMode_(0),
        precBinWidth(BinSpectra::kBinWidth), precBinShift(BinSpectra::kBinShift), // TODO: check if prec masses are binned in same way as peaks
        peakBinWidth(BinSpectra::kBinWidth), peakBinShift(BinSpectra::kBinShift),
        nTotalPeaks(0u), maxCharge(3u), relativeToPrecMz(false),
        intThresh(0.0), truncatePeaks(true) {
      peakCountMatrices.resize(maxCharge);
      specCountVectors.resize(maxCharge);
      peakDistCache_.resize(maxCharge);
    }
    
    void setSmoothingMode(int mode) { smoothingMode_ = mode; }
		void setRelativeToPrecMz(double fracPeakBinWidth);
		inline void setPrecBinWidth(double d) { precBinWidth = d; }
		inline void setPrecBinShift(double d) { precBinShift = d; }
		inline void setPeakBinWidth(double d) { peakBinWidth = d; }
		inline void setPeakBinShift(double d) { peakBinShift = d; }
		
		inline unsigned int getMaxCharge() const { return maxCharge; }
		inline double getPrecBinWidth() const { return precBinWidth; }
		inline double getPrecBinShift() const { return precBinShift; }
		inline double getPeakBinWidth() const { return peakBinWidth; }
		inline double getPeakBinShift() const { return peakBinShift; }
		inline bool isRelativeToPrecMz() const { return relativeToPrecMz; }
		
		inline unsigned int getChargeBin(const unsigned int charge) const { return (std::min)(charge - 1, maxCharge - 1); }
		
		inline PeakCountMatrix& getPeakCountMatrix(const unsigned int charge) { return peakCountMatrices.at(getChargeBin(charge)); }
		inline SpectrumCountVector& getSpecCountVector(const unsigned int charge) { return specCountVectors.at(getChargeBin(charge)); }
		inline unsigned int getNumTotalPeaks() { return nTotalPeaks; }
		
		void addSpectrum(std::vector<MZIntensityPair>& mziPairs, 
		    const double precMz, const unsigned int charge, double precMass, 
		    const unsigned int numQueryPeaks);
		void addSpectrum(std::vector<BinnedMZIntensityPair>& mziPairs, 
		    const double precMz, const unsigned int charge, 
		    const unsigned int numQueryPeaks);
		
		void add(PeakCounts& otherPeakCounts);
		void subtract(PeakCounts& otherPeakCounts);
		double getDistance(const PeakCounts& otherPeakCounts) const;
		
		void generateSinglePeakDistribution(double precMz, unsigned int charge, PeakDistribution& distribution);
		void generatePeakDistribution(double precMz, unsigned int charge, PeakDistribution& distribution,
		                               unsigned int numQueryPeaks);
		void generateRelativePeakDistribution(const unsigned int charge, PeakDistribution& distribution) const;
		
		void print(const std::string& resultBaseFN);
		void print(std::ostream& os, const unsigned int charge = 2u);
		
		void readFromFile(const std::string& peakCountFN);
		void serialize(std::string& peakCountsSerialized);
    void deserialize(std::string& peakCountsSerialized);
		
		static void serializePeakCounts(PeakCounts& peakCounts, std::string& peakCountsSerialized);
    static void deserializePeakCounts(std::string& peakCountsSerialized, PeakCounts& peakCounts);
    static bool peakCountsSerializationUnitTest();
		
  private:
    std::vector<PeakCountMatrix> peakCountMatrices;
    std::vector<SpectrumCountVector> specCountVectors;
    
    std::vector< std::map<unsigned int, PeakDistribution> > peakDistCache_;
    
    int smoothingMode_;
    
    double precBinWidth, precBinShift;
    double peakBinWidth, peakBinShift;
    
    unsigned int nTotalPeaks;
    unsigned int maxCharge;
    bool relativeToPrecMz;
    
    double intThresh;
    bool truncatePeaks;
    
    inline unsigned int getPrecBin(double precMz) const { return BinSpectra::getBin(precMz, precBinWidth, precBinShift); }
    inline unsigned int getPeakBin(double peakMz) const { return BinSpectra::getBin(peakMz, peakBinWidth, peakBinShift); }
    inline double getPrecMz(unsigned int precBin) const { return BinSpectra::getMZ(precBin, precBinWidth, precBinShift); }
    inline unsigned int getRelPeakBin(double peakMz, double precMz, double fracPeakBinWidth) const { 
      return static_cast<unsigned int>(peakMz / precMz / fracPeakBinWidth + 0.5); 
    }
    
		bool equalParams(const PeakCounts& otherPeakCounts) const;
		
		friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
      ar & precBinWidth;
      ar & precBinShift;
      ar & peakBinWidth;
      ar & peakBinShift;
      ar & nTotalPeaks;
      ar & maxCharge;
      ar & truncatePeaks;
      ar & relativeToPrecMz;
      for (unsigned int charge = 0; charge < maxCharge; ++charge) {
        ar & peakCountMatrices[charge];
        ar & specCountVectors[charge];
      }
    }
    
};

#endif // PEAK_COUNTS_H
