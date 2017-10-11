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

#ifndef MARACLUSTER_SPARSEEDGE_H_
#define MARACLUSTER_SPARSEEDGE_H_

#include "ScanId.h"

namespace maracluster {

struct SparseEdge {
  SparseEdge(double _value, ScanId _row, ScanId _col) : 
      value(_value), row(_row), col(_col) {}
  SparseEdge() : value(0.0), row(), col() {}
  
  double value;
  ScanId row, col;
  
  inline bool operator<(const SparseEdge& r2) const {
    return value > r2.value
      || (value == r2.value && row < r2.row)
      || (value == r2.value && row == r2.row && col < r2.col);
  }
  
  bool sameEdge(const SparseEdge& r2) {
    return (row == r2.row && col == r2.col);
  }
};

struct SparseMissingEdge : public SparseEdge {
  SparseMissingEdge(double _value, ScanId _row, ScanId _col, 
             unsigned int _numEdges) : SparseEdge(_value, _row, _col), 
                                       numEdges(_numEdges) {}
  SparseMissingEdge() : SparseEdge(), numEdges(1) {}
  
  unsigned int numEdges;
};

} /* namespace maracluster */

#endif /* MARACLUSTER_SPARSEEDGE_H_ */
