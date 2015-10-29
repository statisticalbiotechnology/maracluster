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
 
#ifndef MASS_CHARGE_CANDIDATE_H
#define MASS_CHARGE_CANDIDATE_H

struct MassChargeCandidate {
  MassChargeCandidate(unsigned int _charge, double _precMz, double _mass) : charge(_charge), precMz(_precMz), mass(_mass) {}
  
  unsigned int charge;
  double precMz;
  double mass;
  
  unsigned int getMassInt() const { return static_cast<unsigned int>(mass); }
  double getMassMod() const { return mass - getMassInt(); }
  
  inline static bool lessMass(const MassChargeCandidate& a, const MassChargeCandidate& b) { return (a.mass < b.mass); }
  inline static bool lessChargeMass(const MassChargeCandidate& a, const MassChargeCandidate& b) { 
      return (a.charge < b.charge) || 
             ((a.charge == b.charge) && (a.mass < b.mass)); 
  }
};

#endif // MASS_CHARGE_CANDIDATE_H
