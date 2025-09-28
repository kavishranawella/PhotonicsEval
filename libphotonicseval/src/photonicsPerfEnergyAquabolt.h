// File: photonicsPerfEnergyAquabolt.h
// PHOTONICSeval Simulator - Performance Energy Models
// Copyright (c) 2024 University of Virginia
// This file is licensed under the MIT License.
// See the LICENSE file in the root of this repository for more details.

#ifndef LAVA_PHOTONICS_PERF_ENERGY_AQUABOLT_H
#define LAVA_PHOTONICS_PERF_ENERGY_AQUABOLT_H

#include "libphotonicseval.h"                // for PhotonicsDeviceEnum, PhotonicsDataType
#include "photonicsParamsDram.h"             // for photonicsParamsDram
#include "photonicsCmd.h"                    // for PhotonicsCmdEnum
#include "photonicsResMgr.h"                 // for photonicsObjInfo
#include "photonicsPerfEnergyBase.h"         // for photonicsPerfEnergyBase


//! @class  photonicsPerfEnergyAquabolt
//! @brief  PHOTONICS performance energy model for bank-level PHOTONICS family
class photonicsPerfEnergyAquabolt : public photonicsPerfEnergyBase
{
public:
  photonicsPerfEnergyAquabolt(const photonicsPerfEnergyModelParams& params) : photonicsPerfEnergyBase(params) {}
  virtual ~photonicsPerfEnergyAquabolt() {}

  virtual photonicseval::perfEnergy getPerfEnergyForFunc1(PhotonicsCmdEnum cmdType, const photonicsObjInfo& objSrc, const photonicsObjInfo& objDest) const override;
  virtual photonicseval::perfEnergy getPerfEnergyForFunc2(PhotonicsCmdEnum cmdType, const photonicsObjInfo& objSrc1, const photonicsObjInfo& objSrc2, const photonicsObjInfo& objDest) const override;
  virtual photonicseval::perfEnergy getPerfEnergyForReduction(PhotonicsCmdEnum cmdType, const photonicsObjInfo& obj, unsigned numPass) const override;
  virtual photonicseval::perfEnergy getPerfEnergyForBroadcast(PhotonicsCmdEnum cmdType, const photonicsObjInfo& obj) const override;
  virtual photonicseval::perfEnergy getPerfEnergyForRotate(PhotonicsCmdEnum cmdType, const photonicsObjInfo& obj) const override;
  
protected:
  unsigned m_aquaboltFPUBitWidth = 16;
  // TODO: Update for Aquabolt
  double m_aquaboltArithmeticEnergy = 0.0000000004992329586; // mJ
};

#endif

