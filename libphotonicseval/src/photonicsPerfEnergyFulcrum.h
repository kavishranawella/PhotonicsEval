// File: photonicsPerfEnergyFulcrum.h
// PHOTONICSeval Simulator - Performance Energy Models
// Copyright (c) 2024 University of Virginia
// This file is licensed under the MIT License.
// See the LICENSE file in the root of this repository for more details.

#ifndef LAVA_PHOTONICS_PERF_ENERGY_FULCRUM_H
#define LAVA_PHOTONICS_PERF_ENERGY_FULCRUM_H

#include "libphotonicseval.h"                // for PhotonicsDeviceEnum, PhotonicsDataType
#include "photonicsParamsDram.h"             // for photonicsParamsDram
#include "photonicsCmd.h"                    // for PhotonicsCmdEnum
#include "photonicsResMgr.h"                 // for photonicsObjInfo
#include "photonicsPerfEnergyBase.h"         // for photonicsPerfEnergyBase


//! @class  photonicsPerfEnergyBitFulcrum
//! @brief  PHOTONICS performance energy model for Fulcrum family
class photonicsPerfEnergyFulcrum : public photonicsPerfEnergyBase
{
public:
  photonicsPerfEnergyFulcrum(const photonicsPerfEnergyModelParams& params) : photonicsPerfEnergyBase(params) {}
  virtual ~photonicsPerfEnergyFulcrum() {}

  virtual photonicseval::perfEnergy getPerfEnergyForFunc1(PhotonicsCmdEnum cmdType, const photonicsObjInfo& objSrc, const photonicsObjInfo& objDest) const override;
  virtual photonicseval::perfEnergy getPerfEnergyForFunc2(PhotonicsCmdEnum cmdType, const photonicsObjInfo& objSrc1, const photonicsObjInfo& objSrc2, const photonicsObjInfo& objDest) const override;
  virtual photonicseval::perfEnergy getPerfEnergyForIter(const photonicsObjInfo& objSrc1, const photonicsObjInfo& objSrc2, const photonicsObjInfo& objDest, int8_t numLoops = 1) const override;
  virtual photonicseval::perfEnergy getPerfEnergyForReduction(PhotonicsCmdEnum cmdType, const photonicsObjInfo& obj, unsigned numPass) const override;
  virtual photonicseval::perfEnergy getPerfEnergyForBroadcast(PhotonicsCmdEnum cmdType, const photonicsObjInfo& obj) const override;
  virtual photonicseval::perfEnergy getPerfEnergyForRotate(PhotonicsCmdEnum cmdType, const photonicsObjInfo& obj) const override;
  virtual photonicseval::perfEnergy getPerfEnergyForPrefixSum(PhotonicsCmdEnum cmdType, const photonicsObjInfo& obj) const override;

protected:
  double m_fulcrumMulLatency = 0.00000609; // 6.09ns
  double m_fulcrumAddLatency = 0.00000120; // 1.20ns
  unsigned m_fulcrumAluBitWidth = 32;
  // Following values are taken from fulcrum paper.
  double m_fulcrumMulEnergy = 0.0000000004992329586; // mJ
  double m_fulcrumAddEnergy = 0.0000000001467846411; // mJ
  double m_fulcrumShiftEnergy = 0.0000000075; // mJ
};

#endif

