// File: photonicsPerfEnergyBankLevel.h
// PHOTONICSeval Simulator - Performance Energy Models
// Copyright (c) 2024 University of Virginia
// This file is licensed under the MIT License.
// See the LICENSE file in the root of this repository for more details.

#ifndef LAVA_PHOTONICS_PERF_ENERGY_BANK_LEVEL_H
#define LAVA_PHOTONICS_PERF_ENERGY_BANK_LEVEL_H

#include "libphotonicseval.h"                // for PhotonicsDeviceEnum, PhotonicsDataType
#include "photonicsParamsDram.h"             // for photonicsParamsDram
#include "photonicsCmd.h"                    // for PhotonicsCmdEnum
#include "photonicsResMgr.h"                 // for photonicsObjInfo
#include "photonicsPerfEnergyBase.h"         // for photonicsPerfEnergyBase


//! @class  photonicsPerfEnergyBankLevel
//! @brief  PHOTONICS performance energy model for bank-level PHOTONICS family
class photonicsPerfEnergyBankLevel : public photonicsPerfEnergyBase
{
public:
  photonicsPerfEnergyBankLevel(const photonicsPerfEnergyModelParams& params) : photonicsPerfEnergyBase(params) {}
  virtual ~photonicsPerfEnergyBankLevel() {}

  virtual photonicseval::perfEnergy getPerfEnergyForFunc1(PhotonicsCmdEnum cmdType, const photonicsObjInfo& objSrc, const photonicsObjInfo& objDest) const override;
  virtual photonicseval::perfEnergy getPerfEnergyForFunc2(PhotonicsCmdEnum cmdType, const photonicsObjInfo& objSrc1, const photonicsObjInfo& objSrc2, const photonicsObjInfo& objDest) const override;
  virtual photonicseval::perfEnergy getPerfEnergyForReduction(PhotonicsCmdEnum cmdType, const photonicsObjInfo& obj, unsigned numPass) const override;
  virtual photonicseval::perfEnergy getPerfEnergyForBroadcast(PhotonicsCmdEnum cmdType, const photonicsObjInfo& obj) const override;
  virtual photonicseval::perfEnergy getPerfEnergyForRotate(PhotonicsCmdEnum cmdType, const photonicsObjInfo& obj) const override;
  virtual photonicseval::perfEnergy getPerfEnergyForPrefixSum(PhotonicsCmdEnum cmdType, const photonicsObjInfo& obj) const override;

protected:
  double m_blimpLatency = m_tCCD_L * m_tCK;
  unsigned m_blimpCoreBitWidth = m_GDLWidth;
  unsigned m_simdUnitCount = m_blimpCoreBitWidth / 32; // 32-bit SIMD unit
  // Following values are taken from fulcrum paper as BLIMP paper does not model energy
  double m_blimpArithmeticEnergy = 0.0000000004992329586 * m_simdUnitCount; // mJ
  double m_blimpLogicalEnergy = 0.0000000001467846411 * m_simdUnitCount; // mJ
};

#endif

