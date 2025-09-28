// File: photonicsPerfEnergyBitSerial.h
// PHOTONICSeval Simulator - Performance Energy Models
// Copyright (c) 2024 University of Virginia
// This file is licensed under the MIT License.
// See the LICENSE file in the root of this repository for more details.

#ifndef LAVA_PHOTONICS_PERF_ENERGY_BIT_SERIAL_H
#define LAVA_PHOTONICS_PERF_ENERGY_BIT_SERIAL_H

#include "libphotonicseval.h"                // for PhotonicsDeviceEnum, PhotonicsDataType
#include "photonicsParamsDram.h"             // for photonicsParamsDram
#include "photonicsCmd.h"                    // for PhotonicsCmdEnum
#include "photonicsResMgr.h"                 // for photonicsObjInfo
#include "photonicsPerfEnergyBase.h"         // for photonicsPerfEnergyBase


//! @class  photonicsPerfEnergyBitSerial
//! @brief  PHOTONICS performance energy model for bit-serial PHOTONICS family
class photonicsPerfEnergyBitSerial : public photonicsPerfEnergyBase
{
public:
  photonicsPerfEnergyBitSerial(const photonicsPerfEnergyModelParams& params) : photonicsPerfEnergyBase(params) {}
  virtual ~photonicsPerfEnergyBitSerial() {}

  virtual photonicseval::perfEnergy getPerfEnergyForFunc1(PhotonicsCmdEnum cmdType, const photonicsObjInfo& objSrc, const photonicsObjInfo& objDest) const override;
  virtual photonicseval::perfEnergy getPerfEnergyForFunc2(PhotonicsCmdEnum cmdType, const photonicsObjInfo& objSrc1, const photonicsObjInfo& objSrc2, const photonicsObjInfo& objDest) const override;
  virtual photonicseval::perfEnergy getPerfEnergyForReduction(PhotonicsCmdEnum cmdType, const photonicsObjInfo& obj, unsigned numPass) const override;
  virtual photonicseval::perfEnergy getPerfEnergyForBroadcast(PhotonicsCmdEnum cmdType, const photonicsObjInfo& obj) const override;
  virtual photonicseval::perfEnergy getPerfEnergyForRotate(PhotonicsCmdEnum cmdType, const photonicsObjInfo& obj) const override;
  virtual photonicseval::perfEnergy getPerfEnergyForPrefixSum(PhotonicsCmdEnum cmdType, const photonicsObjInfo& obj) const override;

protected:
  photonicseval::perfEnergy getPerfEnergyBitSerial(PhotonicsDeviceEnum deviceType, PhotonicsCmdEnum cmdType, unsigned numPass, const photonicsObjInfo& objSrc1, const photonicsObjInfo& objSrc2, const photonicsObjInfo& objDest) const;
  photonicseval::perfEnergy getPerfEnergyTypeConversion(PhotonicsDeviceEnum deviceType, PhotonicsCmdEnum cmdType, const photonicsObjInfo& objSrc, const photonicsObjInfo& objDest) const;

  // Popcount logc Params from DRAM-CAM paper
  double m_pclNsDelay = 0.76; // 64-bit popcount logic ns delay, using LUT no pipeline design
  double m_pclUwPower = 0.03; // 64-bit popcount logic uW power, using LUT no pipeline design
};

#endif

