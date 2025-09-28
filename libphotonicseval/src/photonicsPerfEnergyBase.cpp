// File: photonicsPerfEnergyBase.cc
// PHOTONICSeval Simulator - Performance Energy Models
// Copyright (c) 2024 University of Virginia
// This file is licensed under the MIT License.
// See the LICENSE file in the root of this repository for more details.

#include "photonicsPerfEnergyBase.h"
#include "photonicsCmd.h"
#include "photonicsPerfEnergyBitSerial.h"
#include "photonicsPerfEnergyFulcrum.h"
#include "photonicsPerfEnergyBankLevel.h"
#include "photonicsPerfEnergyAquabolt.h"
#include "photonicsPerfEnergyAim.h"
#include <cstdint>
#include <cstdio>


//! @brief  A factory function to create perf energy model for sim target
std::unique_ptr<photonicsPerfEnergyBase>
photonicsPerfEnergyFactory::createPerfEnergyModel(const photonicsPerfEnergyModelParams& params)
{
  switch (params.getSimTarget()) {
    case PHOTONICS_DEVICE_BITSIMD_V:
    case PHOTONICS_DEVICE_BITSIMD_V_AP:
    case PHOTONICS_DEVICE_BITSIMD_H:
    case PHOTONICS_DEVICE_SIMDRAM:
      printf("PHOTONICS-Info: Created performance energy model for bit-serial PHOTONICS\n");
      return std::make_unique<photonicsPerfEnergyBitSerial>(params);
    case PHOTONICS_DEVICE_FULCRUM:
      printf("PHOTONICS-Info: Created performance energy model for Fulcrum\n");
      return std::make_unique<photonicsPerfEnergyFulcrum>(params);
    case PHOTONICS_DEVICE_BANK_LEVEL:
      printf("PHOTONICS-Info: Created performance energy model for bank-level PHOTONICS\n");
      return std::make_unique<photonicsPerfEnergyBankLevel>(params);
    case PHOTONICS_DEVICE_AQUABOLT:
      printf("PHOTONICS-Info: Created performance energy model for AQUABOLT\n");
      return std::make_unique<photonicsPerfEnergyAquabolt>(params);
    case PHOTONICS_DEVICE_AIM:
      printf("PHOTONICS-Info: Created performance energy model for AiM\n");
      return std::make_unique<photonicsPerfEnergyAim>(params);
    default:
      printf("PHOTONICS-Warning: Created performance energy base model for unrecognized simulation target\n");
  }
  return std::make_unique<photonicsPerfEnergyBase>(params);
}


//! @brief  photonicsPerfEnergyBase ctor
photonicsPerfEnergyBase::photonicsPerfEnergyBase(const photonicsPerfEnergyModelParams& params)
  : m_simTarget(params.getSimTarget()),
    m_numRanks(params.getNumRanks()),
    m_paramsDram(params.getParamsDram())
{
  m_tR = m_paramsDram.getNsRowRead() / m_nano_to_milli;
  m_tW = m_paramsDram.getNsRowWrite() / m_nano_to_milli;
  m_tACT = m_paramsDram.getNsRowActivate() / m_nano_to_milli; // Row activate latency in ms
  m_tPRE = m_paramsDram.getNsRowPrecharge() / m_nano_to_milli; // Row precharge latency in ms
  m_tL = m_paramsDram.getNsTCCD_S() / m_nano_to_milli;
  m_tGDL = m_paramsDram.getNsTCCD_L() / m_nano_to_milli;
  m_eAP = m_paramsDram.getPjActPre() / m_pico_to_milli; // Convert pJ to mJ
  m_eL = m_paramsDram.getPjLogic() / m_pico_to_milli; // Convert pJ to mJ
  m_eR = m_paramsDram.getPjRead() / m_pico_to_milli;
  m_eW = m_paramsDram.getPjWrite() / m_pico_to_milli;
  m_eACT = m_paramsDram.getPjActivate() / m_pico_to_milli; // Convert pJ to mJ
  m_ePRE = m_paramsDram.getPjPrecharge() / m_pico_to_milli; // Convert pJ to mJ
  // m_pBCore = (m_paramsDram.getMwIDD3N() - m_paramsDram.getMwIDD2N()) / 1000.0; // Convert mW to W, so that W * ms = mJ
  m_pBChip = m_paramsDram.getMwIDD3N() / 1000.0; // Convert mW to W, so that W * ms = mJ
  m_GDLWidth = m_paramsDram.getBurstLength() * m_paramsDram.getDeviceWidth();
  m_numChipsPerRank = m_paramsDram.getNumChipsPerRank();
  m_typicalRankBW = m_paramsDram.getTypicalRankBW(); // GB/s
  m_tCK = m_paramsDram.gettCK() / m_nano_to_milli; // Convert ns to ms 
  m_tCCD_S = m_paramsDram.gettCCD_S();
  m_tCCD_L = m_paramsDram.gettCCD_L();
  m_tRCD = m_paramsDram.gettRCD();
  m_tRP = m_paramsDram.gettRP();
  m_tCAS = m_paramsDram.getNsTCAS() / m_nano_to_milli; // Convert ns to ms
  m_tRAS = m_paramsDram.gettRAS();
}

//! @brief  Perf energy model of data transfer between CPU memory and PHOTONICS memory
photonicseval::perfEnergy
photonicsPerfEnergyBase::getPerfEnergyForBytesTransfer(PhotonicsCmdEnum cmdType, uint64_t numBytes) const
{
  //TODO: fine grain perf-energy modeling 
  double mjEnergy = 0.0;
  double msRead = 0.0;
  double msWrite = 0.0;
  double msCompute = 0.0;
  uint64_t mTotalOP = 0;
  double msRuntime = static_cast<double>(numBytes) / (m_typicalRankBW * m_numRanks * 1024 * 1024 * 1024 / 1000);
  switch (cmdType) {
    case PhotonicsCmdEnum::COPY_H2D:
    {
      mjEnergy = m_eW * msRuntime * m_numChipsPerRank * m_numRanks;
      mjEnergy += m_pBChip * m_numChipsPerRank * m_numRanks * msRuntime;
      break;
    }
    case PhotonicsCmdEnum::COPY_D2H:
    {
      mjEnergy = m_eR * msRuntime * m_numChipsPerRank * m_numRanks;
      mjEnergy += m_pBChip * m_numChipsPerRank * m_numRanks * msRuntime;
      break;
    }
    case PhotonicsCmdEnum::COPY_D2D:
    {
      // One row read, one row write within a subarray
      mjEnergy = m_eAP * 2 * msRuntime * m_numChipsPerRank * m_numRanks;
      mjEnergy += m_pBChip * m_numChipsPerRank * m_numRanks * msRuntime;
      break;
    }
    default:
    {
      printf("PHOTONICS-Warning: Perf energy model not available for PHOTONICS command %s\n", photonicsCmd::getName(cmdType, "").c_str());
      break;
    }
  }
  return photonicseval::perfEnergy(msRuntime, mjEnergy, msRead, msWrite, msCompute, mTotalOP);
}

//! @brief  Perf energy model of base class for func1 (placeholder)
photonicseval::perfEnergy
photonicsPerfEnergyBase::getPerfEnergyForFunc1(PhotonicsCmdEnum cmdType, const photonicsObjInfo& objSrc, const photonicsObjInfo& objDest) const
{
  double msRuntime = 1e10;
  double mjEnergy = 999999999.9;
  double msRead = 0.0;
  double msWrite = 0.0;
  double msCompute = 0.0;
  uint64_t mTotalOP = 0;
  return photonicseval::perfEnergy(msRuntime, mjEnergy, msRead, msWrite, msCompute, mTotalOP);
}

//! @brief  Perf energy model of base class for func2 (placeholder)
photonicseval::perfEnergy
photonicsPerfEnergyBase::getPerfEnergyForFunc2(PhotonicsCmdEnum cmdType, const photonicsObjInfo& objSrc1, const photonicsObjInfo& objSrc2, const photonicsObjInfo& objDest) const
{
  double msRuntime = 1e10;
  double mjEnergy = 999999999.9;
  double msRead = 0.0;
  double msWrite = 0.0;
  double msCompute = 0.0;
  uint64_t mTotalOP = 0;
  return photonicseval::perfEnergy(msRuntime, mjEnergy, msRead, msWrite, msCompute, mTotalOP);
}

//! @brief  Perf energy model of base class for iter (placeholder)
photonicseval::perfEnergy
photonicsPerfEnergyBase::getPerfEnergyForIter(const photonicsObjInfo& objSrc1, const photonicsObjInfo& objSrc2, const photonicsObjInfo& objDest, int8_t numLoops) const
{
  double msRuntime = 1e10 * numLoops;
  double mjEnergy = 999999999.9 * numLoops;
  double msRead = 0.0;
  double msWrite = 0.0;
  double msCompute = 0.0;
  uint64_t mTotalOP = 0;
  return photonicseval::perfEnergy(msRuntime, mjEnergy, msRead, msWrite, msCompute, mTotalOP);
}

//! @brief  Perf energy model of base class for reduction sum (placeholder)
photonicseval::perfEnergy
photonicsPerfEnergyBase::getPerfEnergyForReduction(PhotonicsCmdEnum cmdType, const photonicsObjInfo& obj, unsigned numPass) const
{
  double msRuntime = 1e10;
  double mjEnergy = 999999999.9;
  double msRead = 0.0;
  double msWrite = 0.0;
  double msCompute = 0.0;
  uint64_t mTotalOP = 0;
  return photonicseval::perfEnergy(msRuntime, mjEnergy, msRead, msWrite, msCompute, mTotalOP);
}

//! @brief  Perf energy model of base class for broadcast (placeholder)
photonicseval::perfEnergy
photonicsPerfEnergyBase::getPerfEnergyForBroadcast(PhotonicsCmdEnum cmdType, const photonicsObjInfo& obj) const
{
  double msRuntime = 1e10;
  double mjEnergy = 999999999.9;
  double msRead = 0.0;
  double msWrite = 0.0;
  double msCompute = 0.0;
  uint64_t mTotalOP = 0;
  return photonicseval::perfEnergy(msRuntime, mjEnergy, msRead, msWrite, msCompute, mTotalOP);
}

//! @brief  Perf energy model of base class for rotate (placeholder)
photonicseval::perfEnergy
photonicsPerfEnergyBase::getPerfEnergyForRotate(PhotonicsCmdEnum cmdType, const photonicsObjInfo& obj) const
{
  double msRuntime = 1e10;
  double mjEnergy = 999999999.9;
  double msRead = 0.0;
  double msWrite = 0.0;
  double msCompute = 0.0;
  uint64_t mTotalOP = 0;
  return photonicseval::perfEnergy(msRuntime, mjEnergy, msRead, msWrite, msCompute, mTotalOP);
}

//! @brief  Perf energy model of base class for prefixsum
photonicseval::perfEnergy
photonicsPerfEnergyBase::getPerfEnergyForPrefixSum(PhotonicsCmdEnum cmdType, const photonicsObjInfo& obj) const
{
  double msRuntime = 1e10;
  double mjEnergy = 999999999.9;
  double msRead = 0.0;
  double msWrite = 0.0;
  double msCompute = 0.0;
  uint64_t mTotalOP = 0;
  return photonicseval::perfEnergy(msRuntime, mjEnergy, msRead, msWrite, msCompute, mTotalOP);
}

//! @brief  Perf energy model of base class for MAC
photonicseval::perfEnergy
photonicsPerfEnergyBase::getPerfEnergyForMac(PhotonicsCmdEnum cmdType, const photonicsObjInfo& obj) const
{
  double msRuntime = 1e10;
  double mjEnergy = 999999999.9;
  double msRead = 0.0;
  double msWrite = 0.0;
  double msCompute = 0.0;
  uint64_t mTotalOP = 0;
  return photonicseval::perfEnergy(msRuntime, mjEnergy, msRead, msWrite, msCompute, mTotalOP);
}
