// File: photonicsPerfEnergyBase.h
// PHOTONICSeval Simulator - Performance Energy Models
// Copyright (c) 2024 University of Virginia
// This file is licensed under the MIT License.
// See the LICENSE file in the root of this repository for more details.

#ifndef LAVA_PHOTONICS_PERF_ENERGY_BASE_H
#define LAVA_PHOTONICS_PERF_ENERGY_BASE_H

#include "libphotonicseval.h"                // for PhotonicsDeviceEnum, PhotonicsDataType
#include "photonicsParamsDram.h"             // for photonicsParamsDram
#include "photonicsCmd.h"                    // for PhotonicsCmdEnum
#include "photonicsResMgr.h"                 // for photonicsObjInfo
#include <cstdint>
#include <memory>                      // for std::unique_ptr


namespace photonicseval {
  class perfEnergy
  {
    public:
      perfEnergy() : m_msRuntime(0.0), m_mjEnergy(0.0), m_msRead(0.0), m_msWrite(0.0), m_msCompute(0.0), m_totalOp(0) {}
      perfEnergy(double msRuntime, double mjEnergy, double msRead, double msWrite, double msCompute, uint64_t totalOp) : m_msRuntime(msRuntime), m_mjEnergy(mjEnergy), m_msRead(msRead), m_msWrite(msWrite), m_msCompute(msCompute), m_totalOp(totalOp)  {}

      double m_msRuntime;
      double m_mjEnergy;
      double m_msRead;
      double m_msWrite;
      double m_msCompute;
      uint64_t m_totalOp;
  };
}

//! @class  photonicsPerfEnergyModelParams
//! @brief  Parameters for creating perf energy models
class photonicsPerfEnergyModelParams
{
public:
  photonicsPerfEnergyModelParams(PhotonicsDeviceEnum simTarget, unsigned numRanks, const photonicsParamsDram& paramsDram)
    : m_simTarget(simTarget), m_numRanks(numRanks), m_paramsDram(paramsDram) {}
  PhotonicsDeviceEnum getSimTarget() const { return m_simTarget; }
  unsigned getNumRanks() const { return m_numRanks; }
  const photonicsParamsDram& getParamsDram() const { return m_paramsDram; }
private:
  PhotonicsDeviceEnum m_simTarget;
  unsigned m_numRanks;
  const photonicsParamsDram& m_paramsDram;
};

//! @class  photonicsPerfEnergyFactory
//! @brief  PHOTONICS performance energy model factory
class photonicsPerfEnergyBase;
class photonicsPerfEnergyFactory
{
public:
  // photonicsDevice is not fully constructed at this point. Do not call photonicsSim::get() in photonicsPerfEnergyBase ctor.
  static std::unique_ptr<photonicsPerfEnergyBase> createPerfEnergyModel(const photonicsPerfEnergyModelParams& params);
};

//! @class  photonicsPerfEnergyBase
//! @brief  PHOTONICS performance energy model base class
class photonicsPerfEnergyBase
{
public:
  photonicsPerfEnergyBase(const photonicsPerfEnergyModelParams& params);
  virtual ~photonicsPerfEnergyBase() {}

  virtual photonicseval::perfEnergy getPerfEnergyForBytesTransfer(PhotonicsCmdEnum cmdType, uint64_t numBytes) const;
  virtual photonicseval::perfEnergy getPerfEnergyForFunc1(PhotonicsCmdEnum cmdType, const photonicsObjInfo& objSrc, const photonicsObjInfo& objDest) const;
  virtual photonicseval::perfEnergy getPerfEnergyForFunc2(PhotonicsCmdEnum cmdType, const photonicsObjInfo& objSrc1, const photonicsObjInfo& objSrc2, const photonicsObjInfo& objDest) const;
  virtual photonicseval::perfEnergy getPerfEnergyForIter(const photonicsObjInfo& objSrc1, const photonicsObjInfo& objSrc2, const photonicsObjInfo& objDest, int8_t numLoops = 1) const;
  virtual photonicseval::perfEnergy getPerfEnergyForReduction(PhotonicsCmdEnum cmdType, const photonicsObjInfo& obj, unsigned numPass) const;
  virtual photonicseval::perfEnergy getPerfEnergyForBroadcast(PhotonicsCmdEnum cmdType, const photonicsObjInfo& obj) const;
  virtual photonicseval::perfEnergy getPerfEnergyForRotate(PhotonicsCmdEnum cmdType, const photonicsObjInfo& obj) const;
  virtual photonicseval::perfEnergy getPerfEnergyForPrefixSum(PhotonicsCmdEnum cmdType, const photonicsObjInfo& obj) const;
  virtual photonicseval::perfEnergy getPerfEnergyForMac(PhotonicsCmdEnum cmdType, const photonicsObjInfo& obj) const;

protected:
  PhotonicsDeviceEnum m_simTarget;
  unsigned m_numRanks;
  const photonicsParamsDram& m_paramsDram;

  const double m_nano_to_milli = 1000000.0;
  const double m_pico_to_milli = 1000000000.0;
  double m_tR; // Row read latency in ms
  double m_tW; // Row write latency in ms
  double m_tACT; // Row read(ACT) latency in ms
  double m_tPRE; // Row precharge latency in ms
  double m_tL; // Logic operation for bitserial / tCCD in ms
  double m_tGDL; // Fetch data from local row buffer to global row buffer
  double m_tCAS; // CAS time in ms
  int m_GDLWidth; // Number of bits that can be fetched from local to global row buffer.
  int m_numChipsPerRank; // Number of chips per rank
  double m_typicalRankBW; // typical rank data transfer bandwidth in GB/s

  double m_eAP; // Row read(ACT) energy in mJ microjoule
  double m_eL; // Logic energy in mJ microjoule
  double m_eR; // local row buffer to global row buffer
  double m_eW; // global row buffer to local row buffer
  double m_eACT; // Row activate energy in mJ 
  double m_ePRE; // Row precharge energy in mJ
  double m_pBCore; // background power for each core in W
  double m_pBChip; // background power for each core in W
  double m_tCK; // Clock cycle time in ms
  unsigned m_tCCD_S; // Short command delay in cycles
  unsigned m_tCCD_L; // Long command delay in cycles
  unsigned m_tRCD; // RCD in cycles
  unsigned m_tRP; // RP in cycles
  unsigned m_tRAS; // RAS in cycles
};

#endif

