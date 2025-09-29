// File: photonicsDevice.h
// PHOTONICSeval Simulator - PHOTONICS Device
// Copyright (c) 2024 University of Virginia
// This file is licensed under the MIT License.
// See the LICENSE file in the root of this repository for more details.

#ifndef LAVA_PHOTONICS_DEVICE_H
#define LAVA_PHOTONICS_DEVICE_H

#include "libphotonicseval.h"
#include "photonicsSimConfig.h"
#include "photonicsCore.h"
#include "photonicsCmd.h"
#include "photonicsPerfEnergyBase.h"
#ifdef DRAMSIM3_INTEG
#include "cpu.h"
#endif
#include <memory>

class photonicsResMgr;


//! @class  photonicsDevice
//! @brief  PHOTONICS device
class photonicsDevice
{
public:
  photonicsDevice(const photonicsSimConfig& config);
  ~photonicsDevice();

  const photonicsSimConfig& getConfig() const { return m_config; }

  PhotonicsDeviceEnum getDeviceType() const { return m_config.getDeviceType(); }
  PhotonicsDeviceEnum getSimTarget() const { return m_config.getSimTarget(); }
  unsigned getNumRanks() const { return m_config.getNumRanks(); }
  unsigned getNumBankPerRank() const { return m_config.getNumBankPerRank(); }
  unsigned getNumSubarrayPerBank() const { return m_config.getNumSubarrayPerBank(); }
  unsigned getNumRowPerSubarray() const { return m_config.getNumRowPerSubarray(); }
  unsigned getNumColPerSubarray() const { return m_config.getNumColPerSubarray(); }
  unsigned getOnChipBufferSize() const { return m_config.getBufferSize(); }

  unsigned getNumCores() const { return m_numCores; }
  unsigned getNumRows() const { return m_numRows; }
  unsigned getNumCols() const { return m_numCols; }
  unsigned getBufferSize() const { return m_bufferSize; }
  bool isValid() const { return m_isValid; }

  bool isVLayoutDevice() const;
  bool isHLayoutDevice() const;
  bool isHybridLayoutDevice() const;

  PhotonicsObjId photonicsAlloc(PhotonicsAllocEnum allocType, uint64_t numElements, PhotonicsDataType dataType);
  PhotonicsObjId photonicsAllocMat(uint64_t numElements, PhotonicsDataType dataType);
  PhotonicsObjId photonicsAllocAssociated(PhotonicsObjId assocId, PhotonicsDataType dataType);
  PhotonicsObjId photonicsAllocAssociatedSrcVec(PhotonicsObjId assocId, PhotonicsDataType dataType);
  PhotonicsObjId photonicsAllocAssociatedDestVec(PhotonicsObjId assocId, PhotonicsDataType dataType);
  PhotonicsObjId photonicsAllocBuffer(uint32_t numElements, PhotonicsDataType dataType);
  bool photonicsFree(PhotonicsObjId obj, unsigned numCores);
  PhotonicsObjId photonicsCreateRangedRef(PhotonicsObjId refId, uint64_t idxBegin, uint64_t idxEnd);
  PhotonicsObjId photonicsCreateDualContactRef(PhotonicsObjId refId);

  bool photonicsCopyMainToDevice(void* src, PhotonicsObjId dest, uint64_t idxBegin = 0, uint64_t idxEnd = 0);
  bool photonicsCopyDeviceToMain(PhotonicsObjId src, void* dest, uint64_t idxBegin = 0, uint64_t idxEnd = 0);
  bool photonicsCopyMainToDeviceWithType(PhotonicsCopyEnum copyType, void* src, PhotonicsObjId dest, uint64_t idxBegin = 0, uint64_t idxEnd = 0);
  bool photonicsCopyDeviceToMainWithType(PhotonicsCopyEnum copyType, PhotonicsObjId src, void* dest, uint64_t idxBegin = 0, uint64_t idxEnd = 0);
  bool photonicsCopyDeviceToDevice(PhotonicsObjId src, PhotonicsObjId dest, uint64_t idxBegin = 0, uint64_t idxEnd = 0);

  photonicsResMgr* getResMgr() { return m_resMgr.get(); }
  photonicsPerfEnergyBase* getPerfEnergyModel() { return m_perfEnergyModel.get(); }
  photonicsCore& getCore(PhotonicsCoreId coreId) { return m_cores[coreId]; }
  bool executeCmd(std::unique_ptr<photonicsCmd> cmd);

private:
  bool init();
  bool adjustConfigForSimTarget(unsigned& numRanks, unsigned& numBankPerRank, unsigned& numSubarrayPerBank, unsigned& numRows, unsigned& numCols);

  const photonicsSimConfig& m_config;
  unsigned m_numCores = 0;
  unsigned m_numRows = 0;
  unsigned m_numCols = 0;
  unsigned m_bufferSize = 0;
  bool m_isValid = false;
  bool m_isInit = false;
  std::unique_ptr<photonicsResMgr> m_resMgr;
  std::unique_ptr<photonicsPerfEnergyBase> m_perfEnergyModel;
  std::vector<photonicsCore> m_cores;

#ifdef DRAMSIM3_INTEG
  dramsim3::PHOTONICSCPU* m_hostMemory = nullptr;
  dramsim3::PHOTONICSCPU* m_deviceMemory = nullptr;
  dramsim3::Config* m_deviceMemoryConfig = nullptr;
#endif
};

#endif

