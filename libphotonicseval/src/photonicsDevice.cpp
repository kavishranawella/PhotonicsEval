// File: photonicsDevice.cpp
// PHOTONICSeval Simulator - PHOTONICS Device
// Copyright (c) 2024 University of Virginia
// This file is licensed under the MIT License.
// See the LICENSE file in the root of this repository for more details.

#include "photonicsDevice.h"
#include "photonicsResMgr.h"
#include "photonicsSim.h"
#include "libphotonicseval.h"
#include "photonicsUtils.h"
#include <cstdio>
#include <memory>
#include <cassert>
#include <string>


//! @brief  photonicsDevice ctor
photonicsDevice::photonicsDevice(const photonicsSimConfig& config)
  : m_config(config)
{
  init();
}

//! @brief  photonicsDevice dtor
photonicsDevice::~photonicsDevice()
{
}

//! @brief  Adjust config for modeling different simulation target with same inputs
bool
photonicsDevice::adjustConfigForSimTarget(unsigned& numRanks, unsigned& numBankPerRank, unsigned& numSubarrayPerBank, unsigned& numRows, unsigned& numCols)
{
  switch (getSimTarget()) {
  case PHOTONICS_DEVICE_BITSIMD_V:
  case PHOTONICS_DEVICE_BITSIMD_V_NAND:
  case PHOTONICS_DEVICE_BITSIMD_V_MAJ:
  case PHOTONICS_DEVICE_BITSIMD_V_AP:
  case PHOTONICS_DEVICE_DRISA_NOR:
  case PHOTONICS_DEVICE_DRISA_MIXED:
  case PHOTONICS_DEVICE_SIMDRAM:
  case PHOTONICS_DEVICE_BITSIMD_H:
    std::printf("PHOTONICS-Info: Aggregate every two subarrays as a single core\n");
    if (numSubarrayPerBank % 2 != 0) {
      std::printf("PHOTONICS-Error: Please config even number of subarrays in each bank\n");
      return false;
    }
    numRows *= 2;
    numSubarrayPerBank /= 2;
    break;
  
  case PHOTONICS_DEVICE_FULCRUM:
    break;
  case PHOTONICS_DEVICE_BANK_LEVEL:
  case PHOTONICS_DEVICE_AIM:
    std::printf("PHOTONICS-Info: Aggregate all subarrays within a bank as a single core\n");
    numRows *= numSubarrayPerBank;
    numSubarrayPerBank = 1;
    break;
  case PHOTONICS_DEVICE_AQUABOLT:
    std::printf("PHOTONICS-Info: Aggregate all subarrays of two consecutive banks as a single core\n");
    if (numBankPerRank % 2 != 0) {
      std::printf("PHOTONICS-Error: Number of banks must be an even number\n");
      return false;
    }
    numRows *= numSubarrayPerBank*2;
    numSubarrayPerBank = 1;
    numBankPerRank /= 2; 
    break;
  default:
    assert(0);
  }
  return true;
}

//! @brief  If a PHOTONICS device uses vertical data layout
bool
photonicsDevice::isVLayoutDevice() const
{
  return photonicsUtils::getDeviceDataLayout(getSimTarget()) == PhotonicsDataLayout::V;
}

//! @brief  If a PHOTONICS device uses horizontal data layout
bool
photonicsDevice::isHLayoutDevice() const
{
  return photonicsUtils::getDeviceDataLayout(getSimTarget()) == PhotonicsDataLayout::H;
}

//! @brief  If a PHOTONICS device uses hybrid data layout
bool
photonicsDevice::isHybridLayoutDevice() const
{
  return photonicsUtils::getDeviceDataLayout(getSimTarget()) == PhotonicsDataLayout::HYBRID;
}

//! @brief  Init PHOTONICS device
bool
photonicsDevice::init()
{
  assert(!m_isInit);

  // Adjust dimension for simulation, e.g., with subarray aggregation
  unsigned numRanks = getNumRanks();
  unsigned numBankPerRank = getNumBankPerRank();
  unsigned numSubarrayPerBank = getNumSubarrayPerBank();
  unsigned numRows = getNumRowPerSubarray();
  unsigned numCols = getNumColPerSubarray();
  unsigned bufferSize = getOnChipBufferSize();
  if (adjustConfigForSimTarget(numRanks, numBankPerRank, numSubarrayPerBank, numRows, numCols)) {
    m_numCores = numRanks * numBankPerRank * numSubarrayPerBank;
    m_numRows = numRows;
    m_numCols = numCols;
    m_bufferSize = bufferSize;
  } else {
    return false;
  }

#ifdef DRAMSIM3_INTEG
  std::string configFile = m_config.getSimConfigFile();
  //TODO: DRAMSim3 requires an output directory but for our purpose we do not need it so sending empty string
  m_deviceMemory = new dramsim3::PHOTONICSCPU(configFile, "");
  m_deviceMemoryConfig = m_deviceMemory->getMemorySystem()->getConfig();
  u_int64_t rowsPerBank = m_deviceMemoryConfig->rows, columnPerRow = m_deviceMemoryConfig->columns * m_deviceMemoryConfig->device_width;

  // todo: adjust for sim target
  m_numCores = 16;
  m_numRows = rowsPerBank/m_numCores;
  m_numCols = columnPerRow;
#endif

  m_isValid = (m_numCores > 0 && m_numRows > 0 && m_numCols > 0);
  if (m_numCols % 8 != 0) {
    std::printf("PHOTONICS-Error: Number of columns %u is not a multiple of 8\n", m_numCols);
    return false;
  }
  if (!m_isValid) {
    std::printf("PHOTONICS-Error: Incorrect device parameters: %u cores, %u rows, %u columns\n", m_numCores, m_numRows, m_numCols);
    return false;
  }

  m_resMgr = std::make_unique<photonicsResMgr>(this);
  const photonicsParamsDram& paramsDram = photonicsSim::get()->getParamsDram(); // created before photonicsDevice ctor
  photonicsPerfEnergyModelParams params(getSimTarget(), getNumRanks(), paramsDram);
  m_perfEnergyModel = photonicsPerfEnergyFactory::createPerfEnergyModel(params);

  // Disable simulated memory creation for functional simulation
  if (getDeviceType() != PHOTONICS_FUNCTIONAL) {
    m_cores.resize(m_numCores, photonicsCore(m_numRows, m_numCols));
  }

  if (getSimTarget() != PHOTONICS_DEVICE_AIM && m_bufferSize > 0) {
    std::printf("PHOTONICS-Error: Device Does not support On-Chip Buffer\n");
    m_isInit = false;
    m_isValid = false;
    return m_isValid;
  }

  std::printf("PHOTONICS-Info: Created PHOTONICS device with %u cores of %u rows and %u columns.\n", m_numCores, m_numRows, m_numCols);

  m_isInit = true;
  return m_isValid;
}

//! @brief  Alloc a PHOTONICS object
PhotonicsObjId
photonicsDevice::photonicsAlloc(PhotonicsAllocEnum allocType, uint64_t numElements, PhotonicsDataType dataType)
{
  if (allocType == PHOTONICS_ALLOC_AUTO) {
    if (isVLayoutDevice()) {
      allocType = PHOTONICS_ALLOC_V;
    } else if (isHLayoutDevice()) {
      allocType = PHOTONICS_ALLOC_H;
    } else {
      assert(0);
    }
  }
  return m_resMgr->photonicsAlloc(allocType, numElements, dataType);
}

//! @brief  Alloc a PHOTONICS object
PhotonicsObjId
photonicsDevice::photonicsAllocMat(uint64_t numElements, PhotonicsDataType dataType)
{
  return m_resMgr->photonicsAllocMat(numElements, dataType);
}

 //! @brief  Allocate a PHOTONICS buffer
PhotonicsObjId
photonicsDevice::photonicsAllocBuffer(uint32_t numElements, PhotonicsDataType dataType)
{
  if (getSimTarget() != PHOTONICS_DEVICE_AIM) {
    std::printf("PHOTONICS-Error: Device does not support On-Chip Buffer\n");
    return -1;
  }
  return m_resMgr->photonicsAllocBuffer(numElements, dataType);
}

//! @brief  Allocate a PHOTONICS object associated with another PHOTONICS object
PhotonicsObjId
photonicsDevice::photonicsAllocAssociated(PhotonicsObjId assocId, PhotonicsDataType dataType)
{
  return m_resMgr->photonicsAllocAssociated(assocId, dataType);
}

//! @brief  Allocate a PHOTONICS object associated with another PHOTONICS object
PhotonicsObjId
photonicsDevice::photonicsAllocAssociatedSrcVec(PhotonicsObjId assocId, PhotonicsDataType dataType)
{
  return m_resMgr->photonicsAllocAssociatedSrcVec(assocId, dataType);
}

//! @brief  Allocate a PHOTONICS object associated with another PHOTONICS object
PhotonicsObjId
photonicsDevice::photonicsAllocAssociatedDestVec(PhotonicsObjId assocId, PhotonicsDataType dataType)
{
  return m_resMgr->photonicsAllocAssociatedDestVec(assocId, dataType);
}

//! @brief  Free a PHOTONICS object
bool
photonicsDevice::photonicsFree(PhotonicsObjId obj)
{
  return m_resMgr->photonicsFree(obj);
}

//! @brief  Create an obj referencing to a range of an existing obj
PhotonicsObjId
photonicsDevice::photonicsCreateRangedRef(PhotonicsObjId refId, uint64_t idxBegin, uint64_t idxEnd)
{
  return m_resMgr->photonicsCreateRangedRef(refId, idxBegin, idxEnd);
}

//! @brief  Create an obj referencing to negation of an existing obj based on dual-contact memory cells
PhotonicsObjId
photonicsDevice::photonicsCreateDualContactRef(PhotonicsObjId refId)
{
  return m_resMgr->photonicsCreateDualContactRef(refId);
}

//! @brief  Copy data from host to PHOTONICS within a range
bool
photonicsDevice::photonicsCopyMainToDevice(void* src, PhotonicsObjId dest, uint64_t idxBegin, uint64_t idxEnd)
{
  PhotonicsCopyEnum copyType = m_resMgr->isHLayoutObj(dest) ? PHOTONICS_COPY_H : PHOTONICS_COPY_V;
  return photonicsCopyMainToDeviceWithType(copyType, src, dest, idxBegin, idxEnd);
}

//! @brief  Copy data from PHOTONICS to host within a range
bool
photonicsDevice::photonicsCopyDeviceToMain(PhotonicsObjId src, void* dest, uint64_t idxBegin, uint64_t idxEnd)
{
  PhotonicsCopyEnum copyType = m_resMgr->isHLayoutObj(src) ? PHOTONICS_COPY_H : PHOTONICS_COPY_V;
  return photonicsCopyDeviceToMainWithType(copyType, src, dest, idxBegin, idxEnd);
}

//! @brief  Copy data from host to PHOTONICS within a range
bool
photonicsDevice::photonicsCopyMainToDeviceWithType(PhotonicsCopyEnum copyType, void* src, PhotonicsObjId dest, uint64_t idxBegin, uint64_t idxEnd)
{
  std::unique_ptr<photonicsCmd> cmd =
    std::make_unique<photonicsCmdCopy>(PhotonicsCmdEnum::COPY_H2D, copyType, src, dest, idxBegin, idxEnd);
  return executeCmd(std::move(cmd));
}

//! @brief  Copy data from PHOTONICS to host within a range
bool
photonicsDevice::photonicsCopyDeviceToMainWithType(PhotonicsCopyEnum copyType, PhotonicsObjId src, void* dest, uint64_t idxBegin, uint64_t idxEnd)
{
  std::unique_ptr<photonicsCmd> cmd =
    std::make_unique<photonicsCmdCopy>(PhotonicsCmdEnum::COPY_D2H, copyType, src, dest, idxBegin, idxEnd);
  return executeCmd(std::move(cmd));
}

//! @brief  Copy data from PHOTONICS to PHOTONICS within a range
bool
photonicsDevice::photonicsCopyDeviceToDevice(PhotonicsObjId src, PhotonicsObjId dest, uint64_t idxBegin, uint64_t idxEnd)
{
  const photonicsObjInfo& obj = m_resMgr->getObjInfo(src);
  PhotonicsCopyEnum copyType = obj.isVLayout() ? PHOTONICS_COPY_V : PHOTONICS_COPY_H;
  std::unique_ptr<photonicsCmd> cmd =
    std::make_unique<photonicsCmdCopy>(PhotonicsCmdEnum::COPY_D2D, copyType, src, dest, idxBegin, idxEnd);
  return executeCmd(std::move(cmd));
}

//! @brief  Execute a PHOTONICS command
bool
photonicsDevice::executeCmd(std::unique_ptr<photonicsCmd> cmd)
{
  cmd->setDevice(this);
  bool ok = cmd->execute();

  return ok;
}

