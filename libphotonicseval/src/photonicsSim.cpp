// File: photonicsSim.cpp
// PHOTONICSeval Simulator - PHOTONICS Simulator Main Entry
// Copyright (c) 2024 University of Virginia
// This file is licensed under the MIT License.
// See the LICENSE file in the root of this repository for more details.

#include "photonicsSim.h"
#include "photonicsCmd.h"
#include "photonicsCmdFuse.h"
#include "photonicsParamsDram.h"
#include "photonicsStats.h"
#include "photonicsUtils.h"
#include <cstdio>
#include <memory>
#include <algorithm>
#include <string>

// The photonicsSim singleton
photonicsSim* photonicsSim::s_instance = nullptr;

//! @brief  Get or create the photonicsSim singleton
photonicsSim*
photonicsSim::get()
{
  if (!s_instance) {
    s_instance = new photonicsSim();
  }
  return s_instance;
}

//! @brief  Destroy the photonicsSim singleton
void
photonicsSim::destroy()
{
  if (s_instance) {
    delete s_instance;
    s_instance = nullptr;
  }
}

//! @brief  photonicsSim ctor
photonicsSim::photonicsSim()
{
}

//! @brief  photonicsSim dtor
photonicsSim::~photonicsSim()
{
  uninit();
}

//! @brief  Uninitialize photonicsSim member classes
void
photonicsSim::uninit()
{
  m_device.reset();
  m_threadPool.reset();
  m_statsMgr.reset();
  m_paramsDram.reset();
  m_config.uninit();
}

//! @brief  Create a PHOTONICS device
bool
photonicsSim::createDevice(PhotonicsDeviceEnum deviceType, unsigned numRanks, unsigned numBankPerRank, unsigned numSubarrayPerBank, unsigned numRows, unsigned numCols, unsigned bufferSize)
{
  photonicsPerfMon perfMon("createDevice");
  uninit();
  bool success = m_config.init(deviceType, numRanks, numBankPerRank, numSubarrayPerBank, numRows, numCols, bufferSize);
  if (!success) {
    return false;
  }
  return createDeviceCommon();
}

//! @brief  Create a PHOTONICS device from a config file
bool
photonicsSim::createDeviceFromConfig(PhotonicsDeviceEnum deviceType, const char* configFilePath)
{
  photonicsPerfMon perfMon("createDeviceFromConfig");
  uninit();
  bool success = m_config.init(deviceType, configFilePath);
  if (!success) {
    return false;
  }
  return createDeviceCommon();
}

//! @brief  Common code to create a PHOTONICS device
bool
photonicsSim::createDeviceCommon()
{
  // Create memory params, which is needed before creating photonicsDevice
  if (!m_config.getMemConfigFile().empty()) {
    m_paramsDram = photonicsParamsDram::createFromConfig(m_config.getMemConfigFile());
  } else {
    m_paramsDram = photonicsParamsDram::create(m_config.getMemoryProtocol());
  }

  // Create PHOTONICS device
  m_device = std::make_unique<photonicsDevice>(m_config);

  if (!m_device->isValid()) {
    uninit();
    std::printf("PHOTONICS-Error: Failed to create PHOTONICS device of type %s\n", photonicsUtils::photonicsDeviceEnumToStr(m_config.getDeviceType()).c_str());
    return false;
  }

  // Create stats mgr
  m_statsMgr = std::make_unique<photonicsStatsMgr>();

  // Create thread pool
  if (getNumThreads() > 1) {
    m_threadPool = std::make_unique<photonicsUtils::threadPool>(getNumThreads());
  }
  return true;
}

//! @brief  Delete PHOTONICS device
bool
photonicsSim::deleteDevice()
{
  uninit();
  return true;
}

//! @brief  Get device properties
bool
photonicsSim::getDeviceProperties(PhotonicsDeviceProperties* deviceProperties) {
  photonicsPerfMon perfMon("getDeviceProperties");
  if (!m_device) {
    std::printf("PHOTONICS-Error: No PHOTONICS device exists.\n");
    return false;
  }
  deviceProperties->deviceType = m_device->getDeviceType();
  deviceProperties->simTarget = m_device->getSimTarget();
  deviceProperties->numRanks = m_device->getNumRanks();
  deviceProperties->numBankPerRank = m_device->getNumBankPerRank();
  deviceProperties->numSubarrayPerBank = m_device->getNumSubarrayPerBank();
  deviceProperties->numRowPerSubarray = m_device->getNumRowPerSubarray();
  deviceProperties->numColPerSubarray = m_device->getNumColPerSubarray();
  deviceProperties->isHLayoutDevice = m_device->isHLayoutDevice();
  deviceProperties->numPHOTONICSCores = m_device->getNumCores();
  return true;
}

//! @brief  Check if device is valid
bool
photonicsSim::isValidDevice(bool showMsg) const
{
  bool isValid = m_device && m_device->isValid();
  if (!isValid && showMsg) {
    std::printf("PHOTONICS-Error: Invalid PHOTONICS device\n");
  }
  return isValid;
}

//! @brief  Get number of PHOTONICS cores
unsigned
photonicsSim::getNumCores() const
{
  if (m_device && m_device->isValid()) {
    return m_device->getNumCores();
  }
  return 0;
}

//! @brief  Get number of rows per PHOTONICS core
unsigned
photonicsSim::getNumRows() const
{
  if (m_device && m_device->isValid()) {
    return m_device->getNumRows();
  }
  return 0;
}

//! @brief  Get number of columns per PHOTONICS core
unsigned
photonicsSim::getNumCols() const
{
  if (m_device && m_device->isValid()) {
    return m_device->getNumCols();
  }
  return 0;
}

//! @brief  Get number of columns per PHOTONICS core
photonicsPerfEnergyBase*
photonicsSim::getPerfEnergyModel()
{
  if (m_device && m_device->isValid()) {
    return m_device->getPerfEnergyModel();
  }
  return nullptr;
}

//! @brief  Start timer for a PHOTONICS kernel to measure CPU runtime and DRAM refresh
void
photonicsSim::startKernelTimer() const
{
  m_statsMgr->startKernelTimer();
}

//! @brief  End timer for a PHOTONICS kernel to measure CPU runtime and DRAM refresh
void
photonicsSim::endKernelTimer() const
{
  m_statsMgr->endKernelTimer();
}

//! @brief  Show PHOTONICS command stats
void
photonicsSim::showStats() const
{
  m_statsMgr->showStats();
}

//! @brief  Reset PHOTONICS command stats
void
photonicsSim::resetStats() const
{
  m_statsMgr->resetStats();
}

//! @brief  Allocate a PHOTONICS object
PhotonicsObjId
photonicsSim::photonicsAlloc(PhotonicsAllocEnum allocType, uint64_t numElements, PhotonicsDataType dataType)
{
  photonicsPerfMon perfMon("photonicsAlloc");
  if (!isValidDevice()) { return -1; }
  return m_device->photonicsAlloc(allocType, numElements, dataType);
}

//! @brief  Allocate a PHOTONICS object
PhotonicsObjId
photonicsSim::photonicsAllocMat(uint64_t numElements, PhotonicsDataType dataType)
{
  photonicsPerfMon perfMon("photonicsAllocMat");
  if (!isValidDevice()) { return -1; }
  return m_device->photonicsAllocMat(numElements, dataType);
}

//! @brief  Allocate a PHOTONICS object that is associated with an existing ojbect
PhotonicsObjId
photonicsSim::photonicsAllocAssociated(PhotonicsObjId assocId, PhotonicsDataType dataType)
{
  photonicsPerfMon perfMon("photonicsAllocAssociated");
  if (!isValidDevice()) { return -1; }
  return m_device->photonicsAllocAssociated(assocId, dataType);
}

//! @brief  Allocate a PHOTONICS object that is associated with an existing ojbect
PhotonicsObjId
photonicsSim::photonicsAllocAssociatedSrcVec(PhotonicsObjId assocId, PhotonicsDataType dataType)
{
  photonicsPerfMon perfMon("photonicsAllocAssociatedSrcVec");
  if (!isValidDevice()) { return -1; }
  return m_device->photonicsAllocAssociatedSrcVec(assocId, dataType);
}

//! @brief  Allocate a PHOTONICS object that is associated with an existing ojbect
PhotonicsObjId
photonicsSim::photonicsAllocAssociatedDestVec(PhotonicsObjId assocId, PhotonicsDataType dataType)
{
  photonicsPerfMon perfMon("photonicsAllocAssociatedDestVec");
  if (!isValidDevice()) { return -1; }
  return m_device->photonicsAllocAssociatedDestVec(assocId, dataType);
}

PhotonicsObjId
photonicsSim::photonicsAllocBuffer(uint32_t numElements, PhotonicsDataType dataType)
{
  photonicsPerfMon perfMon("photonicsAllocBuffer");
  if (!isValidDevice()) { return -1; }
  return m_device->photonicsAllocBuffer(numElements, dataType);
}

// @brief  Free a PHOTONICS object
bool
photonicsSim::photonicsFreeMat(PhotonicsObjId obj)
{
  photonicsPerfMon perfMon("photonicsFree");
  if (!isValidDevice()) { return false; }
  return m_device->photonicsFree(obj, m_device->getNumCores());
}

// @brief  Free a PHOTONICS object
bool
photonicsSim::photonicsFreeSrcVec(PhotonicsObjId obj)
{
  photonicsPerfMon perfMon("photonicsFree");
  if (!isValidDevice()) { return false; }
  return m_device->photonicsFree(obj, m_device->getNumBankPerRank());
}

// @brief  Free a PHOTONICS object
bool
photonicsSim::photonicsFreeDestVec(PhotonicsObjId obj)
{
  photonicsPerfMon perfMon("photonicsFree");
  if (!isValidDevice()) { return false; }
  return m_device->photonicsFree(obj, m_device->getNumCores());
}

//! @brief  Create an obj referencing to a range of an existing obj
PhotonicsObjId
photonicsSim::photonicsCreateRangedRef(PhotonicsObjId refId, uint64_t idxBegin, uint64_t idxEnd)
{
  photonicsPerfMon perfMon("photonicsCreateRangedRef");
  if (!isValidDevice()) { return -1; }
  return m_device->photonicsCreateRangedRef(refId, idxBegin, idxEnd);
}

//! @brief  Create an obj referencing to negation of an existing obj based on dual-contact memory cells
PhotonicsObjId
photonicsSim::photonicsCreateDualContactRef(PhotonicsObjId refId)
{
  photonicsPerfMon perfMon("photonicsCreateDualContactRef");
  if (!isValidDevice()) { return -1; }
  return m_device->photonicsCreateDualContactRef(refId);
}

// @brief  Copy data from main memory to PHOTONICS device within a range
bool
photonicsSim::photonicsCopyMainToDevice(void* src, PhotonicsObjId dest, uint64_t idxBegin, uint64_t idxEnd)
{
  photonicsPerfMon perfMon("photonicsCopyMainToDevice");
  if (!isValidDevice()) { return false; }
  return m_device->photonicsCopyMainToDevice(src, dest, idxBegin, idxEnd);
}

// @brief  Copy data from PHOTONICS device to main memory within a range
bool
photonicsSim::photonicsCopyDeviceToMain(PhotonicsObjId src, void* dest, uint64_t idxBegin, uint64_t idxEnd)
{
  photonicsPerfMon perfMon("photonicsCopyDeviceToMain");
  if (!isValidDevice()) { return false; }
  return m_device->photonicsCopyDeviceToMain(src, dest, idxBegin, idxEnd);
}

// @brief  Copy data from main memory to PHOTONICS device with type within a range
bool
photonicsSim::photonicsCopyMainToDeviceWithType(PhotonicsCopyEnum copyType, void* src, PhotonicsObjId dest, uint64_t idxBegin, uint64_t idxEnd)
{
  photonicsPerfMon perfMon("photonicsCopyMainToDevice");
  if (!isValidDevice()) { return false; }
  return m_device->photonicsCopyMainToDeviceWithType(copyType, src, dest, idxBegin, idxEnd);
}

// @brief  Copy data from PHOTONICS device to main memory with type within a range
bool
photonicsSim::photonicsCopyDeviceToMainWithType(PhotonicsCopyEnum copyType, PhotonicsObjId src, void* dest, uint64_t idxBegin, uint64_t idxEnd)
{
  photonicsPerfMon perfMon("photonicsCopyDeviceToMain");
  if (!isValidDevice()) { return false; }
  return m_device->photonicsCopyDeviceToMainWithType(copyType, src, dest, idxBegin, idxEnd);
}

// @brief  Copy data from PHOTONICS device to device within a range
bool
photonicsSim::photonicsCopyDeviceToDevice(PhotonicsObjId src, PhotonicsObjId dest, uint64_t idxBegin, uint64_t idxEnd)
{
  photonicsPerfMon perfMon("photonicsCopyDeviceToDevice");
  if (!isValidDevice()) { return false; }
  return m_device->photonicsCopyDeviceToDevice(src, dest, idxBegin, idxEnd);
}

bool photonicsSim::photonicsCopyObjectToObject(PhotonicsObjId src, PhotonicsObjId dest)
{
  photonicsPerfMon perfMon("photonicsCopyObjectToObject");
  if (!isValidDevice()) { return false; }
  std::unique_ptr<photonicsCmd> cmd = std::make_unique<photonicsCmdFunc1>(PhotonicsCmdEnum::COPY_O2O, src, dest);
  return m_device->executeCmd(std::move(cmd));
}

bool
photonicsSim::photonicsConvertType(PhotonicsObjId src, PhotonicsObjId dest)
{
  photonicsPerfMon perfMon("photonicsConvertType");
  if (!isValidDevice()) { return false; }
  std::unique_ptr<photonicsCmd> cmd = std::make_unique<photonicsCmdFunc1>(PhotonicsCmdEnum::CONVERT_TYPE, src, dest);
  return m_device->executeCmd(std::move(cmd));
}

// @brief  Load vector with a scalar value
template <typename T> bool
photonicsSim::photonicsBroadcast(PhotonicsObjId dest, T value)
{
  photonicsPerfMon perfMon("photonicsBroadcast");
  if (!isValidDevice()) { return false; }
  uint64_t signExtBits = photonicsUtils::castTypeToBits(value);
  std::unique_ptr<photonicsCmd> cmd = std::make_unique<photonicsCmdBroadcast>(PhotonicsCmdEnum::BROADCAST, dest, signExtBits);
  return m_device->executeCmd(std::move(cmd));
}

// @brief  PHOTONICS OP: add
bool
photonicsSim::photonicsAdd(PhotonicsObjId src1, PhotonicsObjId src2, PhotonicsObjId dest)
{
  photonicsPerfMon perfMon("photonicsAdd");
  if (!isValidDevice()) { return false; }
  std::unique_ptr<photonicsCmd> cmd = std::make_unique<photonicsCmdFunc2>(PhotonicsCmdEnum::ADD, src1, src2, dest);
  return m_device->executeCmd(std::move(cmd));
}

// @brief  PHOTONICS OP: Matrix-Vector Multiplication
bool
photonicsSim::photonicsMvm(PhotonicsObjId src1, PhotonicsObjId src2, PhotonicsObjId dest)
{
  photonicsPerfMon perfMon("photonicsMvm");
  if (!isValidDevice()) { return false; }
  std::unique_ptr<photonicsCmd> cmd = std::make_unique<photonicsCmdMvm>(src1, src2, dest);
  return m_device->executeCmd(std::move(cmd));
}

// @brief  PHOTONICS OP: iterative loop
bool
photonicsSim::photonicsIter(PhotonicsObjId src1, PhotonicsObjId src2, PhotonicsObjId dest, int8_t numLoops)
{
  photonicsPerfMon perfMon("photonicsIter");
  if (!isValidDevice()) { return false; }
  std::unique_ptr<photonicsCmd> cmd = std::make_unique<photonicsCmdIter>(src1, src2, dest, numLoops);
  return m_device->executeCmd(std::move(cmd));
}

// @brief  PHOTONICS OP: Matrix-Matrix Multiplication
bool
photonicsSim::photonicsMmm(PhotonicsObjId src1, PhotonicsObjId src2, PhotonicsObjId dest)
{
  photonicsPerfMon perfMon("photonicsMmm");
  if (!isValidDevice()) { return false; }
  std::unique_ptr<photonicsCmd> cmd = std::make_unique<photonicsCmdMmm>(src1, src2, dest);
  return m_device->executeCmd(std::move(cmd));
}

// @brief  PHOTONICS OP: sub
bool
photonicsSim::photonicsSub(PhotonicsObjId src1, PhotonicsObjId src2, PhotonicsObjId dest)
{
  photonicsPerfMon perfMon("photonicsSub");
  if (!isValidDevice()) { return false; }
  std::unique_ptr<photonicsCmd> cmd = std::make_unique<photonicsCmdFunc2>(PhotonicsCmdEnum::SUB, src1, src2, dest);
  return m_device->executeCmd(std::move(cmd));
}

// @brief PHOTONICS OP: div
bool
photonicsSim::photonicsDiv(PhotonicsObjId src1, PhotonicsObjId src2, PhotonicsObjId dest)
{
  photonicsPerfMon perfMon("photonicsDiv");
  if (!isValidDevice()) { return false; }
  std::unique_ptr<photonicsCmd> cmd = std::make_unique<photonicsCmdFunc2>(PhotonicsCmdEnum::DIV, src1, src2, dest);
  return m_device->executeCmd(std::move(cmd));
}

// @brief  PHOTONICS OP: abs
bool
photonicsSim::photonicsAbs(PhotonicsObjId src, PhotonicsObjId dest)
{
  photonicsPerfMon perfMon("photonicsAbs");
  if (!isValidDevice()) { return false; }
  std::unique_ptr<photonicsCmd> cmd = std::make_unique<photonicsCmdFunc1>(PhotonicsCmdEnum::ABS, src, dest);
  return m_device->executeCmd(std::move(cmd));
}

// @brief  PHOTONICS OP: mul
bool
photonicsSim::photonicsMul(PhotonicsObjId src1, PhotonicsObjId src2, PhotonicsObjId dest)
{
  photonicsPerfMon perfMon("photonicsMul");
  if (!isValidDevice()) { return false; }
  std::unique_ptr<photonicsCmd> cmd = std::make_unique<photonicsCmdFunc2>(PhotonicsCmdEnum::MUL, src1, src2, dest);
  return m_device->executeCmd(std::move(cmd));
}

// @brief  PHOTONICS OP: not
bool
photonicsSim::photonicsNot(PhotonicsObjId src, PhotonicsObjId dest)
{
  photonicsPerfMon perfMon("photonicsNot");
  if (!isValidDevice()) { return false; }
  std::unique_ptr<photonicsCmd> cmd = std::make_unique<photonicsCmdFunc1>(PhotonicsCmdEnum::NOT, src, dest);
  return m_device->executeCmd(std::move(cmd));
}

// @brief  PHOTONICS OP: and
bool
photonicsSim::photonicsAnd(PhotonicsObjId src1, PhotonicsObjId src2, PhotonicsObjId dest)
{
  photonicsPerfMon perfMon("photonicsAnd");
  if (!isValidDevice()) { return false; }
  std::unique_ptr<photonicsCmd> cmd = std::make_unique<photonicsCmdFunc2>(PhotonicsCmdEnum::AND, src1, src2, dest);
  return m_device->executeCmd(std::move(cmd));
}

// @brief  PHOTONICS OP: or
bool
photonicsSim::photonicsOr(PhotonicsObjId src1, PhotonicsObjId src2, PhotonicsObjId dest)
{
  photonicsPerfMon perfMon("photonicsOr");
  if (!isValidDevice()) { return false; }
  std::unique_ptr<photonicsCmd> cmd = std::make_unique<photonicsCmdFunc2>(PhotonicsCmdEnum::OR, src1, src2, dest);
  return m_device->executeCmd(std::move(cmd));
}

// @brief  PHOTONICS OP: xor
bool
photonicsSim::photonicsXor(PhotonicsObjId src1, PhotonicsObjId src2, PhotonicsObjId dest)
{
  photonicsPerfMon perfMon("photonicsXor");
  if (!isValidDevice()) { return false; }
  std::unique_ptr<photonicsCmd> cmd = std::make_unique<photonicsCmdFunc2>(PhotonicsCmdEnum::XOR, src1, src2, dest);
  return m_device->executeCmd(std::move(cmd));
}

// @brief  PHOTONICS OP: xnor
bool
photonicsSim::photonicsXnor(PhotonicsObjId src1, PhotonicsObjId src2, PhotonicsObjId dest)
{
  photonicsPerfMon perfMon("photonicsXnor");
  if (!isValidDevice()) { return false; }
  std::unique_ptr<photonicsCmd> cmd = std::make_unique<photonicsCmdFunc2>(PhotonicsCmdEnum::XNOR, src1, src2, dest);
  return m_device->executeCmd(std::move(cmd));
}

// @brief  PHOTONICS OP: gt
bool
photonicsSim::photonicsGT(PhotonicsObjId src1, PhotonicsObjId src2, PhotonicsObjId dest)
{
  photonicsPerfMon perfMon("photonicsGT");
  if (!isValidDevice()) { return false; }
  std::unique_ptr<photonicsCmd> cmd = std::make_unique<photonicsCmdFunc2>(PhotonicsCmdEnum::GT, src1, src2, dest);
  return m_device->executeCmd(std::move(cmd));
}

// @brief  PHOTONICS OP: lt
bool
photonicsSim::photonicsLT(PhotonicsObjId src1, PhotonicsObjId src2, PhotonicsObjId dest)
{
  photonicsPerfMon perfMon("photonicsLT");
  if (!isValidDevice()) { return false; }
  std::unique_ptr<photonicsCmd> cmd = std::make_unique<photonicsCmdFunc2>(PhotonicsCmdEnum::LT, src1, src2, dest);
  return m_device->executeCmd(std::move(cmd));
}

// @brief  PHOTONICS OP: eq
bool
photonicsSim::photonicsEQ(PhotonicsObjId src1, PhotonicsObjId src2, PhotonicsObjId dest)
{
  photonicsPerfMon perfMon("photonicsEQ");
  if (!isValidDevice()) { return false; }
  std::unique_ptr<photonicsCmd> cmd = std::make_unique<photonicsCmdFunc2>(PhotonicsCmdEnum::EQ, src1, src2, dest);
  return m_device->executeCmd(std::move(cmd));
}

// @brief  PHOTONICS OP: ne
bool
photonicsSim::photonicsNE(PhotonicsObjId src1, PhotonicsObjId src2, PhotonicsObjId dest)
{
  photonicsPerfMon perfMon("photonicsNE");
  if (!isValidDevice()) { return false; }
  std::unique_ptr<photonicsCmd> cmd = std::make_unique<photonicsCmdFunc2>(PhotonicsCmdEnum::NE, src1, src2, dest);
  return m_device->executeCmd(std::move(cmd));
}

// @brief  PHOTONICS OP: min
bool
photonicsSim::photonicsMin(PhotonicsObjId src1, PhotonicsObjId src2, PhotonicsObjId dest)
{
  photonicsPerfMon perfMon("photonicsMin");
  if (!isValidDevice()) { return false; }
  std::unique_ptr<photonicsCmd> cmd = std::make_unique<photonicsCmdFunc2>(PhotonicsCmdEnum::MIN, src1, src2, dest);
  return m_device->executeCmd(std::move(cmd));
}

// @brief  PHOTONICS OP: max
bool
photonicsSim::photonicsMax(PhotonicsObjId src1, PhotonicsObjId src2, PhotonicsObjId dest)
{
  photonicsPerfMon perfMon("photonicsMax");
  if (!isValidDevice()) { return false; }
  std::unique_ptr<photonicsCmd> cmd = std::make_unique<photonicsCmdFunc2>(PhotonicsCmdEnum::MAX, src1, src2, dest);
  return m_device->executeCmd(std::move(cmd));
}

bool photonicsSim::photonicsAdd(PhotonicsObjId src, PhotonicsObjId dest, uint64_t scalarValue)
{
  photonicsPerfMon perfMon("photonicsAddScalar");
  if (!isValidDevice()) { return false; }
  std::unique_ptr<photonicsCmd> cmd = std::make_unique<photonicsCmdFunc1>(PhotonicsCmdEnum::ADD_SCALAR, src, dest, scalarValue);
  return m_device->executeCmd(std::move(cmd));
}

bool photonicsSim::photonicsSub(PhotonicsObjId src, PhotonicsObjId dest, uint64_t scalarValue)
{
  photonicsPerfMon perfMon("photonicsSubScalar");
  if (!isValidDevice()) { return false; }
  std::unique_ptr<photonicsCmd> cmd = std::make_unique<photonicsCmdFunc1>(PhotonicsCmdEnum::SUB_SCALAR, src, dest, scalarValue);
  return m_device->executeCmd(std::move(cmd));
}

bool photonicsSim::photonicsMul(PhotonicsObjId src, PhotonicsObjId dest, uint64_t scalarValue)
{
  photonicsPerfMon perfMon("photonicsMulScalar");
  if (!isValidDevice()) { return false; }
  std::unique_ptr<photonicsCmd> cmd = std::make_unique<photonicsCmdFunc1>(PhotonicsCmdEnum::MUL_SCALAR, src, dest, scalarValue);
  return m_device->executeCmd(std::move(cmd));
}

bool photonicsSim::photonicsDiv(PhotonicsObjId src, PhotonicsObjId dest, uint64_t scalarValue)
{
  photonicsPerfMon perfMon("photonicsDivScalar");
  if (!isValidDevice()) { return false; }
  std::unique_ptr<photonicsCmd> cmd = std::make_unique<photonicsCmdFunc1>(PhotonicsCmdEnum::DIV_SCALAR, src, dest, scalarValue);
  return m_device->executeCmd(std::move(cmd));
}

bool photonicsSim::photonicsAnd(PhotonicsObjId src, PhotonicsObjId dest, uint64_t scalarValue)
{
  photonicsPerfMon perfMon("photonicsAndScalar");
  if (!isValidDevice()) { return false; }
  std::unique_ptr<photonicsCmd> cmd = std::make_unique<photonicsCmdFunc1>(PhotonicsCmdEnum::AND_SCALAR, src, dest, scalarValue);
  return m_device->executeCmd(std::move(cmd));
}

bool photonicsSim::photonicsOr(PhotonicsObjId src, PhotonicsObjId dest, uint64_t scalarValue)
{
  photonicsPerfMon perfMon("photonicsOrScalar");
  if (!isValidDevice()) { return false; }
  std::unique_ptr<photonicsCmd> cmd = std::make_unique<photonicsCmdFunc1>(PhotonicsCmdEnum::OR_SCALAR, src, dest, scalarValue);
  return m_device->executeCmd(std::move(cmd));
}

bool photonicsSim::photonicsXor(PhotonicsObjId src, PhotonicsObjId dest, uint64_t scalarValue)
{
  photonicsPerfMon perfMon("photonicsXorScalar");
  if (!isValidDevice()) { return false; }
  std::unique_ptr<photonicsCmd> cmd = std::make_unique<photonicsCmdFunc1>(PhotonicsCmdEnum::XOR_SCALAR, src, dest, scalarValue);
  return m_device->executeCmd(std::move(cmd));
}

bool photonicsSim::photonicsXnor(PhotonicsObjId src, PhotonicsObjId dest, uint64_t scalarValue)
{
  photonicsPerfMon perfMon("photonicsXnorScalar");
  if (!isValidDevice()) { return false; }
  std::unique_ptr<photonicsCmd> cmd = std::make_unique<photonicsCmdFunc1>(PhotonicsCmdEnum::XNOR_SCALAR, src, dest, scalarValue);
  return m_device->executeCmd(std::move(cmd));
}

bool photonicsSim::photonicsGT(PhotonicsObjId src, PhotonicsObjId dest, uint64_t scalarValue)
{
  photonicsPerfMon perfMon("photonicsGTScalar");
  if (!isValidDevice()) { return false; }
  std::unique_ptr<photonicsCmd> cmd = std::make_unique<photonicsCmdFunc1>(PhotonicsCmdEnum::GT_SCALAR, src, dest, scalarValue);
  return m_device->executeCmd(std::move(cmd));
}

bool photonicsSim::photonicsLT(PhotonicsObjId src, PhotonicsObjId dest, uint64_t scalarValue)
{
  photonicsPerfMon perfMon("photonicsLTScalar");
  if (!isValidDevice()) { return false; }
  std::unique_ptr<photonicsCmd> cmd = std::make_unique<photonicsCmdFunc1>(PhotonicsCmdEnum::LT_SCALAR, src, dest, scalarValue);
  return m_device->executeCmd(std::move(cmd));
}

bool photonicsSim::photonicsEQ(PhotonicsObjId src, PhotonicsObjId dest, uint64_t scalarValue)
{
  photonicsPerfMon perfMon("photonicsEQScalar");
  if (!isValidDevice()) { return false; }
  std::unique_ptr<photonicsCmd> cmd = std::make_unique<photonicsCmdFunc1>(PhotonicsCmdEnum::EQ_SCALAR, src, dest, scalarValue);
  return m_device->executeCmd(std::move(cmd));
}

bool photonicsSim::photonicsNE(PhotonicsObjId src, PhotonicsObjId dest, uint64_t scalarValue)
{
  photonicsPerfMon perfMon("photonicsNEScalar");
  if (!isValidDevice()) { return false; }
  std::unique_ptr<photonicsCmd> cmd = std::make_unique<photonicsCmdFunc1>(PhotonicsCmdEnum::NE_SCALAR, src, dest, scalarValue);
  return m_device->executeCmd(std::move(cmd));
}

bool photonicsSim::photonicsMin(PhotonicsObjId src, PhotonicsObjId dest, uint64_t scalarValue)
{
  photonicsPerfMon perfMon("photonicsMinScalar");
  if (!isValidDevice()) { return false; }
  std::unique_ptr<photonicsCmd> cmd = std::make_unique<photonicsCmdFunc1>(PhotonicsCmdEnum::MIN_SCALAR, src, dest, scalarValue);
  return m_device->executeCmd(std::move(cmd));
}

bool photonicsSim::photonicsMax(PhotonicsObjId src, PhotonicsObjId dest, uint64_t scalarValue)
{
  photonicsPerfMon perfMon("photonicsMaxScalar");
  if (!isValidDevice()) { return false; }
  std::unique_ptr<photonicsCmd> cmd = std::make_unique<photonicsCmdFunc1>(PhotonicsCmdEnum::MAX_SCALAR, src, dest, scalarValue);
  return m_device->executeCmd(std::move(cmd));
}

bool photonicsSim::photonicsScaledAdd(PhotonicsObjId src1, PhotonicsObjId src2, PhotonicsObjId dest, uint64_t scalarValue) {
  photonicsPerfMon perfMon("photonicsScaledAdd");
  if (!isValidDevice()) { return false; }
  std::unique_ptr<photonicsCmd> cmd = std::make_unique<photonicsCmdFunc2>(PhotonicsCmdEnum::SCALED_ADD, src1, src2, dest, scalarValue);
  return m_device->executeCmd(std::move(cmd));
}

// @brief  PHOTONICS OP: popcount
bool
photonicsSim::photonicsPopCount(PhotonicsObjId src, PhotonicsObjId dest)
{
  photonicsPerfMon perfMon("photonicsPopCount");
  if (!isValidDevice()) { return false; }
  std::unique_ptr<photonicsCmd> cmd = std::make_unique<photonicsCmdFunc1>(PhotonicsCmdEnum::POPCOUNT, src, dest);
  return m_device->executeCmd(std::move(cmd));
}

// @brief  PHOTONICS OP: prefixsum
bool
photonicsSim::photonicsPrefixSum(PhotonicsObjId src, PhotonicsObjId dest)
{
  photonicsPerfMon perfMon("photonicsPrefixSum");
  if (!isValidDevice()) { return false; }
  std::unique_ptr<photonicsCmd> cmd = std::make_unique<photonicsCmdPrefixSum>(PhotonicsCmdEnum::PREFIX_SUM, src, dest);
  return m_device->executeCmd(std::move(cmd));
}

 //! @brief  PHOTONICS OP: multiply-accumulate
bool photonicsSim::photonicsMAC(PhotonicsObjId src1, PhotonicsObjId src2, void* dest)
{
  photonicsPerfMon perfMon("photonicsMAC");
  if (!isValidDevice()) { return false; }
  const PhotonicsDataType dataType = m_device->getResMgr()->getObjInfo(src1).getDataType();
  std::unique_ptr<photonicsCmd> cmd;
  PhotonicsCmdEnum cmdType = PhotonicsCmdEnum::MAC;
  
  switch (dataType) {
    case PhotonicsDataType::PHOTONICS_INT8:
      cmd = std::make_unique<photonicsCmdMAC<int8_t>>(cmdType, src1, src2, dest);
      break;
    case PhotonicsDataType::PHOTONICS_INT16:
      cmd = std::make_unique<photonicsCmdMAC<int16_t>>(cmdType, src1, src2, dest);
      break;
    case PhotonicsDataType::PHOTONICS_INT32:
      cmd = std::make_unique<photonicsCmdMAC<int32_t>>(cmdType, src1, src2, dest);
      break;
    case PhotonicsDataType::PHOTONICS_INT64:
      cmd = std::make_unique<photonicsCmdMAC<int64_t>>(cmdType, src1, src2, dest);
      break;
    case PhotonicsDataType::PHOTONICS_UINT8:
      cmd = std::make_unique<photonicsCmdMAC<uint8_t>>(cmdType, src1, src2, dest);
      break;
    case PhotonicsDataType::PHOTONICS_UINT16:
      cmd = std::make_unique<photonicsCmdMAC<uint16_t>>(cmdType, src1, src2, dest);
      break;
    case PhotonicsDataType::PHOTONICS_UINT32:
      cmd = std::make_unique<photonicsCmdMAC<uint32_t>>(cmdType, src1, src2, dest);
      break;
    case PhotonicsDataType::PHOTONICS_UINT64:
      cmd = std::make_unique<photonicsCmdMAC<uint64_t>>(cmdType, src1, src2, dest);
      break;
    case PhotonicsDataType::PHOTONICS_FP8:
    case PhotonicsDataType::PHOTONICS_FP16:
    case PhotonicsDataType::PHOTONICS_BF16:
    case PhotonicsDataType::PHOTONICS_FP32:
      cmd = std::make_unique<photonicsCmdMAC<float>>(cmdType, src1, src2, dest);
      break;
    default:
      std::printf("PHOTONICS-Error: photonicsRedMin does not support data type %s\n", photonicsUtils::photonicsDataTypeEnumToStr(dataType).c_str());
      return false;
  }
  return m_device->executeCmd(std::move(cmd));
}

//! @brief  Min reduction operation
bool photonicsSim::photonicsRedMin(PhotonicsObjId src, void* min, uint64_t idxBegin, uint64_t idxEnd) {
  std::string tag = (idxBegin != idxEnd && idxBegin < idxEnd) ? "photonicsRedMinRanged" : "photonicsRedMin";
  photonicsPerfMon perfMon(tag);
  if (!isValidDevice()) { return false; }
  if (!min) { return false; }

  // Create the reduction command for Min operation
  const PhotonicsDataType dataType = m_device->getResMgr()->getObjInfo(src).getDataType();
  std::unique_ptr<photonicsCmd> cmd;
  PhotonicsCmdEnum cmdType = (idxBegin < idxEnd && idxEnd > 0) ? PhotonicsCmdEnum::REDMIN_RANGE : PhotonicsCmdEnum::REDMIN;
  switch (dataType) {
    case PhotonicsDataType::PHOTONICS_INT8:
      cmd = std::make_unique<photonicsCmdReduction<int8_t>>(cmdType, src, min, idxBegin, idxEnd);
      break;
    case PhotonicsDataType::PHOTONICS_INT16:
      cmd = std::make_unique<photonicsCmdReduction<int16_t>>(cmdType, src, min, idxBegin, idxEnd);
      break;
    case PhotonicsDataType::PHOTONICS_INT32:
      cmd = std::make_unique<photonicsCmdReduction<int32_t>>(cmdType, src, min, idxBegin, idxEnd);
      break;
    case PhotonicsDataType::PHOTONICS_INT64:
      cmd = std::make_unique<photonicsCmdReduction<int64_t>>(cmdType, src, min, idxBegin, idxEnd);
      break;
    case PhotonicsDataType::PHOTONICS_UINT8:
      cmd = std::make_unique<photonicsCmdReduction<uint8_t>>(cmdType, src, min, idxBegin, idxEnd);
      break;
    case PhotonicsDataType::PHOTONICS_UINT16:
      cmd = std::make_unique<photonicsCmdReduction<uint16_t>>(cmdType, src, min, idxBegin, idxEnd);
      break;
    case PhotonicsDataType::PHOTONICS_UINT32:
      cmd = std::make_unique<photonicsCmdReduction<uint32_t>>(cmdType, src, min, idxBegin, idxEnd);
      break;
    case PhotonicsDataType::PHOTONICS_UINT64:
      cmd = std::make_unique<photonicsCmdReduction<uint64_t>>(cmdType, src, min, idxBegin, idxEnd);
      break;
    case PhotonicsDataType::PHOTONICS_FP8:
    case PhotonicsDataType::PHOTONICS_FP16:
    case PhotonicsDataType::PHOTONICS_BF16:
    case PhotonicsDataType::PHOTONICS_FP32:
      cmd = std::make_unique<photonicsCmdReduction<float>>(cmdType, src, min, idxBegin, idxEnd);
      break;
    default:
      std::printf("PHOTONICS-Error: photonicsRedMin does not support data type %s\n", photonicsUtils::photonicsDataTypeEnumToStr(dataType).c_str());
      return false;
  }
  return m_device->executeCmd(std::move(cmd));
}

//! @brief  Max reduction operation
bool photonicsSim::photonicsRedMax(PhotonicsObjId src, void* max, uint64_t idxBegin, uint64_t idxEnd) {
  std::string tag = (idxBegin != idxEnd && idxBegin < idxEnd) ? "photonicsRedMaxRanged" : "photonicsRedMax";
  photonicsPerfMon perfMon(tag);
  if (!isValidDevice()) { return false; }
  if (!max) { return false; }

  // Create the reduction command for Max operation
  const PhotonicsDataType dataType = m_device->getResMgr()->getObjInfo(src).getDataType();
  std::unique_ptr<photonicsCmd> cmd;
  PhotonicsCmdEnum cmdType = (idxBegin < idxEnd && idxEnd > 0) ? PhotonicsCmdEnum::REDMAX_RANGE : PhotonicsCmdEnum::REDMAX;
  switch (dataType) {
    case PhotonicsDataType::PHOTONICS_INT8:
      cmd = std::make_unique<photonicsCmdReduction<int8_t>>(cmdType, src, max, idxBegin, idxEnd);
      break;
    case PhotonicsDataType::PHOTONICS_INT16:
      cmd = std::make_unique<photonicsCmdReduction<int16_t>>(cmdType, src, max, idxBegin, idxEnd);
      break;
    case PhotonicsDataType::PHOTONICS_INT32:
      cmd = std::make_unique<photonicsCmdReduction<int32_t>>(cmdType, src, max, idxBegin, idxEnd);
      break;
    case PhotonicsDataType::PHOTONICS_INT64:
      cmd = std::make_unique<photonicsCmdReduction<int64_t>>(cmdType, src, max, idxBegin, idxEnd);
      break;
    case PhotonicsDataType::PHOTONICS_UINT8:
      cmd = std::make_unique<photonicsCmdReduction<uint8_t>>(cmdType, src, max, idxBegin, idxEnd);
      break;
    case PhotonicsDataType::PHOTONICS_UINT16:
      cmd = std::make_unique<photonicsCmdReduction<uint16_t>>(cmdType, src, max, idxBegin, idxEnd);
      break;
    case PhotonicsDataType::PHOTONICS_UINT32:
      cmd = std::make_unique<photonicsCmdReduction<uint32_t>>(cmdType, src, max, idxBegin, idxEnd);
      break;
    case PhotonicsDataType::PHOTONICS_UINT64:
      cmd = std::make_unique<photonicsCmdReduction<uint64_t>>(cmdType, src, max, idxBegin, idxEnd);
      break;
    case PhotonicsDataType::PHOTONICS_FP8:
    case PhotonicsDataType::PHOTONICS_FP16:
    case PhotonicsDataType::PHOTONICS_BF16:
    case PhotonicsDataType::PHOTONICS_FP32:
      cmd = std::make_unique<photonicsCmdReduction<float>>(cmdType, src, max, idxBegin, idxEnd);
      break;
    default:
      std::printf("PHOTONICS-Error: photonicsRedMax does not support data type %s\n", photonicsUtils::photonicsDataTypeEnumToStr(dataType).c_str());
      return false;
  }
  return m_device->executeCmd(std::move(cmd));
}

bool
photonicsSim::photonicsRedSum(PhotonicsObjId src, void* sum, uint64_t idxBegin, uint64_t idxEnd)
{
  std::string tag = (idxBegin != idxEnd && idxBegin < idxEnd) ? "photonicsRedSumRanged" : "photonicsRedSum";
  photonicsPerfMon perfMon(tag);
  if (!isValidDevice()) { return false; }
  if (!sum) { return false; }

  const PhotonicsDataType dataType = m_device->getResMgr()->getObjInfo(src).getDataType();
  std::unique_ptr<photonicsCmd> cmd;
  PhotonicsCmdEnum cmdType = (idxBegin < idxEnd && idxEnd > 0) ? PhotonicsCmdEnum::REDSUM_RANGE : PhotonicsCmdEnum::REDSUM;
  switch (dataType) {
    case PhotonicsDataType::PHOTONICS_INT8:
    case PhotonicsDataType::PHOTONICS_INT16:
    case PhotonicsDataType::PHOTONICS_INT32:
    case PhotonicsDataType::PHOTONICS_INT64:
      cmd = std::make_unique<photonicsCmdReduction<int64_t>>(cmdType, src, sum, idxBegin, idxEnd);
      break;
    case PhotonicsDataType::PHOTONICS_BOOL:
    case PhotonicsDataType::PHOTONICS_UINT8:
    case PhotonicsDataType::PHOTONICS_UINT16:
    case PhotonicsDataType::PHOTONICS_UINT32:
    case PhotonicsDataType::PHOTONICS_UINT64:
      cmd = std::make_unique<photonicsCmdReduction<uint64_t>>(cmdType, src, sum, idxBegin, idxEnd);
      break;
    case PhotonicsDataType::PHOTONICS_FP8:
    case PhotonicsDataType::PHOTONICS_FP16:
    case PhotonicsDataType::PHOTONICS_BF16:
    case PhotonicsDataType::PHOTONICS_FP32:
      cmd = std::make_unique<photonicsCmdReduction<float>>(cmdType, src, sum, idxBegin, idxEnd);
      break;
    default:
      std::printf("PHOTONICS-Error: photonicsRedSum does not support data type %s\n", photonicsUtils::photonicsDataTypeEnumToStr(dataType).c_str());
      return false;
  }
  return m_device->executeCmd(std::move(cmd));
}

//! @brief  Extract a bit slice from a data vector. Dest must be BOOL type
bool
photonicsSim::photonicsBitSliceExtract(PhotonicsObjId src, PhotonicsObjId destBool, unsigned bitIdx)
{
  photonicsPerfMon perfMon("photonicsBitSliceExtract");
  if (!isValidDevice()) { return false; }
  std::unique_ptr<photonicsCmd> cmd = std::make_unique<photonicsCmdFunc1>(PhotonicsCmdEnum::BIT_SLICE_EXTRACT, src, destBool, bitIdx);
  return m_device->executeCmd(std::move(cmd));
}

//! @brief  Insert a bit slice to a data vector. Src must be BOOL type
bool
photonicsSim::photonicsBitSliceInsert(PhotonicsObjId srcBool, PhotonicsObjId dest, unsigned bitIdx)
{
  photonicsPerfMon perfMon("photonicsBitSliceInsert");
  if (!isValidDevice()) { return false; }
  std::unique_ptr<photonicsCmd> cmd = std::make_unique<photonicsCmdFunc1>(PhotonicsCmdEnum::BIT_SLICE_INSERT, srcBool, dest, bitIdx);
  return m_device->executeCmd(std::move(cmd));
}

bool
photonicsSim::photonicsCondCopy(PhotonicsObjId condBool, PhotonicsObjId src, PhotonicsObjId dest)
{
  photonicsPerfMon perfMon("photonicsCondCopy");
  if (!isValidDevice()) { return false; }
  std::unique_ptr<photonicsCmd> cmd = std::make_unique<photonicsCmdCond>(PhotonicsCmdEnum::COND_COPY, condBool, src, dest);
  return m_device->executeCmd(std::move(cmd));
}

bool
photonicsSim::photonicsCondBroadcast(PhotonicsObjId condBool, uint64_t scalarBits, PhotonicsObjId dest)
{
  photonicsPerfMon perfMon("photonicsCondBroadcast");
  if (!isValidDevice()) { return false; }
  std::unique_ptr<photonicsCmd> cmd = std::make_unique<photonicsCmdCond>(PhotonicsCmdEnum::COND_BROADCAST, condBool, scalarBits, dest);
  return m_device->executeCmd(std::move(cmd));
}

bool
photonicsSim::photonicsCondSelect(PhotonicsObjId condBool, PhotonicsObjId src1, PhotonicsObjId src2, PhotonicsObjId dest)
{
  photonicsPerfMon perfMon("photonicsCondSelect");
  if (!isValidDevice()) { return false; }
  std::unique_ptr<photonicsCmd> cmd = std::make_unique<photonicsCmdCond>(PhotonicsCmdEnum::COND_SELECT, condBool, src1, src2, dest);
  return m_device->executeCmd(std::move(cmd));
}

bool
photonicsSim::photonicsCondSelectScalar(PhotonicsObjId condBool, PhotonicsObjId src1, uint64_t scalarBits, PhotonicsObjId dest)
{
  photonicsPerfMon perfMon("photonicsCondSelectScalar");
  if (!isValidDevice()) { return false; }
  std::unique_ptr<photonicsCmd> cmd = std::make_unique<photonicsCmdCond>(PhotonicsCmdEnum::COND_SELECT_SCALAR, condBool, src1, scalarBits, dest);
  return m_device->executeCmd(std::move(cmd));
}

bool
photonicsSim::photonicsRotateElementsRight(PhotonicsObjId src)
{
  photonicsPerfMon perfMon("photonicsRotateElementsRight");
  if (!isValidDevice()) { return false; }
  std::unique_ptr<photonicsCmd> cmd = std::make_unique<photonicsCmdRotate>(PhotonicsCmdEnum::ROTATE_ELEM_R, src);
  return m_device->executeCmd(std::move(cmd));
}

bool
photonicsSim::photonicsRotateElementsLeft(PhotonicsObjId src)
{
  photonicsPerfMon perfMon("photonicsRotateElementsLeft");
  if (!isValidDevice()) { return false; }
  std::unique_ptr<photonicsCmd> cmd = std::make_unique<photonicsCmdRotate>(PhotonicsCmdEnum::ROTATE_ELEM_L, src);
  return m_device->executeCmd(std::move(cmd));
}

bool
photonicsSim::photonicsShiftElementsRight(PhotonicsObjId src)
{
  photonicsPerfMon perfMon("photonicsShiftElementsRight");
  if (!isValidDevice()) { return false; }
  std::unique_ptr<photonicsCmd> cmd = std::make_unique<photonicsCmdRotate>(PhotonicsCmdEnum::SHIFT_ELEM_R, src);
  return m_device->executeCmd(std::move(cmd));
}

bool
photonicsSim::photonicsShiftElementsLeft(PhotonicsObjId src)
{
  photonicsPerfMon perfMon("photonicsShiftElementsLeft");
  if (!isValidDevice()) { return false; }
  std::unique_ptr<photonicsCmd> cmd = std::make_unique<photonicsCmdRotate>(PhotonicsCmdEnum::SHIFT_ELEM_L, src);
  return m_device->executeCmd(std::move(cmd));
}

bool
photonicsSim::photonicsShiftBitsRight(PhotonicsObjId src, PhotonicsObjId dest, unsigned shiftAmount)
{
  photonicsPerfMon perfMon("photonicsShiftBitsRight");
  if (!isValidDevice()) { return false; }
  std::unique_ptr<photonicsCmd> cmd = std::make_unique<photonicsCmdFunc1>(PhotonicsCmdEnum::SHIFT_BITS_R, src, dest, shiftAmount);
  return m_device->executeCmd(std::move(cmd));
}

bool
photonicsSim::photonicsShiftBitsLeft(PhotonicsObjId src, PhotonicsObjId dest, unsigned shiftAmount)
{
  photonicsPerfMon perfMon("photonicsShiftBitsLeft");
  if (!isValidDevice()) { return false; }
  std::unique_ptr<photonicsCmd> cmd = std::make_unique<photonicsCmdFunc1>(PhotonicsCmdEnum::SHIFT_BITS_L, src, dest, shiftAmount);
  return m_device->executeCmd(std::move(cmd));
}

bool 
photonicsSim::photonicsAesSbox(PhotonicsObjId src, PhotonicsObjId dest, const std::vector<uint8_t>& lut)
{
  photonicsPerfMon perfMon("photonicsAesSbox");
  if (!isValidDevice()) { return false; }
  std::unique_ptr<photonicsCmd> cmd = std::make_unique<photonicsCmdFunc1>(PhotonicsCmdEnum::AES_SBOX, src, dest, lut);
  return m_device->executeCmd(std::move(cmd));
}

bool 
photonicsSim::photonicsAesInverseSbox(PhotonicsObjId src, PhotonicsObjId dest, const std::vector<uint8_t>& lut)
{
  photonicsPerfMon perfMon("photonicsAesInverseSbox");
  if (!isValidDevice()) { return false; }
  std::unique_ptr<photonicsCmd> cmd = std::make_unique<photonicsCmdFunc1>(PhotonicsCmdEnum::AES_INVERSE_SBOX, src, dest, lut);
  return m_device->executeCmd(std::move(cmd));
}

bool
photonicsSim::photonicsFuse(PhotonicsProg prog)
{
  photonicsPerfMon perfMon("photonicsFuse");
  if (!isValidDevice()) { return false; }
  std::unique_ptr<photonicsCmd> cmd = std::make_unique<photonicsCmdFuse>(prog);
  return m_device->executeCmd(std::move(cmd));
}

bool
photonicsSim::photonicsOpReadRowToSa(PhotonicsObjId objId, unsigned ofst)
{
  photonicsPerfMon perfMon("photonicsOpReadRowToSa");
  if (!isValidDevice()) { return false; }
  std::unique_ptr<photonicsCmd> cmd = std::make_unique<photonicsCmdReadRowToSa>(PhotonicsCmdEnum::ROW_R, objId, ofst);
  return m_device->executeCmd(std::move(cmd));
}

bool
photonicsSim::photonicsOpWriteSaToRow(PhotonicsObjId objId, unsigned ofst)
{
  photonicsPerfMon perfMon("photonicsOpWriteSaToRow");
  if (!isValidDevice()) { return false; }
  std::unique_ptr<photonicsCmd> cmd = std::make_unique<photonicsCmdWriteSaToRow>(PhotonicsCmdEnum::ROW_W, objId, ofst);
  return m_device->executeCmd(std::move(cmd));
}

bool
photonicsSim::photonicsOpTRA(PhotonicsObjId src1, unsigned ofst1, PhotonicsObjId src2, unsigned ofst2, PhotonicsObjId src3, unsigned ofst3)
{
  photonicsPerfMon perfMon("photonicsOpTRA");
  return false;
}

bool
photonicsSim::photonicsOpMove(PhotonicsObjId objId, PhotonicsRowReg src, PhotonicsRowReg dest)
{
  photonicsPerfMon perfMon("photonicsOpMove");
  if (!isValidDevice()) { return false; }
  std::unique_ptr<photonicsCmd> cmd = std::make_unique<photonicsCmdRRegOp>(PhotonicsCmdEnum::RREG_MOV, objId, dest, src);
  return m_device->executeCmd(std::move(cmd));
}

bool
photonicsSim::photonicsOpSet(PhotonicsObjId objId, PhotonicsRowReg dest, bool val)
{
  photonicsPerfMon perfMon("photonicsOpSet");
  if (!isValidDevice()) { return false; }
  std::unique_ptr<photonicsCmd> cmd = std::make_unique<photonicsCmdRRegOp>(PhotonicsCmdEnum::RREG_SET, objId, dest, val);
  return m_device->executeCmd(std::move(cmd));
}

bool
photonicsSim::photonicsOpNot(PhotonicsObjId objId, PhotonicsRowReg src, PhotonicsRowReg dest)
{
  photonicsPerfMon perfMon("photonicsOpNot");
  if (!isValidDevice()) { return false; }
  std::unique_ptr<photonicsCmd> cmd = std::make_unique<photonicsCmdRRegOp>(PhotonicsCmdEnum::RREG_NOT, objId, dest, src);
  return m_device->executeCmd(std::move(cmd));
}

bool
photonicsSim::photonicsOpAnd(PhotonicsObjId objId, PhotonicsRowReg src1, PhotonicsRowReg src2, PhotonicsRowReg dest)
{
  photonicsPerfMon perfMon("photonicsOpAnd");
  if (!isValidDevice()) { return false; }
  std::unique_ptr<photonicsCmd> cmd = std::make_unique<photonicsCmdRRegOp>(PhotonicsCmdEnum::RREG_AND, objId, dest, src1, src2);
  return m_device->executeCmd(std::move(cmd));
}

bool
photonicsSim::photonicsOpOr(PhotonicsObjId objId, PhotonicsRowReg src1, PhotonicsRowReg src2, PhotonicsRowReg dest)
{
  photonicsPerfMon perfMon("photonicsOpOr");
  if (!isValidDevice()) { return false; }
  std::unique_ptr<photonicsCmd> cmd = std::make_unique<photonicsCmdRRegOp>(PhotonicsCmdEnum::RREG_OR, objId, dest, src1, src2);
  return m_device->executeCmd(std::move(cmd));
}

bool
photonicsSim::photonicsOpNand(PhotonicsObjId objId, PhotonicsRowReg src1, PhotonicsRowReg src2, PhotonicsRowReg dest)
{
  photonicsPerfMon perfMon("photonicsOpNand");
  if (!isValidDevice()) { return false; }
  std::unique_ptr<photonicsCmd> cmd = std::make_unique<photonicsCmdRRegOp>(PhotonicsCmdEnum::RREG_NAND, objId, dest, src1, src2);
  return m_device->executeCmd(std::move(cmd));
}

bool
photonicsSim::photonicsOpNor(PhotonicsObjId objId, PhotonicsRowReg src1, PhotonicsRowReg src2, PhotonicsRowReg dest)
{
  photonicsPerfMon perfMon("photonicsOpNor");
  if (!isValidDevice()) { return false; }
  std::unique_ptr<photonicsCmd> cmd = std::make_unique<photonicsCmdRRegOp>(PhotonicsCmdEnum::RREG_NOR, objId, dest, src1, src2);
  return m_device->executeCmd(std::move(cmd));
}

bool
photonicsSim::photonicsOpXor(PhotonicsObjId objId, PhotonicsRowReg src1, PhotonicsRowReg src2, PhotonicsRowReg dest)
{
  photonicsPerfMon perfMon("photonicsOpXor");
  if (!isValidDevice()) { return false; }
  std::unique_ptr<photonicsCmd> cmd = std::make_unique<photonicsCmdRRegOp>(PhotonicsCmdEnum::RREG_XOR, objId, dest, src1, src2);
  return m_device->executeCmd(std::move(cmd));
}

bool
photonicsSim::photonicsOpXnor(PhotonicsObjId objId, PhotonicsRowReg src1, PhotonicsRowReg src2, PhotonicsRowReg dest)
{
  photonicsPerfMon perfMon("photonicsOpXnor");
  if (!isValidDevice()) { return false; }
  std::unique_ptr<photonicsCmd> cmd = std::make_unique<photonicsCmdRRegOp>(PhotonicsCmdEnum::RREG_XNOR, objId, dest, src1, src2);
  return m_device->executeCmd(std::move(cmd));
}

bool
photonicsSim::photonicsOpMaj(PhotonicsObjId objId, PhotonicsRowReg src1, PhotonicsRowReg src2, PhotonicsRowReg src3, PhotonicsRowReg dest)
{
  photonicsPerfMon perfMon("photonicsOpMaj");
  if (!isValidDevice()) { return false; }
  std::unique_ptr<photonicsCmd> cmd = std::make_unique<photonicsCmdRRegOp>(PhotonicsCmdEnum::RREG_MAJ, objId, dest, src1, src2, src3);
  return m_device->executeCmd(std::move(cmd));
}

bool
photonicsSim::photonicsOpSel(PhotonicsObjId objId, PhotonicsRowReg cond, PhotonicsRowReg src1, PhotonicsRowReg src2, PhotonicsRowReg dest)
{
  photonicsPerfMon perfMon("photonicsOpSel");
  if (!isValidDevice()) { return false; }
  std::unique_ptr<photonicsCmd> cmd = std::make_unique<photonicsCmdRRegOp>(PhotonicsCmdEnum::RREG_SEL, objId, dest, cond, src1, src2);
  return m_device->executeCmd(std::move(cmd));
}

bool
photonicsSim::photonicsOpRotateRH(PhotonicsObjId objId, PhotonicsRowReg src)
{
  photonicsPerfMon perfMon("photonicsOpRotateRH");
  if (!isValidDevice()) { return false; }
  std::unique_ptr<photonicsCmd> cmd = std::make_unique<photonicsCmdRRegRotate>(PhotonicsCmdEnum::RREG_ROTATE_R, objId, src);
  return m_device->executeCmd(std::move(cmd));
}

bool
photonicsSim::photonicsOpRotateLH(PhotonicsObjId objId, PhotonicsRowReg src)
{
  photonicsPerfMon perfMon("photonicsOpRotateLH");
  if (!isValidDevice()) { return false; }
  std::unique_ptr<photonicsCmd> cmd = std::make_unique<photonicsCmdRRegRotate>(PhotonicsCmdEnum::RREG_ROTATE_L, objId, src);
  return m_device->executeCmd(std::move(cmd));
}

bool
photonicsSim::photonicsOpAP(int numSrc, va_list args)
{
  photonicsPerfMon perfMon("photonicsOpAP");
  if (!isValidDevice()) { return false; }

  std::vector<std::pair<PhotonicsObjId, unsigned>> srcRows;
  for (int i = 0; i < numSrc; ++i) {
    PhotonicsObjId objId = va_arg(args, PhotonicsObjId);
    unsigned ofst = va_arg(args, unsigned);
    srcRows.emplace_back(objId, ofst);
  }

  std::unique_ptr<photonicsCmd> cmd = std::make_unique<photonicsCmdAnalogAAP>(PhotonicsCmdEnum::ROW_AP, srcRows);
  return m_device->executeCmd(std::move(cmd));
}

bool
photonicsSim::photonicsOpAAP(int numSrc, int numDest, va_list args)
{
  photonicsPerfMon perfMon("photonicsOpAAP");
  if (!isValidDevice()) { return false; }

  std::vector<std::pair<PhotonicsObjId, unsigned>> srcRows;
  for (int i = 0; i < numSrc; ++i) {
    PhotonicsObjId objId = va_arg(args, PhotonicsObjId);
    int ofst = va_arg(args, unsigned);
    srcRows.emplace_back(objId, ofst);
  }
  std::vector<std::pair<PhotonicsObjId, unsigned>> destRows;
  for (int i = 0; i < numDest; ++i) {
    PhotonicsObjId objId = va_arg(args, PhotonicsObjId);
    int ofst = va_arg(args, unsigned);
    destRows.emplace_back(objId, ofst);
  }

  std::unique_ptr<photonicsCmd> cmd = std::make_unique<photonicsCmdAnalogAAP>(PhotonicsCmdEnum::ROW_AAP, srcRows, destRows);
  return m_device->executeCmd(std::move(cmd));
}

// Explicit template instantiations
template bool photonicsSim::photonicsBroadcast<uint64_t>(PhotonicsObjId dest, uint64_t value);
template bool photonicsSim::photonicsBroadcast<int64_t>(PhotonicsObjId dest, int64_t value);
template bool photonicsSim::photonicsBroadcast<float>(PhotonicsObjId dest, float value);

