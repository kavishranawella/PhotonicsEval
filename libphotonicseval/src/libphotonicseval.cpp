// File: libphotonicseval.cpp
// PHOTONICSeval Simulator - Library Interface
// Copyright (c) 2024 University of Virginia
// This file is licensed under the MIT License.
// See the LICENSE file in the root of this repository for more details.

#include "libphotonicseval.h"
#include "photonicsSim.h"
#include "photonicsUtils.h"

//! @brief  Create a PHOTONICS device
PhotonicsStatus
photonicsCreateDevice(PhotonicsDeviceEnum deviceType, unsigned numRanks, unsigned numBankPerRank, unsigned numSubarrayPerBank, unsigned numRows, unsigned numCols, unsigned matrixSize, unsigned bufferSize)
{
  bool ok = photonicsSim::get()->createDevice(deviceType, numRanks, numBankPerRank, numSubarrayPerBank, numRows, numCols, bufferSize, matrixSize);
  return ok ? PHOTONICS_OK : PHOTONICS_ERROR;
}

//! @brief  Create a PHOTONICS device from config file
PhotonicsStatus
photonicsCreateDeviceFromConfig(PhotonicsDeviceEnum deviceType, const char* configFileName)
{
  bool ok = photonicsSim::get()->createDeviceFromConfig(deviceType, configFileName);
  return ok ? PHOTONICS_OK : PHOTONICS_ERROR;
}

//! @brief  Get PHOTONICS device properties
PhotonicsStatus
photonicsGetDeviceProperties(PhotonicsDeviceProperties* deviceProperties)
{
  bool ok = photonicsSim::get()->getDeviceProperties(deviceProperties);
  return ok ? PHOTONICS_OK : PHOTONICS_ERROR;
}

//! @brief  Delete a PHOTONICS device
PhotonicsStatus
photonicsDeleteDevice()
{
  bool ok = photonicsSim::get()->deleteDevice();
  return ok ? PHOTONICS_OK : PHOTONICS_ERROR;
}

PhotonicsStatus
photonicsPrefixSum(PhotonicsObjId src, PhotonicsObjId dest)
{
  bool ok = photonicsSim::get()->photonicsPrefixSum(src, dest);
  return ok ? PHOTONICS_OK : PHOTONICS_ERROR;
}

//! @brief  Start timer for a PHOTONICS kernel to measure CPU runtime and DRAM refresh
void
photonicsStartTimer()
{
  photonicsSim::get()->startKernelTimer();
}

//! @brief  End timer for a PHOTONICS kernel to measure CPU runtime and DRAM refresh
void
photonicsEndTimer()
{
  photonicsSim::get()->endKernelTimer();
}

//! @brief  Show PHOTONICS command stats
void
photonicsShowStats()
{
  photonicsSim::get()->showStats();
}

//! @brief  Reset PHOTONICS command stats
void
photonicsResetStats()
{
  photonicsSim::get()->resetStats();
}

//! @brief  Is analysis mode. Call this after device creation
bool
photonicsIsAnalysisMode()
{
  return photonicsSim::get()->isAnalysisMode();
}

//! @brief  Allocate a PHOTONICS resource
PhotonicsObjId
photonicsAlloc(PhotonicsAllocEnum allocType, uint64_t numElements, PhotonicsDataType dataType)
{
  return photonicsSim::get()->photonicsAlloc(allocType, numElements, dataType);
}

//! @brief  Allocate a PHOTONICS resource
PhotonicsObjId
photonicsAllocMat(uint64_t numElements, PhotonicsDataType dataType)
{
  return photonicsSim::get()->photonicsAllocMat(numElements, dataType);
}

//! @brief  Allocate a PHOTONICS resource, with an associated object as reference
PhotonicsObjId
photonicsAllocAssociated(PhotonicsObjId assocId, PhotonicsDataType dataType)
{
  return photonicsSim::get()->photonicsAllocAssociated(assocId, dataType);
}

//! @brief  Allocate a PHOTONICS resource, with an associated object as reference
PhotonicsObjId
photonicsAllocAssociatedSrcVec(PhotonicsObjId assocId, PhotonicsDataType dataType)
{
  return photonicsSim::get()->photonicsAllocAssociatedSrcVec(assocId, dataType);
}

//! @brief  Allocate a PHOTONICS resource, with an associated object as reference
PhotonicsObjId
photonicsAllocAssociatedDestVec(PhotonicsObjId assocId, PhotonicsDataType dataType)
{
  return photonicsSim::get()->photonicsAllocAssociatedDestVec(assocId, dataType);
}

//! @brief  Allocate a global buffer for broadcasting data to all PHOTONICS cores
PhotonicsObjId
photonicsAllocBuffer(uint32_t numElements, PhotonicsDataType dataType)
{
  return photonicsSim::get()->photonicsAllocBuffer(numElements, dataType);
}

//! @brief  Free a PHOTONICS resource
PhotonicsStatus
photonicsFreeMat(PhotonicsObjId obj)
{
  bool ok = photonicsSim::get()->photonicsFreeMat(obj);
  return ok ? PHOTONICS_OK : PHOTONICS_ERROR;
}

//! @brief  Free a PHOTONICS resource
PhotonicsStatus
photonicsFreeSrcVec(PhotonicsObjId obj)
{
  bool ok = photonicsSim::get()->photonicsFreeSrcVec(obj);
  return ok ? PHOTONICS_OK : PHOTONICS_ERROR;
}

//! @brief  Free a PHOTONICS resource
PhotonicsStatus
photonicsFreeDestVec(PhotonicsObjId obj)
{
  bool ok = photonicsSim::get()->photonicsFreeDestVec(obj);
  return ok ? PHOTONICS_OK : PHOTONICS_ERROR;
}

//! @brief  Create an obj referencing to a range of an existing obj
PhotonicsObjId
photonicsCreateRangedRef(PhotonicsObjId refId, uint64_t idxBegin, uint64_t idxEnd)
{
  return photonicsSim::get()->photonicsCreateRangedRef(refId, idxBegin, idxEnd);
}

//! @brief  Create an obj referencing to negation of an existing obj based on dual-contact memory cells
PhotonicsObjId
photonicsCreateDualContactRef(PhotonicsObjId refId)
{
  return photonicsSim::get()->photonicsCreateDualContactRef(refId);
}

//! @brief  Copy data from main memory to PHOTONICS device for a range of elements within the PHOTONICS object
PhotonicsStatus
photonicsCopyHostToDevice(void* src, PhotonicsObjId dest, uint64_t idxBegin, uint64_t idxEnd)
{
  bool ok = photonicsSim::get()->photonicsCopyMainToDevice(src, dest, idxBegin, idxEnd);
  return ok ? PHOTONICS_OK : PHOTONICS_ERROR;
}

//! @brief  Copy data from main memory to PHOTONICS device for a range of elements within the PHOTONICS object
PhotonicsStatus
photonicsCopyHostToDeviceMat(void* src, PhotonicsObjId dest, uint64_t idxBegin, uint64_t idxEnd)
{
  bool ok = photonicsSim::get()->photonicsCopyMainToDevice(src, dest, idxBegin, idxEnd);
  return ok ? PHOTONICS_OK : PHOTONICS_ERROR;
}

//! @brief  Copy data from main memory to PHOTONICS device for a range of elements within the PHOTONICS object
PhotonicsStatus
photonicsCopyHostToDeviceVec(void* src, PhotonicsObjId dest, uint64_t idxBegin, uint64_t idxEnd)
{
  bool ok = photonicsSim::get()->photonicsCopyMainToDevice(src, dest, idxBegin, idxEnd);
  return ok ? PHOTONICS_OK : PHOTONICS_ERROR;
}

//! @brief  Copy data from PHOTONICS device to main memory for a range of elements within the PHOTONICS object
PhotonicsStatus
photonicsCopyDeviceToHost(PhotonicsObjId src, void* dest, uint64_t idxBegin, uint64_t idxEnd)
{
  bool ok = photonicsSim::get()->photonicsCopyDeviceToMain(src, dest, idxBegin, idxEnd);
  return ok ? PHOTONICS_OK : PHOTONICS_ERROR;
}

//! @brief  Copy data from PHOTONICS device to main memory for a range of elements within the PHOTONICS object
PhotonicsStatus
photonicsCopyDeviceToHostVec(PhotonicsObjId src, void* dest, uint64_t idxBegin, uint64_t idxEnd)
{
  bool ok = photonicsSim::get()->photonicsCopyDeviceToMain(src, dest, idxBegin, idxEnd);
  return ok ? PHOTONICS_OK : PHOTONICS_ERROR;
}

//! @brief  Copy data from main memory to PHOTONICS device with type for a range of elements within the PHOTONICS object
PhotonicsStatus
photonicsCopyHostToDeviceWithType(PhotonicsCopyEnum copyType, void* src, PhotonicsObjId dest, uint64_t idxBegin, uint64_t idxEnd)
{
  bool ok = photonicsSim::get()->photonicsCopyMainToDeviceWithType(copyType, src, dest, idxBegin, idxEnd);
  return ok ? PHOTONICS_OK : PHOTONICS_ERROR;
}

//! @brief  Copy data from PHOTONICS device to main memory with type for a range of elements within the PHOTONICS object
PhotonicsStatus
photonicsCopyDeviceToHostWithType(PhotonicsCopyEnum copyType, PhotonicsObjId src, void* dest, uint64_t idxBegin, uint64_t idxEnd)
{
  bool ok = photonicsSim::get()->photonicsCopyDeviceToMainWithType(copyType, src, dest, idxBegin, idxEnd);
  return ok ? PHOTONICS_OK : PHOTONICS_ERROR;
}

//! @brief  Copy data from PHOTONICS device to device for a range of elements within the PHOTONICS object
PhotonicsStatus
photonicsCopyDeviceToDevice(PhotonicsObjId src, PhotonicsObjId dest, uint64_t idxBegin, uint64_t idxEnd)
{
  bool ok = photonicsSim::get()->photonicsCopyDeviceToDevice(src, dest, idxBegin, idxEnd);
  return ok ? PHOTONICS_OK : PHOTONICS_ERROR;
}

PhotonicsStatus photonicsCopyObjectToObject(PhotonicsObjId src, PhotonicsObjId dest)
{
  bool ok = photonicsSim::get()->photonicsCopyObjectToObject(src, dest);
  return ok ? PHOTONICS_OK : PHOTONICS_ERROR;
}

//! @brief  Convert data type between two associated PHOTONICS objects of different data types
PhotonicsStatus photonicsConvertType(PhotonicsObjId src, PhotonicsObjId dest)
{
  bool ok = photonicsSim::get()->photonicsConvertType(src, dest);
  return ok ? PHOTONICS_OK : PHOTONICS_ERROR;
}

//! @brief  Load vector with a signed int value
PhotonicsStatus
photonicsBroadcastInt(PhotonicsObjId dest, int64_t value)
{
  bool ok = photonicsSim::get()->photonicsBroadcast(dest, value);
  return ok ? PHOTONICS_OK : PHOTONICS_ERROR;
}

//! @brief  Load vector with an unsigned int value
PhotonicsStatus
photonicsBroadcastUInt(PhotonicsObjId dest, uint64_t value)
{
  bool ok = photonicsSim::get()->photonicsBroadcast(dest, value);
  return ok ? PHOTONICS_OK : PHOTONICS_ERROR;
}

//! @brief  Load vector with a float32 value
PhotonicsStatus
photonicsBroadcastFP(PhotonicsObjId dest, float value)
{
  bool ok = photonicsSim::get()->photonicsBroadcast(dest, value);
  return ok ? PHOTONICS_OK : PHOTONICS_ERROR;
}

//! @brief  PHOTONICS add
PhotonicsStatus
photonicsAdd(PhotonicsObjId src1, PhotonicsObjId src2, PhotonicsObjId dest)
{
  bool ok = photonicsSim::get()->photonicsAdd(src1, src2, dest);
  return ok ? PHOTONICS_OK : PHOTONICS_ERROR;
}

//! @brief  PHOTONICS Marix-Vector Multiplication
PhotonicsStatus
photonicsMvm(PhotonicsObjId src1, PhotonicsObjId src2, PhotonicsObjId dest)
{
  bool ok = photonicsSim::get()->photonicsMvm(src1, src2, dest);
  return ok ? PHOTONICS_OK : PHOTONICS_ERROR;
}

//! @brief  PHOTONICS iterative loop
PhotonicsStatus
photonicsIter(PhotonicsObjId src1, PhotonicsObjId src2, PhotonicsObjId dest, int8_t numLoops)
{
  bool ok = photonicsSim::get()->photonicsIter(src1, src2, dest, numLoops);
  return ok ? PHOTONICS_OK : PHOTONICS_ERROR;
}

//! @brief  PHOTONICS Marix-Matrix Multiplication
PhotonicsStatus
photonicsMmm(PhotonicsObjId src1, PhotonicsObjId src2, PhotonicsObjId dest)
{
  bool ok = photonicsSim::get()->photonicsMmm(src1, src2, dest);
  return ok ? PHOTONICS_OK : PHOTONICS_ERROR;
}

//! @brief  PHOTONICS sub
PhotonicsStatus
photonicsSub(PhotonicsObjId src1, PhotonicsObjId src2, PhotonicsObjId dest)
{
  bool ok = photonicsSim::get()->photonicsSub(src1, src2, dest);
  return ok ? PHOTONICS_OK : PHOTONICS_ERROR;
}

//! @brief  PHOTONICS div
PhotonicsStatus
photonicsDiv(PhotonicsObjId src1, PhotonicsObjId src2, PhotonicsObjId dest)
{
  bool ok = photonicsSim::get()->photonicsDiv(src1, src2, dest);
  return ok ? PHOTONICS_OK : PHOTONICS_ERROR;
}

//! @brief  PHOTONICS not
PhotonicsStatus
photonicsNot(PhotonicsObjId src, PhotonicsObjId dest)
{
  bool ok = photonicsSim::get()->photonicsNot(src, dest);;
  return ok ? PHOTONICS_OK : PHOTONICS_ERROR;
}

//! @brief  PHOTONICS or
PhotonicsStatus
photonicsOr(PhotonicsObjId src1, PhotonicsObjId src2, PhotonicsObjId dest)
{
  bool ok = photonicsSim::get()->photonicsOr(src1, src2, dest);
  return ok ? PHOTONICS_OK : PHOTONICS_ERROR;
}

//! @brief  PHOTONICS and
PhotonicsStatus
photonicsAnd(PhotonicsObjId src1, PhotonicsObjId src2, PhotonicsObjId dest)
{
  bool ok = photonicsSim::get()->photonicsAnd(src1, src2, dest);
  return ok ? PHOTONICS_OK : PHOTONICS_ERROR;
}

//! @brief  PHOTONICS xor
PhotonicsStatus
photonicsXor(PhotonicsObjId src1, PhotonicsObjId src2, PhotonicsObjId dest)
{
  bool ok = photonicsSim::get()->photonicsXor(src1, src2, dest);
  return ok ? PHOTONICS_OK : PHOTONICS_ERROR;
}

//! @brief  PHOTONICS xnor
PhotonicsStatus
photonicsXnor(PhotonicsObjId src1, PhotonicsObjId src2, PhotonicsObjId dest)
{
  bool ok = photonicsSim::get()->photonicsXnor(src1, src2, dest);
  return ok ? PHOTONICS_OK : PHOTONICS_ERROR;
}

//! @brief  PHOTONICS abs
PhotonicsStatus
photonicsAbs(PhotonicsObjId src, PhotonicsObjId dest)
{
  bool ok = photonicsSim::get()->photonicsAbs(src, dest);;
  return ok ? PHOTONICS_OK : PHOTONICS_ERROR;
}

//! @brief  PHOTONICS multiplication
PhotonicsStatus
photonicsMul(PhotonicsObjId src1, PhotonicsObjId src2, PhotonicsObjId dest)
{
  bool ok = photonicsSim::get()->photonicsMul(src1, src2, dest);
  return ok ? PHOTONICS_OK : PHOTONICS_ERROR;
}

//! @brief  PHOTONICS GT
PhotonicsStatus
photonicsGT(PhotonicsObjId src1, PhotonicsObjId src2, PhotonicsObjId dest)
{
  bool ok = photonicsSim::get()->photonicsGT(src1, src2, dest);
  return ok ? PHOTONICS_OK : PHOTONICS_ERROR;
}

//! @brief  PHOTONICS LT
PhotonicsStatus
photonicsLT(PhotonicsObjId src1, PhotonicsObjId src2, PhotonicsObjId dest)
{
  bool ok = photonicsSim::get()->photonicsLT(src1, src2, dest);
  return ok ? PHOTONICS_OK : PHOTONICS_ERROR;
}

//! @brief  PHOTONICS EQ
PhotonicsStatus
photonicsEQ(PhotonicsObjId src1, PhotonicsObjId src2, PhotonicsObjId dest)
{
  bool ok = photonicsSim::get()->photonicsEQ(src1, src2, dest);
  return ok ? PHOTONICS_OK : PHOTONICS_ERROR;
}

//! @brief  PHOTONICS NE
PhotonicsStatus
photonicsNE(PhotonicsObjId src1, PhotonicsObjId src2, PhotonicsObjId dest)
{
  bool ok = photonicsSim::get()->photonicsNE(src1, src2, dest);
  return ok ? PHOTONICS_OK : PHOTONICS_ERROR;
}

//! @brief  PHOTONICS Min
PhotonicsStatus
photonicsMin(PhotonicsObjId src1, PhotonicsObjId src2, PhotonicsObjId dest)
{
  bool ok = photonicsSim::get()->photonicsMin(src1, src2, dest);
  return ok ? PHOTONICS_OK : PHOTONICS_ERROR;
}

//! @brief  PHOTONICS Max
PhotonicsStatus
photonicsMax(PhotonicsObjId src1, PhotonicsObjId src2, PhotonicsObjId dest)
{
  bool ok = photonicsSim::get()->photonicsMax(src1, src2, dest);
  return ok ? PHOTONICS_OK : PHOTONICS_ERROR;
}

PhotonicsStatus photonicsAddScalar(PhotonicsObjId src, PhotonicsObjId dest, uint64_t scalarValue)
{
  bool ok = photonicsSim::get()->photonicsAdd(src, dest, scalarValue);
  return ok ? PHOTONICS_OK : PHOTONICS_ERROR;
}

PhotonicsStatus photonicsSubScalar(PhotonicsObjId src, PhotonicsObjId dest, uint64_t scalarValue)
{
  bool ok = photonicsSim::get()->photonicsSub(src, dest, scalarValue);
  return ok ? PHOTONICS_OK : PHOTONICS_ERROR;
}

PhotonicsStatus photonicsMulScalar(PhotonicsObjId src, PhotonicsObjId dest, uint64_t scalarValue)
{
  bool ok = photonicsSim::get()->photonicsMul(src, dest, scalarValue);
  return ok ? PHOTONICS_OK : PHOTONICS_ERROR;
}

PhotonicsStatus photonicsDivScalar(PhotonicsObjId src, PhotonicsObjId dest, uint64_t scalarValue)
{
  bool ok = photonicsSim::get()->photonicsDiv(src, dest, scalarValue);
  return ok ? PHOTONICS_OK : PHOTONICS_ERROR;
}

PhotonicsStatus photonicsAndScalar(PhotonicsObjId src, PhotonicsObjId dest, uint64_t scalarValue)
{
  bool ok = photonicsSim::get()->photonicsAnd(src, dest, scalarValue);
  return ok ? PHOTONICS_OK : PHOTONICS_ERROR;
}

PhotonicsStatus photonicsOrScalar(PhotonicsObjId src, PhotonicsObjId dest, uint64_t scalarValue)
{
  bool ok = photonicsSim::get()->photonicsOr(src, dest, scalarValue);
  return ok ? PHOTONICS_OK : PHOTONICS_ERROR;
}

PhotonicsStatus photonicsXorScalar(PhotonicsObjId src, PhotonicsObjId dest, uint64_t scalarValue)
{
  bool ok = photonicsSim::get()->photonicsXor(src, dest, scalarValue);
  return ok ? PHOTONICS_OK : PHOTONICS_ERROR;
}

PhotonicsStatus photonicsXnorScalar(PhotonicsObjId src, PhotonicsObjId dest, uint64_t scalarValue)
{
  bool ok = photonicsSim::get()->photonicsXnor(src, dest, scalarValue);
  return ok ? PHOTONICS_OK : PHOTONICS_ERROR;
}

PhotonicsStatus photonicsGTScalar(PhotonicsObjId src, PhotonicsObjId dest, uint64_t scalarValue)
{
  bool ok = photonicsSim::get()->photonicsGT(src, dest, scalarValue);
  return ok ? PHOTONICS_OK : PHOTONICS_ERROR;
}

PhotonicsStatus photonicsLTScalar(PhotonicsObjId src, PhotonicsObjId dest, uint64_t scalarValue)
{
  bool ok = photonicsSim::get()->photonicsLT(src, dest, scalarValue);
  return ok ? PHOTONICS_OK : PHOTONICS_ERROR;
}

PhotonicsStatus photonicsEQScalar(PhotonicsObjId src, PhotonicsObjId dest, uint64_t scalarValue)
{
  bool ok = photonicsSim::get()->photonicsEQ(src, dest, scalarValue);
  return ok ? PHOTONICS_OK : PHOTONICS_ERROR;
}

PhotonicsStatus photonicsNEScalar(PhotonicsObjId src, PhotonicsObjId dest, uint64_t scalarValue)
{
  bool ok = photonicsSim::get()->photonicsNE(src, dest, scalarValue);
  return ok ? PHOTONICS_OK : PHOTONICS_ERROR;
}

PhotonicsStatus photonicsMinScalar(PhotonicsObjId src, PhotonicsObjId dest, uint64_t scalarValue)
{
  bool ok = photonicsSim::get()->photonicsMin(src, dest, scalarValue);
  return ok ? PHOTONICS_OK : PHOTONICS_ERROR;
}

PhotonicsStatus photonicsMaxScalar(PhotonicsObjId src, PhotonicsObjId dest, uint64_t scalarValue)
{
  bool ok = photonicsSim::get()->photonicsMax(src, dest, scalarValue);
  return ok ? PHOTONICS_OK : PHOTONICS_ERROR;
}

PhotonicsStatus photonicsScaledAdd(PhotonicsObjId src1, PhotonicsObjId src2, PhotonicsObjId dest, uint64_t scalarValue) 
{
  bool ok = photonicsSim::get()->photonicsScaledAdd(src1, src2, dest, scalarValue);
  return ok ? PHOTONICS_OK : PHOTONICS_ERROR;
}

//! @brief  PHOTONICS Pop Count
PhotonicsStatus
photonicsPopCount(PhotonicsObjId src, PhotonicsObjId dest)
{
  bool ok = photonicsSim::get()->photonicsPopCount(src, dest);
  return ok ? PHOTONICS_OK : PHOTONICS_ERROR;
}

//! @brief  Extract a bit slice from a data vector. Dest must be BOOL type
PhotonicsStatus
photonicsBitSliceExtract(PhotonicsObjId src, PhotonicsObjId destBool, unsigned bitIdx)
{
  bool ok = photonicsSim::get()->photonicsBitSliceExtract(src, destBool, bitIdx);
  return ok ? PHOTONICS_OK : PHOTONICS_ERROR;
}

//! @brief  Insert a bit slice to a data vector. Src must be BOOL type
PhotonicsStatus
photonicsBitSliceInsert(PhotonicsObjId srcBool, PhotonicsObjId dest, unsigned bitIdx)
{
  bool ok = photonicsSim::get()->photonicsBitSliceInsert(srcBool, dest, bitIdx);
  return ok ? PHOTONICS_OK : PHOTONICS_ERROR;
}

//! @brief  Conditional copy: dest[i] = cond ? src[i] : dest[i]
PhotonicsStatus
photonicsCondCopy(PhotonicsObjId condBool, PhotonicsObjId src, PhotonicsObjId dest)
{
  bool ok = photonicsSim::get()->photonicsCondCopy(condBool, src, dest);
  return ok ? PHOTONICS_OK : PHOTONICS_ERROR;
}

//! @brief  Conditional broadcast: dest[i] = cond ? scalar : dest[i]
PhotonicsStatus
photonicsCondBroadcast(PhotonicsObjId condBool, uint64_t scalarBits, PhotonicsObjId dest)
{
  bool ok = photonicsSim::get()->photonicsCondBroadcast(condBool, scalarBits, dest);
  return ok ? PHOTONICS_OK : PHOTONICS_ERROR;
}

//! @brief  Conditional select: dest[i] = cond ? src1[i] : src2[i]
PhotonicsStatus
photonicsCondSelect(PhotonicsObjId condBool, PhotonicsObjId src1, PhotonicsObjId src2, PhotonicsObjId dest)
{
  bool ok = photonicsSim::get()->photonicsCondSelect(condBool, src1, src2, dest);
  return ok ? PHOTONICS_OK : PHOTONICS_ERROR;
}

//! @brief  Conditional select scalar: dest[i] = cond ? src1[i] : scalar
PhotonicsStatus
photonicsCondSelectScalar(PhotonicsObjId condBool, PhotonicsObjId src1, uint64_t scalarBits, PhotonicsObjId dest)
 {
  bool ok = photonicsSim::get()->photonicsCondSelectScalar(condBool, src1, scalarBits, dest);
  return ok ? PHOTONICS_OK : PHOTONICS_ERROR;
 }

//! @brief  AES Sbox: dest[i] = lut[src[i]]
PhotonicsStatus 
photonicsAesSbox(PhotonicsObjId src, PhotonicsObjId dest, const std::vector<uint8_t>& lut)
{
  bool ok = photonicsSim::get()->photonicsAesSbox(src, dest, lut);
  return ok ? PHOTONICS_OK : PHOTONICS_ERROR;
}

//! @brief  AES Sbox: dest[i] = lut[src[i]] (similar to AES sbox, different in perforamance and energy model for the bit-serial architecture)
PhotonicsStatus 
photonicsAesInverseSbox(PhotonicsObjId src, PhotonicsObjId dest, const std::vector<uint8_t>& lut)
{
  bool ok = photonicsSim::get()->photonicsAesInverseSbox(src, dest, lut);
  return ok ? PHOTONICS_OK : PHOTONICS_ERROR;
}

// Implementation of min reduction
PhotonicsStatus photonicsRedMin(PhotonicsObjId src, void* min, uint64_t idxBegin, uint64_t idxEnd) {
    bool ok = photonicsSim::get()->photonicsRedMin(src, min, idxBegin, idxEnd);
    return ok ? PHOTONICS_OK : PHOTONICS_ERROR;
}

// Implementation of max reduction
PhotonicsStatus photonicsRedMax(PhotonicsObjId src, void* max, uint64_t idxBegin, uint64_t idxEnd) {
    bool ok = photonicsSim::get()->photonicsRedMax(src, max, idxBegin, idxEnd);
    return ok ? PHOTONICS_OK : PHOTONICS_ERROR;
}

//! @brief  PHOTONICS MAC operation: dest += src1 * src2
PhotonicsStatus photonicsMAC(PhotonicsObjId src1, PhotonicsObjId src2, void *dest)
{
  bool ok = photonicsSim::get()->photonicsMAC(src1, src2, dest);
  return ok ? PHOTONICS_OK : PHOTONICS_ERROR;
}

//! @brief  PHOTONICS reduction sum for signed int. Result returned to a host variable
PhotonicsStatus
photonicsRedSum(PhotonicsObjId src, void* sum, uint64_t idxBegin, uint64_t idxEnd)
{
  bool ok = photonicsSim::get()->photonicsRedSum(src, sum, idxBegin, idxEnd);
  return ok ? PHOTONICS_OK : PHOTONICS_ERROR;
}

//! @brief  Rotate all elements of an obj by one step to the right
PhotonicsStatus
photonicsRotateElementsRight(PhotonicsObjId src)
{
  bool ok = photonicsSim::get()->photonicsRotateElementsRight(src);
  return ok ? PHOTONICS_OK : PHOTONICS_ERROR;
}

//! @brief  Rotate all elements of an obj by one step to the left
PhotonicsStatus
photonicsRotateElementsLeft(PhotonicsObjId src)
{
  bool ok = photonicsSim::get()->photonicsRotateElementsLeft(src);
  return ok ? PHOTONICS_OK : PHOTONICS_ERROR;
}

//! @brief  Shift elements of an obj by one step to the right and fill zero
PhotonicsStatus
photonicsShiftElementsRight(PhotonicsObjId src)
{
  bool ok = photonicsSim::get()->photonicsShiftElementsRight(src);
  return ok ? PHOTONICS_OK : PHOTONICS_ERROR;
}

//! @brief  Shift elements of an obj by one step to the left and fill zero
PhotonicsStatus
photonicsShiftElementsLeft(PhotonicsObjId src)
{
  bool ok = photonicsSim::get()->photonicsShiftElementsLeft(src);
  return ok ? PHOTONICS_OK : PHOTONICS_ERROR;
}

//! @brief  Shift bits of each elements of an obj by shiftAmount to the right. This currently implements arithmetic shift.
PhotonicsStatus
photonicsShiftBitsRight(PhotonicsObjId src, PhotonicsObjId dest, unsigned shiftAmount)
{
  bool ok = photonicsSim::get()->photonicsShiftBitsRight(src, dest, shiftAmount);
  return ok ? PHOTONICS_OK : PHOTONICS_ERROR;
}

//! @brief  Shift bits of each elements of an obj by shiftAmount to the left.
PhotonicsStatus
photonicsShiftBitsLeft(PhotonicsObjId src, PhotonicsObjId dest, unsigned shiftAmount)
{
  bool ok = photonicsSim::get()->photonicsShiftBitsLeft(src, dest, shiftAmount);
  return ok ? PHOTONICS_OK : PHOTONICS_ERROR;
}

//! @brief  Execute fused PHOTONICS APIs
PhotonicsStatus
photonicsFuse(PhotonicsProg prog)
{
  bool ok = photonicsSim::get()->photonicsFuse(prog);
  return ok ? PHOTONICS_OK : PHOTONICS_ERROR;
}

//! @brief  BitSIMD-V: Read a row to SA
PhotonicsStatus
photonicsOpReadRowToSa(PhotonicsObjId src, unsigned ofst)
{
  bool ok = photonicsSim::get()->photonicsOpReadRowToSa(src, ofst);
  return ok ? PHOTONICS_OK : PHOTONICS_ERROR;
}

//! @brief  BitSIMD-V: Write SA to a row
PhotonicsStatus
photonicsOpWriteSaToRow(PhotonicsObjId src, unsigned ofst)
{
  bool ok = photonicsSim::get()->photonicsOpWriteSaToRow(src, ofst);
  return ok ? PHOTONICS_OK : PHOTONICS_ERROR;
}

//! @brief  BitSIMD-V: Triple row activation to SA
PhotonicsStatus
photonicsOpTRA(PhotonicsObjId src1, unsigned ofst1, PhotonicsObjId src2, unsigned ofst2, PhotonicsObjId src3, unsigned ofst3)
{
  bool ok = photonicsSim::get()->photonicsOpTRA(src1, ofst1, src2, ofst2, src3, ofst3);
  return ok ? PHOTONICS_OK : PHOTONICS_ERROR;
}

//! @brief  BitSIMD-V: Move value between two regs
PhotonicsStatus
photonicsOpMove(PhotonicsObjId objId, PhotonicsRowReg src, PhotonicsRowReg dest)
{
  bool ok = photonicsSim::get()->photonicsOpMove(objId, src, dest);
  return ok ? PHOTONICS_OK : PHOTONICS_ERROR;
}

//! @brief  BitSIMD-V: Set value of a reg
PhotonicsStatus
photonicsOpSet(PhotonicsObjId objId, PhotonicsRowReg src, bool val)
{
  bool ok = photonicsSim::get()->photonicsOpSet(objId, src, val);
  return ok ? PHOTONICS_OK : PHOTONICS_ERROR;
}

//! @brief  BitSIMD-V: Not of a reg
PhotonicsStatus
photonicsOpNot(PhotonicsObjId objId, PhotonicsRowReg src, PhotonicsRowReg dest)
{
  bool ok = photonicsSim::get()->photonicsOpNot(objId, src, dest);
  return ok ? PHOTONICS_OK : PHOTONICS_ERROR;
}

//! @brief  BitSIMD-V: And of two regs
PhotonicsStatus
photonicsOpAnd(PhotonicsObjId objId, PhotonicsRowReg src1, PhotonicsRowReg src2, PhotonicsRowReg dest)
{
  bool ok = photonicsSim::get()->photonicsOpAnd(objId, src1, src2, dest);
  return ok ? PHOTONICS_OK : PHOTONICS_ERROR;
}

//! @brief  BitSIMD-V: Or of two regs
PhotonicsStatus
photonicsOpOr(PhotonicsObjId objId, PhotonicsRowReg src1, PhotonicsRowReg src2, PhotonicsRowReg dest)
{
  bool ok = photonicsSim::get()->photonicsOpOr(objId, src1, src2, dest);
  return ok ? PHOTONICS_OK : PHOTONICS_ERROR;
}

//! @brief  BitSIMD-V: Nand of two regs
PhotonicsStatus
photonicsOpNand(PhotonicsObjId objId, PhotonicsRowReg src1, PhotonicsRowReg src2, PhotonicsRowReg dest)
{
  bool ok = photonicsSim::get()->photonicsOpNand(objId, src1, src2, dest);
  return ok ? PHOTONICS_OK : PHOTONICS_ERROR;
}

//! @brief  BitSIMD-V: Nor of two regs
PhotonicsStatus
photonicsOpNor(PhotonicsObjId objId, PhotonicsRowReg src1, PhotonicsRowReg src2, PhotonicsRowReg dest)
{
  bool ok = photonicsSim::get()->photonicsOpNor(objId, src1, src2, dest);
  return ok ? PHOTONICS_OK : PHOTONICS_ERROR;
}

//! @brief  BitSIMD-V: Xor of two regs
PhotonicsStatus
photonicsOpXor(PhotonicsObjId objId, PhotonicsRowReg src1, PhotonicsRowReg src2, PhotonicsRowReg dest)
{
  bool ok = photonicsSim::get()->photonicsOpXor(objId, src1, src2, dest);
  return ok ? PHOTONICS_OK : PHOTONICS_ERROR;
}

//! @brief  BitSIMD-V: Xnor of two regs
PhotonicsStatus
photonicsOpXnor(PhotonicsObjId objId, PhotonicsRowReg src1, PhotonicsRowReg src2, PhotonicsRowReg dest)
{
  bool ok = photonicsSim::get()->photonicsOpXnor(objId, src1, src2, dest);
  return ok ? PHOTONICS_OK : PHOTONICS_ERROR;
}

//! @brief  BitSIMD-V: Maj of three regs
PhotonicsStatus
photonicsOpMaj(PhotonicsObjId objId, PhotonicsRowReg src1, PhotonicsRowReg src2, PhotonicsRowReg src3, PhotonicsRowReg dest)
{
  bool ok = photonicsSim::get()->photonicsOpMaj(objId, src1, src2, src3, dest);
  return ok ? PHOTONICS_OK : PHOTONICS_ERROR;
}

//! @brief  BitSIMD-V: Conditional selecion: dest = cond ? src1 : src2
PhotonicsStatus
photonicsOpSel(PhotonicsObjId objId, PhotonicsRowReg cond, PhotonicsRowReg src1, PhotonicsRowReg src2, PhotonicsRowReg dest)
{
  bool ok = photonicsSim::get()->photonicsOpSel(objId, cond, src1, src2, dest);
  return ok ? PHOTONICS_OK : PHOTONICS_ERROR;
}

//! @brief  BitSIMD-V: Rotate a reg to the right, using srcId for range
PhotonicsStatus
photonicsOpRotateRH(PhotonicsObjId objId, PhotonicsRowReg src)
{
  bool ok = photonicsSim::get()->photonicsOpRotateRH(objId, src);
  return ok ? PHOTONICS_OK : PHOTONICS_ERROR;
}

//! @brief  BitSIMD-V: Rotate a reg to the left, using srcId for range
PhotonicsStatus
photonicsOpRotateLH(PhotonicsObjId objId, PhotonicsRowReg src)
{
  bool ok = photonicsSim::get()->photonicsOpRotateLH(objId, src);
  return ok ? PHOTONICS_OK : PHOTONICS_ERROR;
}

// @brief  SIMDRAM: AP operation
PhotonicsStatus
photonicsOpAP(int numSrc, ...)
{
  va_list args;
  va_start(args, numSrc);
  bool ok = photonicsSim::get()->photonicsOpAP(numSrc, args);
  va_end(args);
  return ok ? PHOTONICS_OK : PHOTONICS_ERROR;
}

// @brief  SIMDRAM: AAP operation
PhotonicsStatus
photonicsOpAAP(int numSrc, int numDest, ...)
{
  va_list args;
  va_start(args, numDest);
  bool ok = photonicsSim::get()->photonicsOpAAP(numSrc, numDest, args);
  va_end(args);
  return ok ? PHOTONICS_OK : PHOTONICS_ERROR;
}

