// File: libphotonicseval.h
// PHOTONICSeval Simulator - Library Interface
// Copyright (c) 2024 University of Virginia
// This file is licensed under the MIT License.
// See the LICENSE file in the root of this repository for more details.

#ifndef LAVA_LIB_PHOTONICS_EVAL_H
#define LAVA_LIB_PHOTONICS_EVAL_H

#include <cstdint>
#include <cstdarg>
#include <vector>
#include <functional>

//! @brief  Photonics API return status
enum PhotonicsStatus {
  PHOTONICS_ERROR = 0,
  PHOTONICS_OK,
};

//! @brief  Photonics device types
enum PhotonicsDeviceEnum {
  PHOTONICS_DEVICE_NONE = 0,
  PHOTONICS_FUNCTIONAL,
  PHOTONICS_DEVICE_BITSIMD_V,
  PHOTONICS_DEVICE_BITSIMD_V_NAND,
  PHOTONICS_DEVICE_BITSIMD_V_MAJ,
  PHOTONICS_DEVICE_BITSIMD_V_AP,
  PHOTONICS_DEVICE_DRISA_NOR,
  PHOTONICS_DEVICE_DRISA_MIXED,
  PHOTONICS_DEVICE_SIMDRAM,
  PHOTONICS_DEVICE_BITSIMD_H,
  PHOTONICS_DEVICE_FULCRUM,
  PHOTONICS_DEVICE_BANK_LEVEL,
  PHOTONICS_DEVICE_AQUABOLT,
  PHOTONICS_DEVICE_AIM,
};

/**
 * @enum PhotonicsDeviceProtocol
 * @brief Enum representing different memory protocols.
 *
 * @var PHOTONICS_DEVICE_PROTOCOL_DDR
 * Standard DDR protocol. Typically used in general-purpose memory systems.
 *
 * @var PHOTONICS_DEVICE_PROTOCOL_LPDDR
 * Low Power DDR (LPDDR) protocol.
 *
 * @var PHOTONICS_DEVICE_PROTOCOL_HBM
 * High Bandwidth Memory (HBM) protocol.
 * 
 * @var PHOTONICS_DEVICE_PROTOCOL_GDDR
 * Graphics Double Data Rate (GDDR) protocol.
*/
enum PhotonicsDeviceProtocolEnum {
  PHOTONICS_DEVICE_PROTOCOL_DDR = 0,
  PHOTONICS_DEVICE_PROTOCOL_LPDDR,
  PHOTONICS_DEVICE_PROTOCOL_HBM,
  PHOTONICS_DEVICE_PROTOCOL_GDDR,
};

//! @brief  PHOTONICS allocation types
enum PhotonicsAllocEnum {
  PHOTONICS_ALLOC_AUTO = 0, // Auto determine vertical or horizontal layout based on device type
  PHOTONICS_ALLOC_V,        // V layout, multiple regions per core
  PHOTONICS_ALLOC_H,        // H layout, multiple regions per core
  PHOTONICS_ALLOC_V1,       // V layout, at most 1 region per core
  PHOTONICS_ALLOC_H1,       // H layout, at most 1 region per core
};

//! @brief  PHOTONICS data copy types
enum PhotonicsCopyEnum {
  PHOTONICS_COPY_V,
  PHOTONICS_COPY_H,
};

//! @brief  PHOTONICS datatypes
enum PhotonicsDataType {
  PHOTONICS_BOOL = 0,
  PHOTONICS_INT8,
  PHOTONICS_INT16,
  PHOTONICS_INT32,
  PHOTONICS_INT64,
  PHOTONICS_UINT8,
  PHOTONICS_UINT16,
  PHOTONICS_UINT32,
  PHOTONICS_UINT64,
  PHOTONICS_FP32,
  PHOTONICS_FP16,
  PHOTONICS_BF16,
  PHOTONICS_FP8,
};

//! @brief  PHOTONICS device properties
struct PhotonicsDeviceProperties {
  PhotonicsDeviceEnum deviceType = PHOTONICS_DEVICE_NONE;
  PhotonicsDeviceEnum simTarget = PHOTONICS_DEVICE_NONE;
  unsigned numRanks = 0;
  unsigned numBankPerRank = 0;
  unsigned numSubarrayPerBank = 0;
  unsigned numRowPerSubarray = 0;
  unsigned numColPerSubarray = 0;
  unsigned numPHOTONICSCores = 0;
  bool isHLayoutDevice = false;
};

typedef int PhotonicsCoreId;
typedef int PhotonicsObjId;

// PHOTONICSeval simulation
// CPU runtime between start/end timer will be measured for modeling DRAM refresh
void photonicsStartTimer();
void photonicsEndTimer();
void photonicsShowStats();
void photonicsResetStats();
bool photonicsIsAnalysisMode();

// Device creation and deletion
/**
 * @brief Creates and initializes a PHOTONICS (Processing-In-Memory) device with the specified configuration.
 *
 * @param deviceType      The type of PHOTONICS device to create (see PhotonicsDeviceEnum).
 * @param numRanks        Number of ranks in the device.
 * @param numBankPerRank  Number of banks per rank.
 * @param numSubarrayPerBank Number of subarrays per bank.
 * @param numRows         Number of rows in each subarray.
 * @param numCols         Number of columns in each row.
 * @param bufferSize      Optional on-chip buffer size (B) for the device (default is 0). This parameter is only applicable for AiM.
 * @return PhotonicsStatus      Status code indicating success or failure of device creation.
 */
PhotonicsStatus photonicsCreateDevice(PhotonicsDeviceEnum deviceType, unsigned numRanks, unsigned numBankPerRank, unsigned numSubarrayPerBank, unsigned numRows, unsigned numCols, unsigned bufferSize = 0);
PhotonicsStatus photonicsCreateDeviceFromConfig(PhotonicsDeviceEnum deviceType, const char* configFileName);
PhotonicsStatus photonicsGetDeviceProperties(PhotonicsDeviceProperties* deviceProperties);
PhotonicsStatus photonicsDeleteDevice();

// Resource allocation and deletion
PhotonicsObjId photonicsAlloc(PhotonicsAllocEnum allocType, uint64_t numElements, PhotonicsDataType dataType);
PhotonicsObjId photonicsAllocMat(uint64_t numElements, PhotonicsDataType dataType);
PhotonicsObjId photonicsAllocAssociated(PhotonicsObjId assocId, PhotonicsDataType dataType);
PhotonicsObjId photonicsAllocAssociatedSrcVec(PhotonicsObjId assocId, PhotonicsDataType dataType);
PhotonicsObjId photonicsAllocAssociatedDestVec(PhotonicsObjId assocId, PhotonicsDataType dataType);
// Buffer will always be allocated in H layout; Current assumption is buffer is global and shared across all PHOTONICS cores in a chip/device. This assumption is based on AiM.
// The buffer is used for broadcasting data to all PHOTONICS cores in a chip/device.
// Please note that each chip/device will hold the same data in their respective buffers.
// TODO: Support per-core buffers (like UPMEM)
PhotonicsObjId photonicsAllocBuffer(uint32_t numElements, PhotonicsDataType dataType);
PhotonicsStatus photonicsFreeMat(PhotonicsObjId obj);
PhotonicsStatus photonicsFreeSrcVec(PhotonicsObjId obj);
PhotonicsStatus photonicsFreeDestVec(PhotonicsObjId obj);

// Data transfer
// Note: idxBegin and idxEnd specify the range of indexes to be processed by the PHOTONICS.
// The size of the host-side vector should match the size of this range on the PHOTONICS side.
// If the default values for idxBegin and idxEnd are used, the entire range of the PHOTONICS object will be considered.
// For PHOTONICS_BOOL type, please use std::vector<uint8_t> instead of std::vector<bool> as host data.
PhotonicsStatus photonicsCopyHostToDevice(void* src, PhotonicsObjId dest, uint64_t idxBegin = 0, uint64_t idxEnd = 0);
PhotonicsStatus photonicsCopyDeviceToHost(PhotonicsObjId src, void* dest, uint64_t idxBegin = 0, uint64_t idxEnd = 0);
PhotonicsStatus photonicsCopyHostToDeviceMat(void* src, PhotonicsObjId dest, uint64_t idxBegin = 0, uint64_t idxEnd = 0);
PhotonicsStatus photonicsCopyHostToDeviceVec(void* src, PhotonicsObjId dest, uint64_t idxBegin = 0, uint64_t idxEnd = 0);
PhotonicsStatus photonicsCopyDeviceToHostVec(PhotonicsObjId src, void* dest, uint64_t idxBegin = 0, uint64_t idxEnd = 0);
PhotonicsStatus photonicsCopyDeviceToDevice(PhotonicsObjId src, PhotonicsObjId dest, uint64_t idxBegin = 0, uint64_t idxEnd = 0);
PhotonicsStatus photonicsCopyObjectToObject(PhotonicsObjId src, PhotonicsObjId dest);
PhotonicsStatus photonicsConvertType(PhotonicsObjId src, PhotonicsObjId dest);

// Logic and Arithmetic Operation
// Mixed data type extensions:
// - photonicsAdd/photonicsSub: If src1 is an integer vector, src2 can be a Boolean vector for accumulation purposes.
PhotonicsStatus photonicsIter(PhotonicsObjId src1, PhotonicsObjId src2, PhotonicsObjId dest, int8_t numLoops = 1);
PhotonicsStatus photonicsAdd(PhotonicsObjId src1, PhotonicsObjId src2, PhotonicsObjId dest);
PhotonicsStatus photonicsSub(PhotonicsObjId src1, PhotonicsObjId src2, PhotonicsObjId dest);
PhotonicsStatus photonicsMul(PhotonicsObjId src1, PhotonicsObjId src2, PhotonicsObjId dest);
PhotonicsStatus photonicsDiv(PhotonicsObjId src1, PhotonicsObjId src2, PhotonicsObjId dest);
PhotonicsStatus photonicsAbs(PhotonicsObjId src, PhotonicsObjId dest);
PhotonicsStatus photonicsNot(PhotonicsObjId src, PhotonicsObjId dest);
PhotonicsStatus photonicsAnd(PhotonicsObjId src1, PhotonicsObjId src2, PhotonicsObjId dest);
PhotonicsStatus photonicsOr(PhotonicsObjId src1, PhotonicsObjId src2, PhotonicsObjId dest);
PhotonicsStatus photonicsXor(PhotonicsObjId src1, PhotonicsObjId src2, PhotonicsObjId dest);
PhotonicsStatus photonicsXnor(PhotonicsObjId src1, PhotonicsObjId src2, PhotonicsObjId dest);
PhotonicsStatus photonicsMin(PhotonicsObjId src1, PhotonicsObjId src2, PhotonicsObjId dest);
PhotonicsStatus photonicsMax(PhotonicsObjId src1, PhotonicsObjId src2, PhotonicsObjId dest);
PhotonicsStatus photonicsAddScalar(PhotonicsObjId src, PhotonicsObjId dest, uint64_t scalarValue);
PhotonicsStatus photonicsSubScalar(PhotonicsObjId src, PhotonicsObjId dest, uint64_t scalarValue);
PhotonicsStatus photonicsMulScalar(PhotonicsObjId src, PhotonicsObjId dest, uint64_t scalarValue);
PhotonicsStatus photonicsDivScalar(PhotonicsObjId src, PhotonicsObjId dest, uint64_t scalarValue);
PhotonicsStatus photonicsAndScalar(PhotonicsObjId src, PhotonicsObjId dest, uint64_t scalarValue);
PhotonicsStatus photonicsOrScalar(PhotonicsObjId src, PhotonicsObjId dest, uint64_t scalarValue);
PhotonicsStatus photonicsXorScalar(PhotonicsObjId src, PhotonicsObjId dest, uint64_t scalarValue);
PhotonicsStatus photonicsXnorScalar(PhotonicsObjId src, PhotonicsObjId dest, uint64_t scalarValue);
PhotonicsStatus photonicsMinScalar(PhotonicsObjId src, PhotonicsObjId dest, uint64_t scalarValue);
PhotonicsStatus photonicsMaxScalar(PhotonicsObjId src, PhotonicsObjId dest, uint64_t scalarValue);

// Relational operations - Dest object is BOOL type
PhotonicsStatus photonicsGT(PhotonicsObjId src1, PhotonicsObjId src2, PhotonicsObjId destBool);
PhotonicsStatus photonicsLT(PhotonicsObjId src1, PhotonicsObjId src2, PhotonicsObjId destBool);
PhotonicsStatus photonicsEQ(PhotonicsObjId src1, PhotonicsObjId src2, PhotonicsObjId destBool);
PhotonicsStatus photonicsNE(PhotonicsObjId src1, PhotonicsObjId src2, PhotonicsObjId destBool);
PhotonicsStatus photonicsGTScalar(PhotonicsObjId src, PhotonicsObjId destBool, uint64_t scalarValue);
PhotonicsStatus photonicsLTScalar(PhotonicsObjId src, PhotonicsObjId destBool, uint64_t scalarValue);
PhotonicsStatus photonicsEQScalar(PhotonicsObjId src, PhotonicsObjId destBool, uint64_t scalarValue);
PhotonicsStatus photonicsNEScalar(PhotonicsObjId src, PhotonicsObjId destBool, uint64_t scalarValue);

// multiply src1 with scalarValue and add the multiplication result with src2. Save the result to dest. 
PhotonicsStatus photonicsScaledAdd(PhotonicsObjId src1, PhotonicsObjId src2, PhotonicsObjId dest, uint64_t scalarValue);
PhotonicsStatus photonicsPopCount(PhotonicsObjId src, PhotonicsObjId dest);

// Only supported by bit-parallel PHOTONICS
PhotonicsStatus photonicsPrefixSum(PhotonicsObjId src, PhotonicsObjId dest);

// MAC operation: dest += src1 * src2
// Note: src2 is a global buffer that holds a vector of values to be multiplied with src1.
// Note: dest must be of the same data type as src1 and src2; Size of dest must be equal to the total number of PHOTONICS cores in the device.
// Note: The MAC operation is performed in parallel across all PHOTONICS cores, and each PHOTONICS core writes its local MAC value to the specific id of the dest.
// Note: User needs to ensure that dest vector is of size equal to the total number of PHOTONICS cores in the device, and contains `0` or any value that the user wants it to have as initial values.
PhotonicsStatus photonicsMAC(PhotonicsObjId src1, PhotonicsObjId src2, void* dest);

// Note: Reduction sum range is [idxBegin, idxEnd)
PhotonicsStatus photonicsRedSum(PhotonicsObjId src, void* sum, uint64_t idxBegin = 0, uint64_t idxEnd = 0);
// Min/Max Reduction APIs
PhotonicsStatus photonicsRedMin(PhotonicsObjId src, void* min, uint64_t idxBegin = 0, uint64_t idxEnd = 0);
PhotonicsStatus photonicsRedMax(PhotonicsObjId src, void* max, uint64_t idxBegin = 0, uint64_t idxEnd = 0);

// Bit slice operations
PhotonicsStatus photonicsBitSliceExtract(PhotonicsObjId src, PhotonicsObjId destBool, unsigned bitIdx);
PhotonicsStatus photonicsBitSliceInsert(PhotonicsObjId srcBool, PhotonicsObjId dest, unsigned bitIdx);
// Conditional operations
// photonicsCondCopy:         dest[i] = cond ? src[i] : dest[i]
// photonicsCondBroadcast:    dest[i] = cond ? scalar : dest[i]
// photonicsCondSelect:       dest[i] = cond ? src1[i] : src2[i]
// photonicsCondSelectScalar: dest[i] = cond ? src[i] : scalar
PhotonicsStatus photonicsCondCopy(PhotonicsObjId condBool, PhotonicsObjId src, PhotonicsObjId dest);
PhotonicsStatus photonicsCondBroadcast(PhotonicsObjId condBool, uint64_t scalarBits, PhotonicsObjId dest);
PhotonicsStatus photonicsCondSelect(PhotonicsObjId condBool, PhotonicsObjId src1, PhotonicsObjId src2, PhotonicsObjId dest);
PhotonicsStatus photonicsCondSelectScalar(PhotonicsObjId condBool, PhotonicsObjId src1, uint64_t scalarBits, PhotonicsObjId dest);

PhotonicsStatus photonicsBroadcastInt(PhotonicsObjId dest, int64_t value);
PhotonicsStatus photonicsBroadcastUInt(PhotonicsObjId dest, uint64_t value);
PhotonicsStatus photonicsBroadcastFP(PhotonicsObjId dest, float value);
PhotonicsStatus photonicsRotateElementsRight(PhotonicsObjId src);
PhotonicsStatus photonicsRotateElementsLeft(PhotonicsObjId src);
PhotonicsStatus photonicsShiftElementsRight(PhotonicsObjId src);
PhotonicsStatus photonicsShiftElementsLeft(PhotonicsObjId src);
PhotonicsStatus photonicsShiftBitsRight(PhotonicsObjId src, PhotonicsObjId dest, unsigned shiftAmount);
PhotonicsStatus photonicsShiftBitsLeft(PhotonicsObjId src, PhotonicsObjId dest, unsigned shiftAmount);

// AES sbox and inverse-box APIs
// Note: AES S-box and inverse S-box are treated separately because their bit-serial performance models differ.
// However, it is the user's responsibility to provide the appropriate LUT to ensure correct functionality.
// The function photonicsAesInverseSbox expects an inverse S-box LUT as its input.
PhotonicsStatus photonicsAesSbox(PhotonicsObjId src, PhotonicsObjId dest, const std::vector<uint8_t>& lut); 
PhotonicsStatus photonicsAesInverseSbox(PhotonicsObjId src, PhotonicsObjId dest, const std::vector<uint8_t>& lut); 

////////////////////////////////////////////////////////////////////////////////
// Experimental Feature: PHOTONICS API Fusion                                       //
////////////////////////////////////////////////////////////////////////////////
struct PhotonicsProg {
  template <typename... Args>
  void add(PhotonicsStatus(*api)(Args...), Args... args) {
    m_apis.push_back([=]() { return api(args...); });
  }
  std::vector<std::function<PhotonicsStatus()>> m_apis;
};
PhotonicsStatus photonicsFuse(PhotonicsProg prog);

////////////////////////////////////////////////////////////////////////////////
// Warning: Avoid using below customized APIs for functional simulation       //
//          Some are PHOTONICS architecture dependent, some are in progress         //
////////////////////////////////////////////////////////////////////////////////

// Data copy APIs that supports data transposition between V/H layout
PhotonicsStatus photonicsCopyHostToDeviceWithType(PhotonicsCopyEnum copyType, void* src, PhotonicsObjId dest, uint64_t idxBegin = 0, uint64_t idxEnd = 0);
PhotonicsStatus photonicsCopyDeviceToHostWithType(PhotonicsCopyEnum copyType, PhotonicsObjId src, void* dest, uint64_t idxBegin = 0, uint64_t idxEnd = 0);

// Dual contact reference: Create a new PhotonicsObjId that references to the negation of the original PhotonicsObjId
// Do not use a dual contact reference PhotonicsObjId as refId
PhotonicsObjId photonicsCreateDualContactRef(PhotonicsObjId refId);

// Ranged reference: Create a new PhotonicsObjId that references to a range of the original PhotonicsObjId
// This is not available for now
PhotonicsObjId photonicsCreateRangedRef(PhotonicsObjId refId, uint64_t idxBegin, uint64_t idxEnd);


////////////////////////////////////////////////////////////////////////////////
// Warning: Do not use below micro-ops level APIs for functional simulation   //
////////////////////////////////////////////////////////////////////////////////

// BitSIMD micro ops
// Note: Below APIs are for low-level micro-ops programming but not for functional simulation
// BitSIMD-V: Row-wide bit registers per subarray
enum PhotonicsRowReg {
  PHOTONICS_RREG_NONE = 0,
  PHOTONICS_RREG_SA,
  PHOTONICS_RREG_R1,
  PHOTONICS_RREG_R2,
  PHOTONICS_RREG_R3,
  PHOTONICS_RREG_R4,
  PHOTONICS_RREG_R5,
  PHOTONICS_RREG_R6,
  PHOTONICS_RREG_R7,
  PHOTONICS_RREG_R8,
  PHOTONICS_RREG_R9,
  PHOTONICS_RREG_R10,
  PHOTONICS_RREG_R11,
  PHOTONICS_RREG_R12,
  PHOTONICS_RREG_R13,
  PHOTONICS_RREG_R14,
  PHOTONICS_RREG_R15,
  PHOTONICS_RREG_R16,
  PHOTONICS_RREG_MAX
};

// BitSIMD-V micro ops
PhotonicsStatus photonicsOpReadRowToSa(PhotonicsObjId src, unsigned ofst);
PhotonicsStatus photonicsOpWriteSaToRow(PhotonicsObjId src, unsigned ofst);
PhotonicsStatus photonicsOpTRA(PhotonicsObjId src1, unsigned ofst1, PhotonicsObjId src2, unsigned ofst2, PhotonicsObjId src3, unsigned ofst3);
PhotonicsStatus photonicsOpMove(PhotonicsObjId objId, PhotonicsRowReg src, PhotonicsRowReg dest);
PhotonicsStatus photonicsOpSet(PhotonicsObjId objId, PhotonicsRowReg src, bool val);
PhotonicsStatus photonicsOpNot(PhotonicsObjId objId, PhotonicsRowReg src, PhotonicsRowReg dest);
PhotonicsStatus photonicsOpAnd(PhotonicsObjId objId, PhotonicsRowReg src1, PhotonicsRowReg src2, PhotonicsRowReg dest);
PhotonicsStatus photonicsOpOr(PhotonicsObjId objId, PhotonicsRowReg src1, PhotonicsRowReg src2, PhotonicsRowReg dest);
PhotonicsStatus photonicsOpNand(PhotonicsObjId objId, PhotonicsRowReg src1, PhotonicsRowReg src2, PhotonicsRowReg dest);
PhotonicsStatus photonicsOpNor(PhotonicsObjId objId, PhotonicsRowReg src1, PhotonicsRowReg src2, PhotonicsRowReg dest);
PhotonicsStatus photonicsOpXor(PhotonicsObjId objId, PhotonicsRowReg src1, PhotonicsRowReg src2, PhotonicsRowReg dest);
PhotonicsStatus photonicsOpXnor(PhotonicsObjId objId, PhotonicsRowReg src1, PhotonicsRowReg src2, PhotonicsRowReg dest);
PhotonicsStatus photonicsOpMaj(PhotonicsObjId objId, PhotonicsRowReg src1, PhotonicsRowReg src2, PhotonicsRowReg src3, PhotonicsRowReg dest);
PhotonicsStatus photonicsOpSel(PhotonicsObjId objId, PhotonicsRowReg cond, PhotonicsRowReg src1, PhotonicsRowReg src2, PhotonicsRowReg dest);
PhotonicsStatus photonicsOpRotateRH(PhotonicsObjId objId, PhotonicsRowReg src);
PhotonicsStatus photonicsOpRotateLH(PhotonicsObjId objId, PhotonicsRowReg src);

// SIMDRAM micro ops
// AP:
//   - Functionality: {srcRows} = MAJ(srcRows)
//   - Action: Activate srcRows simultaneously, followed by a precharge
//   - Example: photonicsOpAP(3, T0, 0, T1, 0, T2, 0) // T0, T1, T2 = MAJ(T0, T1, T2)
// AAP:
//   - Functionality: {srcRows, destRows} = MAJ(srcRows)
//   - Action: Activate srcRows simultaneously, copy result to all destRows, followed by a precharge
//   - Example: photonicsOpAAP(1, 2, DCC0N, 0, T0, 0, T3, 0) // T0, T3 = DCC0N
// Requirements:
//   - numSrc must be odd (1 or 3) to perform MAJ operation
//   - Number of var args must be 2*numSrc for AP and 2*(numDest+numSrc) for AAP
//   - Var args must be a list of (PhotonicsObjId, unsigned ofst) pairs
PhotonicsStatus photonicsOpAP(int numSrc, ...);
PhotonicsStatus photonicsOpAAP(int numSrc, int numDest, ...);

#endif

