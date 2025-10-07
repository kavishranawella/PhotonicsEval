// File: photonicsSim.h
// PHOTONICSeval Simulator - PHOTONICS Simulator Main Entry
// Copyright (c) 2024 University of Virginia
// This file is licensed under the MIT License.
// See the LICENSE file in the root of this repository for more details.

#ifndef LAVA_PHOTONICS_SIM_H
#define LAVA_PHOTONICS_SIM_H

#include "libphotonicseval.h"
#include "photonicsSimConfig.h"
#include "photonicsDevice.h"
#include "photonicsParamsDram.h"
#include "photonicsPerfEnergyBase.h"
#include "photonicsStats.h"
#include <cstdarg>
#include <memory>


//! @class  photonicsSim
//! @brief  PHOTONICS simulator singleton class
class photonicsSim
{
public:
  static photonicsSim* get();
  static void destroy();

  // Device creation and deletion
  bool createDevice(PhotonicsDeviceEnum deviceType, unsigned numRanks, unsigned numBankPerRank, unsigned numSubarrayPerBank, unsigned numRows, unsigned numCols, unsigned bufferSize);
  bool createDeviceFromConfig(PhotonicsDeviceEnum deviceType, const char* configFileName);
  bool getDeviceProperties(PhotonicsDeviceProperties* deviceProperties);
  bool deleteDevice();
  bool isValidDevice(bool showMsg = true) const;

  // From photonicsSimConfig
  const photonicsSimConfig& getConfig() const { return m_config; }
  PhotonicsDeviceEnum getDeviceType() const { return m_config.getDeviceType(); }
  PhotonicsDeviceEnum getSimTarget() const { return m_config.getSimTarget(); }
  unsigned getNumRanks() const { return m_config.getNumRanks(); }
  unsigned getNumBankPerRank() const { return m_config.getNumBankPerRank(); }
  unsigned getNumSubarrayPerBank() const { return m_config.getNumSubarrayPerBank(); }
  unsigned getNumRowPerSubarray() const { return m_config.getNumRowPerSubarray(); }
  unsigned getNumColPerSubarray() const { return m_config.getNumColPerSubarray(); }
  bool isAnalysisMode() const { return m_config.isAnalysisMode(); }
  unsigned getNumThreads() const { return m_config.getNumThreads(); }
  bool isDebug(photonicsSimConfig::photonicsDebugFlags flag) const { return m_config.getDebug() & flag; }

  unsigned getNumCores() const;
  unsigned getNumRows() const;
  unsigned getNumCols() const;

  void startKernelTimer() const;
  void endKernelTimer() const;
  void showStats() const;
  void resetStats() const;
  photonicsStatsMgr* getStatsMgr() { return m_statsMgr.get(); }
  const photonicsParamsDram& getParamsDram() const { assert(m_paramsDram); return *m_paramsDram; }
  photonicsPerfEnergyBase* getPerfEnergyModel();

  photonicsUtils::threadPool* getThreadPool() { return m_threadPool.get(); }

  // Resource allocation and deletion
  PhotonicsObjId photonicsAlloc(PhotonicsAllocEnum allocType, uint64_t numElements, PhotonicsDataType dataType);
  PhotonicsObjId photonicsAllocMat(uint64_t numElements, PhotonicsDataType dataType);
  PhotonicsObjId photonicsAllocAssociated(PhotonicsObjId assocId, PhotonicsDataType dataType);
  PhotonicsObjId photonicsAllocAssociatedSrcVec(PhotonicsObjId assocId, PhotonicsDataType dataType);
  PhotonicsObjId photonicsAllocAssociatedDestVec(PhotonicsObjId assocId, PhotonicsDataType dataType);
  PhotonicsObjId photonicsAllocBuffer(uint32_t numElements, PhotonicsDataType dataType);
  bool photonicsFreeMat(PhotonicsObjId obj);
  bool photonicsFreeSrcVec(PhotonicsObjId obj);
  bool photonicsFreeDestVec(PhotonicsObjId obj);
  PhotonicsObjId photonicsCreateRangedRef(PhotonicsObjId refId, uint64_t idxBegin, uint64_t idxEnd);
  PhotonicsObjId photonicsCreateDualContactRef(PhotonicsObjId refId);

  // Data transfer
  bool photonicsCopyMainToDevice(void* src, PhotonicsObjId dest, uint64_t idxBegin = 0, uint64_t idxEnd = 0);
  bool photonicsCopyDeviceToMain(PhotonicsObjId src, void* dest, uint64_t idxBegin = 0, uint64_t idxEnd = 0);
  bool photonicsCopyMainToDeviceWithType(PhotonicsCopyEnum copyType, void* src, PhotonicsObjId dest, uint64_t idxBegin = 0, uint64_t idxEnd = 0);
  bool photonicsCopyDeviceToMainWithType(PhotonicsCopyEnum copyType, PhotonicsObjId src, void* dest, uint64_t idxBegin = 0, uint64_t idxEnd = 0);
  bool photonicsCopyDeviceToDevice(PhotonicsObjId src, PhotonicsObjId dest, uint64_t idxBegin = 0, uint64_t idxEnd = 0);
  bool photonicsCopyObjectToObject(PhotonicsObjId src, PhotonicsObjId dest);
  bool photonicsConvertType(PhotonicsObjId src, PhotonicsObjId dest);

  // Computation
  bool photonicsMvm(PhotonicsObjId src1, PhotonicsObjId src2, PhotonicsObjId dest);
  bool photonicsIter(PhotonicsObjId src1, PhotonicsObjId src2, PhotonicsObjId dest, int8_t numLoops = 1);
  bool photonicsMmm(PhotonicsObjId src1, PhotonicsObjId src2, PhotonicsObjId dest);
  bool photonicsAdd(PhotonicsObjId src1, PhotonicsObjId src2, PhotonicsObjId dest);
  bool photonicsSub(PhotonicsObjId src1, PhotonicsObjId src2, PhotonicsObjId dest);
  bool photonicsDiv(PhotonicsObjId src1, PhotonicsObjId src2, PhotonicsObjId dest);
  bool photonicsAbs(PhotonicsObjId src, PhotonicsObjId dest);
  bool photonicsMul(PhotonicsObjId src1, PhotonicsObjId src2, PhotonicsObjId dest);
  bool photonicsNot(PhotonicsObjId src, PhotonicsObjId dest);
  bool photonicsOr(PhotonicsObjId src1, PhotonicsObjId src2, PhotonicsObjId dest);
  bool photonicsAnd(PhotonicsObjId src1, PhotonicsObjId src2, PhotonicsObjId dest);
  bool photonicsXor(PhotonicsObjId src1, PhotonicsObjId src2, PhotonicsObjId dest);
  bool photonicsXnor(PhotonicsObjId src1, PhotonicsObjId src2, PhotonicsObjId dest);
  bool photonicsGT(PhotonicsObjId src1, PhotonicsObjId src2, PhotonicsObjId dest);
  bool photonicsLT(PhotonicsObjId src1, PhotonicsObjId src2, PhotonicsObjId dest);
  bool photonicsEQ(PhotonicsObjId src1, PhotonicsObjId src2, PhotonicsObjId dest);
  bool photonicsNE(PhotonicsObjId src1, PhotonicsObjId src2, PhotonicsObjId dest);
  bool photonicsMin(PhotonicsObjId src1, PhotonicsObjId src2, PhotonicsObjId dest);
  bool photonicsMax(PhotonicsObjId src1, PhotonicsObjId src2, PhotonicsObjId dest);
  bool photonicsAdd(PhotonicsObjId src, PhotonicsObjId dest, uint64_t scalarValue);
  bool photonicsSub(PhotonicsObjId src, PhotonicsObjId dest, uint64_t scalarValue);
  bool photonicsMul(PhotonicsObjId src, PhotonicsObjId dest, uint64_t scalarValue);
  bool photonicsDiv(PhotonicsObjId src, PhotonicsObjId dest, uint64_t scalarValue);
  bool photonicsAnd(PhotonicsObjId src, PhotonicsObjId dest, uint64_t scalarValue);
  bool photonicsOr(PhotonicsObjId src, PhotonicsObjId dest, uint64_t scalarValue);
  bool photonicsXor(PhotonicsObjId src, PhotonicsObjId dest, uint64_t scalarValue);
  bool photonicsXnor(PhotonicsObjId src, PhotonicsObjId dest, uint64_t scalarValue);
  bool photonicsGT(PhotonicsObjId src, PhotonicsObjId dest, uint64_t scalarValue);
  bool photonicsLT(PhotonicsObjId src, PhotonicsObjId dest, uint64_t scalarValue);
  bool photonicsEQ(PhotonicsObjId src, PhotonicsObjId dest, uint64_t scalarValue);
  bool photonicsNE(PhotonicsObjId src, PhotonicsObjId dest, uint64_t scalarValue);
  bool photonicsMin(PhotonicsObjId src, PhotonicsObjId dest, uint64_t scalarValue);
  bool photonicsMax(PhotonicsObjId src, PhotonicsObjId dest, uint64_t scalarValue);
  bool photonicsScaledAdd(PhotonicsObjId src1, PhotonicsObjId src2, PhotonicsObjId dest, uint64_t scalarValue);
  bool photonicsPopCount(PhotonicsObjId src, PhotonicsObjId dest);
  bool photonicsRedSum(PhotonicsObjId src, void* sum, uint64_t idxBegin = 0, uint64_t idxEnd = 0);
  bool photonicsRedMin(PhotonicsObjId src, void* min, uint64_t idxBegin = 0, uint64_t idxEnd = 0);
  bool photonicsRedMax(PhotonicsObjId src, void* max, uint64_t idxBegin = 0, uint64_t idxEnd = 0);
  bool photonicsBitSliceExtract(PhotonicsObjId src, PhotonicsObjId destBool, unsigned bitIdx);
  bool photonicsBitSliceInsert(PhotonicsObjId srcBool, PhotonicsObjId dest, unsigned bitIdx);
  bool photonicsCondCopy(PhotonicsObjId condBool, PhotonicsObjId src, PhotonicsObjId dest);
  bool photonicsCondBroadcast(PhotonicsObjId condBool, uint64_t scalarBits, PhotonicsObjId dest);
  bool photonicsCondSelect(PhotonicsObjId condBool, PhotonicsObjId src1, PhotonicsObjId src2, PhotonicsObjId dest);
  bool photonicsCondSelectScalar(PhotonicsObjId condBool, PhotonicsObjId src1, uint64_t scalarBits, PhotonicsObjId dest);
  template <typename T> bool photonicsBroadcast(PhotonicsObjId dest, T value);
  bool photonicsRotateElementsRight(PhotonicsObjId src);
  bool photonicsRotateElementsLeft(PhotonicsObjId src);
  bool photonicsShiftElementsRight(PhotonicsObjId src);
  bool photonicsShiftElementsLeft(PhotonicsObjId src);
  bool photonicsShiftBitsRight(PhotonicsObjId src, PhotonicsObjId dest, unsigned shiftAmount);
  bool photonicsShiftBitsLeft(PhotonicsObjId src, PhotonicsObjId dest, unsigned shiftAmount);
  bool photonicsAesSbox(PhotonicsObjId src, PhotonicsObjId dest, const std::vector<uint8_t>& lut); 
  bool photonicsAesInverseSbox(PhotonicsObjId src, PhotonicsObjId dest, const std::vector<uint8_t>& lut); 
  bool photonicsPrefixSum(PhotonicsObjId src, PhotonicsObjId dest);
  bool photonicsMAC(PhotonicsObjId src1, PhotonicsObjId src2, void* dest);

  // PHOTONICS API Fusion
  bool photonicsFuse(PhotonicsProg prog);

  // BitSIMD-V micro ops
  bool photonicsOpReadRowToSa(PhotonicsObjId src, unsigned ofst);
  bool photonicsOpWriteSaToRow(PhotonicsObjId src, unsigned ofst);
  bool photonicsOpTRA(PhotonicsObjId src1, unsigned ofst1, PhotonicsObjId src2, unsigned ofst2, PhotonicsObjId src3, unsigned ofst3);
  bool photonicsOpMove(PhotonicsObjId objId, PhotonicsRowReg src, PhotonicsRowReg dest);
  bool photonicsOpSet(PhotonicsObjId objId, PhotonicsRowReg dest, bool val);
  bool photonicsOpNot(PhotonicsObjId objId, PhotonicsRowReg src, PhotonicsRowReg dest);
  bool photonicsOpAnd(PhotonicsObjId objId, PhotonicsRowReg src1, PhotonicsRowReg src2, PhotonicsRowReg dest);
  bool photonicsOpOr(PhotonicsObjId objId, PhotonicsRowReg src1, PhotonicsRowReg src2, PhotonicsRowReg dest);
  bool photonicsOpNand(PhotonicsObjId objId, PhotonicsRowReg src1, PhotonicsRowReg src2, PhotonicsRowReg dest);
  bool photonicsOpNor(PhotonicsObjId objId, PhotonicsRowReg src1, PhotonicsRowReg src2, PhotonicsRowReg dest);
  bool photonicsOpXor(PhotonicsObjId objId, PhotonicsRowReg src1, PhotonicsRowReg src2, PhotonicsRowReg dest);
  bool photonicsOpXnor(PhotonicsObjId objId, PhotonicsRowReg src1, PhotonicsRowReg src2, PhotonicsRowReg dest);
  bool photonicsOpMaj(PhotonicsObjId objId, PhotonicsRowReg src1, PhotonicsRowReg src2, PhotonicsRowReg src3, PhotonicsRowReg dest);
  bool photonicsOpSel(PhotonicsObjId objId, PhotonicsRowReg cond, PhotonicsRowReg src1, PhotonicsRowReg src2, PhotonicsRowReg dest);
  bool photonicsOpRotateRH(PhotonicsObjId objId, PhotonicsRowReg src);
  bool photonicsOpRotateLH(PhotonicsObjId objId, PhotonicsRowReg src);

  // SIMDRAM micro ops
  bool photonicsOpAP(int numSrc, va_list args);
  bool photonicsOpAAP(int numSrc, int numDest, va_list args);

private:
  photonicsSim();
  ~photonicsSim();
  photonicsSim(const photonicsSim&) = delete;
  photonicsSim operator=(const photonicsSim&) = delete;
  bool createDeviceCommon();
  void uninit();

  static photonicsSim* s_instance;
  photonicsSimConfig m_config;

  // support one device for now
  std::unique_ptr<photonicsDevice> m_device;
  std::unique_ptr<photonicsParamsDram> m_paramsDram;
  std::unique_ptr<photonicsStatsMgr> m_statsMgr;
  std::unique_ptr<photonicsUtils::threadPool> m_threadPool;

};

#endif

