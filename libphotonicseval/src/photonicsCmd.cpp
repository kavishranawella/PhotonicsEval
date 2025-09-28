// File: photonicsCmd.cpp
// PHOTONICSeval Simulator - PHOTONICS Commands
// Copyright (c) 2024 University of Virginia
// This file is licensed under the MIT License.
// See the LICENSE file in the root of this repository for more details.

#include "photonicsCmd.h"          // for photonicsCmd
#include "photonicsSim.h"          // for photonicsSim
#include "photonicsSimConfig.h"    // for photonicsSimConfig
#include "photonicsDevice.h"       // for photonicsDevice
#include "photonicsCore.h"         // for photonicsCore
#include "photonicsResMgr.h"       // for photonicsResMgr
#include "libphotonicseval.h"      // for PhotonicsObjId
#include <cstdio>
#include <cmath>
#include <unordered_map>
#include <unordered_set>
#include <climits>
#include <cinttypes>         // for PRIu64, PRIx64

//! @brief  Get PHOTONICS command name from command type enum
std::string
photonicsCmd::getName(PhotonicsCmdEnum cmdType, const std::string& suffix)
{
  static const std::unordered_map<PhotonicsCmdEnum, std::string> cmdNames = {
    { PhotonicsCmdEnum::NOOP, "noop" },
    { PhotonicsCmdEnum::COPY_H2D, "copy_h2d" },
    { PhotonicsCmdEnum::COPY_D2H, "copy_d2h" },
    { PhotonicsCmdEnum::COPY_D2D, "copy_d2d" },
    { PhotonicsCmdEnum::COPY_O2O, "copy_o2o" },
    { PhotonicsCmdEnum::ABS, "abs" },
    { PhotonicsCmdEnum::POPCOUNT, "popcount" },
    { PhotonicsCmdEnum::SHIFT_BITS_R, "shift_bits_r" },
    { PhotonicsCmdEnum::SHIFT_BITS_L, "shift_bits_l" },
    { PhotonicsCmdEnum::BROADCAST, "broadcast" },
    { PhotonicsCmdEnum::ADD, "add" },
    { PhotonicsCmdEnum::SUB, "sub" },
    { PhotonicsCmdEnum::MUL, "mul" },
    { PhotonicsCmdEnum::SCALED_ADD, "scaled_add" },
    { PhotonicsCmdEnum::DIV, "div" },
    { PhotonicsCmdEnum::NOT, "not" },
    { PhotonicsCmdEnum::AND, "and" },
    { PhotonicsCmdEnum::OR, "or" },
    { PhotonicsCmdEnum::XOR, "xor" },
    { PhotonicsCmdEnum::XNOR, "xnor" },
    { PhotonicsCmdEnum::GT, "gt" },
    { PhotonicsCmdEnum::LT, "lt" },
    { PhotonicsCmdEnum::EQ, "eq" },
    { PhotonicsCmdEnum::NE, "ne" },
    { PhotonicsCmdEnum::MIN, "min" },
    { PhotonicsCmdEnum::MAX, "max" },
    { PhotonicsCmdEnum::ADD_SCALAR, "add_scalar" },
    { PhotonicsCmdEnum::SUB_SCALAR, "sub_scalar" },
    { PhotonicsCmdEnum::MUL_SCALAR, "mul_scalar" },
    { PhotonicsCmdEnum::DIV_SCALAR, "div_scalar" },
    { PhotonicsCmdEnum::AND_SCALAR, "and_scalar" },
    { PhotonicsCmdEnum::OR_SCALAR, "or_scalar" },
    { PhotonicsCmdEnum::XOR_SCALAR, "xor_scalar" },
    { PhotonicsCmdEnum::XNOR_SCALAR, "xnor_scalar" },
    { PhotonicsCmdEnum::GT_SCALAR, "gt_scalar" },
    { PhotonicsCmdEnum::LT_SCALAR, "lt_scalar" },
    { PhotonicsCmdEnum::EQ_SCALAR, "eq_scalar" },
    { PhotonicsCmdEnum::NE_SCALAR, "ne_scalar" },
    { PhotonicsCmdEnum::MIN_SCALAR, "min_scalar" },
    { PhotonicsCmdEnum::MAX_SCALAR, "max_scalar" },
    { PhotonicsCmdEnum::CONVERT_TYPE, "convert_type" },
    { PhotonicsCmdEnum::BIT_SLICE_EXTRACT, "bit_slice_extract" },
    { PhotonicsCmdEnum::BIT_SLICE_INSERT, "bit_slice_insert" },
    { PhotonicsCmdEnum::COND_COPY, "cond_copy" },
    { PhotonicsCmdEnum::COND_BROADCAST, "cond_broadcast" },
    { PhotonicsCmdEnum::COND_SELECT, "cond_select" },
    { PhotonicsCmdEnum::COND_SELECT_SCALAR, "cond_select_scalar" },
    { PhotonicsCmdEnum::AES_SBOX, "aes_sbox" },
    { PhotonicsCmdEnum::AES_INVERSE_SBOX, "aes_inverse_sbox" },
    { PhotonicsCmdEnum::PREFIX_SUM, "prefix_sum"},
    { PhotonicsCmdEnum::REDSUM, "redsum" },
    { PhotonicsCmdEnum::REDSUM_RANGE, "redsum_range" },
    { PhotonicsCmdEnum::REDMIN, "redmin" },
    { PhotonicsCmdEnum::REDMIN_RANGE, "redmin_range" },
    { PhotonicsCmdEnum::REDMAX, "redmax" },
    { PhotonicsCmdEnum::REDMAX_RANGE, "redmax_range" },
    { PhotonicsCmdEnum::MAC, "mac" },
    { PhotonicsCmdEnum::ROTATE_ELEM_R, "rotate_elem_r" },
    { PhotonicsCmdEnum::ROTATE_ELEM_L, "rotate_elem_l" },
    { PhotonicsCmdEnum::SHIFT_ELEM_R, "shift_elem_r" },
    { PhotonicsCmdEnum::SHIFT_ELEM_L, "shift_elem_l" },
    { PhotonicsCmdEnum::ROW_R, "row_r" },
    { PhotonicsCmdEnum::ROW_W, "row_w" },
    { PhotonicsCmdEnum::RREG_MOV, "rreg.mov" },
    { PhotonicsCmdEnum::RREG_SET, "rreg.set" },
    { PhotonicsCmdEnum::RREG_NOT, "rreg.not" },
    { PhotonicsCmdEnum::RREG_AND, "rreg.and" },
    { PhotonicsCmdEnum::RREG_OR, "rreg.or" },
    { PhotonicsCmdEnum::RREG_NAND, "rreg.nand" },
    { PhotonicsCmdEnum::RREG_NOR, "rreg.nor" },
    { PhotonicsCmdEnum::RREG_XOR, "rreg.xor" },
    { PhotonicsCmdEnum::RREG_XNOR, "rreg.xnor" },
    { PhotonicsCmdEnum::RREG_MAJ, "rreg.maj" },
    { PhotonicsCmdEnum::RREG_SEL, "rreg.sel" },
    { PhotonicsCmdEnum::RREG_ROTATE_R, "rreg.rotate_r" },
    { PhotonicsCmdEnum::RREG_ROTATE_L, "rreg.rotate_l" },
    { PhotonicsCmdEnum::ROW_AP, "row_ap" },
    { PhotonicsCmdEnum::ROW_AAP, "row_aap" },
  };
  auto it = cmdNames.find(cmdType);
  return it != cmdNames.end() ? it->second + suffix : "unknown";
}

//! @brief  photonicsCmd constructor
photonicsCmd::photonicsCmd(PhotonicsCmdEnum cmdType)
  : m_cmdType(cmdType)
{
  m_debugCmds = photonicsSim::get()->isDebug(photonicsSimConfig::DEBUG_CMDS);
}

//! @brief  Check if an obj ID is valid
bool
photonicsCmd::isValidObjId(photonicsResMgr* resMgr, PhotonicsObjId objId) const
{
  if (!resMgr->isValidObjId(objId)) {
    std::printf("PHOTONICS-Error: Invalid object id %d\n", objId);
    return false;
  }
  return true;
}

//! @brief  Check if two objects are associated
bool
photonicsCmd::isAssociated(const photonicsObjInfo& obj1, const photonicsObjInfo& obj2) const
{
  if (obj1.getAssocObjId() != obj2.getAssocObjId()) {
    std::printf("PHOTONICS-Error: Object id %d and %d are not associated\n", obj1.getObjId(), obj2.getObjId());
    return false;
  }
  return true;
}

//! @brief  Check if two objects have compatible type.
bool
photonicsCmd::isCompatibleType(const photonicsObjInfo& obj1, const photonicsObjInfo& obj2) const
{
  // TODO: Type conversion eg. 32-bit and 64-bit should be compatible
  if (obj1.getDataType() != obj2.getDataType()) {
    std::printf("PHOTONICS-Error: Type mismatch between object %d and %d\n", obj1.getObjId(), obj2.getObjId());
    return false;
  }
  return true;
}

//! @brief  Check if src type can be converted to dest type.
bool
photonicsCmd::isConvertibleType(const photonicsObjInfo& src, const photonicsObjInfo& dest) const
{
  // TODO: Type conversion
  if (src.getDataType() != dest.getDataType()) {
    std::printf("PHOTONICS-Error: Cannot convert from %s to %s\n",
        photonicsUtils::photonicsDataTypeEnumToStr(src.getDataType()).c_str(),
        photonicsUtils::photonicsDataTypeEnumToStr(dest.getDataType()).c_str());
    return false;
  }
  return true;
}

//! @brief  Process all regions in MT used by derived classes
bool
photonicsCmd::computeAllRegions(unsigned numRegions)
{
  // skip PHOTONICS computation in analysis mode
  if (photonicsSim::get()->isAnalysisMode()) {
    return true;
  }
  if (photonicsSim::get()->getNumThreads() > 1) { // MT
    std::vector<photonicsUtils::threadWorker*> workers;
    for (unsigned i = 0; i < numRegions; ++i) {
      workers.push_back(new regionWorker(this, i));
    }
    photonicsSim::get()->getThreadPool()->doWork(workers);
    for (unsigned i = 0; i < numRegions; ++i) {
      delete workers[i];
    }
  } else { // single thread
    for (unsigned i = 0; i < numRegions; ++i) {
      computeRegion(i);
    }
  }
  return true;
}


//! @brief  PHOTONICS Data Copy
bool
photonicsCmdCopy::execute()
{
  if (!sanityCheck()) {
    return false;
  }

  // for non-functional simulation, sync src data from simulated memory
  if (photonicsSim::get()->getDeviceType() != PHOTONICS_FUNCTIONAL) {
    if (m_cmdType == PhotonicsCmdEnum::COPY_D2H || m_cmdType == PhotonicsCmdEnum::COPY_D2D) {
      photonicsObjInfo &objSrc = m_device->getResMgr()->getObjInfo(m_src);
      objSrc.syncFromSimulatedMem();
    }
  }

  if (!photonicsSim::get()->isAnalysisMode()) {
    if (m_cmdType == PhotonicsCmdEnum::COPY_H2D) {
      photonicsObjInfo &objDest = m_device->getResMgr()->getObjInfo(m_dest);
      objDest.copyFromHost(m_ptr, m_idxBegin, m_idxEnd);
    } else if (m_cmdType == PhotonicsCmdEnum::COPY_D2H) {
      const photonicsObjInfo &objSrc = m_device->getResMgr()->getObjInfo(m_src);
      objSrc.copyToHost(m_ptr, m_idxBegin, m_idxEnd);
    } else if (m_cmdType == PhotonicsCmdEnum::COPY_D2D) {
      const photonicsObjInfo &objSrc = m_device->getResMgr()->getObjInfo(m_src);
      photonicsObjInfo &objDest = m_device->getResMgr()->getObjInfo(m_dest);
      objSrc.copyToObj(objDest, m_idxBegin, m_idxEnd);
    } else {
      assert(0);
    }
  }

  // for non-functional simulation, sync dest data to simulated memory
  if (photonicsSim::get()->getDeviceType() != PHOTONICS_FUNCTIONAL) {
    if (m_cmdType == PhotonicsCmdEnum::COPY_H2D || m_cmdType == PhotonicsCmdEnum::COPY_D2D) {
      const photonicsObjInfo &objDest = m_device->getResMgr()->getObjInfo(m_dest);
      objDest.syncToSimulatedMem();
    }
  }

  updateStats();
  return true;
}

//! @brief  PHOTONICS Data Copy - sanity check
bool
photonicsCmdCopy::sanityCheck() const
{
  photonicsResMgr* resMgr = m_device->getResMgr();
  uint64_t numElements = 0;
  switch (m_cmdType) {
  case PhotonicsCmdEnum::COPY_H2D:
  {
    if (!m_ptr) {
      std::printf("PHOTONICS-Error: Invalid null pointer as copy source\n");
      return false;
    }
    if (!resMgr->isValidObjId(m_dest)) {
      std::printf("PHOTONICS-Error: Invalid PHOTONICS object ID %d as copy destination\n", m_dest);
      return false;
    }
    const photonicsObjInfo &objDest = m_device->getResMgr()->getObjInfo(m_dest);
    numElements = objDest.getNumElements();
    break;
  }
  case PhotonicsCmdEnum::COPY_D2H:
  {
    if (!resMgr->isValidObjId(m_src)) {
      std::printf("PHOTONICS-Error: Invalid PHOTONICS object ID %d as copy source\n", m_src);
      return false;
    }
    if (!m_ptr) {
      std::printf("PHOTONICS-Error: Invalid null pointer as copy destination\n");
      return false;
    }
    const photonicsObjInfo &objSrc = m_device->getResMgr()->getObjInfo(m_src);
    numElements = objSrc.getNumElements();
    break;
  }
  case PhotonicsCmdEnum::COPY_D2D:
  {
    if (!resMgr->isValidObjId(m_src)) {
      std::printf("PHOTONICS-Error: Invalid PHOTONICS object ID %d as copy source\n", m_src);
      return false;
    }
    if (!resMgr->isValidObjId(m_dest)) {
      std::printf("PHOTONICS-Error: Invalid PHOTONICS object ID %d as copy destination\n", m_dest);
      return false;
    }
    const photonicsObjInfo &objSrc = resMgr->getObjInfo(m_src);
    const photonicsObjInfo &objDest = resMgr->getObjInfo(m_dest);
    if (!isAssociated(objSrc, objDest)) {
      std::printf("PHOTONICS-Error: PHOTONICS object IDs %d and %d are not associated for device-to-device copying\n", m_src, m_dest);
      return false;
    }
    numElements = objSrc.getNumElements();
    break;
  }
  default:
    assert(0);
  }
  if (!m_copyFullRange) {
    if (m_idxBegin > numElements) {
      std::printf("PHOTONICS-Error: The beginning of the copy range for PHOTONICS object ID %d is greater than the number of elements\n", m_dest);
      return false;
    }
    if (m_idxEnd > numElements) {
      std::printf("PHOTONICS-Error: The end of the copy range for PHOTONICS object ID %d is greater than the number of elements\n", m_dest);
      return false;
    }
    if (m_idxEnd < m_idxBegin) {
      std::printf("PHOTONICS-Error: The end of the copy range for PHOTONICS object ID %d is less than its beginning\n", m_dest);
      return false;
    }
  }

  return true;
}

//! @brief  PHOTONICS Data Copy - update stats
bool
photonicsCmdCopy::updateStats() const
{
   if (m_cmdType == PhotonicsCmdEnum::COPY_H2D) {
    const photonicsObjInfo &objDest = m_device->getResMgr()->getObjInfo(m_dest);
    uint64_t numElements = objDest.getNumElements();
    if (!m_copyFullRange) {
      numElements = m_idxEnd - m_idxBegin;
    }
    unsigned bitsPerElement = 8;
    photonicseval::perfEnergy mPerfEnergy = photonicsSim::get()->getPerfEnergyModel()->getPerfEnergyForBytesTransfer(m_cmdType, numElements * bitsPerElement / 8);
    photonicsSim::get()->getStatsMgr()->recordCopyMainToDevice(numElements * bitsPerElement, mPerfEnergy);

    if (m_debugCmds) {
      std::printf("PHOTONICS-Cmd: Copied %" PRIu64 " elements of %u bits from host to PHOTONICS obj %d\n",
                  numElements, bitsPerElement, m_dest);
    }
  } else if (m_cmdType == PhotonicsCmdEnum::COPY_D2H) {
    const photonicsObjInfo &objSrc = m_device->getResMgr()->getObjInfo(m_src);
    uint64_t numElements = 8;
    if (!m_copyFullRange) {
      numElements = m_idxEnd - m_idxBegin;
    }
    unsigned bitsPerElement = objSrc.getBitsPerElement(PhotonicsBitWidth::ACTUAL);
    photonicseval::perfEnergy mPerfEnergy = photonicsSim::get()->getPerfEnergyModel()->getPerfEnergyForBytesTransfer(m_cmdType, numElements * bitsPerElement / 8);
    photonicsSim::get()->getStatsMgr()->recordCopyDeviceToMain(numElements * bitsPerElement, mPerfEnergy);

    if (m_debugCmds) {
      std::printf("PHOTONICS-Cmd: Copied %" PRIu64 " elements of %u bits from PHOTONICS obj %d to host\n",
                  numElements, bitsPerElement, m_src);
    }
  } else if (m_cmdType == PhotonicsCmdEnum::COPY_D2D) {
    const photonicsObjInfo &objSrc = m_device->getResMgr()->getObjInfo(m_src);
    uint64_t numElements = objSrc.getNumElements();
    if (!m_copyFullRange) {
      numElements = m_idxEnd - m_idxBegin;
    }
    unsigned bitsPerElement = objSrc.getBitsPerElement(PhotonicsBitWidth::ACTUAL);
    photonicseval::perfEnergy mPerfEnergy = photonicsSim::get()->getPerfEnergyModel()->getPerfEnergyForBytesTransfer(m_cmdType, numElements * bitsPerElement / 8);
    photonicsSim::get()->getStatsMgr()->recordCopyDeviceToDevice(numElements * bitsPerElement, mPerfEnergy);

    if (m_debugCmds) {
      std::printf("PHOTONICS-Cmd: Copied %" PRIu64 " elements of %u bits from PHOTONICS obj %d to PHOTONICS obj %d\n",
                  numElements, bitsPerElement, m_src, m_dest);
    }
  } else {
    assert(0);
  }
  return true;
}


//! @brief  PHOTONICS CMD: Functional 1-operand
bool
photonicsCmdFunc1::execute()
{
  if (m_debugCmds) {
    std::printf("PHOTONICS-Cmd: %s (obj id %d -> %d)\n", getName().c_str(), m_src, m_dest);
  }

  if (!sanityCheck()) {
    return false;
  }

  if (photonicsSim::get()->getDeviceType() != PHOTONICS_FUNCTIONAL) {
    photonicsObjInfo &objSrc = m_device->getResMgr()->getObjInfo(m_src);
    objSrc.syncFromSimulatedMem();
    if (m_cmdType == PhotonicsCmdEnum::BIT_SLICE_INSERT) {  // require dest data to be synced
      photonicsObjInfo &objDest = m_device->getResMgr()->getObjInfo(m_dest);
      objDest.syncFromSimulatedMem();
    }
  }

  const photonicsObjInfo& objSrc = m_device->getResMgr()->getObjInfo(m_src);
  unsigned numRegions = objSrc.getRegions().size();
  computeAllRegions(numRegions);

  if (photonicsSim::get()->getDeviceType() != PHOTONICS_FUNCTIONAL) {
    const photonicsObjInfo &objDest = m_device->getResMgr()->getObjInfo(m_dest);
    objDest.syncToSimulatedMem();
  }

  updateStats();
  return true;
}

//! @brief  PHOTONICS CMD: Functional 1-operand - sanity check
bool
photonicsCmdFunc1::sanityCheck() const
{
  photonicsResMgr* resMgr = m_device->getResMgr();
  if (!isValidObjId(resMgr, m_src) || !isValidObjId(resMgr, m_dest)) {
    return false;
  }
  const photonicsObjInfo& objSrc = resMgr->getObjInfo(m_src);
  const photonicsObjInfo& objDest = resMgr->getObjInfo(m_dest);
  if (!isAssociated(objSrc, objDest)) {
    return false;
  }
  if (objSrc.getDataType() == PHOTONICS_BOOL) {
    switch (m_cmdType) {
      case PhotonicsCmdEnum::NOT:
      case PhotonicsCmdEnum::CONVERT_TYPE:
      case PhotonicsCmdEnum::BIT_SLICE_EXTRACT:
      case PhotonicsCmdEnum::BIT_SLICE_INSERT:
      case PhotonicsCmdEnum::COPY_O2O:
        break;
      default:
        std::printf("PHOTONICS-Error: PHOTONICS command %s does not support PHOTONICS_BOOL type\n", getName().c_str());
        return false;
    }
  }
  // Define command specific data type rules
  switch (m_cmdType) {
    case PhotonicsCmdEnum::CONVERT_TYPE:
      break;
    case PhotonicsCmdEnum::GT_SCALAR:
    case PhotonicsCmdEnum::LT_SCALAR:
    case PhotonicsCmdEnum::EQ_SCALAR:
    case PhotonicsCmdEnum::NE_SCALAR:
      if (objDest.getDataType() != PHOTONICS_BOOL) {
        std::printf("PHOTONICS-Error: PHOTONICS command %s destination operand must be PHOTONICS_BOOL type\n", getName().c_str());
        return false;
      }
      break;
    case PhotonicsCmdEnum::BIT_SLICE_EXTRACT: // src, destBool, bitIdx
      if (objDest.getDataType() != PHOTONICS_BOOL) {
        std::printf("PHOTONICS-Error: PHOTONICS command %s destination operand must be PHOTONICS_BOOL type\n", getName().c_str());
        return false;
      }
      if (m_scalarValue >= objSrc.getBitsPerElement(PhotonicsBitWidth::SIM)) {
        std::printf("PHOTONICS-Error: PHOTONICS command %s bit index %" PRIu64 " out of range of %s type\n", getName().c_str(),
                    m_scalarValue, photonicsUtils::photonicsDataTypeEnumToStr(objSrc.getDataType()).c_str());
        return false;
      }
      break;
    case PhotonicsCmdEnum::BIT_SLICE_INSERT: // srcBool, dest, bitIdx
      if (objSrc.getDataType() != PHOTONICS_BOOL) {
        std::printf("PHOTONICS-Error: PHOTONICS command %s source operand must be PHOTONICS_BOOL type\n", getName().c_str());
        return false;
      }
      if (m_scalarValue >= objDest.getBitsPerElement(PhotonicsBitWidth::SIM)) {
        std::printf("PHOTONICS-Error: PHOTONICS command %s bit index %" PRIu64 " out of range of %s type\n", getName().c_str(),
                    m_scalarValue, photonicsUtils::photonicsDataTypeEnumToStr(objDest.getDataType()).c_str());
        return false;
      }
      break;
    case PhotonicsCmdEnum::AES_SBOX:
    case PhotonicsCmdEnum::AES_INVERSE_SBOX:
      if (objSrc.getDataType() != PHOTONICS_UINT8) {
        return false;
      }
      if (objDest.getDataType() != PHOTONICS_UINT8) {
        return false;
      }
      if (m_lut.size() != 256) {
        return false;
      }
      break;
    default:
      if (objSrc.getDataType() != objDest.getDataType()) {
        std::printf("PHOTONICS-Error: PHOTONICS command %s does not support mixed data type\n", getName().c_str());
        return false;
      }
  }
  return true;
}

//! @brief  PHOTONICS CMD: Functional 1-operand - compute region
bool
photonicsCmdFunc1::computeRegion(unsigned index)
{
  const photonicsObjInfo& objSrc = m_device->getResMgr()->getObjInfo(m_src);
  photonicsObjInfo& objDest = m_device->getResMgr()->getObjInfo(m_dest);

  PhotonicsDataType dataType = objSrc.getDataType();
  unsigned bitsPerElementSrc = objSrc.getBitsPerElement(PhotonicsBitWidth::SIM);
  const photonicsRegion& srcRegion = objSrc.getRegions()[index];

  // perform the computation
  uint64_t elemIdxBegin = srcRegion.getElemIdxBegin();
  unsigned numElementsInRegion = srcRegion.getNumElemInRegion();
  for (unsigned j = 0; j < numElementsInRegion; ++j) {
    uint64_t elemIdx = elemIdxBegin + j;
    if (m_cmdType == PhotonicsCmdEnum::CONVERT_TYPE) {
      convertType(objSrc, objDest, elemIdx);
      continue;
    } else if (m_cmdType == PhotonicsCmdEnum::BIT_SLICE_EXTRACT) {
      bitSliceExtract(objSrc, objDest, m_scalarValue, elemIdx);
      continue;
    } else if (m_cmdType == PhotonicsCmdEnum::BIT_SLICE_INSERT) {
      bitSliceInsert(objSrc, objDest, m_scalarValue, elemIdx);
      continue;
    }
    if (photonicsUtils::isSigned(dataType)) {
      int64_t signedOperand = objSrc.getElementBits(elemIdx);
      int64_t result = 0;
      if(!computeResult(signedOperand, m_cmdType, (int64_t)m_scalarValue, result, bitsPerElementSrc)) return false;
      objDest.setElement(elemIdx, result);
    } else if (photonicsUtils::isUnsigned(dataType)) {
      uint64_t unsignedOperand = objSrc.getElementBits(elemIdx);
      uint64_t result = 0;
      if(!computeResult(unsignedOperand, m_cmdType, m_scalarValue, result, bitsPerElementSrc)) return false;
      objDest.setElement(elemIdx, result);
    } else if (photonicsUtils::isFP(dataType)) {
      uint64_t bits = objSrc.getElementBits(elemIdx);
      float floatOperand = photonicsUtils::castBitsToType<float>(bits);
      float result = 0.0;
      if(!computeResultFP(floatOperand, m_cmdType, photonicsUtils::castBitsToType<float>(m_scalarValue), result)) return false;
      if (objDest.getDataType() == PHOTONICS_BOOL) {
        bool resultBool = result > 0;
        objDest.setElement(elemIdx, resultBool);
      } else {
        objDest.setElement(elemIdx, result);
      }
    } else {
      assert(0); // todo: data type
    }
  }
  return true;
}

//! @brief  PHOTONICS CMD: Functional 1-operand - compute region - convert data type
bool
photonicsCmdFunc1::convertType(const photonicsObjInfo& objSrc, photonicsObjInfo& objDest, uint64_t elemIdx) const
{
  PhotonicsDataType dataTypeSrc = objSrc.getDataType();
  PhotonicsDataType dataTypeDest = objDest.getDataType();
  if (photonicsUtils::isSigned(dataTypeSrc)) {
    int64_t signedVal = objSrc.getElementBits(elemIdx);
    if (photonicsUtils::isSigned(dataTypeDest)) {
      int64_t result = signedVal;
      objDest.setElement(elemIdx, result);
    } else if (photonicsUtils::isUnsigned(dataTypeDest)) {
      uint64_t result = static_cast<uint64_t>(signedVal);
      objDest.setElement(elemIdx, result);
    } else if (photonicsUtils::isFP(dataTypeDest)) {
      assert(0); // todo
    }
  } else if (photonicsUtils::isUnsigned(dataTypeSrc)) {
    uint64_t unsignedVal = objSrc.getElementBits(elemIdx);
    if (photonicsUtils::isSigned(dataTypeDest)) {
      int64_t result = static_cast<int64_t>(unsignedVal);
      objDest.setElement(elemIdx, result);
    } else if (photonicsUtils::isUnsigned(dataTypeDest)) {
      uint64_t result = unsignedVal;
      objDest.setElement(elemIdx, result);
    } else if (photonicsUtils::isFP(dataTypeDest)) {
      assert(0); // todo
    }
  } else if (photonicsUtils::isFP(dataTypeSrc)) {
    assert(0); // todo
  }
  return true;
}

//! @brief  PHOTONICS CMD: Functional 1-operand - compute region - bit slice extract
bool
photonicsCmdFunc1::bitSliceExtract(const photonicsObjInfo& objSrc, photonicsObjInfo& objDestBool, uint64_t bitIdx, uint64_t elemIdx) const
{
  uint64_t src = objSrc.getElementBits(elemIdx);
  uint64_t result = (src >> bitIdx) & 1L;
  objDestBool.setElement(elemIdx, result);
  return true;
}

//! @brief  PHOTONICS CMD: Functional 1-operand - compute region - bit slice insert
bool
photonicsCmdFunc1::bitSliceInsert(const photonicsObjInfo& objSrcBool, photonicsObjInfo& objDest, uint64_t bitIdx, uint64_t elemIdx) const
{
  uint64_t src = objSrcBool.getElementBits(elemIdx);
  uint64_t dest = objDest.getElementBits(elemIdx);
  uint64_t result = (dest & ~(1L << bitIdx)) | (src << bitIdx);
  objDest.setElement(elemIdx, result);
  return true;
}

//! @brief  PHOTONICS CMD: Functional 1-operand - update stats
bool
photonicsCmdFunc1::updateStats() const
{
  // Special handling: Use dest for performance energy calculation of bit-slice insert
  bool useDestAsSrc = (m_cmdType == PhotonicsCmdEnum::BIT_SLICE_INSERT);
  const photonicsObjInfo& objSrc = (useDestAsSrc? m_device->getResMgr()->getObjInfo(m_dest) : m_device->getResMgr()->getObjInfo(m_src));
  const photonicsObjInfo& objDest = m_device->getResMgr()->getObjInfo(m_dest);
  PhotonicsDataType dataType = objSrc.getDataType();
  bool isVLayout = objSrc.isVLayout();

  photonicseval::perfEnergy mPerfEnergy = photonicsSim::get()->getPerfEnergyModel()->getPerfEnergyForFunc1(m_cmdType, objSrc, objDest);
  photonicsSim::get()->getStatsMgr()->recordCmd(getName(dataType, isVLayout), mPerfEnergy);
  return true;
}

//! @brief  PHOTONICS CMD: Functional 2-operand
bool
photonicsCmdFunc2::execute()
{
  if (m_debugCmds) {
    std::printf("PHOTONICS-Cmd: %s (obj id %d - %d -> %d)\n", getName().c_str(), m_src1, m_src2, m_dest);
  }

  if (!sanityCheck()) {
    return false;
  }

  if (photonicsSim::get()->getDeviceType() != PHOTONICS_FUNCTIONAL) {
    photonicsObjInfo &objSrc1 = m_device->getResMgr()->getObjInfo(m_src1);
    photonicsObjInfo &objSrc2 = m_device->getResMgr()->getObjInfo(m_src2);
    objSrc1.syncFromSimulatedMem();
    objSrc2.syncFromSimulatedMem();
  }

  const photonicsObjInfo& objSrc1 = m_device->getResMgr()->getObjInfo(m_src1);
  unsigned numRegions = objSrc1.getRegions().size();
  computeAllRegions(numRegions);

  if (photonicsSim::get()->getDeviceType() != PHOTONICS_FUNCTIONAL) {
    const photonicsObjInfo &objDest = m_device->getResMgr()->getObjInfo(m_dest);
    objDest.syncToSimulatedMem();
  }

  updateStats();
  return true;
}

//! @brief  PHOTONICS CMD: Functional 2-operand - sanity check
bool
photonicsCmdFunc2::sanityCheck() const
{
  photonicsResMgr* resMgr = m_device->getResMgr();
  if (!isValidObjId(resMgr, m_src1) || !isValidObjId(resMgr, m_src2) || !isValidObjId(resMgr, m_dest)) {
    return false;
  }
  const photonicsObjInfo& objSrc1 = resMgr->getObjInfo(m_src1);
  const photonicsObjInfo& objSrc2 = resMgr->getObjInfo(m_src2);
  const photonicsObjInfo& objDest = resMgr->getObjInfo(m_dest);
  if (!isAssociated(objSrc1, objSrc2) || !isAssociated(objSrc1, objDest)) {
    return false;
  }
  // Define command specific data type rules
  bool isBoolSrc1Allowed = false;
  bool isBoolSrc2Allowed = false;
  bool isBoolDestRequired = false;
  bool isSrc1Src2SameType = true;
  bool isSrc1DestSameType = true;
  switch (m_cmdType) {
    case PhotonicsCmdEnum::AND:
    case PhotonicsCmdEnum::OR:
    case PhotonicsCmdEnum::XOR:
    case PhotonicsCmdEnum::XNOR:
      isBoolSrc1Allowed = true;
      isBoolSrc2Allowed = true;
      break;
    case PhotonicsCmdEnum::GT:
    case PhotonicsCmdEnum::LT:
    case PhotonicsCmdEnum::EQ:
    case PhotonicsCmdEnum::NE:
      isBoolDestRequired = true;
      isSrc1DestSameType = false;
      break;
    case PhotonicsCmdEnum::ADD:
    case PhotonicsCmdEnum::SUB:
      isBoolSrc2Allowed = true;
      isSrc1Src2SameType = false;
      if (m_cmdType == PhotonicsCmdEnum::ADD && objSrc2.getDataType() == PHOTONICS_BOOL) {
        isBoolSrc1Allowed = true;  // support photonicsAdd bool + bool = int
        isSrc1DestSameType = false;
      }
      // extra checks
      if ((photonicsUtils::isFP(objSrc1.getDataType()) || photonicsUtils::isFP(objDest.getDataType())) && objSrc2.getDataType() == PHOTONICS_BOOL) {
        std::printf("PHOTONICS-Error: PHOTONICS command %s does not support mixed FP and PHOTONICS_BOOL types\n", getName().c_str());
        return false;
      }
      if (objSrc1.getDataType() != objSrc2.getDataType() && objSrc2.getDataType() != PHOTONICS_BOOL) {
        std::printf("PHOTONICS-Error: PHOTONICS command %s can only support mixed data types if src2 is of PHOTONICS_BOOL type\n", getName().c_str());
        return false;
      }
      break;
    default:
      ; // pass
  }
  if (!isBoolSrc1Allowed && objSrc1.getDataType() == PHOTONICS_BOOL) {
    std::printf("PHOTONICS-Error: PHOTONICS command %s src1 cannot be of PHOTONICS_BOOL type\n", getName().c_str());
    return false;
  }
  if (!isBoolSrc2Allowed && objSrc2.getDataType() == PHOTONICS_BOOL) {
    std::printf("PHOTONICS-Error: PHOTONICS command %s src2 cannot be of PHOTONICS_BOOL type\n", getName().c_str());
    return false;
  }
  if (isBoolDestRequired && objDest.getDataType() != PHOTONICS_BOOL) {
    std::printf("PHOTONICS-Error: PHOTONICS command %s dest must be of PHOTONICS_BOOL type\n", getName().c_str());
    return false;
  }
  if (isSrc1Src2SameType && objSrc1.getDataType() != objSrc2.getDataType()) {
    std::printf("PHOTONICS-Error: PHOTONICS command %s src1 and src2 must be of same data type\n", getName().c_str());
    return false;
  }
  if (isSrc1DestSameType && objSrc1.getDataType() != objDest.getDataType()) {
    std::printf("PHOTONICS-Error: PHOTONICS command %s src1 and dest must be of same data type\n", getName().c_str());
    return false;
  }

  return true;
}

//! @brief  PHOTONICS CMD: Functional 2-operand - compute region
bool
photonicsCmdFunc2::computeRegion(unsigned index)
{
  const photonicsObjInfo& objSrc1 = m_device->getResMgr()->getObjInfo(m_src1);
  const photonicsObjInfo& objSrc2 = m_device->getResMgr()->getObjInfo(m_src2);
  photonicsObjInfo& objDest = m_device->getResMgr()->getObjInfo(m_dest);

  PhotonicsDataType dataType = objSrc1.getDataType();

  const photonicsRegion& src1Region = objSrc1.getRegions()[index];

  // perform the computation
  uint64_t elemIdxBegin = src1Region.getElemIdxBegin();
  unsigned numElementsInRegion = src1Region.getNumElemInRegion();
  for (unsigned j = 0; j < numElementsInRegion; ++j) {
    uint64_t elemIdx = elemIdxBegin + j;
    if (photonicsUtils::isSigned(dataType)) {
      uint64_t operandBits1 = objSrc1.getElementBits(elemIdx);
      uint64_t operandBits2 = objSrc2.getElementBits(elemIdx);
      int64_t operand1 = photonicsUtils::signExt(operandBits1, dataType);
      int64_t operand2 = photonicsUtils::signExt(operandBits2, dataType);
      int64_t result = 0;
      if(!computeResult(operand1, operand2, m_cmdType, (int64_t)m_scalarValue, result)) return false;
      objDest.setElement(elemIdx, result);
    } else if (photonicsUtils::isUnsigned(dataType)) {
      uint64_t unsignedOperand1 = objSrc1.getElementBits(elemIdx);
      uint64_t unsignedOperand2 = objSrc2.getElementBits(elemIdx);
      uint64_t result = 0;
      if(!computeResult(unsignedOperand1, unsignedOperand2, m_cmdType, m_scalarValue, result)) return false;
      objDest.setElement(elemIdx, result);
    } else if (photonicsUtils::isFP(dataType)) {
      uint64_t operandBits1 = objSrc1.getElementBits(elemIdx);
      uint64_t operandBits2 = objSrc2.getElementBits(elemIdx);
      float floatOperand1 = photonicsUtils::castBitsToType<float>(operandBits1);
      float floatOperand2 = photonicsUtils::castBitsToType<float>(operandBits2);
      float result = 0.0;
      if(!computeResultFP(floatOperand1, floatOperand2, m_cmdType, photonicsUtils::castBitsToType<float>(m_scalarValue), result)) return false;
      if (objDest.getDataType() == PHOTONICS_BOOL) {
        bool resultBool = result > 0;
        objDest.setElement(elemIdx, resultBool);
      } else {
        objDest.setElement(elemIdx, result);
      }
    } else {
      assert(0); // todo: data type
    }
  }
  return true;
}

//! @brief  PHOTONICS CMD: Functional 2-operand - update stats
bool
photonicsCmdFunc2::updateStats() const
{
  const photonicsObjInfo& objSrc1 = m_device->getResMgr()->getObjInfo(m_src1);
  const photonicsObjInfo& objSrc2 = m_device->getResMgr()->getObjInfo(m_src2);
  const photonicsObjInfo& objDest = m_device->getResMgr()->getObjInfo(m_dest);
  PhotonicsDataType dataType = objSrc1.getDataType();
  bool isVLayout = objSrc1.isVLayout();

  photonicseval::perfEnergy mPerfEnergy = photonicsSim::get()->getPerfEnergyModel()->getPerfEnergyForFunc2(m_cmdType, objSrc1, objSrc2, objDest);
  photonicsSim::get()->getStatsMgr()->recordCmd(getName(dataType, isVLayout), mPerfEnergy);
  return true;
}

//! @brief  PHOTONICS CMD: Iterative loop
bool
photonicsCmdIter::execute()
{
  if (m_debugCmds) {
    std::printf("PHOTONICS-Cmd: %s (obj id %d * %d -> %d) for %i loops.\n", getName().c_str(), m_src1, m_src2, m_dest, m_numLoops);
  }

  if (!sanityCheck()) {
    return false;
  }

  m_numVChips = m_device->getNumRanks();
  m_numHChips = m_device->getNumBankPerRank();
  m_numVectorsPerCore = m_device->getNumSubarrayPerBank();
  photonicsObjInfo &objSrc1 = m_device->getResMgr()->getObjInfo(m_src1);
  photonicsObjInfo &objSrc2 = m_device->getResMgr()->getObjInfo(m_src2);
  photonicsObjInfo &objDest = m_device->getResMgr()->getObjInfo(m_dest);

  if (photonicsSim::get()->getDeviceType() != PHOTONICS_FUNCTIONAL) {
    objSrc1.syncFromSimulatedMem();
    objSrc2.syncFromSimulatedMem();
  }

  const std::vector<photonicsRegion>& src2Regions = objSrc2.getRegions();
  const std::vector<photonicsRegion>& destRegions = objDest.getRegions();
  unsigned numRegions = destRegions.size();

  for (int i = 0; i < m_numLoops; i++) {
    computeAllRegions(numRegions);
    for (unsigned j = 0; j < m_numVChips; ++j) {
      float sumFinal = 0.0;
      uint64_t idx;
      for (unsigned k = 0; k < m_numHChips; ++k) {
        idx = destRegions[j*m_numHChips+k].getElemIdxBegin();
        uint32_t valBits = objDest.getElementBits(idx);
        float sumPartial = photonicsUtils::castBitsToType<float>(valBits);
        sumFinal += sumPartial;
      }
      idx = destRegions[j*m_numHChips].getElemIdxBegin();
      objDest.setElement(idx, sumFinal);
      idx = src2Regions[j].getElemIdxBegin();
      objSrc2.setElement(idx, sumFinal);
    }
  }

  if (photonicsSim::get()->getDeviceType() != PHOTONICS_FUNCTIONAL) {
    objDest.syncToSimulatedMem();
  }

  updateStats();
  return true;
}

//! @brief  PHOTONICS CMD: Iterative loop - sanity check
bool
photonicsCmdIter::sanityCheck() const
{
  photonicsResMgr* resMgr = m_device->getResMgr();
  if (!isValidObjId(resMgr, m_src1) || !isValidObjId(resMgr, m_src2) || !isValidObjId(resMgr, m_dest)) {
    return false;
  }
  const photonicsObjInfo& objSrc1 = resMgr->getObjInfo(m_src1);
  const photonicsObjInfo& objSrc2 = resMgr->getObjInfo(m_src2);
  const photonicsObjInfo& objDest = resMgr->getObjInfo(m_dest);
  if (!isAssociated(objSrc1, objSrc2) || !isAssociated(objSrc1, objDest)) {
    return false;
  }

  return true;
}

//! @brief  PHOTONICS CMD: Iterative loop - compute region
bool
photonicsCmdIter::computeRegion(unsigned index)
{
  const photonicsObjInfo& objSrc1 = m_device->getResMgr()->getObjInfo(m_src1);
  const photonicsObjInfo& objSrc2 = m_device->getResMgr()->getObjInfo(m_src2);
  photonicsObjInfo& objDest = m_device->getResMgr()->getObjInfo(m_dest);

  const photonicsRegion& src1Region = objSrc1.getRegions()[index];
  const photonicsRegion& src2Region = objSrc2.getRegions()[index%m_numHChips];
  const photonicsRegion& destRegion = objDest.getRegions()[index];

  // perform the computation
  uint64_t src1IdxBegin = src1Region.getElemIdxBegin();
  uint64_t src2IdxBegin = src2Region.getElemIdxBegin();
  uint64_t destIdxBegin = destRegion.getElemIdxBegin();
  unsigned numElementsInRegion = src1Region.getNumElemInRegion();
  float result = 0.0;
  for (unsigned j = 0; j < numElementsInRegion; ++j) {
    uint64_t src1ElemIdx = src1IdxBegin + j;
    uint64_t src2ElemIdx = src2IdxBegin + j;
    uint32_t operandBits1 = objSrc1.getElementBits(src1ElemIdx);
    uint32_t operandBits2 = objSrc2.getElementBits(src2ElemIdx);
    float floatOperand1 = photonicsUtils::castBitsToType<float>(operandBits1);
    float floatOperand2 = photonicsUtils::castBitsToType<float>(operandBits2);
    result = result + floatOperand1 * floatOperand2;
  }
  objDest.setElement(destIdxBegin, result);
  return true;
}

//! @brief  PHOTONICS CMD: Iterative loop - update stats
bool
photonicsCmdIter::updateStats() const
{
  const photonicsObjInfo& objSrc1 = m_device->getResMgr()->getObjInfo(m_src1);
  const photonicsObjInfo& objSrc2 = m_device->getResMgr()->getObjInfo(m_src2);
  const photonicsObjInfo& objDest = m_device->getResMgr()->getObjInfo(m_dest);
  PhotonicsDataType dataType = objSrc1.getDataType();
  bool isVLayout = objSrc1.isVLayout();

  photonicseval::perfEnergy mPerfEnergy = photonicsSim::get()->getPerfEnergyModel()->getPerfEnergyForIter(objSrc1, objSrc2, objDest, m_numLoops);
  photonicsSim::get()->getStatsMgr()->recordCmd(getName(dataType, isVLayout), mPerfEnergy);
  return true;
}

//! @brief  PHOTONICS CMD: Conditional Operations
bool
photonicsCmdCond::execute()
{
  if (m_debugCmds) {
    std::printf("PHOTONICS-Cmd: %s (obj ids: bool %d, src1 %d, src2 %d, dest %d, scalar 0x%" PRIx64 ")\n",
        getName().c_str(), m_condBool, m_src1, m_src2, m_dest, m_scalarBits);
  }

  if (!sanityCheck()) {
    return false;
  }

  if (photonicsSim::get()->getDeviceType() != PHOTONICS_FUNCTIONAL) {
    photonicsObjInfo &objBool = m_device->getResMgr()->getObjInfo(m_condBool);
    objBool.syncFromSimulatedMem();
    if (m_cmdType == PhotonicsCmdEnum::COND_COPY || m_cmdType == PhotonicsCmdEnum::COND_SELECT || m_cmdType == PhotonicsCmdEnum::COND_SELECT_SCALAR) {
      photonicsObjInfo &objSrc1 = m_device->getResMgr()->getObjInfo(m_src1);
      objSrc1.syncFromSimulatedMem();
    }
    if (m_cmdType == PhotonicsCmdEnum::COND_SELECT) {
      photonicsObjInfo &objSrc2 = m_device->getResMgr()->getObjInfo(m_src2);
      objSrc2.syncFromSimulatedMem();
    }
    if (m_cmdType == PhotonicsCmdEnum::COND_COPY || m_cmdType == PhotonicsCmdEnum::COND_BROADCAST) {  // require dest data to be synced
      photonicsObjInfo &objDest = m_device->getResMgr()->getObjInfo(m_dest);
      objDest.syncFromSimulatedMem();
    }
  }

  const photonicsObjInfo& objDest = m_device->getResMgr()->getObjInfo(m_dest);
  unsigned numRegions = objDest.getRegions().size();
  computeAllRegions(numRegions);

  if (photonicsSim::get()->getDeviceType() != PHOTONICS_FUNCTIONAL) {
    const photonicsObjInfo &objDest = m_device->getResMgr()->getObjInfo(m_dest);
    objDest.syncToSimulatedMem();
  }

  updateStats();
  return true;
}

//! @brief  PHOTONICS CMD: Conditional Operations - sanity check
bool
photonicsCmdCond::sanityCheck() const
{
  photonicsResMgr* resMgr = m_device->getResMgr();

  // common checks
  if (!isValidObjId(resMgr, m_condBool) || !isValidObjId(resMgr, m_dest)) {
    return false;
  }
  const photonicsObjInfo& objBool = resMgr->getObjInfo(m_condBool);
  if (objBool.getDataType() != PHOTONICS_BOOL) {
    std::printf("PHOTONICS-Error: PHOTONICS command %s condition must be PHOTONICS_BOOL type\n", getName().c_str());
    return false;
  }
  const photonicsObjInfo& objDest = resMgr->getObjInfo(m_dest);
  if (!isAssociated(objBool, objDest)) {
    return false;
  }

  // command specific checks
  if (m_cmdType == PhotonicsCmdEnum::COND_COPY || m_cmdType == PhotonicsCmdEnum::COND_SELECT || m_cmdType == PhotonicsCmdEnum::COND_SELECT_SCALAR) {
    if (!isValidObjId(resMgr, m_src1)) {
      return false;
    }
    const photonicsObjInfo& objSrc1 = resMgr->getObjInfo(m_src1);
    if (!isAssociated(objSrc1, objBool)) {
      return false;
    }
    if (objSrc1.getDataType() != objDest.getDataType()) {
      std::printf("PHOTONICS-Error: PHOTONICS command %s does not support mixed data type\n", getName().c_str());
      return false;
    }
  }
  if (m_cmdType == PhotonicsCmdEnum::COND_SELECT) {
    if (!isValidObjId(resMgr, m_src2)) {
      return false;
    }
    const photonicsObjInfo& objSrc2 = resMgr->getObjInfo(m_src2);
    if (!isAssociated(objSrc2, objBool)) {
      return false;
    }
    if (objSrc2.getDataType() != objDest.getDataType()) {
      std::printf("PHOTONICS-Error: PHOTONICS command %s does not support mixed data type\n", getName().c_str());
      return false;
    }
  }
  return true;
}

//! @brief  PHOTONICS CMD: Conditional Operations - compute region
bool
photonicsCmdCond::computeRegion(unsigned index)
{
  const photonicsObjInfo& objBool = m_device->getResMgr()->getObjInfo(m_condBool);
  photonicsObjInfo& objDest = m_device->getResMgr()->getObjInfo(m_dest);

  // perform the computation
  const photonicsRegion& destRegion = objDest.getRegions()[index];
  uint64_t elemIdxBegin = destRegion.getElemIdxBegin();
  unsigned numElementsInRegion = destRegion.getNumElemInRegion();
  switch (m_cmdType) {
    case PhotonicsCmdEnum::COND_COPY: {
      const photonicsObjInfo& objSrc1 = m_device->getResMgr()->getObjInfo(m_src1);
      for (unsigned j = 0; j < numElementsInRegion; ++j) {
        uint64_t elemIdx = elemIdxBegin + j;
        uint64_t bitsBool = objBool.getElementBits(elemIdx);
        uint64_t bitsSrc1 = objSrc1.getElementBits(elemIdx);
        uint64_t bitsDest = objDest.getElementBits(elemIdx);
        uint64_t bitsResult = bitsBool ? bitsSrc1 : bitsDest;
        objDest.setElement(elemIdx, bitsResult);
      }
      break;
    }
    case PhotonicsCmdEnum::COND_BROADCAST: {
      for (unsigned j = 0; j < numElementsInRegion; ++j) {
        uint64_t elemIdx = elemIdxBegin + j;
        uint64_t bitsBool = objBool.getElementBits(elemIdx);
        uint64_t bitsDest = objDest.getElementBits(elemIdx);
        uint64_t bitsResult = bitsBool ? m_scalarBits : bitsDest;
        objDest.setElement(elemIdx, bitsResult);
      }
      break;
    }
    case PhotonicsCmdEnum::COND_SELECT: {
      const photonicsObjInfo& objSrc1 = m_device->getResMgr()->getObjInfo(m_src1);
      const photonicsObjInfo& objSrc2 = m_device->getResMgr()->getObjInfo(m_src2);
      for (unsigned j = 0; j < numElementsInRegion; ++j) {
        uint64_t elemIdx = elemIdxBegin + j;
        uint64_t bitsBool = objBool.getElementBits(elemIdx);
        uint64_t bitsSrc1 = objSrc1.getElementBits(elemIdx);
        uint64_t bitsSrc2 = objSrc2.getElementBits(elemIdx);
        uint64_t bitsResult = bitsBool ? bitsSrc1 : bitsSrc2;
        objDest.setElement(elemIdx, bitsResult);
      }
      break;
    }
    case PhotonicsCmdEnum::COND_SELECT_SCALAR: {
      const photonicsObjInfo& objSrc1 = m_device->getResMgr()->getObjInfo(m_src1);
      for (unsigned j = 0; j < numElementsInRegion; ++j) {
        uint64_t elemIdx = elemIdxBegin + j;
        uint64_t bitsBool = objBool.getElementBits(elemIdx);
        uint64_t bitsSrc1 = objSrc1.getElementBits(elemIdx);
        uint64_t bitsResult = bitsBool ? bitsSrc1 : m_scalarBits;
        objDest.setElement(elemIdx, bitsResult);
      }
      break;
    }
    default:
      assert(0);
  }
  return true;
}

//! @brief  PHOTONICS CMD: Conditional Operations - update stats
bool
photonicsCmdCond::updateStats() const
{
  const photonicsObjInfo& objDest = m_device->getResMgr()->getObjInfo(m_dest);
  PhotonicsDataType dataType = objDest.getDataType();
  bool isVLayout = objDest.isVLayout();

  // Reuse func2 to calculate performance and energy
  photonicseval::perfEnergy mPerfEnergy = photonicsSim::get()->getPerfEnergyModel()->getPerfEnergyForFunc2(m_cmdType, objDest, objDest, objDest);
  photonicsSim::get()->getStatsMgr()->recordCmd(getName(dataType, isVLayout), mPerfEnergy);
  return true;
}
 
//! @brief  PHOTONICS CMD: redsum non-ranged/ranged - sanity check
template <typename T> bool
photonicsCmdReduction<T>::sanityCheck() const
{
  photonicsResMgr* resMgr = m_device->getResMgr();
  if (!isValidObjId(resMgr, m_src) || !m_result) {
    return false;
  }

  uint64_t numElements = m_device->getResMgr()->getObjInfo(m_src).getNumElements();
  if (m_idxBegin > numElements) {
    std::printf("PHOTONICS-Error: The beginning of the reduction range for PHOTONICS object ID %d is greater than the number of elements\n", m_src);
    return false;
  }
  if (m_idxEnd < m_idxBegin) {
    std::printf("PHOTONICS-Error: The end index of the reduction range for PHOTONICS object ID %d is less than the start index\n", m_src);
    return false;
  }
  return true;
}

template <typename T> bool
photonicsCmdReduction<T>::execute()
{
  if (m_debugCmds) {
    std::printf("PHOTONICS-Cmd: %s (obj id %d)\n", getName().c_str(), m_src);
  }

  if (!sanityCheck()) {
    return false;
  }

  photonicsObjInfo &objSrc = m_device->getResMgr()->getObjInfo(m_src);
  if (photonicsSim::get()->getDeviceType() != PHOTONICS_FUNCTIONAL) {
    objSrc.syncFromSimulatedMem();
  }

  unsigned numRegions = objSrc.getRegions().size();

  // prepare per-region storage
  //reduction
  for (unsigned i = 0; i < numRegions; ++i) {
    if (m_cmdType == PhotonicsCmdEnum::REDSUM || m_cmdType == PhotonicsCmdEnum::REDSUM_RANGE) {
      m_regionResult.resize(numRegions, 0);
    } else if (m_cmdType == PhotonicsCmdEnum::REDMIN || m_cmdType == PhotonicsCmdEnum::REDMIN_RANGE) {
      m_regionResult.resize(numRegions, std::numeric_limits<T>::max());
    } else if (m_cmdType == PhotonicsCmdEnum::REDMAX || m_cmdType == PhotonicsCmdEnum::REDMAX_RANGE) {
      m_regionResult.resize(numRegions, std::numeric_limits<T>::lowest());
    }
  }

  computeAllRegions(numRegions);
  
  //reduction
  for (unsigned i = 0; i < numRegions; ++i) {
    if (m_cmdType == PhotonicsCmdEnum::REDSUM || m_cmdType == PhotonicsCmdEnum::REDSUM_RANGE) {
      if (std::is_integral_v<T> && std::is_signed_v<T>)
      {
        *static_cast<int64_t *>(m_result) += static_cast<int64_t>(m_regionResult[i]);
      }
      else if (std::is_integral_v<T> && std::is_unsigned_v<T>)
      {
        *static_cast<uint64_t *>(m_result) += static_cast<uint64_t>(m_regionResult[i]);
      }
      else
      {
        *static_cast<float *>(m_result) += static_cast<float>(m_regionResult[i]);
      }
    }
    else if (m_cmdType == PhotonicsCmdEnum::REDMIN || m_cmdType == PhotonicsCmdEnum::REDMIN_RANGE)
    {
      *static_cast<T *>(m_result) = *static_cast<T *>(m_result) > static_cast<T>(m_regionResult[i]) ? static_cast<T>(m_regionResult[i]) : *static_cast<T *>(m_result);
    }
    else if (m_cmdType == PhotonicsCmdEnum::REDMAX || m_cmdType == PhotonicsCmdEnum::REDMAX_RANGE)
    {
      *static_cast<T *>(m_result) = *static_cast<T *>(m_result) < static_cast<T>(m_regionResult[i]) ? static_cast<T>(m_regionResult[i]) : *static_cast<T *>(m_result);
    }
  }

  updateStats();
  return true;
}

template <typename T> bool
photonicsCmdReduction<T>::computeRegion(unsigned index)
{
  const photonicsObjInfo& objSrc = m_device->getResMgr()->getObjInfo(m_src);
  const photonicsRegion& srcRegion = objSrc.getRegions()[index];
  PhotonicsDataType dataType = objSrc.getDataType();
  unsigned numElementsInRegion = srcRegion.getNumElemInRegion();
  uint64_t currIdx = srcRegion.getElemIdxBegin();

  for (unsigned j = 0; j < numElementsInRegion && currIdx < m_idxEnd; ++j) {
    if (currIdx >= m_idxBegin) {
      uint64_t operandBits = objSrc.getElementBits(currIdx);
      bool isFP = photonicsUtils::isFP(dataType);
      bool isSigned = photonicsUtils::isSigned(dataType);
      if (!isFP)
      {
        T integerOperand;
        if (isSigned)
        {
          integerOperand = photonicsUtils::signExt(operandBits, dataType);
        }
        else
        {
          integerOperand = static_cast<T>(operandBits);
        }

        if (m_cmdType == PhotonicsCmdEnum::REDSUM || m_cmdType == PhotonicsCmdEnum::REDSUM_RANGE)
        {
          m_regionResult[index] += integerOperand;
        }
        else if (m_cmdType == PhotonicsCmdEnum::REDMIN || m_cmdType == PhotonicsCmdEnum::REDMIN_RANGE)
        {
          m_regionResult[index] = m_regionResult[index] > integerOperand ? integerOperand : m_regionResult[index];
        }
        else if (m_cmdType == PhotonicsCmdEnum::REDMAX || m_cmdType == PhotonicsCmdEnum::REDMAX_RANGE)
        {
          m_regionResult[index] = m_regionResult[index] < integerOperand ? integerOperand : m_regionResult[index];
        }
      }
      else if (isFP)
      {
        float floatOperand = photonicsUtils::castBitsToType<float>(operandBits);
        if (m_cmdType == PhotonicsCmdEnum::REDSUM || m_cmdType == PhotonicsCmdEnum::REDSUM_RANGE)
        {
          m_regionResult[index] += floatOperand;
        }
        else if (m_cmdType == PhotonicsCmdEnum::REDMIN || m_cmdType == PhotonicsCmdEnum::REDMIN_RANGE)
        {
          m_regionResult[index] = m_regionResult[index] > floatOperand ? floatOperand : m_regionResult[index];
        }
        else if (m_cmdType == PhotonicsCmdEnum::REDMAX || m_cmdType == PhotonicsCmdEnum::REDMAX_RANGE)
        {
          m_regionResult[index] = m_regionResult[index] < floatOperand ? floatOperand : m_regionResult[index];
        }
      }
      else
      {
        assert(0); // Unexpected data type
      }
    }
    currIdx += 1;
  }
  return true;
}

template <typename T> bool
photonicsCmdReduction<T>::updateStats() const
{
  const photonicsObjInfo& objSrc = m_device->getResMgr()->getObjInfo(m_src);
  PhotonicsDataType dataType = objSrc.getDataType();
  bool isVLayout = objSrc.isVLayout();

  unsigned numPass = 0;
  if (m_cmdType == PhotonicsCmdEnum::REDSUM_RANGE || m_cmdType == PhotonicsCmdEnum::REDMIN_RANGE || m_cmdType == PhotonicsCmdEnum::REDMAX_RANGE) {
    // determine numPass for ranged reduction
    std::unordered_map<PhotonicsCoreId, unsigned> activeRegionPerCore;
    uint64_t index = 0;
    for (const auto& region : objSrc.getRegions()) {
      PhotonicsCoreId coreId = region.getCoreId();
      unsigned numElementsInRegion = region.getNumElemInRegion();
      bool isActive = index < m_idxEnd && index + numElementsInRegion - 1 >= m_idxBegin;
      if (isActive) {
        activeRegionPerCore[coreId]++;
      }
      index += numElementsInRegion;
    }
    for (const auto& [coreId, count] : activeRegionPerCore) {
      if (numPass < count) {
        numPass = count;
      }
    }
  } else {
    numPass = objSrc.getMaxNumRegionsPerCore();
  }

  photonicseval::perfEnergy mPerfEnergy = photonicsSim::get()->getPerfEnergyModel()->getPerfEnergyForReduction(m_cmdType, objSrc, numPass);
  photonicsSim::get()->getStatsMgr()->recordCmd(getName(dataType, isVLayout), mPerfEnergy);
  return true;
}

//! @brief  PHOTONICS CMD: broadcast a value to all elements
bool
photonicsCmdBroadcast::execute()
{
  if (m_debugCmds) {
    std::printf("PHOTONICS-Cmd: %s (obj id %d value %" PRIu64 ")\n", getName().c_str(), m_dest, m_signExtBits);
  }

  if (!sanityCheck()) {
    return false;
  }

  const photonicsObjInfo& objDest = m_device->getResMgr()->getObjInfo(m_dest);
  unsigned numRegions = objDest.getRegions().size();
  computeAllRegions(numRegions);

  if (photonicsSim::get()->getDeviceType() != PHOTONICS_FUNCTIONAL) {
    const photonicsObjInfo &objDest = m_device->getResMgr()->getObjInfo(m_dest);
    objDest.syncToSimulatedMem();
  }

  updateStats();
  return true;
}

//! @brief  PHOTONICS CMD: broadcast a value to all elements - sanity check
bool
photonicsCmdBroadcast::sanityCheck() const
{
  photonicsResMgr* resMgr = m_device->getResMgr();
  if (!isValidObjId(resMgr, m_dest)) {
    return false;
  }
  return true;
}

//! @brief  PHOTONICS CMD: broadcast a value to all elements - compute region
bool
photonicsCmdBroadcast::computeRegion(unsigned index)
{
  photonicsObjInfo& objDest = m_device->getResMgr()->getObjInfo(m_dest);
  const photonicsRegion& destRegion = objDest.getRegions()[index];

  uint64_t elemIdxBegin = destRegion.getElemIdxBegin();
  unsigned numElementsInRegion = destRegion.getNumElemInRegion();

  for (unsigned j = 0; j < numElementsInRegion; ++j) {
    objDest.setElement(elemIdxBegin + j, m_signExtBits);
  }
  return true;
}

//! @brief  PHOTONICS CMD: broadcast a value to all elements - update stats
bool
photonicsCmdBroadcast::updateStats() const
{
  const photonicsObjInfo& objDest = m_device->getResMgr()->getObjInfo(m_dest);
  PhotonicsDataType dataType = objDest.getDataType();
  bool isVLayout = objDest.isVLayout();

  photonicseval::perfEnergy mPerfEnergy = photonicsSim::get()->getPerfEnergyModel()->getPerfEnergyForBroadcast(m_cmdType, objDest);
  photonicsSim::get()->getStatsMgr()->recordCmd(getName(dataType, isVLayout), mPerfEnergy);
  return true;
}


//! @brief  PHOTONICS CMD: rotate right/left
bool
photonicsCmdRotate::execute()
{
  if (m_debugCmds) {
    std::printf("PHOTONICS-Cmd: %s (obj id %d)\n", getName().c_str(), m_src);
  }

  if (!sanityCheck()) {
    return false;
  }

  if (photonicsSim::get()->getDeviceType() != PHOTONICS_FUNCTIONAL) {
    photonicsObjInfo &objSrc = m_device->getResMgr()->getObjInfo(m_src);
    objSrc.syncFromSimulatedMem();
  }

  photonicsObjInfo& objSrc = m_device->getResMgr()->getObjInfo(m_src);
  unsigned numRegions = objSrc.getRegions().size();
  m_regionBoundary.resize(numRegions, 0);

  computeAllRegions(numRegions);

  // handle region boundaries
  if (m_cmdType == PhotonicsCmdEnum::ROTATE_ELEM_R || m_cmdType == PhotonicsCmdEnum::SHIFT_ELEM_R) {
    for (unsigned i = 0; i < numRegions; ++i) {
      const photonicsRegion &srcRegion = objSrc.getRegions()[i];
      uint64_t elemIdxBegin = srcRegion.getElemIdxBegin();
      uint64_t val = 0;
      if (i == 0 && m_cmdType == PhotonicsCmdEnum::ROTATE_ELEM_R) {
        val = m_regionBoundary[numRegions - 1];
      } else if (i > 0) {
        val = m_regionBoundary[i - 1];
      }
      objSrc.setElement(elemIdxBegin, val);
    }
  } else if (m_cmdType == PhotonicsCmdEnum::ROTATE_ELEM_L || m_cmdType == PhotonicsCmdEnum::SHIFT_ELEM_L) {
    for (unsigned i = 0; i < numRegions; ++i) {
      const photonicsRegion &srcRegion = objSrc.getRegions()[i];
      unsigned numElementsInRegion = srcRegion.getNumElemInRegion();
      uint64_t elemIdxBegin = srcRegion.getElemIdxBegin();
      uint64_t val = 0;
      if (i == numRegions - 1 && m_cmdType == PhotonicsCmdEnum::ROTATE_ELEM_L) {
        val = m_regionBoundary[0];
      } else if (i < numRegions - 1) {
        val = m_regionBoundary[i + 1];
      }
      objSrc.setElement(elemIdxBegin + numElementsInRegion - 1, val);
    }
  } else {
    assert(0);
  }

  if (photonicsSim::get()->getDeviceType() != PHOTONICS_FUNCTIONAL) {
    const photonicsObjInfo &objSrc = m_device->getResMgr()->getObjInfo(m_src);
    objSrc.syncToSimulatedMem();
  }

  updateStats();
  return true;
}

//! @brief  PHOTONICS CMD: rotate right/left - sanity check
bool
photonicsCmdRotate::sanityCheck() const
{
  photonicsResMgr* resMgr = m_device->getResMgr();
  if (!isValidObjId(resMgr, m_src)) {
    return false;
  }
  return true;
}

//! @brief  PHOTONICS CMD: rotate right/left - compute region
bool
photonicsCmdRotate::computeRegion(unsigned index)
{
  photonicsObjInfo& objSrc = m_device->getResMgr()->getObjInfo(m_src);

  const photonicsRegion& srcRegion = objSrc.getRegions()[index];

  // read out values
  uint64_t elemIdxBegin = srcRegion.getElemIdxBegin();
  unsigned numElementsInRegion = srcRegion.getNumElemInRegion();
  std::vector<uint64_t> regionVector(numElementsInRegion);
  for (unsigned j = 0; j < numElementsInRegion; ++j) {
    regionVector[j] = objSrc.getElementBits(elemIdxBegin + j);
  }

  // perform rotation
  if (m_cmdType == PhotonicsCmdEnum::ROTATE_ELEM_R || m_cmdType == PhotonicsCmdEnum::SHIFT_ELEM_R) {
    m_regionBoundary[index] = regionVector[numElementsInRegion - 1];
    uint64_t carry = 0;
    for (unsigned j = 0; j < numElementsInRegion; ++j) {
      uint64_t temp = regionVector[j];
      regionVector[j] = carry;
      carry = temp;
    }
  } else if (m_cmdType == PhotonicsCmdEnum::ROTATE_ELEM_L || m_cmdType == PhotonicsCmdEnum::SHIFT_ELEM_L) {
    m_regionBoundary[index] = regionVector[0];
    uint64_t carry = 0;
    for (int j = numElementsInRegion - 1; j >= 0; --j) {
      uint64_t temp = regionVector[j];
      regionVector[j] = carry;
      carry = temp;
    }
  } else {
    assert(0);
  }

  // write back values
  for (unsigned j = 0; j < numElementsInRegion; ++j) {
    objSrc.setElement(elemIdxBegin + j, regionVector[j]);
  }
  return true;
}

//! @brief  PHOTONICS CMD: rotate right/left - update stats
bool
photonicsCmdRotate::updateStats() const
{
  const photonicsObjInfo& objSrc = m_device->getResMgr()->getObjInfo(m_src);
  PhotonicsDataType dataType = objSrc.getDataType();
  bool isVLayout = objSrc.isVLayout();

  photonicseval::perfEnergy mPerfEnergy = photonicsSim::get()->getPerfEnergyModel()->getPerfEnergyForRotate(m_cmdType, objSrc);
  photonicsSim::get()->getStatsMgr()->recordCmd(getName(dataType, isVLayout), mPerfEnergy);
  return true;
}

//! @brief  PHOTONICS CMD: prefix sum - sanity check
bool
photonicsCmdPrefixSum::sanityCheck() const
{
  photonicsResMgr* resMgr = m_device->getResMgr();
  if (!isValidObjId(resMgr, m_src) || !isValidObjId(resMgr, m_dst)) {
    return false;
  }

  if (!(isAssociated(resMgr->getObjInfo(m_src), resMgr->getObjInfo(m_dst)))) {
    return false;
  }
  return true;
}

bool
photonicsCmdPrefixSum::execute()
{
  if (m_debugCmds) {
    std::printf("PHOTONICS-Cmd: %s (obj id %d)\n", getName().c_str(), m_src);
  }

  if (!sanityCheck()) {
    return false;
  }

  photonicsObjInfo &objSrc = m_device->getResMgr()->getObjInfo(m_src);
  if (photonicsSim::get()->getDeviceType() != PHOTONICS_FUNCTIONAL) {
    objSrc.syncFromSimulatedMem();
  }

  unsigned numRegions = objSrc.getRegions().size();
  computeAllRegions(numRegions);
  updateStats();
  return true;
}

bool
photonicsCmdPrefixSum::computeRegion(unsigned index)
{
  //TODO: Make it parallel
  if (index > 0) {
    return true;
  }
  const photonicsObjInfo& objSrc = m_device->getResMgr()->getObjInfo(m_src);
  photonicsObjInfo& objDst = m_device->getResMgr()->getObjInfo(m_dst);
  PhotonicsDataType dataType = objSrc.getDataType();

  for (uint64_t j = 0; j < objSrc.getNumElements(); ++j) {
      if (photonicsUtils::isSigned(dataType)) {
      uint64_t srcOperandBits = objSrc.getElementBits(j);
      if (j == 0) objDst.setElement(j, photonicsUtils::signExt(srcOperandBits, dataType));
      else
      {
        uint64_t dstOperandBits = objDst.getElementBits(j-1);
        int64_t operand1 = photonicsUtils::signExt(srcOperandBits, dataType);
        int64_t operand2 = photonicsUtils::signExt(dstOperandBits, dataType);
        int64_t result = operand1 + operand2;
        objDst.setElement(j, result);
      }
    } else if (photonicsUtils::isUnsigned(dataType)) {
      uint64_t unsignedOperand1 = objSrc.getElementBits(j);
      if (j == 0) objDst.setElement(j, unsignedOperand1);
      else {
        uint64_t unsignedOperand2 = objDst.getElementBits(j-1);
        uint64_t result = unsignedOperand1 + unsignedOperand2;
        objDst.setElement(j, result);
      }
    } else if (photonicsUtils::isFP(dataType)) {
      uint64_t operandBits1 = objSrc.getElementBits(j);
      float floatOperand1 = photonicsUtils::castBitsToType<float>(operandBits1);
      if (j == 0) objDst.setElement(j, floatOperand1);
      else {
        uint64_t operandBits2 = objDst.getElementBits(j-1);
        float floatOperand2 = photonicsUtils::castBitsToType<float>(operandBits2);
        float result = floatOperand1 + floatOperand2;
        objDst.setElement(j, result);
      }
    } else {
      assert(0); // todo: data type
    }
  }
  return true;
}

bool
photonicsCmdPrefixSum::updateStats() const
{
  const photonicsObjInfo& objSrc = m_device->getResMgr()->getObjInfo(m_src);
  PhotonicsDataType dataType = objSrc.getDataType();
  bool isVLayout = objSrc.isVLayout();
  photonicseval::perfEnergy mPerfEnergy = photonicsSim::get()->getPerfEnergyModel()->getPerfEnergyForPrefixSum(m_cmdType, objSrc);
  photonicsSim::get()->getStatsMgr()->recordCmd(getName(dataType, isVLayout), mPerfEnergy);
  return true;
}

//! @brief  PHOTONICS CMD: MAC - sanity check
template <typename T> bool
photonicsCmdMAC<T>::sanityCheck() const
{
  photonicsResMgr* resMgr = m_device->getResMgr();
  if (m_device->getSimTarget() != PHOTONICS_DEVICE_AIM) {
    std::printf("PHOTONICS-Error: PHOTONICS CMD %s is only supported on AiM.\n", getName().c_str());
    return false;
  }

  if (!isValidObjId(resMgr, m_src1) || !isValidObjId(resMgr, m_src2)) {
    return false;
  }

  if (!m_device->getResMgr()->getObjInfo(m_src2).isBuffer()) {
    std::printf("PHOTONICS-Error: PHOTONICS CMD %s requires source object %d to be buffers\n", getName().c_str(), m_src2);
    return false;
  }

  if (m_device->getResMgr()->getObjInfo(m_src1).getDataType() != m_device->getResMgr()->getObjInfo(m_src2).getDataType()) {
    std::printf("PHOTONICS-Error: PHOTONICS command %s does not support mixed data types.\n", getName().c_str());
    return false;
  }

  return true;
}

template <typename T> bool
photonicsCmdMAC<T>::execute()
{
  if (m_debugCmds) {
    std::printf("PHOTONICS-Cmd: %s (obj id %d and obj id %d)\n", getName().c_str(), m_src1, m_src2);
  }

  if (!sanityCheck()) {
    return false;
  }

  photonicsObjInfo &objSrc1 = m_device->getResMgr()->getObjInfo(m_src1);
  photonicsObjInfo &objSrc2 = m_device->getResMgr()->getObjInfo(m_src2);
  if (photonicsSim::get()->getDeviceType() != PHOTONICS_FUNCTIONAL) {
    objSrc1.syncFromSimulatedMem();
    objSrc2.syncFromSimulatedMem();
  }

  unsigned numRegions = objSrc1.getRegions().size();
  for (unsigned i = 0; i < numRegions; ++i) {
    m_regionResult.resize(numRegions, 0);
  }
  computeAllRegions(numRegions);
  
  //reduction
  for (unsigned i = 0; i < numRegions; ++i) {
    if (std::is_integral_v<T> && std::is_signed_v<T>)
    {
      switch (objSrc1.getDataType())
      {
      case PHOTONICS_INT8:
        static_cast<int8_t *>(m_dest)[objSrc1.getRegions()[i].getCoreId()] += static_cast<int8_t>(m_regionResult[i]);
        break;
      case PHOTONICS_INT16:
        static_cast<int16_t *>(m_dest)[objSrc1.getRegions()[i].getCoreId()] += static_cast<int16_t>(m_regionResult[i]);
        break;
      case PHOTONICS_INT32:
        static_cast<int32_t *>(m_dest)[objSrc1.getRegions()[i].getCoreId()] += static_cast<int32_t>(m_regionResult[i]);
        break;
      case PHOTONICS_INT64:
        static_cast<int64_t *>(m_dest)[objSrc1.getRegions()[i].getCoreId()] += static_cast<int64_t>(m_regionResult[i]);
        break;
      default:
        break;
      }
    }
    else if (std::is_integral_v<T> && std::is_unsigned_v<T>)
    {
      switch (objSrc1.getDataType())
      {
      case PHOTONICS_UINT8:
        static_cast<int8_t *>(m_dest)[objSrc1.getRegions()[i].getCoreId()] += static_cast<int8_t>(m_regionResult[i]);
        break;
      case PHOTONICS_UINT16:
        static_cast<int16_t *>(m_dest)[objSrc1.getRegions()[i].getCoreId()] += static_cast<int16_t>(m_regionResult[i]);
        break;
      case PHOTONICS_UINT32:
        static_cast<int32_t *>(m_dest)[objSrc1.getRegions()[i].getCoreId()] += static_cast<int32_t>(m_regionResult[i]);
        break;
      case PHOTONICS_UINT64:
        static_cast<int64_t *>(m_dest)[objSrc1.getRegions()[i].getCoreId()] += static_cast<int64_t>(m_regionResult[i]);
        break;
      default:
        break;
      }
    }
    else
    {
      static_cast<float *>(m_dest)[objSrc1.getRegions()[i].getCoreId()] += static_cast<float>(m_regionResult[i]);
    }
  }
  updateStats();
  return true;
}

template <typename T> bool
photonicsCmdMAC<T>::computeRegion(unsigned index)
{
  //TODO: Make it parallel
  const photonicsObjInfo& objSrc1 = m_device->getResMgr()->getObjInfo(m_src1);
  const photonicsObjInfo& objSrc2 = m_device->getResMgr()->getObjInfo(m_src2);

  const photonicsRegion& src1Region = objSrc1.getRegions()[index];
  PhotonicsDataType dataType = objSrc1.getDataType();
  uint64_t elemIdxBegin = src1Region.getElemIdxBegin();
  unsigned numElementsInRegion = src1Region.getNumElemInRegion();

  for (uint64_t j = 0; j < numElementsInRegion; ++j) {
    uint64_t elemIdx = elemIdxBegin + j;
    if (photonicsUtils::isSigned(dataType)) {
      uint64_t operandBits1 = objSrc1.getElementBits(elemIdx);
      uint64_t operandBits2 = objSrc2.getElementBits(j);
      int64_t operand1 = photonicsUtils::signExt(operandBits1, dataType);
      int64_t operand2 = photonicsUtils::signExt(operandBits2, dataType);
      m_regionResult[index] += operand1 * operand2;
    } else if (photonicsUtils::isUnsigned(dataType)) {
      uint64_t unsignedOperand1 = objSrc1.getElementBits(elemIdx);
      uint64_t unsignedOperand2 = objSrc2.getElementBits(elemIdx);
      m_regionResult[index] += unsignedOperand1 * unsignedOperand2;
    } else if (photonicsUtils::isFP(dataType)) {
      uint64_t operandBits1 = objSrc1.getElementBits(elemIdx);
      uint64_t operandBits2 = objSrc2.getElementBits(elemIdx);
      float floatOperand1 = photonicsUtils::castBitsToType<float>(operandBits1);
      float floatOperand2 = photonicsUtils::castBitsToType<float>(operandBits2);
      m_regionResult[index] += floatOperand1 * floatOperand2;
    } else {
      assert(0); // todo: data type
    }
  }
  return true;
}

template <typename T> bool
photonicsCmdMAC<T>::updateStats() const
{
  const photonicsObjInfo& objSrc = m_device->getResMgr()->getObjInfo(m_src1);
  PhotonicsDataType dataType = objSrc.getDataType();
  bool isVLayout = objSrc.isVLayout();
  photonicseval::perfEnergy mPerfEnergy = photonicsSim::get()->getPerfEnergyModel()->getPerfEnergyForMac(m_cmdType, objSrc);
  photonicsSim::get()->getStatsMgr()->recordCmd(getName(dataType, isVLayout), mPerfEnergy);
  return true;
}

//! @brief  Photonics CMD: BitSIMD-V: Read a row to SA
bool
photonicsCmdReadRowToSa::execute()
{
  if (m_debugCmds) {
    std::printf("PHOTONICS-MicroOp: BitSIMD-V ReadRowToSa (obj id %d ofst %u)\n", m_objId, m_ofst);
  }

  photonicsResMgr* resMgr = m_device->getResMgr();
  const photonicsObjInfo& objSrc = resMgr->getObjInfo(m_objId);
  for (unsigned i = 0; i < objSrc.getRegions().size(); ++i) {
    const photonicsRegion& srcRegion = objSrc.getRegions()[i];
    if (m_ofst >= srcRegion.getNumAllocRows()) {
      std::printf("PHOTONICS-Error: Row offset %u out of range [0, %u)\n", m_ofst, srcRegion.getNumAllocRows());
      return false;
    }
    PhotonicsCoreId coreId = srcRegion.getCoreId();
    m_device->getCore(coreId).readRow(srcRegion.getRowIdx() + m_ofst);
  }

  // Update stats
  photonicseval::perfEnergy prfEnrgy;
  photonicsSim::get()->getStatsMgr()->recordCmd(getName(), prfEnrgy);
  return true;
}

//! @brief  Photonics CMD: BitSIMD-V: Write SA to a row
bool
photonicsCmdWriteSaToRow::execute()
{
  if (m_debugCmds) {
    std::printf("PHOTONICS-MicroOp: BitSIMD-V WriteSaToRow (obj id %d ofst %u)\n", m_objId, m_ofst);
  }

  photonicsResMgr* resMgr = m_device->getResMgr();
  const photonicsObjInfo& objSrc = resMgr->getObjInfo(m_objId);
  for (unsigned i = 0; i < objSrc.getRegions().size(); ++i) {
    const photonicsRegion& srcRegion = objSrc.getRegions()[i];
    if (m_ofst >= srcRegion.getNumAllocRows()) {
      std::printf("PHOTONICS-Error: Row offset %u out of range [0, %u)\n", m_ofst, srcRegion.getNumAllocRows());
      return false;
    }
    PhotonicsCoreId coreId = srcRegion.getCoreId();
    m_device->getCore(coreId).writeRow(srcRegion.getRowIdx() + m_ofst);
  }

  // Update stats
  photonicseval::perfEnergy prfEnrgy;
  photonicsSim::get()->getStatsMgr()->recordCmd(getName(), prfEnrgy);
  return true;
}

//! @brief  Photonics CMD: BitSIMD-V: Row reg operations
bool
photonicsCmdRRegOp::execute()
{
  if (m_debugCmds) {
    std::printf("PHOTONICS-MicroOp: BitSIMD-V %s (obj-id %d dest-reg %d src-reg %d %d %d val %d)\n",
                getName().c_str(), m_objId, m_dest, m_src1, m_src2, m_src3, m_val);
  }

  photonicsResMgr* resMgr = m_device->getResMgr();
  const photonicsObjInfo& refObj = resMgr->getObjInfo(m_objId);
  for (unsigned i = 0; i < refObj.getRegions().size(); ++i) {
    const photonicsRegion& refRegion = refObj.getRegions()[i];
    PhotonicsCoreId coreId = refRegion.getCoreId();
    for (unsigned j = 0; j < m_device->getNumCols(); j++) {
      switch (m_cmdType) {
      case PhotonicsCmdEnum::RREG_MOV:
      {
        m_device->getCore(coreId).getRowReg(m_dest)[j] = m_device->getCore(coreId).getRowReg(m_src1)[j];
        break;
      }
      case PhotonicsCmdEnum::RREG_SET:
      {
        m_device->getCore(coreId).getRowReg(m_dest)[j] = m_val;
        break;
      }
      case PhotonicsCmdEnum::RREG_NOT:
      {
        bool src = m_device->getCore(coreId).getRowReg(m_src1)[j];
        m_device->getCore(coreId).getRowReg(m_dest)[j] = !src;
        break;
      }
      case PhotonicsCmdEnum::RREG_AND:
      {
        bool src1 = m_device->getCore(coreId).getRowReg(m_src1)[j];
        bool src2 = m_device->getCore(coreId).getRowReg(m_src2)[j];
        m_device->getCore(coreId).getRowReg(m_dest)[j] = (src1 & src2);
        break;
      }
      case PhotonicsCmdEnum::RREG_OR:
      {
        bool src1 = m_device->getCore(coreId).getRowReg(m_src1)[j];
        bool src2 = m_device->getCore(coreId).getRowReg(m_src2)[j];
        m_device->getCore(coreId).getRowReg(m_dest)[j] = src1 | src2;
        break;
      }
      case PhotonicsCmdEnum::RREG_NAND:
      {
        bool src1 = m_device->getCore(coreId).getRowReg(m_src1)[j];
        bool src2 = m_device->getCore(coreId).getRowReg(m_src2)[j];
        m_device->getCore(coreId).getRowReg(m_dest)[j] = !(src1 & src2);
        break;
      }
      case PhotonicsCmdEnum::RREG_NOR:
      {
        bool src1 = m_device->getCore(coreId).getRowReg(m_src1)[j];
        bool src2 = m_device->getCore(coreId).getRowReg(m_src2)[j];
        m_device->getCore(coreId).getRowReg(m_dest)[j] = !(src1 | src2);
        break;
      }
      case PhotonicsCmdEnum::RREG_XOR:
      {
        bool src1 = m_device->getCore(coreId).getRowReg(m_src1)[j];
        bool src2 = m_device->getCore(coreId).getRowReg(m_src2)[j];
        m_device->getCore(coreId).getRowReg(m_dest)[j] = src1 ^ src2;
        break;
      }
      case PhotonicsCmdEnum::RREG_XNOR:
      {
        bool src1 = m_device->getCore(coreId).getRowReg(m_src1)[j];
        bool src2 = m_device->getCore(coreId).getRowReg(m_src2)[j];
        m_device->getCore(coreId).getRowReg(m_dest)[j] = !(src1 ^ src2);
        break;
      }
      case PhotonicsCmdEnum::RREG_MAJ:
      {
        bool src1 = m_device->getCore(coreId).getRowReg(m_src1)[j];
        bool src2 = m_device->getCore(coreId).getRowReg(m_src2)[j];
        bool src3 = m_device->getCore(coreId).getRowReg(m_src3)[j];
        m_device->getCore(coreId).getRowReg(m_dest)[j] =
            ((src1 & src2) || (src1 & src3) || (src2 & src3));
        break;
      }
      case PhotonicsCmdEnum::RREG_SEL:
      {
        bool cond = m_device->getCore(coreId).getRowReg(m_src1)[j];
        bool src2 = m_device->getCore(coreId).getRowReg(m_src2)[j];
        bool src3 = m_device->getCore(coreId).getRowReg(m_src3)[j];
        m_device->getCore(coreId).getRowReg(m_dest)[j] = (cond ? src2 : src3);
        break;
      }
      default:
        std::printf("PHOTONICS-Error: Unexpected cmd type %d\n", static_cast<int>(m_cmdType));
        assert(0);
      }
    }
  }

  // Update stats
  photonicseval::perfEnergy prfEnrgy;
  photonicsSim::get()->getStatsMgr()->recordCmd(getName(), prfEnrgy);
  return true;
}


//! @brief  Photonics CMD: BitSIMD-V: row reg rotate right/left by one step
bool
photonicsCmdRRegRotate::execute()
{
  if (m_debugCmds) {
    std::printf("PHOTONICS-MicroOp: BitSIMD-V %s (obj-id %d src-reg %d)\n", getName().c_str(), m_objId, m_dest);
  }

  photonicsResMgr* resMgr = m_device->getResMgr();
  const photonicsObjInfo& objSrc = resMgr->getObjInfo(m_objId);
  if (m_cmdType == PhotonicsCmdEnum::RREG_ROTATE_R) {  // Right Rotate
    bool prevVal = 0;
    for (unsigned i = 0; i < objSrc.getRegions().size(); ++i) {
      const photonicsRegion &srcRegion = objSrc.getRegions()[i];
      PhotonicsCoreId coreId = srcRegion.getCoreId();
      for (unsigned j = 0; j < srcRegion.getNumAllocCols(); ++j) {
        unsigned colIdx = srcRegion.getColIdx() + j;
        bool tmp = m_device->getCore(coreId).getRowReg(m_dest)[colIdx];
        m_device->getCore(coreId).getRowReg(m_dest)[colIdx] = prevVal;
        prevVal = tmp;
      }
    }
    // write the last val to the first place
    const photonicsRegion &firstRegion = objSrc.getRegions().front();
    PhotonicsCoreId firstCoreId = firstRegion.getCoreId();
    unsigned firstColIdx = firstRegion.getColIdx();
    m_device->getCore(firstCoreId).getRowReg(m_dest)[firstColIdx] = prevVal;
  } else if (m_cmdType == PhotonicsCmdEnum::RREG_ROTATE_L) {  // Left Rotate
    bool prevVal = 0;
    for (unsigned i = objSrc.getRegions().size(); i > 0; --i) {
      const photonicsRegion &srcRegion = objSrc.getRegions()[i - 1];
      PhotonicsCoreId coreId = srcRegion.getCoreId();
      for (unsigned j = srcRegion.getNumAllocCols(); j > 0; --j) {
        unsigned colIdx = srcRegion.getColIdx() + j - 1;
        bool tmp = m_device->getCore(coreId).getRowReg(m_dest)[colIdx];
        m_device->getCore(coreId).getRowReg(m_dest)[colIdx] = prevVal;
        prevVal = tmp;
      }
    }
    // write the first val to the last place
    const photonicsRegion &lastRegion = objSrc.getRegions().back();
    PhotonicsCoreId lastCoreId = lastRegion.getCoreId();
    unsigned lastColIdx = lastRegion.getColIdx() + lastRegion.getNumAllocCols() - 1;
    m_device->getCore(lastCoreId).getRowReg(m_dest)[lastColIdx] = prevVal;
  }

  // Update stats
  photonicseval::perfEnergy prfEnrgy;
  photonicsSim::get()->getStatsMgr()->recordCmd(getName(), prfEnrgy);
  return true;
}

//! @brief  Photonics CMD: SIMDRAM: Analog based multi-row AP and AAP
bool
photonicsCmdAnalogAAP::execute()
{
  if (m_debugCmds) {
    printDebugInfo();
  }

  if (m_srcRows.empty()) {
    return false;
  }

  photonicsResMgr* resMgr = m_device->getResMgr();
  const photonicsObjInfo& objSrc = resMgr->getObjInfo(m_srcRows[0].first);

  // 1st activate: compute majority
  std::unordered_set<unsigned> visitedRows;
  for (unsigned i = 0; i < objSrc.getRegions().size(); ++i) {
    const photonicsRegion& srcRegion = objSrc.getRegions()[i];
    PhotonicsCoreId coreId = srcRegion.getCoreId();
    photonicsCore &core = m_device->getCore(coreId);

    std::vector<std::pair<unsigned, bool>> rowIdxs;
    for (const auto& objOfst : m_srcRows) {
      if (!isValidObjId(resMgr, objOfst.first)) {
        return false;
      }
      const photonicsObjInfo& obj = resMgr->getObjInfo(objOfst.first);
      if (!isAssociated(objSrc, obj)) {
        return false;
      }
      unsigned ofst = objOfst.second;
      unsigned idx = obj.getRegions()[i].getRowIdx() + ofst;
      bool isDCCN = obj.isDualContactRef();
      rowIdxs.emplace_back(idx, isDCCN);
      if (i == 0) { // sanity check
        if (visitedRows.find(idx) == visitedRows.end()) {
          visitedRows.insert(idx);
        } else {
          std::printf("PHOTONICS-Error: Cannot access same src row multiple times during AP/AAP\n");
          return false;
        }
      }
    }
    core.readMultiRows(rowIdxs);
  }

  // 2nd activate: write multiple rows
  if (!m_destRows.empty()) {
    for (unsigned i = 0; i < objSrc.getRegions().size(); ++i) {
      const photonicsRegion& srcRegion = objSrc.getRegions()[i];
      PhotonicsCoreId coreId = srcRegion.getCoreId();
      photonicsCore &core = m_device->getCore(coreId);

      std::vector<std::pair<unsigned, bool>> rowIdxs;
      for (const auto& objOfst : m_destRows) {
        if (!isValidObjId(resMgr, objOfst.first)) {
          return false;
        }
        const photonicsObjInfo& obj = resMgr->getObjInfo(objOfst.first);
        if (!isAssociated(objSrc, obj)) {
          return false;
        }
        unsigned ofst = objOfst.second;
        unsigned idx = obj.getRegions()[i].getRowIdx() + ofst;
        bool isDCCN = obj.isDualContactRef();
        rowIdxs.emplace_back(idx, isDCCN);
        if (i == 0) { // sanity check
          if (visitedRows.find(idx) == visitedRows.end()) {
            visitedRows.insert(idx);
          } else {
            std::printf("PHOTONICS-Error: Cannot access same src/dest row multiple times during AP/AAP\n");
            return false;
          }
        }
      }
      core.writeMultiRows(rowIdxs);
    }
  }

  // Update stats
  std::string cmdName = getName();
  cmdName += "@" + std::to_string(m_srcRows.size()) + "," + std::to_string(m_destRows.size());
  photonicseval::perfEnergy prfEnrgy;
  photonicsSim::get()->getStatsMgr()->recordCmd(cmdName, prfEnrgy);
  return true;
}

//! @brief  Photonics CMD: SIMDRAM: AP/AAP debug info
void
photonicsCmdAnalogAAP::printDebugInfo() const
{
  std::string msg;
  for (const auto &kv : m_srcRows) {
    msg += " " + std::to_string(kv.first) + "[" + std::to_string(kv.second) + "]";
  }
  if (!m_destRows.empty()) {
    msg += " ->";
  }
  for (const auto &kv : m_destRows) {
    msg += " " + std::to_string(kv.first) + "[" + std::to_string(kv.second) + "]";
  }
  std::printf("PHOTONICS-MicroOp: %s (#src = %lu, #dest = %lu, rows =%s)\n",
              getName().c_str(), m_srcRows.size(), m_destRows.size(), msg.c_str());
}

template class photonicsCmdReduction<int8_t>;
template class photonicsCmdReduction<int16_t>;
template class photonicsCmdReduction<int32_t>;
template class photonicsCmdReduction<int64_t>;
template class photonicsCmdReduction<uint8_t>;
template class photonicsCmdReduction<uint16_t>;
template class photonicsCmdReduction<uint32_t>;
template class photonicsCmdReduction<uint64_t>;
template class photonicsCmdReduction<float>;

template class photonicsCmdMAC<int8_t>;
template class photonicsCmdMAC<int16_t>;
template class photonicsCmdMAC<int32_t>;
template class photonicsCmdMAC<int64_t>;
template class photonicsCmdMAC<uint8_t>;
template class photonicsCmdMAC<uint16_t>;
template class photonicsCmdMAC<uint32_t>;
template class photonicsCmdMAC<uint64_t>;
template class photonicsCmdMAC<float>;
