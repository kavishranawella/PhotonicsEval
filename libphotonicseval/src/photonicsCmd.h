// File: photonicsCmd.h
// PHOTONICSeval Simulator - PHOTONICS Commands
// Copyright (c) 2024 University of Virginia
// This file is licensed under the MIT License.
// See the LICENSE file in the root of this repository for more details.

#ifndef LAVA_PHOTONICS_CMD_H
#define LAVA_PHOTONICS_CMD_H

#include "libphotonicseval.h"      // for PhotonicsDataType, PhotonicsObjId
#include "photonicsResMgr.h"       // for photonicsResMgr, photonicsObjInfo
#include "photonicsCore.h"         // for photonicsCore
#include "photonicsUtils.h"        // for photonicsDataTypeEnumToStr, threadWorker
#include <vector>            // for vector
#include <string>            // for string
#include <climits>            // for numeric_limits
#include <cassert>           // for assert
#include <bitset>            // for bitset
#include <variant>

class photonicsDevice;


enum class PhotonicsCmdEnum {
  NOOP = 0,
  COPY_H2D,
  COPY_D2H,
  COPY_D2D,
  COPY_O2O, // This copies data between two associated memory objects. Hence, will be treated as PHOTONICS command not data copy
  // Functional 1-operand
  ABS,
  POPCOUNT,
  SHIFT_BITS_R,
  SHIFT_BITS_L,
  ADD_SCALAR,
  SUB_SCALAR,
  MUL_SCALAR,
  DIV_SCALAR,
  AND_SCALAR,
  OR_SCALAR,
  XOR_SCALAR,
  XNOR_SCALAR,
  GT_SCALAR,
  LT_SCALAR,
  EQ_SCALAR,
  NE_SCALAR,
  MIN_SCALAR,
  MAX_SCALAR,
  CONVERT_TYPE,
  BIT_SLICE_EXTRACT,
  BIT_SLICE_INSERT,
  // Functional 2-operand
  ADD,
  SUB,
  MUL,
  SCALED_ADD,
  DIV,
  NOT,
  AND,
  OR,
  XOR,
  XNOR,
  GT,
  LT,
  EQ,
  NE,
  MIN,
  MAX,
  // Conditional operations
  COND_COPY,
  COND_BROADCAST,
  COND_SELECT,
  COND_SELECT_SCALAR,
  // Functional special
  REDSUM,
  REDSUM_RANGE,
  REDMIN,
  REDMIN_RANGE,
  REDMAX,
  REDMAX_RANGE,
  BROADCAST,
  ROTATE_ELEM_R,
  ROTATE_ELEM_L,
  SHIFT_ELEM_R,
  SHIFT_ELEM_L,
  AES_SBOX,
  AES_INVERSE_SBOX,
  PREFIX_SUM,
  MAC,

  // BitSIMD v-layout commands
  ROW_R,
  ROW_W,
  RREG_MOV,
  RREG_SET,
  RREG_NOT,
  RREG_AND,
  RREG_OR,
  RREG_NAND,
  RREG_NOR,
  RREG_XOR,
  RREG_XNOR,
  RREG_MAJ,
  RREG_SEL,
  RREG_ROTATE_R,
  RREG_ROTATE_L,
  // SIMDRAM
  ROW_AP,
  ROW_AAP,
  // PHOTONICS
  ITER,
};


//! @class  photonicsCmd
//! @brief  Photonics command base class
class photonicsCmd
{
public:
  photonicsCmd(PhotonicsCmdEnum cmdType);
  virtual ~photonicsCmd() {}

  void setDevice(photonicsDevice* device) { m_device = device; }
  virtual bool execute() = 0;

  std::string getName() const {
    return getName(m_cmdType, "");
  }
  std::string getName(PhotonicsDataType dataType, bool isVLayout) const {
    std::string suffix = "." + photonicsUtils::photonicsDataTypeEnumToStr(dataType);
    suffix += isVLayout ? ".v" : ".h";
    return getName(m_cmdType, suffix);
  }
  static std::string getName(PhotonicsCmdEnum cmdType, const std::string& suffix);

protected:
  bool isValidObjId(photonicsResMgr* resMgr, PhotonicsObjId objId) const;
  bool isAssociated(const photonicsObjInfo& obj1, const photonicsObjInfo& obj2) const;
  bool isCompatibleType(const photonicsObjInfo& obj1, const photonicsObjInfo& obj2) const;
  bool isConvertibleType(const photonicsObjInfo& src, const photonicsObjInfo& dest) const;

  unsigned getNumElementsInRegion(const photonicsRegion& region, unsigned bitsPerElement) const;

  virtual bool sanityCheck() const { return false; }
  virtual bool computeRegion(unsigned index) { return false; }
  virtual bool updateStats() const { return false; }
  bool computeAllRegions(unsigned numRegions);

  //! @brief  Utility: Get bits of an element from a region. The bits are stored as uint64_t without sign extension
  inline uint64_t getBits(const photonicsCore& core, bool isVLayout, unsigned rowLoc, unsigned colLoc, unsigned numBits) const
  {
    return isVLayout ? core.getBitsV(rowLoc, colLoc, numBits) : core.getBitsH(rowLoc, colLoc, numBits);
  }

  //! @brief  Utility: Set bits of an element to a region
  inline void setBits(photonicsCore& core, bool isVLayout, unsigned rowLoc, unsigned colLoc, uint64_t bits, unsigned numBits) const
  {
    if (isVLayout) {
      core.setBitsV(rowLoc, colLoc, bits, numBits);
    } else {
      core.setBitsH(rowLoc, colLoc, bits, numBits);
    }
  }

  PhotonicsCmdEnum m_cmdType;
  photonicsDevice* m_device = nullptr;
  bool m_debugCmds;

  //! @class  photonicsCmd::regionWorker
  //! @brief  Thread worker to process regions in parallel
  class regionWorker : public photonicsUtils::threadWorker {
  public:
    regionWorker(photonicsCmd* cmd, unsigned regionIdx) : m_cmd(cmd), m_regionIdx(regionIdx) {}
    virtual ~regionWorker() {}
    virtual void execute() {
      m_cmd->computeRegion(m_regionIdx);
    }
  private:
    photonicsCmd* m_cmd = nullptr;
    const unsigned m_regionIdx = 0;
  };
};

//! @class  photonicsCmdDataTransfer
//! @brief  Data transfer. Not tracked as a regular Photonics CMD
class photonicsCmdCopy : public photonicsCmd
{
public:
  photonicsCmdCopy(PhotonicsCmdEnum cmdType, PhotonicsCopyEnum copyType, void* src, PhotonicsObjId dest, uint64_t idxBegin = 0, uint64_t idxEnd = 0)
    : photonicsCmd(PhotonicsCmdEnum::COPY_H2D), m_copyType(copyType), m_ptr(src), m_dest(dest), m_idxBegin(idxBegin), m_idxEnd(idxEnd), m_copyFullRange(idxEnd == 0ULL) {}
  photonicsCmdCopy(PhotonicsCmdEnum cmdType, PhotonicsCopyEnum copyType, PhotonicsObjId src, void* dest, uint64_t idxBegin = 0, uint64_t idxEnd = 0)
    : photonicsCmd(PhotonicsCmdEnum::COPY_D2H), m_copyType(copyType), m_ptr(dest), m_src(src), m_idxBegin(idxBegin), m_idxEnd(idxEnd), m_copyFullRange(idxEnd == 0ULL) {}
  photonicsCmdCopy(PhotonicsCmdEnum cmdType, PhotonicsCopyEnum copyType, PhotonicsObjId src, PhotonicsObjId dest, uint64_t idxBegin = 0, uint64_t idxEnd = 0)
    : photonicsCmd(PhotonicsCmdEnum::COPY_D2D), m_copyType(copyType), m_src(src), m_dest(dest), m_idxBegin(idxBegin), m_idxEnd(idxEnd), m_copyFullRange(idxEnd == 0ULL) {}

  virtual ~photonicsCmdCopy() {}
  virtual bool execute() override;
  virtual bool sanityCheck() const override;
  virtual bool updateStats() const override;
protected:
  PhotonicsCopyEnum m_copyType;
  void* m_ptr = nullptr;
  PhotonicsObjId m_src = -1;
  PhotonicsObjId m_dest = -1;
  uint64_t m_idxBegin = 0;
  uint64_t m_idxEnd = 0; 
  bool m_copyFullRange = false;
};

//! @class  photonicsCmdFunc1
//! @brief  Photonics CMD: Functional 1-operand
class photonicsCmdFunc1 : public photonicsCmd
{
public:
  photonicsCmdFunc1(PhotonicsCmdEnum cmdType, PhotonicsObjId src, PhotonicsObjId dest, uint64_t scalarValue = 0)
    : photonicsCmd(cmdType), m_src(src), m_dest(dest), m_scalarValue(scalarValue) {}
  photonicsCmdFunc1(PhotonicsCmdEnum cmdType, PhotonicsObjId src, PhotonicsObjId dest, const std::vector<uint8_t>& lut)
    : photonicsCmd(cmdType), m_src(src), m_dest(dest), m_lut(lut) {}
  virtual ~photonicsCmdFunc1() {}
  virtual bool execute() override;
  virtual bool sanityCheck() const override;
  virtual bool computeRegion(unsigned index) override;
  virtual bool updateStats() const override;
protected:
  PhotonicsObjId m_src;
  PhotonicsObjId m_dest;
  uint64_t m_scalarValue;
  std::vector<uint8_t> m_lut; 
private:
  template<typename T>
  inline bool computeResult(T operand, PhotonicsCmdEnum cmdType, T scalarValue, T& result, int bitsPerElementSrc) {
    result = operand;
    switch (cmdType) {
    case PhotonicsCmdEnum::COPY_O2O: result = operand; break;
    case PhotonicsCmdEnum::ADD_SCALAR: result += scalarValue; break;
    case PhotonicsCmdEnum::SUB_SCALAR: result -= scalarValue; break;
    case PhotonicsCmdEnum::MUL_SCALAR: result *= scalarValue; break;
    case PhotonicsCmdEnum::DIV_SCALAR:
        if (scalarValue == 0) {
            std::printf("PHOTONICS-Error: Division by zero\n");
            return false;
        }
        result /= scalarValue;
        break;
    case PhotonicsCmdEnum::NOT: result = ~operand; break;
    case PhotonicsCmdEnum::AND_SCALAR: result &= scalarValue; break;
    case PhotonicsCmdEnum::OR_SCALAR: result |= scalarValue; break;
    case PhotonicsCmdEnum::XOR_SCALAR: result ^= scalarValue; break;
    case PhotonicsCmdEnum::XNOR_SCALAR: result = ~(operand ^ scalarValue); break;
    case PhotonicsCmdEnum::GT_SCALAR: result = (operand > scalarValue) ? 1 : 0; break;
    case PhotonicsCmdEnum::LT_SCALAR: result = (operand < scalarValue) ? 1 : 0; break;
    case PhotonicsCmdEnum::EQ_SCALAR: result = (operand == scalarValue) ? 1 : 0; break;
    case PhotonicsCmdEnum::NE_SCALAR: result = (operand != scalarValue) ? 1 : 0; break;
    case PhotonicsCmdEnum::MIN_SCALAR: result = std::min(operand, scalarValue); break;
    case PhotonicsCmdEnum::MAX_SCALAR: result = std::max(operand, scalarValue); break;
    case PhotonicsCmdEnum::POPCOUNT:
        switch (bitsPerElementSrc) {
        case 8: result = std::bitset<8>(operand).count(); break;
        case 16: result = std::bitset<16>(operand).count(); break;
        case 32: result = std::bitset<32>(operand).count(); break;
        case 64: result = std::bitset<64>(operand).count(); break;
        default:
            std::printf("PHOTONICS-Error: Unsupported bits per element %u\n", bitsPerElementSrc);
            return false;
        }
        break;
    case PhotonicsCmdEnum::SHIFT_BITS_R: result >>= static_cast<uint64_t>(scalarValue); break;
    case PhotonicsCmdEnum::SHIFT_BITS_L: result <<= static_cast<uint64_t>(scalarValue); break;
    case PhotonicsCmdEnum::ABS:
    {
        if (std::is_signed<T>::value) {
          result = (operand < 0) ? -operand : operand;
        } else {
          result = operand;
        }
        break;
    }
    case PhotonicsCmdEnum::AES_SBOX:
    case PhotonicsCmdEnum::AES_INVERSE_SBOX:
      result = m_lut[operand]; 
      break;
    default:
        std::printf("PHOTONICS-Error: Unexpected cmd type %d\n", static_cast<int>(cmdType));
        assert(0);
    }
    return true;
  }

  template<typename T>
  inline bool computeResultFP(T operand, PhotonicsCmdEnum cmdType, T scalerValue, T& result) {
    result = operand;
    switch (cmdType) {
    case PhotonicsCmdEnum::COPY_O2O: result = operand; break;
    case PhotonicsCmdEnum::ADD_SCALAR: result += scalerValue; break;
    case PhotonicsCmdEnum::SUB_SCALAR: result -= scalerValue; break;
    case PhotonicsCmdEnum::MUL_SCALAR: result *= scalerValue; break;
    case PhotonicsCmdEnum::DIV_SCALAR:
        if (scalerValue == 0) {
            std::printf("PHOTONICS-Error: Division by zero\n");
            return false;
        }
        result /= scalerValue;
        break;
    case PhotonicsCmdEnum::GT_SCALAR: result = (operand > scalerValue) ? 1 : 0; break;
    case PhotonicsCmdEnum::LT_SCALAR: result = (operand < scalerValue) ? 1 : 0; break;
    case PhotonicsCmdEnum::EQ_SCALAR: result = (operand == scalerValue) ? 1 : 0; break;
    case PhotonicsCmdEnum::NE_SCALAR: result = (operand != scalerValue) ? 1 : 0; break;
    case PhotonicsCmdEnum::MIN_SCALAR: result = std::min(operand, scalerValue); break;
    case PhotonicsCmdEnum::MAX_SCALAR: result = std::max(operand, scalerValue); break;
    case PhotonicsCmdEnum::ABS:
    {
        if (std::is_signed<T>::value) {
          result = (operand < 0) ? -operand : operand;
        } else {
          result = operand;
        }
        break;
    }
    case PhotonicsCmdEnum::AND_SCALAR:
    case PhotonicsCmdEnum::OR_SCALAR:
    case PhotonicsCmdEnum::XOR_SCALAR:
    case PhotonicsCmdEnum::XNOR_SCALAR:
    case PhotonicsCmdEnum::POPCOUNT:
    case PhotonicsCmdEnum::SHIFT_BITS_R:
    case PhotonicsCmdEnum::SHIFT_BITS_L:
        std::printf("PHOTONICS-Error: Cannot perform bitwise operation on floating point values.\n");
        return false;
    default:
        std::printf("PHOTONICS-Error: Unexpected cmd type %d\n", static_cast<int>(cmdType));
        assert(0);
    }
    return true;
  }

  bool convertType(const photonicsObjInfo& objSrc, photonicsObjInfo& objDest, uint64_t elemIdx) const;
  bool bitSliceExtract(const photonicsObjInfo& objSrc, photonicsObjInfo& objDestBool, uint64_t bitIdx, uint64_t elemIdx) const;
  bool bitSliceInsert(const photonicsObjInfo& objSrcBool, photonicsObjInfo& objDest, uint64_t bitIdx, uint64_t elemIdx) const;
};

//! @class  photonicsCmdFunc2
//! @brief  Photonics CMD: Functional 2-operand
class photonicsCmdFunc2 : public photonicsCmd
{
public:
  photonicsCmdFunc2(PhotonicsCmdEnum cmdType, PhotonicsObjId src1, PhotonicsObjId src2, PhotonicsObjId dest)
    : photonicsCmd(cmdType), m_src1(src1), m_src2(src2), m_dest(dest) {}
  photonicsCmdFunc2(PhotonicsCmdEnum cmdType, PhotonicsObjId src1, PhotonicsObjId src2, PhotonicsObjId dest, uint64_t scalarValue)
    : photonicsCmd(cmdType), m_src1(src1), m_src2(src2), m_dest(dest), m_scalarValue(scalarValue) {}
  virtual ~photonicsCmdFunc2() {}
  virtual bool execute() override;
  virtual bool sanityCheck() const override;
  virtual bool computeRegion(unsigned index) override;
  virtual bool updateStats() const override;
protected:
  PhotonicsObjId m_src1;
  PhotonicsObjId m_src2;
  PhotonicsObjId m_dest;
  uint64_t m_scalarValue;
private:
  template<typename T>
  inline bool computeResult(T operand1, T operand2, PhotonicsCmdEnum cmdType, T scalarValue, T& result) {
    switch (cmdType) {
    case PhotonicsCmdEnum::ADD: result = operand1 + operand2; break;
    case PhotonicsCmdEnum::SUB: result = operand1 - operand2; break;
    case PhotonicsCmdEnum::MUL: result = operand1 * operand2; break;
    case PhotonicsCmdEnum::DIV:
        if (operand2 == 0) {
            std::printf("PHOTONICS-Error: Division by zero\n");
            return false;
        }
        result = operand1 / operand2;
        break;
    case PhotonicsCmdEnum::AND: result = operand1 & operand2; break;
    case PhotonicsCmdEnum::OR: result = operand1 | operand2; break;
    case PhotonicsCmdEnum::XOR: result = operand1 ^ operand2; break;
    case PhotonicsCmdEnum::XNOR: result = ~(operand1 ^ operand2); break;
    case PhotonicsCmdEnum::GT: result = operand1 > operand2 ? 1 : 0; break;
    case PhotonicsCmdEnum::LT: result = operand1 < operand2 ? 1 : 0; break;
    case PhotonicsCmdEnum::EQ: result = operand1 == operand2 ? 1 : 0; break;
    case PhotonicsCmdEnum::NE: result = operand1 != operand2 ? 1 : 0; break;
    case PhotonicsCmdEnum::MIN: result = (operand1 < operand2) ? operand1 : operand2; break;
    case PhotonicsCmdEnum::MAX: result = (operand1 > operand2) ? operand1 : operand2; break;
    case PhotonicsCmdEnum::SCALED_ADD: result = (operand1 * scalarValue) + operand2; break;
    default:
        std::printf("PHOTONICS-Error: Unexpected cmd type %d\n", static_cast<int>(m_cmdType));
          assert(0);
    }
    return true;
  }

  template<typename T>
  inline bool computeResultFP(T operand1, T operand2, PhotonicsCmdEnum cmdType, T scalarValue, T& result) {
    switch (cmdType) {
    case PhotonicsCmdEnum::ADD: result = operand1 + operand2; break;
    case PhotonicsCmdEnum::SUB: result = operand1 - operand2; break;
    case PhotonicsCmdEnum::MUL: result = operand1 * operand2; break;
    case PhotonicsCmdEnum::DIV:
        if (operand2 == 0) {
            std::printf("PHOTONICS-Error: Division by zero\n");
            return false;
        }
        result = operand1 / operand2;
        break;
    case PhotonicsCmdEnum::GT: result = operand1 > operand2 ? 1 : 0; break;
    case PhotonicsCmdEnum::LT: result = operand1 < operand2 ? 1 : 0; break;
    case PhotonicsCmdEnum::EQ: result = operand1 == operand2 ? 1 : 0; break;
    case PhotonicsCmdEnum::NE: result = operand1 != operand2 ? 1 : 0; break;
    case PhotonicsCmdEnum::MIN: result = (operand1 < operand2) ? operand1 : operand2; break;
    case PhotonicsCmdEnum::MAX: result = (operand1 > operand2) ? operand1 : operand2; break;
    case PhotonicsCmdEnum::SCALED_ADD: result = (operand1 * scalarValue) + operand2; break;
    case PhotonicsCmdEnum::AND:
    case PhotonicsCmdEnum::OR:
    case PhotonicsCmdEnum::XOR:
    case PhotonicsCmdEnum::XNOR:
        std::printf("PHOTONICS-Error: Cannot perform bitwise operation on floating point values.\n");
        return false;
    default:
        std::printf("PHOTONICS-Error: Unexpected cmd type %d\n", static_cast<int>(m_cmdType));
          assert(0);
    }
    return true;
  }
};

//! @class  photonicsCmdIter
//! @brief  Photonics CMD: Iterative loop
class photonicsCmdIter : public photonicsCmd
{
public:
  photonicsCmdIter(PhotonicsObjId src1, PhotonicsObjId src2, PhotonicsObjId dest, int8_t numLoops = 1)
    : photonicsCmd((PhotonicsCmdEnum::ITER)), m_src1(src1), m_src2(src2), m_dest(dest), m_numLoops(numLoops) {}
  virtual ~photonicsCmdIter() {}
  virtual bool execute() override;
  virtual bool sanityCheck() const override;
  virtual bool computeRegion(unsigned index) override;
  virtual bool updateStats() const override;
protected:
  PhotonicsObjId m_src1;
  PhotonicsObjId m_src2;
  PhotonicsObjId m_dest;
  int8_t m_numLoops;
  unsigned m_numVChips;
  unsigned m_numHChips;
  unsigned m_numVectorsPerCore;
private:

};

//! @class  photonicsCmdCond
//! @brief  Photonics CMD: Conditional operations using BOOL as the first operand
//!   COND_COPY:          dest[i] = cond ? src[i] : dest[i]
//!   COND_BROADCAST:     dest[i] = cond ? scalar : dest[i]
//!   COND_SELECT:        dest[i] = cond ? src1[i] : src2[i]
//!   COND_SELECT_SCALAR: dest[i] = cond ? src[1] : scalar
class photonicsCmdCond : public photonicsCmd
{
public:
  photonicsCmdCond(PhotonicsCmdEnum cmdType, PhotonicsObjId condBool, PhotonicsObjId src1, PhotonicsObjId dest)
    : photonicsCmd(cmdType), m_condBool(condBool), m_src1(src1), m_dest(dest)
  {
    assert(cmdType == PhotonicsCmdEnum::COND_COPY);
  }
  photonicsCmdCond(PhotonicsCmdEnum cmdType, PhotonicsObjId condBool, uint64_t scalarBits, PhotonicsObjId dest)
    : photonicsCmd(cmdType), m_condBool(condBool), m_scalarBits(scalarBits), m_dest(dest)
  {
    assert(cmdType == PhotonicsCmdEnum::COND_BROADCAST);
  }
  photonicsCmdCond(PhotonicsCmdEnum cmdType, PhotonicsObjId condBool, PhotonicsObjId src1, PhotonicsObjId src2, PhotonicsObjId dest)
    : photonicsCmd(cmdType), m_condBool(condBool), m_src1(src1), m_src2(src2), m_dest(dest)
  {
    assert(cmdType == PhotonicsCmdEnum::COND_SELECT);
  }
  photonicsCmdCond(PhotonicsCmdEnum cmdType, PhotonicsObjId condBool, PhotonicsObjId src1, uint64_t scalarBits, PhotonicsObjId dest)
  : photonicsCmd(cmdType), m_condBool(condBool), m_src1(src1), m_scalarBits(scalarBits), m_dest(dest)
  {
    assert(cmdType == PhotonicsCmdEnum::COND_SELECT_SCALAR);
  }
  virtual ~photonicsCmdCond() {}
  virtual bool execute() override;
  virtual bool sanityCheck() const override;
  virtual bool computeRegion(unsigned index) override;
  virtual bool updateStats() const override;
protected:
  PhotonicsObjId m_condBool;
  PhotonicsObjId m_src1 = -1;
  PhotonicsObjId m_src2 = -1;
  uint64_t m_scalarBits = 0;
  PhotonicsObjId m_dest;
};

//! @class  photonicsCmdReduction
//! @brief  Photonics CMD: Reduction non-ranged/ranged
template <typename T>
class photonicsCmdReduction : public photonicsCmd
{
public:
  photonicsCmdReduction(PhotonicsCmdEnum cmdType, PhotonicsObjId src, void* result)
    : photonicsCmd(cmdType), m_src(src), m_result(result)
  {
    assert(cmdType == PhotonicsCmdEnum::REDSUM || cmdType == PhotonicsCmdEnum::REDMIN || cmdType == PhotonicsCmdEnum::REDMAX);
  }
  photonicsCmdReduction(PhotonicsCmdEnum cmdType, PhotonicsObjId src, void* result, uint64_t idxBegin, uint64_t idxEnd)
    : photonicsCmd(cmdType), m_src(src), m_result(result), m_idxBegin(idxBegin)
  {
    assert(cmdType == PhotonicsCmdEnum::REDSUM || cmdType == PhotonicsCmdEnum::REDMIN || cmdType == PhotonicsCmdEnum::REDMAX || cmdType == PhotonicsCmdEnum::REDSUM_RANGE || cmdType == PhotonicsCmdEnum::REDMIN_RANGE || cmdType == PhotonicsCmdEnum::REDMAX_RANGE);
    if (idxEnd) m_idxEnd = idxEnd;
  }
  virtual ~photonicsCmdReduction() {}
  virtual bool execute() override;
  virtual bool sanityCheck() const override;
  virtual bool computeRegion(unsigned index) override;
  virtual bool updateStats() const override;
protected:
  PhotonicsObjId m_src;
  void* m_result;
  std::vector<T> m_regionResult;
  uint64_t m_idxBegin = 0;
  uint64_t m_idxEnd = std::numeric_limits<uint64_t>::max();
};

//! @class  photonicsCmdPrefixSum
//! @brief  Photonics CMD: PrefixSum
class photonicsCmdPrefixSum : public photonicsCmd
{
public:
  photonicsCmdPrefixSum(PhotonicsCmdEnum cmdType, PhotonicsObjId src, PhotonicsObjId dest)
    : photonicsCmd(cmdType), m_src(src), m_dst(dest)
  {
    assert(cmdType == PhotonicsCmdEnum::PREFIX_SUM);
  }
  virtual ~photonicsCmdPrefixSum() {}
  virtual bool execute() override;
  virtual bool sanityCheck() const override;
  virtual bool computeRegion(unsigned index) override;
  virtual bool updateStats() const override;
protected:
  PhotonicsObjId m_src, m_dst;
};

//! @class  photonicsCmdMAC
//! @brief  Photonics CMD: Multiply-Accumulate
template <typename T>
class photonicsCmdMAC : public photonicsCmd
{
public:
  photonicsCmdMAC(PhotonicsCmdEnum cmdType, PhotonicsObjId src1, PhotonicsObjId src2, void* dest)
    : photonicsCmd(cmdType), m_src1(src1), m_src2(src2), m_dest(dest)
  {
    assert(cmdType == PhotonicsCmdEnum::MAC);
  }
  virtual ~photonicsCmdMAC() {}
  virtual bool execute() override;
  virtual bool sanityCheck() const override;
  virtual bool computeRegion(unsigned index) override;
  virtual bool updateStats() const override;
protected:
  std::vector<T> m_regionResult;
  PhotonicsObjId m_src1, m_src2;
  void* m_dest; // Pointer to the destination buffer where MAC results will be stored
};


//! @class  photonicsCmdBroadcast
//! @brief  Photonics CMD: Broadcast a value to all elements
class photonicsCmdBroadcast : public photonicsCmd
{
public:
  photonicsCmdBroadcast(PhotonicsCmdEnum cmdType, PhotonicsObjId dest, uint64_t signExtBits)
    : photonicsCmd(cmdType), m_dest(dest), m_signExtBits(signExtBits)
  {
    assert(cmdType == PhotonicsCmdEnum::BROADCAST);
  }
  virtual ~photonicsCmdBroadcast() {}
  virtual bool execute() override;
  virtual bool sanityCheck() const override;
  virtual bool computeRegion(unsigned index) override;
  virtual bool updateStats() const override;
protected:
  PhotonicsObjId m_dest;
  uint64_t m_signExtBits;
};

//! @class  photonicsCmdRotate
//! @brief  Photonics CMD: rotate/shift elements right/left
class photonicsCmdRotate : public photonicsCmd
{
public:
  photonicsCmdRotate(PhotonicsCmdEnum cmdType, PhotonicsObjId src)
    : photonicsCmd(cmdType), m_src(src)
  {
    assert(cmdType == PhotonicsCmdEnum::ROTATE_ELEM_R || cmdType == PhotonicsCmdEnum::ROTATE_ELEM_L ||
           cmdType == PhotonicsCmdEnum::SHIFT_ELEM_R || cmdType == PhotonicsCmdEnum::SHIFT_ELEM_L);
  }
  virtual ~photonicsCmdRotate() {}
  virtual bool execute() override;
  virtual bool sanityCheck() const override;
  virtual bool computeRegion(unsigned index) override;
  virtual bool updateStats() const override;
protected:
  PhotonicsObjId m_src;
  std::vector<uint64_t> m_regionBoundary;
};

//! @class  photonicsCmdReadRowToSa
//! @brief  Photonics CMD: BitSIMD-V: Read a row to SA
class photonicsCmdReadRowToSa : public photonicsCmd
{
public:
  photonicsCmdReadRowToSa(PhotonicsCmdEnum cmdType, PhotonicsObjId objId, unsigned ofst)
    : photonicsCmd(cmdType), m_objId(objId), m_ofst(ofst) {}
  virtual ~photonicsCmdReadRowToSa() {}
  virtual bool execute() override;
protected:
  PhotonicsObjId m_objId;
  unsigned m_ofst;
};

//! @class  photonicsCmdWriteSaToRow
//! @brief  Photonics CMD: BitSIMD-V: Write SA to a row
class photonicsCmdWriteSaToRow : public photonicsCmd
{
public:
  photonicsCmdWriteSaToRow(PhotonicsCmdEnum cmdType, PhotonicsObjId objId, unsigned ofst)
    : photonicsCmd(cmdType), m_objId(objId), m_ofst(ofst) {}
  virtual ~photonicsCmdWriteSaToRow() {}
  virtual bool execute() override;
protected:
  PhotonicsObjId m_objId;
  unsigned m_ofst;
};

//! @class  photonicsCmdRRegOp : public photonicsCmd
//! @brief  Photonics CMD: BitSIMD-V: Row reg operations
class photonicsCmdRRegOp : public photonicsCmd
{
public:
  photonicsCmdRRegOp(PhotonicsCmdEnum cmdType, PhotonicsObjId objId, PhotonicsRowReg dest, bool val)
    : photonicsCmd(cmdType), m_objId(objId), m_dest(dest), m_val(val)
  {
    assert(cmdType == PhotonicsCmdEnum::RREG_SET);
  }
  photonicsCmdRRegOp(PhotonicsCmdEnum cmdType, PhotonicsObjId objId, PhotonicsRowReg dest, PhotonicsRowReg src1)
    : photonicsCmd(cmdType), m_objId(objId), m_dest(dest), m_src1(src1)
  {
    assert(cmdType == PhotonicsCmdEnum::RREG_MOV || cmdType == PhotonicsCmdEnum::RREG_NOT);
  }
  photonicsCmdRRegOp(PhotonicsCmdEnum cmdType, PhotonicsObjId objId, PhotonicsRowReg dest, PhotonicsRowReg src1, PhotonicsRowReg src2)
    : photonicsCmd(cmdType), m_objId(objId), m_dest(dest), m_src1(src1), m_src2(src2)
  {
  }
  photonicsCmdRRegOp(PhotonicsCmdEnum cmdType, PhotonicsObjId objId, PhotonicsRowReg dest, PhotonicsRowReg src1, PhotonicsRowReg src2, PhotonicsRowReg src3)
    : photonicsCmd(cmdType), m_objId(objId), m_dest(dest), m_src1(src1), m_src2(src2), m_src3(src3)
  {
    assert(cmdType == PhotonicsCmdEnum::RREG_MAJ || cmdType == PhotonicsCmdEnum::RREG_SEL);
  }
  virtual ~photonicsCmdRRegOp() {}
  virtual bool execute() override;
protected:
  PhotonicsObjId m_objId;
  PhotonicsRowReg m_dest;
  bool m_val = 0;
  PhotonicsRowReg m_src1 = PHOTONICS_RREG_NONE;
  PhotonicsRowReg m_src2 = PHOTONICS_RREG_NONE;
  PhotonicsRowReg m_src3 = PHOTONICS_RREG_NONE;
};

//! @class  photonicsCmdRRegRotate
//! @brief  Photonics CMD: BitSIMD-V: row reg rotate right by one step
class photonicsCmdRRegRotate : public photonicsCmd
{
public:
  photonicsCmdRRegRotate(PhotonicsCmdEnum cmdType, PhotonicsObjId objId, PhotonicsRowReg dest)
    : photonicsCmd(cmdType), m_objId(objId), m_dest(dest) {}
  virtual ~photonicsCmdRRegRotate() {}
  virtual bool execute() override;
protected:
  PhotonicsObjId m_objId;
  PhotonicsRowReg m_dest;
};

//! @class  photonicsCmdAnalogAAP
//! @brief  Photonics CMD: SIMDRAM: Analog based multi-row AP (activate-precharge) or AAP (activate-activate-precharge)
class photonicsCmdAnalogAAP : public photonicsCmd
{
public:
  photonicsCmdAnalogAAP(PhotonicsCmdEnum cmdType,
                  const std::vector<std::pair<PhotonicsObjId, unsigned>>& srcRows,
                  const std::vector<std::pair<PhotonicsObjId, unsigned>>& destRows = {})
    : photonicsCmd(cmdType), m_srcRows(srcRows), m_destRows(destRows)
  {
    assert(cmdType == PhotonicsCmdEnum::ROW_AP || cmdType == PhotonicsCmdEnum::ROW_AAP);
  }
  virtual ~photonicsCmdAnalogAAP() {}
  virtual bool execute() override;
protected:
  void printDebugInfo() const;
  std::vector<std::pair<PhotonicsObjId, unsigned>> m_srcRows;
  std::vector<std::pair<PhotonicsObjId, unsigned>> m_destRows;
};

#endif

