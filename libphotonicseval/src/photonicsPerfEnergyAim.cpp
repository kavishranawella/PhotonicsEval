// File: photonicsPerfEnergyAim.cc
// PHOTONICSeval Simulator - Performance Energy Models
// Copyright (c) 2024 University of Virginia
// This file is licensed under the MIT License.
// See the LICENSE file in the root of this repository for more details.

#include "photonicsPerfEnergyAim.h"
#include "photonicsCmd.h"
#include <cmath>
#include <cstdio>

// AiM adds a SIMD Multiplier and a Reduction Tree in each bank.
// The supported instructions are: MAC.
// For simplicity, the SIMD lane width is assumed to be determined by the GDL width of the HBM/DDR memory.
// NOTE: The energy model is approximated. 

//! @brief  Perf energy model of aim PHOTONICS for func1
photonicseval::perfEnergy
photonicsPerfEnergyAim::getPerfEnergyForFunc1(PhotonicsCmdEnum cmdType, const photonicsObjInfo& obj, const photonicsObjInfo& objDest) const
{
  double msRuntime = 0.0;
  double mjEnergy = 0.0;
  double msRead = 0.0;
  double msWrite = 0.0;
  double msCompute = 0.0;
  uint64_t totalOp = 0;
  switch (cmdType)
  {
    // Refer to AiM Paper (Table 2, Figure 5). OP Format: GRF = BANK +/* SRF
    case PhotonicsCmdEnum::ADD_SCALAR:
    case PhotonicsCmdEnum::MUL_SCALAR:
    case PhotonicsCmdEnum::AES_SBOX:
    case PhotonicsCmdEnum::AES_INVERSE_SBOX:
    case PhotonicsCmdEnum::POPCOUNT:
    case PhotonicsCmdEnum::ABS:
    case PhotonicsCmdEnum::SUB_SCALAR:
    case PhotonicsCmdEnum::DIV_SCALAR:
    case PhotonicsCmdEnum::AND_SCALAR:
    case PhotonicsCmdEnum::OR_SCALAR:
    case PhotonicsCmdEnum::XOR_SCALAR:
    case PhotonicsCmdEnum::XNOR_SCALAR:
    case PhotonicsCmdEnum::GT_SCALAR:
    case PhotonicsCmdEnum::LT_SCALAR:
    case PhotonicsCmdEnum::EQ_SCALAR:
    case PhotonicsCmdEnum::NE_SCALAR:
    case PhotonicsCmdEnum::MIN_SCALAR:
    case PhotonicsCmdEnum::MAX_SCALAR:
    case PhotonicsCmdEnum::SHIFT_BITS_L:
    case PhotonicsCmdEnum::SHIFT_BITS_R:
    default:
      printf("PHOTONICS-Warning: Perf energy model not available for PHOTONICS command %s\n", photonicsCmd::getName(cmdType, "").c_str());
      break;
  }

  return photonicseval::perfEnergy(msRuntime, mjEnergy, msRead, msWrite, msCompute, totalOp);
}

//! @brief  Perf energy model of aim for func2
photonicseval::perfEnergy
photonicsPerfEnergyAim::getPerfEnergyForFunc2(PhotonicsCmdEnum cmdType, const photonicsObjInfo& obj, const photonicsObjInfo& objSrc2, const photonicsObjInfo& objDest) const
{
  double msRuntime = 0.0;
  double mjEnergy = 0.0;
  double msRead = 0.0;
  double msWrite = 0.0;
  double msCompute = 0.0;
  uint64_t totalOp = 0;
  switch (cmdType)
  {
    // Refer to Aquabolt Paper (Table 2, Figure 5). OP Format: GRF = BANK +/* GRF
    case PhotonicsCmdEnum::ADD:
    case PhotonicsCmdEnum::MUL:
    case PhotonicsCmdEnum::SCALED_ADD:
    case PhotonicsCmdEnum::DIV:
    case PhotonicsCmdEnum::SUB:
    case PhotonicsCmdEnum::AND:
    case PhotonicsCmdEnum::OR:
    case PhotonicsCmdEnum::XOR:
    case PhotonicsCmdEnum::XNOR:
    case PhotonicsCmdEnum::GT:
    case PhotonicsCmdEnum::LT:
    case PhotonicsCmdEnum::EQ:
    case PhotonicsCmdEnum::NE:
    case PhotonicsCmdEnum::MIN:
    case PhotonicsCmdEnum::MAX:
    default:
      printf("PHOTONICS-Warning: Unsupported for AiM: %s\n", photonicsCmd::getName(cmdType, "").c_str());
      break;
  }

  return photonicseval::perfEnergy(msRuntime, mjEnergy, msRead, msWrite, msCompute, totalOp);
}

//! @brief  Perf energy model of aim PHOTONICS for reduction sum
photonicseval::perfEnergy
photonicsPerfEnergyAim::getPerfEnergyForReduction(PhotonicsCmdEnum cmdType, const photonicsObjInfo& obj, unsigned numPass) const
{
  double msRuntime = 0.0;
  double mjEnergy = 0.0;
  double msRead = 0.0;
  double msWrite = 0.0;
  double msCompute = 0.0;
  uint64_t totalOp = 0;

  switch (cmdType) {
    case PhotonicsCmdEnum::REDSUM:
    case PhotonicsCmdEnum::REDSUM_RANGE:
    case PhotonicsCmdEnum::REDMIN:
    case PhotonicsCmdEnum::REDMIN_RANGE:
    case PhotonicsCmdEnum::REDMAX:
    case PhotonicsCmdEnum::REDMAX_RANGE:
    default:
      printf("PHOTONICS-Warning: Unsupported for AiM: %s\n", photonicsCmd::getName(cmdType, "").c_str());
      break;
  }
  return photonicseval::perfEnergy(msRuntime, mjEnergy, msRead, msWrite, msCompute, totalOp);
}

//! @brief  Perf energy model of aim for broadcast
photonicseval::perfEnergy
photonicsPerfEnergyAim::getPerfEnergyForBroadcast(PhotonicsCmdEnum cmdType, const photonicsObjInfo& obj) const
{
  double msRuntime = 0.0;
  double mjEnergy = 0.0;
  double msRead = 0.0;
  double msWrite = 0.0;
  double msCompute = 0.0;
  uint64_t totalOp = 0;

  return photonicseval::perfEnergy(msRuntime, mjEnergy, msRead, msWrite, msCompute, totalOp);
}

//! @brief  Perf energy model of aim for rotate
photonicseval::perfEnergy
photonicsPerfEnergyAim::getPerfEnergyForRotate(PhotonicsCmdEnum cmdType, const photonicsObjInfo& obj) const
{
  double msRuntime = 0.0;
  double mjEnergy = 0.0;
  double msRead = 0.0;
  double msWrite = 0.0;
  double msCompute = 0.0;
  uint64_t totalOp = 0;
  printf("PHOTONICS-Warning: Unsupported for AiM: %s\n", photonicsCmd::getName(cmdType, "").c_str());

  return photonicseval::perfEnergy(msRuntime, mjEnergy, msRead, msWrite, msCompute, totalOp);
}

photonicseval::perfEnergy photonicsPerfEnergyAim::getPerfEnergyForMac(PhotonicsCmdEnum cmdType, const photonicsObjInfo &obj) const
{
  // NumPass is always 1 for MAC operation in AiM. User really needs to make sure that this holds true.
  // Buffer read time is `tCAS - m_tGDL` based on following reasoning:
  // 1. tCAS = cycles required to data available at the I/O interface after a read command.
  // 2. m_tGDL = cycles required for two consecutive read commands to the same bank.
  // Hence, the time to read data from the global AiM buffer to the bank interface is `tCAS - m_tGDL`.
  // User may wonder why buffer read time is not multiplied by number of banks per chip. This is because according the AiM paper mentions that buffer is n-way fanout to n banks in the same chip.
  // AiM paper mentions accumulation reduction tree requires 4 cycles after the multiplier. Hence, the compute time for accumulation is `4 * tCK`.
  // TODO: Energy model
  double msRuntime = 0.0;
  double mjEnergy = 0.0;
  double msRead = 0.0;
  double msWrite = 0.0;
  double msCompute = 0.0;
  uint64_t totalOp = 0;
  unsigned bitsPerElement = obj.getBitsPerElement(PhotonicsBitWidth::ACTUAL);
  unsigned maxElementsPerRegion = obj.getMaxElementsPerRegion();
  unsigned numCore = obj.getNumCoreAvailable();
  unsigned elementsPerCore = std::ceil(obj.getNumElements() * 1.0 / numCore);
  unsigned gdlItr = std::ceil(elementsPerCore * bitsPerElement * 1.0 / m_GDLWidth);
  unsigned numBankPerChip = numCore / m_numChipsPerRank;

  photonicseval::perfEnergy perfEnergyBT = getPerfEnergyForBytesTransfer(PhotonicsCmdEnum::COPY_D2H, (bitsPerElement * numCore) / 8);
  
  msRead = m_tACT + m_tPRE + (m_tCAS - m_tGDL) * gdlItr;
  msWrite = perfEnergyBT.m_msRuntime;
  msCompute = (gdlItr * m_tGDL + 4 * m_tCK * gdlItr);
  msRuntime = msRead + msWrite + msCompute;
  mjEnergy = ((m_eACT + m_ePRE) + (maxElementsPerRegion * m_aquaboltArithmeticEnergy)) * numCore;
  mjEnergy += m_eR * numBankPerChip * m_numRanks * gdlItr; // Energy for reading data from local row buffer to global row buffer
  mjEnergy += perfEnergyBT.m_mjEnergy;
  mjEnergy += m_pBChip * m_numChipsPerRank * m_numRanks * msRuntime;
  totalOp = obj.getNumElements() * 2;
  return photonicseval::perfEnergy(msRuntime, mjEnergy, msRead, msWrite, msCompute, totalOp);
}
