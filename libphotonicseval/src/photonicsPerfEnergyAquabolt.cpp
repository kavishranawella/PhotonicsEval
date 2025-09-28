// File: photonicsPerfEnergyAquabolt.cc
// PHOTONICSeval Simulator - Performance Energy Models
// Copyright (c) 2024 University of Virginia
// This file is licensed under the MIT License.
// See the LICENSE file in the root of this repository for more details.

#include "photonicsPerfEnergyAquabolt.h"
#include "photonicsCmd.h"
#include <cstdio>
#include <cmath>

// Aquabolt adds a SIMD FPU shared between two banks, with only one bank accessing it at a time.
// The supported FPU instructions are: ADD, MUL, MAC, and RELU. However, RELU is currently not implemented in the simulator.
// This model assumes that each FPU operation (ADD, MUL, MAC, or RELU) takes `tCCD_L` cycles to execute.
// Additionally, for simplicity, the SIMD lane width is assumed to be determined by the GDL width of the HBM/DDR memory. 
// This analytical model has been validated against the Aquabolt for vector addition and multiplication using a 100M-element vector of 16-bit integers. 
// The model demonstrates a 1.5x speedup compared to the original Aquabolt.
// NOTE: The energy model is approximated. 

//! @brief  Perf energy model of aquabolt PHOTONICS for func1
photonicseval::perfEnergy
photonicsPerfEnergyAquabolt::getPerfEnergyForFunc1(PhotonicsCmdEnum cmdType, const photonicsObjInfo& obj, const photonicsObjInfo& objDest) const
{
  double msRuntime = 0.0;
  double mjEnergy = 0.0;
  double msRead = 0.0;
  double msWrite = 0.0;
  double msCompute = 0.0;
  unsigned numPass = obj.getMaxNumRegionsPerCore();
  unsigned bitsPerElement = obj.getBitsPerElement(PhotonicsBitWidth::ACTUAL);
  unsigned numCores = obj.getNumCoreAvailable();
  unsigned maxElementsPerRegion = obj.getMaxElementsPerRegion();
  unsigned numberOfOperationPerElement = std::ceil(bitsPerElement * 1.0 / m_aquaboltFPUBitWidth);
  unsigned elementsPerCore = std::ceil(obj.getNumElements() * 1.0 / numCores);
  unsigned minElementPerRegion = elementsPerCore > maxElementsPerRegion ? elementsPerCore - (maxElementsPerRegion * (numPass - 1)) : elementsPerCore;
  unsigned maxGDLItr = std::ceil(maxElementsPerRegion * bitsPerElement * 1.0 / m_GDLWidth);
  unsigned minGDLItr = std::ceil(minElementPerRegion * bitsPerElement * 1.0 / m_GDLWidth);
  double aquaboltCoreCycle = m_tGDL;
  unsigned numActPre = std::ceil(maxElementsPerRegion * bitsPerElement * 1.0 / (16 * 256));
  uint64_t totalOp = 0;
  unsigned numBankPerChip = numCores / m_numChipsPerRank;

  switch (cmdType)
  {
    // Refer to Aquabolt Paper (Table 2, Figure 5). OP Format: GRF = BANK +/* SRF
    // Aquabolt has 16 16-bit vector registers (GRF) per PHOTONICS core.
    // As a result, depending on the bitsPerElement and columns per bank row, same row may be opened multiple times -- this is calculated as numActPre.
    case PhotonicsCmdEnum::ADD_SCALAR:
    case PhotonicsCmdEnum::MUL_SCALAR:
    { 
      msRead = (m_tACT + m_tPRE) * numPass * numActPre;
      msWrite = (m_tACT + m_tPRE) * numPass * numActPre;
      msCompute = (minGDLItr * aquaboltCoreCycle * numberOfOperationPerElement) + ((maxGDLItr * aquaboltCoreCycle * numberOfOperationPerElement) * (numPass - 1));
      msRuntime = msRead + msWrite + msCompute;
      mjEnergy = ((m_eACT + m_ePRE) * numActPre * 2 + (maxElementsPerRegion * m_aquaboltArithmeticEnergy * numberOfOperationPerElement)) * numCores * (numPass - 1);
      mjEnergy += ((m_eACT + m_ePRE) * numActPre * 2 +  (minElementPerRegion * m_aquaboltArithmeticEnergy * numberOfOperationPerElement)) * numCores;
      mjEnergy += (m_eR * maxGDLItr * (numPass-1) * numBankPerChip * m_numRanks + (m_eR * minGDLItr * numBankPerChip * m_numRanks));
      mjEnergy += (m_eW * maxGDLItr * (numPass-1) * numBankPerChip * m_numRanks + (m_eW * minGDLItr * numBankPerChip * m_numRanks));
      mjEnergy += m_pBChip * m_numChipsPerRank * m_numRanks * msRuntime;
      totalOp = obj.getNumElements();
      break;
    }
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

//! @brief  Perf energy model of aquabolt PHOTONICS for func2
photonicseval::perfEnergy
photonicsPerfEnergyAquabolt::getPerfEnergyForFunc2(PhotonicsCmdEnum cmdType, const photonicsObjInfo& obj, const photonicsObjInfo& objSrc2, const photonicsObjInfo& objDest) const
{
  double msRuntime = 0.0;
  double mjEnergy = 0.0;
  double msRead = 0.0;
  double msWrite = 0.0;
  double msCompute = 0.0;
  unsigned numPass = obj.getMaxNumRegionsPerCore();
  unsigned bitsPerElement = obj.getBitsPerElement(PhotonicsBitWidth::ACTUAL);
  unsigned numCoresUsed = obj.getNumCoreAvailable();
  unsigned maxElementsPerRegion = obj.getMaxElementsPerRegion();
  unsigned elementsPerCore = std::ceil(obj.getNumElements() * 1.0 / numCoresUsed);
  unsigned minElementPerRegion = elementsPerCore > maxElementsPerRegion ? elementsPerCore - (maxElementsPerRegion * (numPass - 1)) : elementsPerCore;
  unsigned maxGDLItr = std::ceil(maxElementsPerRegion * bitsPerElement * 1.0 / m_GDLWidth);
  unsigned minGDLItr = std::ceil(minElementPerRegion * bitsPerElement * 1.0 / m_GDLWidth);
  double aquaboltCoreCycle = m_tGDL;
  unsigned numActPre = std::ceil(maxElementsPerRegion * bitsPerElement * 1.0 / (8 * 256));
  uint64_t totalOp = 0;
  unsigned numBankPerChip = numCoresUsed / m_numChipsPerRank;
  switch (cmdType)
  {
    // Refer to Aquabolt Paper (Table 2, Figure 5). OP Format: GRF = BANK +/* GRF
    case PhotonicsCmdEnum::ADD:
    case PhotonicsCmdEnum::MUL:
    {
      unsigned numberOfOperationPerElement = std::ceil(bitsPerElement * 1.0 / m_aquaboltFPUBitWidth);
      msRead = (2 * (m_tACT + m_tPRE) * numPass * numActPre) + (maxGDLItr * m_tGDL * (numPass - 1)) + (minGDLItr * m_tGDL);
      msWrite = ((m_tACT + m_tPRE) * numPass * numActPre) + (maxGDLItr * m_tGDL * (numPass - 1)) + (minGDLItr * m_tGDL);
      msCompute = (maxGDLItr * numberOfOperationPerElement * aquaboltCoreCycle) * (numPass - 1);
      msCompute += (minGDLItr * numberOfOperationPerElement * aquaboltCoreCycle);
      msRuntime = msRead + msWrite + msCompute;
      mjEnergy = (((m_eACT + m_ePRE) * 3 * numActPre) + ((maxElementsPerRegion * m_aquaboltArithmeticEnergy * numberOfOperationPerElement))) * numCoresUsed * (numPass - 1);
      mjEnergy += (((m_eACT + m_ePRE) * 3 * numActPre) + ((minElementPerRegion * m_aquaboltArithmeticEnergy * numberOfOperationPerElement))) * numCoresUsed;
      mjEnergy += (m_eR * maxGDLItr * 2 * (numPass-1) * numBankPerChip * m_numRanks + (m_eR * 2 * minGDLItr * numBankPerChip * m_numRanks));
      mjEnergy += (m_eW * maxGDLItr * (numPass-1) * numBankPerChip * m_numRanks + (m_eW * minGDLItr * numBankPerChip * m_numRanks));
      mjEnergy += m_pBChip * m_numChipsPerRank * m_numRanks * msRuntime;
      totalOp = obj.getNumElements();
      break;
    }
    case PhotonicsCmdEnum::SCALED_ADD:
    {
      /**
       * Performs a multiply-add operation on rows in DRAM.
       *
       * This command executes the following steps:
       * 1. Multiply the elements of a source row by a scalar value.
       * 2. Add the result of the multiplication to the elements of another row.
       * 3. Write the final result back to a row in DRAM.
       *
       * Performance Optimizations:
       * - While performing the multiplication, the next row to be added can be fetched without any additional overhead.
       * - During the addition, the next row to be multiplied can be fetched concurrently.
       *
       * As a result, only one read operation is necessary for the entire pass.
      */
      // OP Format: GRF = BANK * SRF; GRF = BANK + GRF 
      unsigned numberOfOperationPerElement = std::ceil(bitsPerElement * 1.0 / m_aquaboltFPUBitWidth) * 2; // multiplying by 2 as one addition and one multiplication is needed
      msRead = (m_tACT + m_tPRE) * numPass * numActPre * 2;
      msWrite = (m_tACT + m_tPRE) * numPass * numActPre + (maxGDLItr * m_tGDL * (numPass - 1)) + (minGDLItr * m_tGDL);
      msCompute = (maxGDLItr * aquaboltCoreCycle * numberOfOperationPerElement) * (numPass - 1);
      msCompute += (minGDLItr * aquaboltCoreCycle * numberOfOperationPerElement);
      msRuntime = msRead + msWrite + msCompute;
      mjEnergy = (((m_eACT + m_ePRE) * 3 * numActPre) + ((maxElementsPerRegion * m_aquaboltArithmeticEnergy * numberOfOperationPerElement))) * numCoresUsed * (numPass - 1);
      mjEnergy += (((m_eACT + m_ePRE) * 3 * numActPre) + ((minElementPerRegion * m_aquaboltArithmeticEnergy * numberOfOperationPerElement))) * numCoresUsed;
      mjEnergy += (m_eR * maxGDLItr * 2 * (numPass-1) * numBankPerChip * m_numRanks + (m_eR * 2 * minGDLItr * numBankPerChip * m_numRanks));
      mjEnergy += (m_eW * maxGDLItr * (numPass-1) * numBankPerChip * m_numRanks + (m_eW * minGDLItr * numBankPerChip * m_numRanks));
      mjEnergy += m_pBChip * m_numChipsPerRank * m_numRanks * msRuntime;
      totalOp = obj.getNumElements() * 2;
      break;
    }
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
      printf("PHOTONICS-Warning: Unsupported for Aquabolt: %s\n", photonicsCmd::getName(cmdType, "").c_str());
      break;
  }

  return photonicseval::perfEnergy(msRuntime, mjEnergy, msRead, msWrite, msCompute, totalOp);
}

//! @brief  Perf energy model of aquabolt PHOTONICS for reduction sum
photonicseval::perfEnergy
photonicsPerfEnergyAquabolt::getPerfEnergyForReduction(PhotonicsCmdEnum cmdType, const photonicsObjInfo& obj, unsigned numPass) const
{
  double msRuntime = 0.0;
  double mjEnergy = 0.0;
  double msRead = 0.0;
  double msWrite = 0.0;
  double msCompute = 0.0;
  unsigned bitsPerElement = obj.getBitsPerElement(PhotonicsBitWidth::ACTUAL);
  unsigned maxElementsPerRegion = obj.getMaxElementsPerRegion();
  unsigned numCore = obj.getNumCoreAvailable();
  double cpuTDP = 200; // W; AMD EPYC 9124 16 core
  unsigned elementsPerCore = std::ceil(obj.getNumElements() * 1.0 / numCore);
  unsigned minElementPerRegion = elementsPerCore > maxElementsPerRegion ? elementsPerCore - (maxElementsPerRegion * (numPass - 1)) : elementsPerCore;
  unsigned maxGDLItr = std::ceil(maxElementsPerRegion * bitsPerElement * 1.0 / m_GDLWidth);
  unsigned minGDLItr = std::ceil(minElementPerRegion * bitsPerElement * 1.0 / m_GDLWidth);
  unsigned numberOfOperationPerElement = std::ceil(bitsPerElement * 1.0 / m_aquaboltFPUBitWidth);
  double aquaboltCoreCycle = m_tGDL;
  uint64_t totalOp = 0;

  switch (cmdType) {
    case PhotonicsCmdEnum::REDSUM:
    case PhotonicsCmdEnum::REDSUM_RANGE:
    {
      msRead = (m_tR * numPass) + m_tGDL * numPass;
      msCompute = (maxGDLItr * aquaboltCoreCycle * numberOfOperationPerElement) * (numPass - 1);
      msCompute += (m_tR + (minGDLItr * aquaboltCoreCycle * numberOfOperationPerElement));
      msRuntime = msRead + msWrite + msCompute;
      // Refer to fulcrum documentation
      mjEnergy = (m_eAP + ((m_eR * maxGDLItr) + (maxElementsPerRegion * m_aquaboltArithmeticEnergy * numberOfOperationPerElement))) * numPass * numCore;
      // reduction for all regions
      double aggregateMs = static_cast<double>(numCore) / (3200000 * 16);
      msRuntime += aggregateMs;
      mjEnergy += aggregateMs * cpuTDP;
      mjEnergy += m_pBChip * m_numChipsPerRank * m_numRanks * msRuntime;
      totalOp = obj.getNumElements();
      break;
    }
    case PhotonicsCmdEnum::REDMIN:
    case PhotonicsCmdEnum::REDMIN_RANGE:
    case PhotonicsCmdEnum::REDMAX:
    case PhotonicsCmdEnum::REDMAX_RANGE:
    default:
      printf("PHOTONICS-Warning: Unsupported for Aquabolt: %s\n", photonicsCmd::getName(cmdType, "").c_str());
      break;
  }
  return photonicseval::perfEnergy(msRuntime, mjEnergy, msRead, msWrite, msCompute, totalOp);
}

//! @brief  Perf energy model of aquabolt PHOTONICS for broadcast
photonicseval::perfEnergy
photonicsPerfEnergyAquabolt::getPerfEnergyForBroadcast(PhotonicsCmdEnum cmdType, const photonicsObjInfo& obj) const
{
  double msRuntime = 0.0;
  double mjEnergy = 0.0;
  double msRead = 0.0;
  double msWrite = 0.0;
  double msCompute = 0.0;
  unsigned numPass = obj.getMaxNumRegionsPerCore();
  unsigned bitsPerElement = obj.getBitsPerElement(PhotonicsBitWidth::ACTUAL);
  unsigned maxElementsPerRegion = obj.getMaxElementsPerRegion();
  unsigned numCore = obj.getNumCoreAvailable();
  uint64_t totalOp = 0;

  unsigned elementsPerCore = std::ceil(obj.getNumElements() * 1.0 / numCore);
  unsigned minElementPerRegion = elementsPerCore > maxElementsPerRegion ? elementsPerCore - (maxElementsPerRegion * (numPass - 1)) : elementsPerCore;
  unsigned maxGDLItr = std::ceil(maxElementsPerRegion * bitsPerElement * 1.0 / m_GDLWidth);
  unsigned minGDLItr = std::ceil(minElementPerRegion * bitsPerElement * 1.0 / m_GDLWidth);
  msWrite = (m_tW + maxGDLItr * m_tGDL) * (numPass - 1);
  msWrite += (m_tW + minGDLItr * m_tGDL);
  msRuntime = msRead + msWrite + msCompute;
  mjEnergy = (m_eAP + m_eR * maxGDLItr) * (numPass - 1) * numCore;
  mjEnergy += (m_eAP + m_eR * minGDLItr) * numCore;
  mjEnergy += m_pBChip * m_numChipsPerRank * m_numRanks * msRuntime;

  return photonicseval::perfEnergy(msRuntime, mjEnergy, msRead, msWrite, msCompute, totalOp);
}

//! @brief  Perf energy model of aquabolt PHOTONICS for rotate
photonicseval::perfEnergy
photonicsPerfEnergyAquabolt::getPerfEnergyForRotate(PhotonicsCmdEnum cmdType, const photonicsObjInfo& obj) const
{
  double msRuntime = 0.0;
  double mjEnergy = 0.0;
  double msRead = 0.0;
  double msWrite = 0.0;
  double msCompute = 0.0;
  uint64_t totalOp = 0;
  printf("PHOTONICS-Warning: Unsupported for Aquabolt: %s\n", photonicsCmd::getName(cmdType, "").c_str());

  return photonicseval::perfEnergy(msRuntime, mjEnergy, msRead, msWrite, msCompute, totalOp);
}

