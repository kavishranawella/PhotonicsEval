// File: photonicsPerfEnergyFulcrum.cc
// PHOTONICSeval Simulator - Performance Energy Models
// Copyright (c) 2024 University of Virginia
// This file is licensed under the MIT License.
// See the LICENSE file in the root of this repository for more details.

#include "photonicsPerfEnergyFulcrum.h"
#include "photonicsCmd.h"
#include <cmath>
#include <cstdint>
#include <cstdio>


//! @brief  Perf energy model of Fulcrum for func1
photonicseval::perfEnergy
photonicsPerfEnergyFulcrum::getPerfEnergyForFunc1(PhotonicsCmdEnum cmdType, const photonicsObjInfo& obj, const photonicsObjInfo& objDest) const
{

  // Fulcrum utilizes three walkers: two for input operands and one for the output operand.
  // For instructions that operate on a single operand, the next operand is fetched by the walker.
  // Consequently, only one row read operation is required in this case.
  // Additionally, using the walker-renaming technique (refer to the Fulcrum paper for details),
  // the write operation is also pipelined. Thus, only one row write operation is needed.

  double msRuntime = 0.0;
  double mjEnergy = 0.0;
  double msRead = 0.0;
  double msWrite = 0.0;
  double msALU = 0.0;
  uint64_t totalOp = 0;
  unsigned numPass = obj.getMaxNumRegionsPerCore();
  unsigned bitsPerElement = obj.getBitsPerElement(PhotonicsBitWidth::ACTUAL);
  if (cmdType == PhotonicsCmdEnum::CONVERT_TYPE) {
    // for type conversion, ALU parallelism is determined by the wider data type
    bitsPerElement = std::max(bitsPerElement, objDest.getBitsPerElement(PhotonicsBitWidth::ACTUAL));
  }
  unsigned numCores =  obj.isLoadBalanced() ? obj.getNumCoreAvailable() : obj.getNumCoresUsed();
  unsigned maxElementsPerRegion = obj.getMaxElementsPerRegion();
  double numberOfALUOperationPerElement = ((double)bitsPerElement / m_fulcrumAluBitWidth);
  unsigned minElementPerRegion = obj.isLoadBalanced() ? (std::ceil(obj.getNumElements() * 1.0 / numCores) - (maxElementsPerRegion * (numPass - 1))) : maxElementsPerRegion;
  switch (cmdType)
  {
    case PhotonicsCmdEnum::COPY_O2O:
    {
      msRead = m_tR * numPass;
      msWrite = m_tW * numPass;
      msRuntime = msRead + msWrite + msALU;
      mjEnergy = numPass * numCores * m_eAP * 2;
      mjEnergy += m_pBChip * m_numChipsPerRank * m_numRanks * msRuntime;
      break;
    }
    case PhotonicsCmdEnum::POPCOUNT:
    {
      double msPopCount = (m_fulcrumAddLatency * 11 +  m_fulcrumMulLatency); // 4 shifts, 4 ands, 3 add/sub, 1 mul
      msRead = m_tR;
      msWrite = m_tW;
      msALU = ((maxElementsPerRegion * msPopCount * numberOfALUOperationPerElement) * (numPass - 1)) + (minElementPerRegion * msPopCount * numberOfALUOperationPerElement);
      msRuntime = msRead + msWrite + msALU;
      double energyArithmetic = (((maxElementsPerRegion - 1) * 2 *  m_fulcrumShiftEnergy * numberOfALUOperationPerElement) + (maxElementsPerRegion * m_fulcrumMulEnergy * numberOfALUOperationPerElement)) * (numPass - 1);
      energyArithmetic += (((minElementPerRegion - 1) * 2 *  m_fulcrumShiftEnergy * numberOfALUOperationPerElement) + (minElementPerRegion * m_fulcrumMulEnergy * numberOfALUOperationPerElement));
      double energyLogical = ((m_eAP * 2) + (((maxElementsPerRegion - 1) * 2 *  m_fulcrumShiftEnergy) + (maxElementsPerRegion * m_fulcrumAddEnergy * 11))) * (numPass - 1);
      energyLogical += (m_eAP * 2) + (((minElementPerRegion - 1) * 2 *  m_fulcrumShiftEnergy) + (minElementPerRegion * m_fulcrumAddEnergy * 11));
      mjEnergy = (energyArithmetic + energyLogical) * numCores ;
      mjEnergy += m_pBChip * m_numChipsPerRank * m_numRanks * msRuntime;
      totalOp = obj.getNumElements() * 12;
      break;
    }
    case PhotonicsCmdEnum::BIT_SLICE_EXTRACT:
    case PhotonicsCmdEnum::BIT_SLICE_INSERT:
    {
      if (cmdType == PhotonicsCmdEnum::BIT_SLICE_EXTRACT) {
        // Assume one ALU cycle to do this for now
        // numberOfALUOperationPerElement *= 2; // 1 shift, 1 and
      } else if (cmdType == PhotonicsCmdEnum::BIT_SLICE_INSERT) {
        // Assume one ALU cycle to do this for now
        // numberOfALUOperationPerElement *= 5; // 2 shifts, 1 not, 1 and, 1 or
      }
      msRead = m_tR;
      msWrite = m_tW;
      msALU = ((maxElementsPerRegion * m_fulcrumAddLatency * numberOfALUOperationPerElement) * (numPass - 1)) + (minElementPerRegion * m_fulcrumAddLatency * numberOfALUOperationPerElement);
      msRuntime = msRead + msWrite + msALU;
      double energyLogical = ((m_eAP * 2) + ((maxElementsPerRegion - 1) * 2 *  m_fulcrumShiftEnergy) + (maxElementsPerRegion * m_fulcrumAddEnergy * numberOfALUOperationPerElement)) * (numPass - 1);
      energyLogical += ((m_eAP * 2) + (minElementPerRegion - 1) * 2 *  m_fulcrumShiftEnergy) + (minElementPerRegion * m_fulcrumAddEnergy * numberOfALUOperationPerElement);
      mjEnergy = energyLogical * numCores;
      mjEnergy += m_pBChip * m_numChipsPerRank * m_numRanks * msRuntime;
      totalOp = obj.getNumElements();
      break;
    }
    case PhotonicsCmdEnum::MUL_SCALAR:
    case PhotonicsCmdEnum::DIV_SCALAR:
    {
      msRead = m_tR;
      msWrite = m_tW;
      msALU = ((maxElementsPerRegion * m_fulcrumMulLatency * numberOfALUOperationPerElement) * (numPass - 1)) + (minElementPerRegion * m_fulcrumMulLatency * numberOfALUOperationPerElement);
      msRuntime = msRead + msWrite + msALU;
      mjEnergy = (numPass - 1) * numCores * ((m_eAP * 2) + ((maxElementsPerRegion - 1) * 2 *  m_fulcrumShiftEnergy) + (maxElementsPerRegion * m_fulcrumMulEnergy * numberOfALUOperationPerElement));
      mjEnergy += numCores * ((m_eAP * 2) + ((minElementPerRegion - 1) * 2 *  m_fulcrumShiftEnergy) + (minElementPerRegion * m_fulcrumMulEnergy * numberOfALUOperationPerElement));
      mjEnergy += m_pBChip * m_numChipsPerRank * m_numRanks * msRuntime;
      totalOp = obj.getNumElements();
      break;
    }
    case PhotonicsCmdEnum::ABS:
    case PhotonicsCmdEnum::CONVERT_TYPE:
    case PhotonicsCmdEnum::ADD_SCALAR:
    case PhotonicsCmdEnum::SUB_SCALAR:
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
    {
      msRead = m_tR;
      msWrite = m_tW;
      msALU = ((maxElementsPerRegion * m_fulcrumAddLatency * numberOfALUOperationPerElement) * (numPass - 1)) + (minElementPerRegion * m_fulcrumAddLatency * numberOfALUOperationPerElement);
      msRuntime = msRead + msWrite + msALU;
      mjEnergy = (numPass - 1) * numCores * ((m_eAP * 2) + ((maxElementsPerRegion - 1) * 2 *  m_fulcrumShiftEnergy) + (maxElementsPerRegion * m_fulcrumAddEnergy * numberOfALUOperationPerElement));
      mjEnergy = numCores * ((m_eAP * 2) + ((minElementPerRegion - 1) * 2 *  m_fulcrumShiftEnergy) + (minElementPerRegion * m_fulcrumAddEnergy * numberOfALUOperationPerElement));
      mjEnergy += m_pBChip * m_numChipsPerRank * m_numRanks * msRuntime; 
      totalOp = obj.getNumElements();
      break;
    }
    case PhotonicsCmdEnum::AES_SBOX:
    case PhotonicsCmdEnum::AES_INVERSE_SBOX:
      msRuntime = 1e10;
      mjEnergy = 999999999.9;
      printf("PHOTONICS-Warning: Perf energy model not available for PHOTONICS command %s\n", photonicsCmd::getName(cmdType, "").c_str());
      break;
    default:
      printf("PHOTONICS-Warning: Perf energy model not available for PHOTONICS command %s\n", photonicsCmd::getName(cmdType, "").c_str());
      break;
  }

  return photonicseval::perfEnergy(msRuntime, mjEnergy, msRead, msWrite, msALU, totalOp);
}

//! @brief  Perf energy model of Fulcrum for func2
photonicseval::perfEnergy
photonicsPerfEnergyFulcrum::getPerfEnergyForFunc2(PhotonicsCmdEnum cmdType, const photonicsObjInfo& obj, const photonicsObjInfo& objSrc2, const photonicsObjInfo& objDest) const
{
  double msRuntime = 0.0;
  double mjEnergy = 0.0;
  double msRead = 0.0;
  double msWrite = 0.0;
  double msALU = 0.0;
  uint64_t totalOp = 0;
  unsigned numPass = obj.getMaxNumRegionsPerCore();
  unsigned bitsPerElement = obj.getBitsPerElement(PhotonicsBitWidth::ACTUAL);
  unsigned numCoresUsed = obj.isLoadBalanced() ? obj.getNumCoreAvailable() : obj.getNumCoresUsed();
  unsigned maxElementsPerRegion = obj.getMaxElementsPerRegion();
  unsigned minElementPerRegion = obj.isLoadBalanced() ? (std::ceil(obj.getNumElements() * 1.0 / obj.getNumCoreAvailable()) - (maxElementsPerRegion * (numPass - 1))) : maxElementsPerRegion;
  double numberOfALUOperationPerElement = ((double)bitsPerElement / m_fulcrumAluBitWidth);
  switch (cmdType)
  {
    case PhotonicsCmdEnum::MUL:
    case PhotonicsCmdEnum::DIV:
    {
      msRead = 2 * m_tR * numPass;
      msWrite = m_tW * numPass;
      msALU = (maxElementsPerRegion * numberOfALUOperationPerElement * m_fulcrumMulLatency * (numPass - 1)) +  (minElementPerRegion * numberOfALUOperationPerElement * m_fulcrumMulLatency);
      msRuntime = msRead + msWrite + msALU;
      mjEnergy = numCoresUsed * (numPass - 1) * ((m_eAP * 3) + ((maxElementsPerRegion - 1) * 3 *  m_fulcrumShiftEnergy) + (maxElementsPerRegion * m_fulcrumMulEnergy * numberOfALUOperationPerElement));
      mjEnergy += numCoresUsed * ((m_eAP * 3) + ((minElementPerRegion - 1) * 3 *  m_fulcrumShiftEnergy) + ((minElementPerRegion) * m_fulcrumMulEnergy * numberOfALUOperationPerElement));
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
       * - Total execution time for one region of multiplication and addition >>>> reading/writing three DRAM rows as a result using walker renaming, row write is also pipelined
       *
       * As a result, only one read operation and one write operation is necessary for the entire pass.
      */
      msRead = m_tR;
      msWrite = m_tW;
      msALU = (maxElementsPerRegion * numberOfALUOperationPerElement * (m_fulcrumAddLatency + m_fulcrumMulLatency ) * (numPass - 1)) +  (minElementPerRegion * numberOfALUOperationPerElement * (m_fulcrumAddLatency + m_fulcrumMulLatency ));
      msRuntime = msRead + msWrite + msALU;
      mjEnergy = numCoresUsed * (numPass - 1) * ((m_eAP * 2) + ((maxElementsPerRegion - 1) * 2 *  m_fulcrumShiftEnergy * 2) + (maxElementsPerRegion * (m_fulcrumAddEnergy + m_fulcrumMulEnergy) * numberOfALUOperationPerElement));
      mjEnergy += numCoresUsed * ((m_eAP * 2) + ((minElementPerRegion - 1) * 2 *  m_fulcrumShiftEnergy * 2) + (minElementPerRegion * (m_fulcrumAddEnergy + m_fulcrumMulEnergy) * numberOfALUOperationPerElement));
      mjEnergy += m_pBChip * m_numChipsPerRank * m_numRanks * msRuntime;
      totalOp = obj.getNumElements() * 2;
      break;
    }
    case PhotonicsCmdEnum::ADD:
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
    case PhotonicsCmdEnum::COND_BROADCAST: // read from bool and dest, write to dest
    {
      msRead = 2 * m_tR * numPass;
      msWrite = m_tW * numPass;
      msALU = (maxElementsPerRegion * numberOfALUOperationPerElement * m_fulcrumAddLatency * (numPass - 1)) +  (minElementPerRegion * numberOfALUOperationPerElement * m_fulcrumAddLatency);
      msRuntime = msRead + msWrite + msALU;
      mjEnergy = numCoresUsed * (numPass - 1) * ((m_eAP * 3) + ((maxElementsPerRegion - 1) * 3 *  m_fulcrumShiftEnergy) + (maxElementsPerRegion * m_fulcrumAddEnergy * numberOfALUOperationPerElement));
      mjEnergy += numCoresUsed * ((m_eAP * 3) + ((minElementPerRegion - 1) * 3 *  m_fulcrumShiftEnergy) + (minElementPerRegion * m_fulcrumAddEnergy * numberOfALUOperationPerElement));
      mjEnergy += m_pBChip * m_numChipsPerRank * m_numRanks * msRuntime;
      totalOp = obj.getNumElements();
      break;
    }
    default:
      printf("PHOTONICS-Warning: Perf energy model not available for PHOTONICS command %s\n", photonicsCmd::getName(cmdType, "").c_str());
      break;
  } 
  return photonicseval::perfEnergy(msRuntime, mjEnergy, msRead, msWrite, msALU, totalOp);
}

//! @brief  Perf energy model of Fulcrum for func2
photonicseval::perfEnergy
photonicsPerfEnergyFulcrum::getPerfEnergyForIter(const photonicsObjInfo& obj, const photonicsObjInfo& objSrc2, const photonicsObjInfo& objDest, int8_t numLoops) const
{
  double msRuntime = 0.0;
  double mjEnergy = 0.0;
  double msRead = 0.0;
  double msWrite = 0.0;
  double msALU = 0.0;
  uint64_t totalOp = 0;
  unsigned numPass = obj.getMaxNumRegionsPerCore();
  unsigned bitsPerElement = obj.getBitsPerElement(PhotonicsBitWidth::ACTUAL);
  unsigned numCoresUsed = obj.isLoadBalanced() ? obj.getNumCoreAvailable() : obj.getNumCoresUsed();
  unsigned maxElementsPerRegion = obj.getMaxElementsPerRegion();
  unsigned minElementPerRegion = obj.isLoadBalanced() ? (std::ceil(obj.getNumElements() * 1.0 / obj.getNumCoreAvailable()) - (maxElementsPerRegion * (numPass - 1))) : maxElementsPerRegion;
  double numberOfALUOperationPerElement = ((double)bitsPerElement / m_fulcrumAluBitWidth);

  msRead = 0;
  msWrite = 0;
  msALU = 0.00005 * numLoops;
  msRuntime = msRead + msWrite + msALU;
  mjEnergy = numCoresUsed * (numPass - 1) * ((m_eAP * 2) + ((maxElementsPerRegion - 1) * 2 *  m_fulcrumShiftEnergy * 2) + (maxElementsPerRegion * (m_fulcrumAddEnergy + m_fulcrumMulEnergy) * numberOfALUOperationPerElement));
  mjEnergy += numCoresUsed * ((m_eAP * 2) + ((minElementPerRegion - 1) * 2 *  m_fulcrumShiftEnergy * 2) + (minElementPerRegion * (m_fulcrumAddEnergy + m_fulcrumMulEnergy) * numberOfALUOperationPerElement));
  mjEnergy += m_pBChip * m_numChipsPerRank * m_numRanks * msRuntime;
  totalOp = obj.getNumElements() * 2;

  return photonicseval::perfEnergy(msRuntime, mjEnergy, msRead, msWrite, msALU, totalOp);
}

//! @brief  Perf energy model of Fulcrum for reduction sum
photonicseval::perfEnergy
photonicsPerfEnergyFulcrum::getPerfEnergyForReduction(PhotonicsCmdEnum cmdType, const photonicsObjInfo& obj, unsigned numPass) const
{
  double msRuntime = 0.0;
  double mjEnergy = 0.0;
  double msRead = 0.0;
  double msWrite = 0.0;
  double msCompute = 0.0;
  unsigned bitsPerElement = obj.getBitsPerElement(PhotonicsBitWidth::ACTUAL);
  unsigned maxElementsPerRegion = obj.getMaxElementsPerRegion();
  unsigned numCore = obj.isLoadBalanced() ? obj.getNumCoreAvailable() : obj.getNumCoresUsed();
  unsigned minElementPerRegion = obj.isLoadBalanced() ? (std::ceil(obj.getNumElements() * 1.0 / numCore) - (maxElementsPerRegion * (numPass - 1))) : maxElementsPerRegion;
  double cpuTDP = 225; // W; AMD EPYC 7742 64 core
  uint64_t totalOp = 0;

  switch (cmdType)
  {
  case PhotonicsCmdEnum::REDSUM:
  case PhotonicsCmdEnum::REDSUM_RANGE:
  case PhotonicsCmdEnum::REDMIN:
  case PhotonicsCmdEnum::REDMIN_RANGE:
  case PhotonicsCmdEnum::REDMAX:
  case PhotonicsCmdEnum::REDMAX_RANGE:
  {
    // read a row to walker, then reduce in serial
    double numberOfOperationPerElement = ((double)bitsPerElement / m_fulcrumAluBitWidth);
    // TODO: This needs to be flexible
    double aggregateMs = static_cast<double>(obj.getNumCoresUsed()) / 2300000;
    
    msRead = m_tR;
    msWrite = 0;
    msCompute = aggregateMs + (maxElementsPerRegion * m_fulcrumAddLatency * numberOfOperationPerElement * (numPass  - 1)) + (minElementPerRegion * m_fulcrumAddLatency * numberOfOperationPerElement);
    msRuntime = msRead + msWrite + msCompute;
    mjEnergy = (numPass - 1) * numCore * (m_eAP + ((maxElementsPerRegion - 1) *  m_fulcrumShiftEnergy) + (maxElementsPerRegion * m_fulcrumAddEnergy * numberOfOperationPerElement));
    mjEnergy += numCore * (m_eAP + ((minElementPerRegion - 1) *  m_fulcrumShiftEnergy) + (minElementPerRegion * m_fulcrumAddEnergy * numberOfOperationPerElement));
    mjEnergy += aggregateMs * cpuTDP;
    mjEnergy += m_pBChip * m_numChipsPerRank * m_numRanks * msRuntime;
    totalOp = obj.getNumElements();
    break;
  }
  default:
    printf("PHOTONICS-Warning: Perf energy model not available for PHOTONICS command %s\n", photonicsCmd::getName(cmdType, "").c_str());
    break;
  }
  return photonicseval::perfEnergy(msRuntime, mjEnergy, msRead, msWrite, msCompute, totalOp);
}

//! @brief  Perf energy model of Fulcrum for broadcast
photonicseval::perfEnergy
photonicsPerfEnergyFulcrum::getPerfEnergyForBroadcast(PhotonicsCmdEnum cmdType, const photonicsObjInfo& obj) const
{
  double msRuntime = 0.0;
  double mjEnergy = 0.0;
  double msRead = 0.0;
  double msWrite = 0.0;
  double msCompute = 0.0;
  unsigned numPass = obj.getMaxNumRegionsPerCore();
  unsigned bitsPerElement = obj.getBitsPerElement(PhotonicsBitWidth::ACTUAL);
  unsigned maxElementsPerRegion = obj.getMaxElementsPerRegion();
  unsigned numCore = obj.isLoadBalanced() ? obj.getNumCoreAvailable() : obj.getNumCoresUsed();
  unsigned minElementPerRegion = obj.isLoadBalanced() ? (std::ceil(obj.getNumElements() * 1.0 / numCore) - (maxElementsPerRegion * (numPass - 1))) : maxElementsPerRegion;
  uint64_t totalOp = 0;
  // assume taking 1 ALU latency to write an element
  double numberOfOperationPerElement = ((double)bitsPerElement / m_fulcrumAluBitWidth);
  msWrite = m_tW * numPass;
  msCompute = (m_fulcrumAddLatency * maxElementsPerRegion * numberOfOperationPerElement * (numPass - 1)) + (m_fulcrumAddLatency * minElementPerRegion * numberOfOperationPerElement);
  msRuntime = msRead + msWrite + msCompute;
  mjEnergy = (numPass - 1) * numCore * (m_eAP + ((maxElementsPerRegion - 1) *  m_fulcrumShiftEnergy) + ((maxElementsPerRegion) * m_fulcrumAddEnergy * numberOfOperationPerElement));
  mjEnergy += numCore * (m_eAP + ((minElementPerRegion - 1) *  m_fulcrumShiftEnergy) + (minElementPerRegion * m_fulcrumAddEnergy * numberOfOperationPerElement));
  mjEnergy += m_pBChip * m_numChipsPerRank * m_numRanks * msRuntime;
  totalOp = 0;

  return photonicseval::perfEnergy(msRuntime, mjEnergy, msRead, msWrite, msCompute, totalOp);
}

//! @brief  Perf energy model of Fulcrum for rotate
photonicseval::perfEnergy
photonicsPerfEnergyFulcrum::getPerfEnergyForRotate(PhotonicsCmdEnum cmdType, const photonicsObjInfo& obj) const
{
  double msRuntime = 0.0;
  double mjEnergy = 0.0;
  double msRead = 0.0;
  double msWrite = 0.0;
  double msCompute = 0.0;
  unsigned numPass = obj.getMaxNumRegionsPerCore();
  unsigned bitsPerElement = obj.getBitsPerElement(PhotonicsBitWidth::ACTUAL);
  unsigned numRegions = obj.getRegions().size();
  uint64_t totalOp = 0;
  // boundary handling - assume two times copying between device and host for boundary elements
  photonicseval::perfEnergy perfEnergyBT = getPerfEnergyForBytesTransfer(PhotonicsCmdEnum::COPY_D2H, numRegions * bitsPerElement / 8);

  // rotate within subarray:
  // For every bit: Read row to SA; move SA to R1; Shift R1 by N steps; Move R1 to SA; Write SA to row
  // TODO: separate bank level and GDL
  // TODO: energy unimplemented
  msRead = m_tR * numPass;
  msCompute = (bitsPerElement + 2) * m_tL * numPass;
  msWrite = m_tW * numPass;
  msRuntime = msRead + msWrite + msCompute;
  mjEnergy = (m_eAP + (bitsPerElement + 2) * m_eL) * numPass;
  msRuntime += 2 * perfEnergyBT.m_msRuntime;
  mjEnergy += 2 * perfEnergyBT.m_mjEnergy;
  printf("PHOTONICS-Warning: Perf energy model is not precise for PHOTONICS command %s\n", photonicsCmd::getName(cmdType, "").c_str());

  return photonicseval::perfEnergy(msRuntime, mjEnergy, msRead, msWrite, msCompute, totalOp);
}

//! @brief  Perf energy model of Fulcrum for prefix sum
photonicseval::perfEnergy
photonicsPerfEnergyFulcrum::getPerfEnergyForPrefixSum(PhotonicsCmdEnum cmdType, const photonicsObjInfo& obj) const
{
  double msRuntime = 0.0;
  double mjEnergy = 0.0;
  double msRead = 0.0;
  double msWrite = 0.0;
  double msCompute = 0.0;
  unsigned bitsPerElement = obj.getBitsPerElement(PhotonicsBitWidth::ACTUAL);
  unsigned maxElementsPerRegion = obj.getMaxElementsPerRegion();
  unsigned numCore = obj.isLoadBalanced() ? obj.getNumCoreAvailable() : obj.getNumCoresUsed();
  unsigned numPass = obj.getMaxNumRegionsPerCore();
  unsigned minElementPerRegion = obj.isLoadBalanced() ? (std::ceil(obj.getNumElements() * 1.0 / numCore) - (maxElementsPerRegion * (numPass - 1))) : maxElementsPerRegion;
  double cpuTDP = 225; // W; AMD EPYC 7742 64 core
  uint64_t totalOp = 0;

  switch (cmdType)
  {
  case PhotonicsCmdEnum::PREFIX_SUM:
  {

    /**
     * Performs prefix sum: dstVec[i] = dstVec[i-1] + srcVec[i]
     *
     * Execution Steps:
     * 1. Each subarray performs a local prefix sum on its portion of the data.
     * 2. The host CPU fetches the final value from each subarray using `n` DRAM READ. Here, `n = number of subarrays`.
     * 3. The host CPU aggregates these values (i.e., computes the prefix sum across subarrays).
     * 4. The host CPU writes the aggregated values back to DRAM using `n` DRAM WRITE.
     * 5. Each subarray updates its elements using the received value to complete the final prefix sum.
     *
     * Performance Model:
     * - While performing addition, the next row can be
     * fetched concurrently. As a result, `msRead = 2 * m_tR` (multiplied by two because, two prefix sum iterations are
     * required).
     * - `aggregateMs` models the time for host-side aggregation.
     * - `hostRW` accounts for host read/write overhead, including DRAM tR, tW,
     * and GDL delays.
     *
    */

    // read a row to walker, then add in serial
    double numberOfOperationPerElement = ((double)bitsPerElement / m_fulcrumAluBitWidth);
    
    // TODO: This needs to be flexible
    double aggregateMs = static_cast<double>(obj.getNumCoresUsed()) / 2300000;
    double hostRW = (obj.getNumCoresUsed() * 1.0 / m_numChipsPerRank) * (m_tR + m_tW + (m_tGDL * 2));

    msRead = 2 * m_tR;
    msWrite = 2 * m_tW;
    msCompute = aggregateMs + hostRW + (maxElementsPerRegion * m_fulcrumAddLatency * numberOfOperationPerElement * (numPass  - 1) * 2) + (minElementPerRegion * m_fulcrumAddLatency * numberOfOperationPerElement) * 2;
    msRuntime = msRead + msWrite + msCompute;
    mjEnergy = 2 * (numPass - 1) * numCore * (m_eAP + ((maxElementsPerRegion - 1) *  m_fulcrumShiftEnergy) + (maxElementsPerRegion * m_fulcrumAddEnergy * numberOfOperationPerElement * 2));
    mjEnergy += 2 * numCore * (m_eAP + ((minElementPerRegion - 1) *  m_fulcrumShiftEnergy) + (minElementPerRegion * m_fulcrumAddEnergy * numberOfOperationPerElement));
    mjEnergy += aggregateMs * cpuTDP + ((obj.getNumCoresUsed() * 1.0 / m_numChipsPerRank) * ((2 * m_eAP)  + m_eR + m_eW));
    mjEnergy += m_pBChip * m_numChipsPerRank * m_numRanks * msRuntime;
    totalOp = obj.getNumElements() * 2;
    break;
  }
  default:
    printf("PHOTONICS-Warning: Perf energy model not available for PHOTONICS command %s\n", photonicsCmd::getName(cmdType, "").c_str());
    break;
  }
  return photonicseval::perfEnergy(msRuntime, mjEnergy, msRead, msWrite, msCompute, totalOp);
}