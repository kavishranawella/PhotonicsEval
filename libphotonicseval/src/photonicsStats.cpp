// File: photonicsStats.cpp
// PHOTONICSeval Simulator - Stats
// Copyright (c) 2024 University of Virginia
// This file is licensed under the MIT License.
// See the LICENSE file in the root of this repository for more details.

#include "photonicsStats.h"
#include "photonicsSim.h"
#include "photonicsUtils.h"
#include <chrono>            // for chrono
#include <cstdint>           // for uint64_t
#include <cstdio>            // for printf
#include <iomanip>           // for setw, fixed, setprecision


//! @brief  Show PHOTONICS stats
void
photonicsStatsMgr::showStats() const
{
  std::printf("----------------------------------------\n");
  if (photonicsSim::get()->isDebug(photonicsSimConfig::DEBUG_API_CALLS)) {
    showApiStats();
  }
  showDeviceParams();
  showCopyStats();
  showCmdStats();
  std::printf("----------------------------------------\n");
}

//! @brief  Show API stats
void
photonicsStatsMgr::showApiStats() const
{
  std::printf("Simulator API Stats:\n");
  std::printf(" %30s : %10s %14s\n", "PHOTONICS-API", "CNT", "Elapsed(ms)");
  int totCalls = 0;
  int totCallsDevice = 0;
  int totCallsAlloc = 0;
  int totCallsCopy = 0;
  int totCallsCompute = 0;
  double msTotalElapsed = 0.0;
  double msTotalElapsedDevice = 0.0;
  double msTotalElapsedAlloc = 0.0;
  double msTotalElapsedCopy = 0.0;
  double msTotalElapsedCompute = 0.0;
  for (const auto& it : m_msElapsed) {
    std::printf(" %30s : %10d %14f\n", it.first.c_str(), it.second.first, it.second.second);
    totCalls += it.second.first;
    msTotalElapsed += it.second.second;
    if (it.first.find("createDevice") == 0) {
      totCallsDevice += it.second.first;
      msTotalElapsedDevice += it.second.second;
    } else if (it.first.find("photonicsAllocMat") == 0 || it.first.find("photonicsAlloc") == 0 || it.first.find("photonicsFree") == 0) {
      totCallsAlloc += it.second.first;
      msTotalElapsedAlloc += it.second.second;
    } else if (it.first.find("photonicsCopy") == 0) {
      totCallsCopy += it.second.first;
      msTotalElapsedCopy += it.second.second;
    } else {
      totCallsCompute += it.second.first;
      msTotalElapsedCompute += it.second.second;
    }
  }
  std::printf(" %30s : %10d %14f\n", "TOTAL ---------", totCalls, msTotalElapsed);
  std::printf(" %30s : %10d %14f\n", "TOTAL (Device )", totCallsDevice, msTotalElapsedDevice);
  std::printf(" %30s : %10d %14f\n", "TOTAL ( Alloc )", totCallsAlloc, msTotalElapsedAlloc);
  std::printf(" %30s : %10d %14f\n", "TOTAL (  Copy )", totCallsCopy, msTotalElapsedCopy);
  std::printf(" %30s : %10d %14f\n", "TOTAL (Compute)", totCallsCompute, msTotalElapsedCompute);
}

//! @brief  Show PHOTONICS device params
void
photonicsStatsMgr::showDeviceParams() const
{
  const photonicsParamsDram& paramsDram = photonicsSim::get()->getParamsDram();
  std::printf("Photonics Params:\n");
  std::printf(" %30s : %s\n", "Photonics Device Type Enum",
              photonicsUtils::photonicsDeviceEnumToStr(photonicsSim::get()->getDeviceType()).c_str());
  std::printf(" %30s : %s\n", "Photonics Simulation Target",
              photonicsUtils::photonicsDeviceEnumToStr(photonicsSim::get()->getSimTarget()).c_str());
  std::printf(" %30s : %u, %u, %u, %u\n", "VCores, HCores, VectorsPerCore, ElementsPerVector",
              photonicsSim::get()->getNumRanks(),
              photonicsSim::get()->getNumBankPerRank(),
              photonicsSim::get()->getNumSubarrayPerBank(),
              photonicsSim::get()->getNumColPerSubarray()/32);
  std::printf(" %30s : %u\n", "Number of Photonics Cores", photonicsSim::get()->getNumPhotonicCores());
  std::printf(" %30s : %u\n", "Number of Rows per Core", photonicsSim::get()->getNumSubarrayPerBank());
  std::printf(" %30s : %u\n", "Number of Cols per Core", photonicsSim::get()->getNumCols()/32);
  std::printf(" %30s : %u\n", "Matrix Size", photonicsSim::get()->getMatrixSize());
  std::printf(" %30s : %f GB/s\n", "Memory BW", paramsDram.getTypicalRankBW());
}

//! @brief  Show data copy stats
void
photonicsStatsMgr::showCopyStats() const
{
  std::printf("Data Copy Stats:\n");
  uint64_t bytesCopiedMainToDevice = m_bitsCopiedMainToDevice / 8;
  uint64_t bytesCopiedDeviceToMain = m_bitsCopiedDeviceToMain / 8;
  uint64_t bytesCopiedDeviceToDevice = m_bitsCopiedDeviceToDevice / 8;
  uint64_t totalBytes = bytesCopiedMainToDevice + bytesCopiedDeviceToMain;
  double totalMsRuntime = m_elapsedTimeCopiedMainToDevice + m_elapsedTimeCopiedDeviceToMain + m_elapsedTimeCopiedDeviceToDevice;
  double totalMjEnergy = m_mJCopiedMainToDevice + m_mJCopiedDeviceToMain + m_mJCopiedDeviceToDevice;
  std::printf(" %45s : %llu bytes\n", "Host to Device", (unsigned long long)bytesCopiedMainToDevice);
  std::printf(" %45s : %llu bytes\n", "Device to Host", (unsigned long long)bytesCopiedDeviceToMain);
  std::printf(" %45s : %llu bytes\n", "Device to Device", (unsigned long long)bytesCopiedDeviceToDevice);
  std::printf(" %45s : %llu bytes %14.6f ms Estimated Runtime %14.6f mj Estimated Energy\n", "TOTAL ---------", (unsigned long long)totalBytes, totalMsRuntime, totalMjEnergy);
}

//! @brief  Show PHOTONICS cmd and perf stats
void
photonicsStatsMgr::showCmdStats() const
{
  std::printf("PHOTONICS Command Stats:\n");
  std::printf(" %44s : %10s %14s %14s %14s %7s %7s %7s\n", "PHOTONICS-CMD", "CNT", "Runtime(ms)", "Energy(mJ)", "GOPS/W", "%R", "%W", "%L");
  int totalCmd = 0;
  double totalMsRuntime = 0.0;
  double totalMjEnergy = 0.0;
  double totalMsRead = 0.0;
  double totalMsWrite = 0.0;
  double totalMsCompute = 0.0;
  uint64_t totalOp = 0;
  for (const auto& it : m_cmdPerf) {
    double cmdRuntime = it.second.second.m_msRuntime;
    double percentRead = cmdRuntime == 0.0 ? 0.0 : (it.second.second.m_msRead * 100 / cmdRuntime);
    double percentWrite = cmdRuntime == 0.0 ? 0.0 : (it.second.second.m_msWrite * 100 / cmdRuntime);
    double percentCompute = cmdRuntime == 0.0 ? 0.0 : (it.second.second.m_msCompute * 100 / cmdRuntime);
    double cmdEnergy = it.second.second.m_mjEnergy;
    double perfWatt = cmdEnergy == 0.0 ? 0.0 : (it.second.second.m_totalOp * 1.0 / cmdEnergy * 1e-6);
    std::printf(" %44s : %10d %14f %14f %14f %7.2f %7.2f %7.2f\n", it.first.c_str(), it.second.first, it.second.second.m_msRuntime, it.second.second.m_mjEnergy, perfWatt, percentRead, percentWrite, percentCompute);
    totalCmd += it.second.first;
    totalMsRuntime += it.second.second.m_msRuntime;
    totalMjEnergy += it.second.second.m_mjEnergy;
    totalMsRead += it.second.first * percentRead;
    totalMsWrite += it.second.first * percentWrite;
    totalMsCompute += it.second.first * percentCompute;
    totalOp += it.second.second.m_totalOp;
  }
  std::printf(" %44s : %10d %14f %14f %14f %7.2f %7.2f %7.2f\n", "TOTAL ---------", totalCmd, totalMsRuntime, totalMjEnergy, (totalOp * 1.0 / totalMjEnergy * 1e-6), (totalMsRead / totalCmd), (totalMsWrite / totalCmd), (totalMsCompute / totalCmd) );
  // analyze micro-ops
  int numR = 0;
  int numW = 0;
  int numL = 0;
  int numActivate = 0;
  int numPrecharge = 0;
  for (const auto& it : m_cmdPerf) {
    if (it.first == "row_r") {
      numR += it.second.first;
      numActivate += it.second.first;
      numPrecharge += it.second.first;
    } else if (it.first == "row_w") {
      numW += it.second.first;
      numActivate += it.second.first;
      numPrecharge += it.second.first;
    } else if (it.first.find("rreg.") == 0) {
      numL += it.second.first;
    }
  }
  if (numR > 0 || numW > 0 || numL > 0) {
    std::printf(" %30s : %d, %d, %d\n", "Num Read, Write, Logic", numR, numW, numL);
    std::printf(" %30s : %d, %d\n", "Num Activate, Precharge", numActivate, numPrecharge);
  }
}

//! @brief  Reset PHOTONICS stats
void
photonicsStatsMgr::resetStats()
{
  m_cmdPerf.clear();
  m_msElapsed.clear();
  m_bitsCopiedMainToDevice = 0;
  m_bitsCopiedDeviceToMain = 0;
  m_bitsCopiedDeviceToDevice = 0;
}

//! @brief  Record estimated runtime and energy of a PHOTONICS command
void
photonicsStatsMgr::recordCmd(const std::string& cmdName, photonicseval::perfEnergy mPerfEnergy)
{
  auto& item = m_cmdPerf[cmdName];
  item.first++;
  item.second.m_msRuntime += mPerfEnergy.m_msRuntime;
  m_curApiMsEstRuntime += mPerfEnergy.m_msRuntime;
  item.second.m_mjEnergy += mPerfEnergy.m_mjEnergy;
  item.second.m_msRead += mPerfEnergy.m_msRead;
  item.second.m_msWrite += mPerfEnergy.m_msWrite;
  item.second.m_msCompute += mPerfEnergy.m_msCompute;
  item.second.m_totalOp += mPerfEnergy.m_totalOp;
}

//! @brief  Record estimated runtime and energy of data copy
void
photonicsStatsMgr::recordCopyMainToDevice(uint64_t numBits, photonicseval::perfEnergy mPerfEnergy)
{
  m_bitsCopiedMainToDevice += numBits;
  m_elapsedTimeCopiedMainToDevice += mPerfEnergy.m_msRuntime;
  m_curApiMsEstRuntime += mPerfEnergy.m_msRuntime;
  m_mJCopiedMainToDevice += mPerfEnergy.m_mjEnergy;
}

//! @brief  Record estimated runtime and energy of data copy
void
photonicsStatsMgr::recordCopyDeviceToMain(uint64_t numBits, photonicseval::perfEnergy mPerfEnergy)
{
  m_bitsCopiedDeviceToMain += numBits;
  m_elapsedTimeCopiedDeviceToMain += mPerfEnergy.m_msRuntime;
  m_curApiMsEstRuntime += mPerfEnergy.m_msRuntime;
  m_mJCopiedDeviceToMain += mPerfEnergy.m_mjEnergy;
}

//! @brief  Record estimated runtime and energy of data copy
void
photonicsStatsMgr::recordCopyDeviceToDevice(uint64_t numBits, photonicseval::perfEnergy mPerfEnergy)
{
  m_bitsCopiedDeviceToDevice += numBits;
  m_elapsedTimeCopiedDeviceToDevice += mPerfEnergy.m_msRuntime;
  m_curApiMsEstRuntime += mPerfEnergy.m_msRuntime;
  m_mJCopiedDeviceToDevice += mPerfEnergy.m_mjEnergy;
}

//! @brief  Preprocessing at the beginning of a PHOTONICS API scope
void
photonicsStatsMgr::photonicsApiScopeStart()
{
  // Restart for current PHOTONICS API call
  m_curApiMsEstRuntime = 0.0;
}

//! @brief  Postprocessing at the end of a PHOTONICS API scope
void
photonicsStatsMgr::photonicsApiScopeEnd(const std::string& tag, double elapsed)
{
  // Record API stats
  auto& item = m_msElapsed[tag];
  item.first++;
  item.second += elapsed;

  // Update kernel stats
  if (m_isKernelTimerOn) {
    m_kernelMsElapsedSim += elapsed;
    m_kernelMsEstRuntime += m_curApiMsEstRuntime;
  }
}

//! @brief  Start timer for a PHOTONICS kernel to measure CPU runtime and DRAM refresh
void
photonicsStatsMgr::startKernelTimer()
{
  if (m_isKernelTimerOn) {
    std::printf("PHOTONICS-Warning: Kernel timer has already started\n");
    return;
  }
  std::printf("PHOTONICS-Info: Start kernel timer.\n");
  m_isKernelTimerOn = true;
  m_kernelStart = std::chrono::high_resolution_clock::now();
}

//! @brief  End timer for a PHOTONICS kernel to measure CPU runtime and DRAM refresh
void
photonicsStatsMgr::endKernelTimer()
{
  if (!m_isKernelTimerOn) {
    std::printf("PHOTONICS-Warning: Kernel timer has not started\n");
    return;
  }
  auto now = std::chrono::high_resolution_clock::now();
  double kernelMsElapsedTotal = std::chrono::duration<double, std::milli>(now - m_kernelStart).count();
  double kernelMsElapsedCpu = kernelMsElapsedTotal - m_kernelMsElapsedSim;
  std::printf("PHOTONICS-Info: End kernel timer. Runtime = %14f ms, CPU = %14f ms, PHOTONICS = %14f ms\n",
      kernelMsElapsedCpu + m_kernelMsEstRuntime, kernelMsElapsedCpu, m_kernelMsEstRuntime);
  m_kernelStart = std::chrono::high_resolution_clock::time_point(); // reset
  m_isKernelTimerOn = false;
}

//! @brief photonicsPerfMon ctor
photonicsPerfMon::photonicsPerfMon(const std::string& tag)
{
  m_startTime = std::chrono::high_resolution_clock::now();
  m_tag = tag;
  // assumption: photonicsPerfMon is not nested
  if (photonicsSim::get()->getStatsMgr()) {
    photonicsSim::get()->getStatsMgr()->photonicsApiScopeStart();
  }
}

//! @brief photonicsPerfMon dtor
photonicsPerfMon::~photonicsPerfMon()
{
  auto now = std::chrono::high_resolution_clock::now();
  double elapsed = std::chrono::duration<double, std::milli>(now - m_startTime).count();
  if (photonicsSim::get()->getStatsMgr()) {
    photonicsSim::get()->getStatsMgr()->photonicsApiScopeEnd(m_tag, elapsed);
  }
}

