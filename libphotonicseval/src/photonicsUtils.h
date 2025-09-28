// File: photonicsUtils.h
// PHOTONICSeval Simulator - Utilities
// Copyright (c) 2024 University of Virginia
// This file is licensed under the MIT License.
// See the LICENSE file in the root of this repository for more details.

#ifndef LAVA_PHOTONICS_UTILS_H
#define LAVA_PHOTONICS_UTILS_H

#include "libphotonicseval.h"
#include <string>
#include <queue>
#include <vector>
#include <thread>
#include <atomic>
#include <mutex>
#include <condition_variable>
#include <algorithm>
#include <cctype>
#include <locale>
#include <unordered_map>
#include <type_traits>
#include <cstring>
#include <cstdint>


//! @enum   PhotonicsBitWidth
//! @brief  Bit width definitions of PHOTONICS data types under different usage scenarios
enum class PhotonicsBitWidth
{
  ACTUAL = 0,  // bit width of a data type in real hardware
  HOST,        // bit width of a data type used by host for data transfer in PHOTONICSeval
  SIM,         // bit width of a data type used for functional compution in PHOTONICSeval
  PADDED,      // bit width of a data element with association and padding in PHOTONICSeval
};

//! @enum   PhotonicsDataLayout
//! @brief  PHOTONICS data layout definitions
enum class PhotonicsDataLayout
{
  UNKNOWN = 0,  // unknown
  H,            // horizontal data layout
  V,            // vertical data layout
  HYBRID,       // hybrid data layout
};

namespace photonicsUtils
{
  std::string photonicsStatusEnumToStr(PhotonicsStatus status);
  std::string photonicsDeviceEnumToStr(PhotonicsDeviceEnum deviceType);
  PhotonicsDeviceEnum strToPhotonicsDeviceEnum(const std::string& deviceTypeStr);
  std::string photonicsAllocEnumToStr(PhotonicsAllocEnum allocType);
  std::string photonicsCopyEnumToStr(PhotonicsCopyEnum copyType);
  std::string photonicsDataTypeEnumToStr(PhotonicsDataType dataType);
  unsigned getNumBitsOfDataType(PhotonicsDataType dataType, PhotonicsBitWidth bitWidthType);
  bool isSigned(PhotonicsDataType dataType);
  bool isUnsigned(PhotonicsDataType dataType);
  bool isFP(PhotonicsDataType dataType);
  std::string photonicsProtocolEnumToStr(PhotonicsDeviceProtocolEnum protocol);
  PhotonicsDataLayout getDeviceDataLayout(PhotonicsDeviceEnum deviceType);

  // Convert raw bits into sign-extended bits based on PHOTONICS data type.
  // Input: Raw bits represented as uint64_t
  // Output: Sign-extended bits represented as uint64_t
  inline uint64_t signExt(uint64_t bits, PhotonicsDataType dataType) {
    switch (dataType) {
      case PHOTONICS_INT8: return static_cast<uint64_t>(static_cast<int64_t>(static_cast<int8_t>(bits)));
      case PHOTONICS_INT16: return static_cast<uint64_t>(static_cast<int64_t>(static_cast<int16_t>(bits)));
      case PHOTONICS_INT32: return static_cast<uint64_t>(static_cast<int64_t>(static_cast<int32_t>(bits)));
      case PHOTONICS_INT64: return static_cast<uint64_t>(static_cast<int64_t>(static_cast<int64_t>(bits)));
      default: break; // no-op
    }
    return bits;
  }

  // Convert sign-extended bits into specific C++ type.
  // Input: Sign-extended bits represented as uint64_t
  // Output: A value in C++ data type T
  template <typename T> T castBitsToType(uint64_t signExtBits) {
    T val;
    std::memcpy(&val, &signExtBits, sizeof(T));
    return val;
  }

  // Convert specific type into sign-extended bits.
  // Input: A value in C++ data type T
  // Output: sign-extended bits represented as uint64_t
  template <typename T> uint64_t castTypeToBits(T val) {
    uint64_t signExtBits = 0;
    if constexpr (std::is_integral<T>::value && std::is_signed<T>::value) {
      signExtBits = static_cast<uint64_t>(static_cast<int64_t>(val)); // sign ext
    } else {
      std::memcpy(&signExtBits, &val, sizeof(T)); // zero padding
    }
    return signExtBits;
  }

  // Service APIs for file system, config files, env vars
  std::string& ltrim(std::string& s);
  std::string& rtrim(std::string& s);
  std::string& trim(std::string& s);
  bool readFileContent(const char* fileName, std::string& fileContent);
  std::string getParam(const std::unordered_map<std::string, std::string>& params, const std::string& key);
  std::string getOptionalParam(const std::unordered_map<std::string, std::string>& params, const std::string& key, bool& returnStatus);
  std::string removeAfterSemicolon(const std::string &input);

  std::string getDirectoryPath(const std::string& filePath);
  bool getEnvVar(const std::string &varName, std::string &varValue);
  bool convertStringToUnsigned(const std::string& str, unsigned& retVal);
  std::unordered_map<std::string, std::string> readParamsFromConfigFile(const std::string& configFilePath);
  std::unordered_map<std::string, std::string> readParamsFromEnvVars(const std::vector<std::string>& envVarNames);

  const std::unordered_map<PhotonicsDeviceEnum, std::string> enumToStrMap = {
      {PHOTONICS_DEVICE_NONE, "PHOTONICS_DEVICE_NONE"},
      {PHOTONICS_FUNCTIONAL, "PHOTONICS_FUNCTIONAL"},
      {PHOTONICS_DEVICE_BITSIMD_V, "PHOTONICS_DEVICE_BITSIMD_V"},
      {PHOTONICS_DEVICE_BITSIMD_V_NAND, "PHOTONICS_DEVICE_BITSIMD_V_NAND"},
      {PHOTONICS_DEVICE_BITSIMD_V_MAJ, "PHOTONICS_DEVICE_BITSIMD_V_MAJ"},
      {PHOTONICS_DEVICE_BITSIMD_V_AP, "PHOTONICS_DEVICE_BITSIMD_V_AP"},
      {PHOTONICS_DEVICE_DRISA_NOR, "PHOTONICS_DEVICE_DRISA_NOR"},
      {PHOTONICS_DEVICE_DRISA_MIXED, "PHOTONICS_DEVICE_DRISA_MIXED"},
      {PHOTONICS_DEVICE_SIMDRAM, "PHOTONICS_DEVICE_SIMDRAM"},
      {PHOTONICS_DEVICE_BITSIMD_H, "PHOTONICS_DEVICE_BITSIMD_H"},
      {PHOTONICS_DEVICE_FULCRUM, "PHOTONICS_DEVICE_FULCRUM"},
      {PHOTONICS_DEVICE_BANK_LEVEL, "PHOTONICS_DEVICE_BANK_LEVEL"},
      {PHOTONICS_DEVICE_AQUABOLT, "PHOTONICS_DEVICE_AQUABOLT"},
      {PHOTONICS_DEVICE_AIM, "PHOTONICS_DEVICE_AIM"}
  };

  const std::unordered_map<std::string, PhotonicsDeviceEnum> strToEnumMap = {
      {"PHOTONICS_DEVICE_NONE", PHOTONICS_DEVICE_NONE},
      {"PHOTONICS_FUNCTIONAL", PHOTONICS_FUNCTIONAL},
      {"PHOTONICS_DEVICE_BITSIMD_V", PHOTONICS_DEVICE_BITSIMD_V},
      {"PHOTONICS_DEVICE_BITSIMD_V_NAND", PHOTONICS_DEVICE_BITSIMD_V_NAND},
      {"PHOTONICS_DEVICE_BITSIMD_V_MAJ", PHOTONICS_DEVICE_BITSIMD_V_MAJ},
      {"PHOTONICS_DEVICE_BITSIMD_V_AP", PHOTONICS_DEVICE_BITSIMD_V_AP},
      {"PHOTONICS_DEVICE_DRISA_NOR", PHOTONICS_DEVICE_DRISA_NOR},
      {"PHOTONICS_DEVICE_DRISA_MIXED", PHOTONICS_DEVICE_DRISA_MIXED},
      {"PHOTONICS_DEVICE_SIMDRAM", PHOTONICS_DEVICE_SIMDRAM},
      {"PHOTONICS_DEVICE_BITSIMD_H", PHOTONICS_DEVICE_BITSIMD_H},
      {"PHOTONICS_DEVICE_FULCRUM", PHOTONICS_DEVICE_FULCRUM},
      {"PHOTONICS_DEVICE_BANK_LEVEL", PHOTONICS_DEVICE_BANK_LEVEL},
      {"PHOTONICS_DEVICE_AQUABOLT", PHOTONICS_DEVICE_AQUABOLT},
      {"PHOTONICS_DEVICE_AIM", PHOTONICS_DEVICE_AIM}
  };

  //! @class  threadWorker
  //! @brief  Thread worker base class
  class threadWorker {
  public:
    threadWorker() {}
    virtual ~threadWorker() {}
    virtual void execute() = 0;
  };

  //! @class  threadPool
  //! @brief  Thread pool that runs multiple workers in threads
  class threadPool {
  public:
    threadPool(size_t numThreads);
    ~threadPool();
    void doWork(const std::vector<photonicsUtils::threadWorker*>& workers);
  private:
    void workerThread();

    std::vector<std::thread> m_threads;
    std::queue<threadWorker*> m_workers;
    std::mutex m_mutex;
    std::condition_variable m_cond;
    bool m_terminate;
    std::atomic<size_t> m_workersRemaining;
  };

}

#endif

