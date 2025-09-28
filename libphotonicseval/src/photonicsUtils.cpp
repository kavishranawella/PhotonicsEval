// File: photonicsUtils.cc
// PHOTONICSeval Simulator - Utilities
// Copyright (c) 2024 University of Virginia
// This file is licensed under the MIT License.
// See the LICENSE file in the root of this repository for more details.

#include "photonicsUtils.h"
#include "libphotonicseval.h"
#include <fstream>
#include <cstdio>
#include <sstream>
#include <unordered_map>
#include <string>
#include <filesystem>
#include <cstdlib>
#include <cassert>
#include <stdexcept>


//! @brief  Convert PhotonicsStatus enum to string
std::string
photonicsUtils::photonicsStatusEnumToStr(PhotonicsStatus status)
{
  switch (status) {
  case PHOTONICS_ERROR: return "ERROR";
  case PHOTONICS_OK: return "OK";
  }
  return "Unknown";
}

//! @brief  Convert PhotonicsDeviceEnum to string
std::string 
photonicsUtils::photonicsDeviceEnumToStr(PhotonicsDeviceEnum deviceType) {
  auto it = enumToStrMap.find(deviceType);
  if (it != enumToStrMap.end()) {
    return it->second;
  }
  return "Unknown";
}

//! @brief  Convert string to PhotonicsDeviceEnum
PhotonicsDeviceEnum
photonicsUtils::strToPhotonicsDeviceEnum(const std::string& deviceTypeStr) {
  auto it = strToEnumMap.find(deviceTypeStr);
  if (it != strToEnumMap.end()) {
    return it->second;
  }
  return PHOTONICS_DEVICE_NONE;
}

//! @brief  Convert PhotonicsAllocEnum to string
std::string
photonicsUtils::photonicsAllocEnumToStr(PhotonicsAllocEnum allocType)
{
  switch (allocType) {
  case PHOTONICS_ALLOC_AUTO: return "PHOTONICS_ALLOC_AUTO";
  case PHOTONICS_ALLOC_V: return "PHOTONICS_ALLOC_V";
  case PHOTONICS_ALLOC_H: return "PHOTONICS_ALLOC_H";
  case PHOTONICS_ALLOC_V1: return "PHOTONICS_ALLOC_V1";
  case PHOTONICS_ALLOC_H1: return "PHOTONICS_ALLOC_H1";
  }
  return "Unknown";
}

//! @brief  Convert PhotonicsCopyEnum to string
std::string
photonicsUtils::photonicsCopyEnumToStr(PhotonicsCopyEnum copyType)
{
  switch (copyType) {
  case PHOTONICS_COPY_V: return "PHOTONICS_COPY_V";
  case PHOTONICS_COPY_H: return "PHOTONICS_COPY_H";
  }
  return "Unknown";
}

//! @brief  Convert PhotonicsDataType enum to string
std::string
photonicsUtils::photonicsDataTypeEnumToStr(PhotonicsDataType dataType)
{
  switch (dataType) {
  case PHOTONICS_BOOL: return "bool";
  case PHOTONICS_INT8: return "int8";
  case PHOTONICS_INT16: return "int16";
  case PHOTONICS_INT32: return "int32";
  case PHOTONICS_INT64: return "int64";
  case PHOTONICS_UINT8: return "uint8";
  case PHOTONICS_UINT16: return "uint16";
  case PHOTONICS_UINT32: return "uint32";
  case PHOTONICS_UINT64: return "uint64";
  case PHOTONICS_FP32: return "fp32";
  case PHOTONICS_FP16: return "fp16";
  case PHOTONICS_BF16: return "bf16";
  case PHOTONICS_FP8: return "fp8";
  }
  return "Unknown";
}

namespace photonicsUtils {
  //! @brief  Static definitions of bits of data types (see PhotonicsBitWidth)
  //! Notes:
  //! - BOOL: PHOTONICSeval requires host data to store one bool value per byte
  //! - FP16/BF16/FP8: PHOTONICSeval uses FP32 for functional simulation
  static const std::unordered_map<PhotonicsDataType, std::unordered_map<PhotonicsBitWidth, unsigned>> s_bitsOfDataType = {
    { PHOTONICS_BOOL, {{PhotonicsBitWidth::ACTUAL, 1}, {PhotonicsBitWidth::SIM, 1}, {PhotonicsBitWidth::HOST, 8}} },
    { PHOTONICS_INT8, {{PhotonicsBitWidth::ACTUAL, 8}, {PhotonicsBitWidth::SIM, 8}, {PhotonicsBitWidth::HOST, 8}} },
    { PHOTONICS_INT16, {{PhotonicsBitWidth::ACTUAL, 16}, {PhotonicsBitWidth::SIM, 16}, {PhotonicsBitWidth::HOST, 16}} },
    { PHOTONICS_INT32, {{PhotonicsBitWidth::ACTUAL, 32}, {PhotonicsBitWidth::SIM, 32}, {PhotonicsBitWidth::HOST, 32}} },
    { PHOTONICS_INT64, {{PhotonicsBitWidth::ACTUAL, 64}, {PhotonicsBitWidth::SIM, 64}, {PhotonicsBitWidth::HOST, 64}} },
    { PHOTONICS_UINT8, {{PhotonicsBitWidth::ACTUAL, 8}, {PhotonicsBitWidth::SIM, 8}, {PhotonicsBitWidth::HOST, 8}} },
    { PHOTONICS_UINT16, {{PhotonicsBitWidth::ACTUAL, 16}, {PhotonicsBitWidth::SIM, 16}, {PhotonicsBitWidth::HOST, 16}} },
    { PHOTONICS_UINT32, {{PhotonicsBitWidth::ACTUAL, 32}, {PhotonicsBitWidth::SIM, 32}, {PhotonicsBitWidth::HOST, 32}} },
    { PHOTONICS_UINT64, {{PhotonicsBitWidth::ACTUAL, 64}, {PhotonicsBitWidth::SIM, 64}, {PhotonicsBitWidth::HOST, 64}} },
    { PHOTONICS_FP32, {{PhotonicsBitWidth::ACTUAL, 32}, {PhotonicsBitWidth::SIM, 32}, {PhotonicsBitWidth::HOST, 32}} },
    { PHOTONICS_FP16, {{PhotonicsBitWidth::ACTUAL, 16}, {PhotonicsBitWidth::SIM, 32}, {PhotonicsBitWidth::HOST, 32}} },
    { PHOTONICS_BF16, {{PhotonicsBitWidth::ACTUAL, 16}, {PhotonicsBitWidth::SIM, 32}, {PhotonicsBitWidth::HOST, 32}} },
    { PHOTONICS_FP8, {{PhotonicsBitWidth::ACTUAL, 8}, {PhotonicsBitWidth::SIM, 32}, {PhotonicsBitWidth::HOST, 32}} },
  };
}

//! @brief  Get number of bits of a PHOTONICS data type
unsigned
photonicsUtils::getNumBitsOfDataType(PhotonicsDataType dataType, PhotonicsBitWidth bitWidthType)
{
  if (bitWidthType == PhotonicsBitWidth::ACTUAL || bitWidthType == PhotonicsBitWidth::SIM || bitWidthType == PhotonicsBitWidth::HOST) {
    auto it = photonicsUtils::s_bitsOfDataType.find(dataType);
    return it != s_bitsOfDataType.end() ? it->second.at(bitWidthType) : 0;
  }
  return 0;
}

//! @brief  Check if a PHOTONICS data type is signed integer
bool
photonicsUtils::isSigned(PhotonicsDataType dataType)
{
  return dataType == PHOTONICS_INT8 || dataType == PHOTONICS_INT16 || dataType == PHOTONICS_INT32 || dataType == PHOTONICS_INT64;
}

//! @brief  Check if a PHOTONICS data type is unsigned integer
bool
photonicsUtils::isUnsigned(PhotonicsDataType dataType)
{
  return dataType == PHOTONICS_BOOL || dataType == PHOTONICS_UINT8 || dataType == PHOTONICS_UINT16 || dataType == PHOTONICS_UINT32 || dataType == PHOTONICS_UINT64;
}

//! @brief  Check if a PHOTONICS data type is floating point
bool
photonicsUtils::isFP(PhotonicsDataType dataType)
{
  return dataType == PHOTONICS_FP32 || dataType == PHOTONICS_FP16 || dataType == PHOTONICS_BF16 || dataType == PHOTONICS_FP8;
}

//! @brief  Convert PhotonicsDeviceProtocolEnum to string
std::string
photonicsUtils::photonicsProtocolEnumToStr(PhotonicsDeviceProtocolEnum protocol)
{
  switch (protocol) {
    case PHOTONICS_DEVICE_PROTOCOL_DDR: return "DDR";
    case PHOTONICS_DEVICE_PROTOCOL_LPDDR: return "LPDDR";
    case PHOTONICS_DEVICE_PROTOCOL_HBM: return "HBM";
    case PHOTONICS_DEVICE_PROTOCOL_GDDR: return "GDDR";
  }
  return "Unknown";
}

//! @brief  Get device data layout
PhotonicsDataLayout
photonicsUtils::getDeviceDataLayout(PhotonicsDeviceEnum deviceType)
{
  switch (deviceType) {
    case PHOTONICS_DEVICE_BITSIMD_V: return PhotonicsDataLayout::V;
    case PHOTONICS_DEVICE_BITSIMD_V_NAND: return PhotonicsDataLayout::V;
    case PHOTONICS_DEVICE_BITSIMD_V_MAJ: return PhotonicsDataLayout::V;
    case PHOTONICS_DEVICE_BITSIMD_V_AP: return PhotonicsDataLayout::V;
    case PHOTONICS_DEVICE_DRISA_NOR: return PhotonicsDataLayout::V;
    case PHOTONICS_DEVICE_DRISA_MIXED: return PhotonicsDataLayout::V;
    case PHOTONICS_DEVICE_SIMDRAM: return PhotonicsDataLayout::V;
    case PHOTONICS_DEVICE_BITSIMD_H: return PhotonicsDataLayout::H;
    case PHOTONICS_DEVICE_FULCRUM: return PhotonicsDataLayout::H;
    case PHOTONICS_DEVICE_BANK_LEVEL: return PhotonicsDataLayout::H;
    case PHOTONICS_DEVICE_AQUABOLT: return PhotonicsDataLayout::H;
    case PHOTONICS_DEVICE_AIM: return PhotonicsDataLayout::H;
    case PHOTONICS_FUNCTIONAL:
    case PHOTONICS_DEVICE_NONE: return PhotonicsDataLayout::UNKNOWN;
  }
  return PhotonicsDataLayout::UNKNOWN;
}

//! @brief  Thread pool ctor
photonicsUtils::threadPool::threadPool(size_t numThreads)
  : m_terminate(false),
    m_workersRemaining(0)
{
  // reserve one thread for main program
  for (size_t i = 1; i < numThreads; ++i) {
    m_threads.emplace_back([this] { workerThread(); });
  }
  std::printf("PHOTONICS-Info: Created thread pool with %lu threads.\n", m_threads.size());
}

//! @brief  Thread pool dtor
photonicsUtils::threadPool::~threadPool()
{
  {
    std::unique_lock<std::mutex> lock(m_mutex);
    m_terminate = true;
  }
  m_cond.notify_all();
  for (auto& thread : m_threads) {
    if (thread.joinable()) {
      thread.join();
    }
  }
}

//! @brief  Entry to process workers in MT
void
photonicsUtils::threadPool::doWork(const std::vector<photonicsUtils::threadWorker*>& workers)
{
  {
    std::unique_lock<std::mutex> lock(m_mutex);
    for (auto& worker : workers) {
      m_workers.push(worker);
    }
    m_workersRemaining = workers.size();
  }
  m_cond.notify_all();

  // Wait for all workers to be done
  std::unique_lock<std::mutex> lock(m_mutex);
  m_cond.wait(lock, [this] { return m_workersRemaining == 0; });
}

//! @brief  Worker thread that process workers
void
photonicsUtils::threadPool::workerThread()
{
  while (true) {
    threadWorker* worker;
    {
      std::unique_lock<std::mutex> lock(m_mutex);
      m_cond.wait(lock, [this] { return m_terminate || !m_workers.empty(); });
      if (m_terminate && m_workers.empty()) {
        return;
      }
      worker = m_workers.front();
      m_workers.pop();
    }
    worker->execute();
    {
      std::unique_lock<std::mutex> lock(m_mutex);
      --m_workersRemaining;
    }
    m_cond.notify_all();
  }
}

//! @brief Helper function to trim from the start (left) of the string
std::string&
photonicsUtils::ltrim(std::string& s) {
  s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch) {
    return !std::isspace(ch);
  }));
  return s;
 }

//! @brief Helper function to trim from the end (right) of the string
std::string&
photonicsUtils::rtrim(std::string& s) {
  s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) {
    return !std::isspace(ch);
  }).base(), s.end());
  return s;
}

//! @brief Function to trim from both ends
std::string&
photonicsUtils::trim(std::string& s) {
  return ltrim(rtrim(s));
}


//! @brief Reads the content of a file and stores it in a provided std::string reference
bool
photonicsUtils::readFileContent(const char* fileName, std::string& fileContent) {
    std::ifstream fileStream(fileName);

    if (!fileStream.is_open()) {
        printf("PHOTONICS-Error: Could not open the file: %s\n", fileName);
        return false;
    }

    std::stringstream buffer;
    buffer << fileStream.rdbuf();
    fileContent = buffer.str();
    return true;
}

//! @brief Retrieves the value of the specified environment variable.
bool
photonicsUtils::getEnvVar(const std::string &name, std::string &value) {
    const char* evnVal = std::getenv(name.c_str());
    if (evnVal == nullptr) {
        return false;
    } else {
        value = evnVal;
        return true;
    }
}

//! @brief Returns the values of each parameter in the config files.
std::string
photonicsUtils::getParam(const std::unordered_map<std::string, std::string>& params, const std::string& key) {
  auto it = params.find(key);
  if (it == params.end()) {
    throw std::invalid_argument(key);
  }
  return it->second;
}

//! @brief Returns the values of each parameter in the config files. Return empty string and false return status if the parameter value is not found
std::string 
photonicsUtils::getOptionalParam(const std::unordered_map<std::string, std::string>& params, const std::string& key, bool& returnStatus) {
  returnStatus = false;  
  auto it = params.find(key);
  if (it == params.end()) {
    return "";
  }
  returnStatus = true;
  return it->second;
} 

//! @brief Returns a substring from the beginning of the input string up to the first ';' character, or the entire string if ';' is not found
std::string 
photonicsUtils::removeAfterSemicolon(const std::string &input) {
    size_t pos = input.find(';');
    if (pos != std::string::npos) {
        return input.substr(0, pos);
    }
    return input;
}

//! @brief Returns the directory path of the input file.
std::string
photonicsUtils::getDirectoryPath(const std::string& filePath) {
    std::filesystem::path path(filePath);
    return path.parent_path().string() + "/";
}

//! @brief Convert a string to unsigned int. Return false and 0 if invalid
bool
photonicsUtils::convertStringToUnsigned(const std::string& str, unsigned& retVal)
{
  try {
    unsigned long val = std::stoul(str);
    if (val > std::numeric_limits<unsigned int>::max()) { // out of range
      throw std::out_of_range("Value exceeds unsigned int range");
    }
    retVal = static_cast<unsigned int>(val);
  } catch (const std::exception &e) {
    retVal = 0;
    return false;
  }
  return true;
}

//! @brief Given a config file path, read all parameters
std::unordered_map<std::string, std::string>
photonicsUtils::readParamsFromConfigFile(const std::string& configFilePath)
{
  std::unordered_map<std::string, std::string> params;
  if (configFilePath.empty()) {
    return params;
  }
  std::string contents;
  bool success = photonicsUtils::readFileContent(configFilePath.c_str(), contents);
  if (!success) {
    std::printf("PHOTONICS-Error: Failed to read config file %s\n", configFilePath.c_str());
    return params;
  }
  std::istringstream iss(contents);
  std::string line;
  while (std::getline(iss, line)) {
    line = photonicsUtils::removeAfterSemicolon(line);
    if (line.empty() || line[0] == '[') {
      continue;
    }
    size_t equalPos = line.find('=');
    if (equalPos != std::string::npos) {
      std::string key = line.substr(0, equalPos);
      std::string value = line.substr(equalPos + 1);
      params[photonicsUtils::trim(key)] = photonicsUtils::trim(value);
    }
  }
  return params;
}

//! @brief Read environment variables, given a list of env var names
std::unordered_map<std::string, std::string>
photonicsUtils::readParamsFromEnvVars(const std::vector<std::string>& envVarNames)
{
  std::unordered_map<std::string, std::string> params;
  for (const auto& envVar : envVarNames) {
    std::string val;
    bool hasEnv = photonicsUtils::getEnvVar(envVar, val);
    val = photonicsUtils::trim(val);
    if (hasEnv && !val.empty()) {
      params[envVar] = val;
    }
  }
  return params;
}

