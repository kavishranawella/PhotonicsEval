// File: photonicsSimConfig.h
// PHOTONICSeval Simulator - PHOTONICS Simulator Configurations
// Copyright (c) 2024 University of Virginia
// This file is licensed under the MIT License.
// See the LICENSE file in the root of this repository for more details.

#ifndef LAVA_PHOTONICS_SIM_CONFIG_H
#define LAVA_PHOTONICS_SIM_CONFIG_H

#include "libphotonicseval.h"
#include <string>
#include <unordered_map>
#include <vector>


//! @class  photonicsSimConfig
//! @brief  PHOTONICS simulator configurations
//! This class manages global static PHOTONICSeval configurations from config files, env vars, or API parameters,
//! which are for controlling PHOTONICSeval library behavior once it is linked with PHOTONICS application executables.
//!
//! Supported configuration file parameters:
//!   memory_config_file = <ini-file>            // memory config file, e.g., DDR4_8Gb_x16_3200.ini
//!   simulation_target = <PhotonicsDeviceEnum>        // simulation target, e.g., PHOTONICS_DEVICE_BITSIMD_V
//!   num_ranks = <int>                          // number of ranks
//!   num_bank_per_rank = <int>                  // number of banks per rank
//!   num_subarray_per_bank = <int>              // number of subarrays per bank
//!   num_row_per_subarray = <int>               // number of rows per subarray
//!   num_col_per_subarray = <int>               // number of columns per subarray
//!   max_num_threads = <int>                    // maximum number of threads used by simulation
//!   should_load_balance = <0|1>                // distribute data evenly among all cores
//!
//! Supported environment variables:
//!   PHOTONICSEVAL_SIM_CONFIG <abs-path/cfg-file>     // PHOTONICSeval config file, e.g., abs-path/PHOTONICSeval_BitSimdV.cfg
//!   PHOTONICSEVAL_MEM_CONFIG <abs-path/ini-file>     // memory config file, e.g., DDR4_8Gb_x16_3200.ini
//!   PHOTONICSEVAL_SIM_TARGET <PhotonicsDeviceEnum>         // simulation target, e.g., PHOTONICS_DEVICE_BITSIMD_V
//!   PHOTONICSEVAL_NUM_RANKS <int>                    // number of ranks
//!   PHOTONICSEVAL_NUM_BANK_PER_RANK <int>            // number of banks per rank
//!   PHOTONICSEVAL_NUM_SUBARRAY_PER_BANK <int>        // number of subarrays per bank
//!   PHOTONICSEVAL_NUM_ROW_PER_SUBARRAY <int>         // number of rows per subarray
//!   PHOTONICSEVAL_NUM_COL_PER_SUBARRAY <int>         // number of columns per subarray
//!   PHOTONICSEVAL_MAX_NUM_THREADS <int>              // maximum number of threads used by simulation
//!   PHOTONICSEVAL_ANALYSIS_MODE <0|1>                // PHOTONICSeval analysis mode
//!   PHOTONICSEVAL_DEBUG <int>                        // PHOTONICSeval debug flags (see enum photonicsDebugFlags)
//!   PHOTONICSEVAL_LOAD_BALANCE <0|1>                 // distribute data evenly among all cores
//!
//! Precedence rules (highest to lowest priority):
//! * Config file: Either from -c command-line argument or from PHOTONICSEVAL_SIM_CONFIG
//! * Environment variables
//! * Parameters from photonicsCreateDevice API or C++ macros
//!
//! About config file paths:
//! * Simulator config file
//!   - If passing it through -c command line argument, it's already a valid absolute or relative path
//!   - If using the PHOTONICSEVAL_SIM_CONFIG env variable, its value needs to be a valid absolute path
//! * Memory config file
//!   - If it is in the same directory as the simulator config file, specifying a file name is sufficient
//!   - Otherwise, use absolute path
//!
//! How to add a new configuration parameter:
//! * Define a key string below as m_cfgVar* and/or m_envVar*, and update the documentation above
//! * Define a member variable, and update the uninit and show function
//! * Add a private derive function (can reuse deriveMiscEnvVars if it's env only)
//! * Add a public getter function
//! * Add env var to readEnvVars function
//!
class photonicsSimConfig
{
public:
  photonicsSimConfig() { reset(); }
  ~photonicsSimConfig() {}

  // Update PHOTONICSeval simulation configuration parameters at device creation
  bool init(PhotonicsDeviceEnum deviceType, unsigned numRanks, unsigned numBankPerRank, unsigned numSubarrayPerBank, 
      unsigned numRowPerSubarray, unsigned numColPerSubarray, unsigned bufferSize, unsigned matrixSize);
  bool init(PhotonicsDeviceEnum deviceType, const std::string& configFilePath);
  void uninit() { reset(); }
  bool isInit() const { return m_isInit; }
  void show() const;

  // Getters
  const std::string& getSimConfigFile() const { return m_simConfigFile; }
  const std::string& getMemConfigFile() const { return m_memConfigFile; }
  PhotonicsDeviceEnum getDeviceType() const { return m_deviceType; }
  PhotonicsDeviceEnum getSimTarget() const { return m_simTarget; }
  PhotonicsDeviceProtocolEnum getMemoryProtocol() const { return m_memoryProtocol; }
  unsigned getNumRanks() const { return m_numRanks; }
  unsigned getNumBankPerRank() const { return m_numBankPerRank; }
  unsigned getNumSubarrayPerBank() const { return m_numSubarrayPerBank; }
  unsigned getNumRowPerSubarray() const { return m_numRowPerSubarray; }
  unsigned getNumColPerSubarray() const { return m_numColPerSubarray; }
  unsigned getNumThreads() const { return m_numThreads; }
  unsigned getBufferSize() const { return m_bufferSize; }
  unsigned getMatrixSize() const { return m_matrixSize; }
  bool isAnalysisMode() const { return m_analysisMode; }
  unsigned getDebug() const { return m_debug; }
  bool isLoadBalanced() const { return m_loadBalanced; }

  enum photonicsDebugFlags
  {
    DEBUG_PARAMS      = 0x0001,
    DEBUG_API_CALLS   = 0x0002,
    DEBUG_CMDS        = 0x0004,
    DEBUG_ALLOC       = 0x0008,
    DEBUG_PERF        = 0x0010,
  };

private:
  bool deriveConfig(PhotonicsDeviceEnum deviceType,
      const std::string& configFilePath = "",
      unsigned numRanks = 0,
      unsigned numBankPerRank = 0,
      unsigned numSubarrayPerBank = 0,
      unsigned numRowPerSubarray = 0,
      unsigned numColPerSubarray = 0,
      unsigned bufferSize = 0,
      unsigned matrixSize = 0);

  bool deriveDebug();
  std::unordered_map<std::string, std::string> readEnvVars() const;
  bool deriveSimConfigFile(const std::string& configFilePath);
  std::unordered_map<std::string, std::string> readSimConfigFileParams() const;
  bool deriveDeviceType(PhotonicsDeviceEnum deviceType);
  bool deriveSimTarget();
  bool deriveMemConfigFile();
  bool deriveDimension(const std::string& envVar, const std::string& cfgVar, const unsigned apiVal, const unsigned defVal, unsigned& retVal);
  bool deriveDimensions(unsigned numRanks, unsigned numBankPerRank, unsigned numSubarrayPerBank, unsigned numRowPerSubarray, unsigned numColPerSubarray, unsigned bufferSize, unsigned matrixSize);
  bool deriveNumThreads();
  bool deriveMiscEnvVars();
  bool deriveLoadBalance();

  bool parseConfigFromFile(const std::string& config, unsigned& numRanks, unsigned& numBankPerRank, unsigned& numSubarrayPerBank, unsigned& numRows, unsigned& numCols);

  // Config file parameters
  inline static const std::string m_cfgVarMemConfig = "memory_config_file";
  inline static const std::string m_cfgVarSimTarget = "simulation_target";
  inline static const std::string m_cfgVarNumRanks = "num_ranks";
  inline static const std::string m_cfgVarNumBankPerRank = "num_bank_per_rank";
  inline static const std::string m_cfgVarNumSubarrayPerBank = "num_subarray_per_bank";
  inline static const std::string m_cfgVarNumRowPerSubarray = "num_row_per_subarray";
  inline static const std::string m_cfgVarNumColPerSubarray = "num_col_per_subarray";
  inline static const std::string m_cfgVarMaxNumThreads = "max_num_threads";
  inline static const std::string m_cfgVarLoadBalance = "should_load_balance";
  inline static const std::string m_cfgVarBufferSize = "buffer_size";
  inline static const std::string m_cfgVarMatrixSize = "matrix_size";

  // Environment variables
  inline static const std::string m_envVarSimConfig = "PHOTONICSEVAL_SIM_CONFIG";
  inline static const std::string m_envVarMemConfig = "PHOTONICSEVAL_MEM_CONFIG";
  inline static const std::string m_envVarSimTarget = "PHOTONICSEVAL_SIM_TARGET";
  inline static const std::string m_envVarNumRanks = "PHOTONICSEVAL_NUM_RANKS";
  inline static const std::string m_envVarNumBankPerRank = "PHOTONICSEVAL_NUM_BANK_PER_RANK";
  inline static const std::string m_envVarNumSubarrayPerBank = "PHOTONICSEVAL_NUM_SUBARRAY_PER_BANK";
  inline static const std::string m_envVarNumRowPerSubarray = "PHOTONICSEVAL_NUM_ROW_PER_SUBARRAY";
  inline static const std::string m_envVarNumColPerSubarray = "PHOTONICSEVAL_NUM_COL_PER_SUBARRAY";
  inline static const std::string m_envVarBufferSize = "PHOTONICSEVAL_BUFFER_SIZE";
  inline static const std::string m_envVarMatrixSize = "PHOTONICSEVAL_MATRIX_SIZE";
  inline static const std::string m_envVarMaxNumThreads = "PHOTONICSEVAL_MAX_NUM_THREADS";
  inline static const std::string m_envVarAnalysisMode = "PHOTONICSEVAL_ANALYSIS_MODE";
  inline static const std::string m_envVarDebug = "PHOTONICSEVAL_DEBUG";
  inline static const std::string m_envVarLoadBalance = "PHOTONICSEVAL_LOAD_BALANCE";

  // Add env vars to this list for readEnvVars
  inline static const std::vector<std::string> m_envVarList = {
    m_envVarSimConfig,
    m_envVarMemConfig,
    m_envVarSimTarget,
    m_envVarNumRanks,
    m_envVarNumBankPerRank,
    m_envVarNumSubarrayPerBank,
    m_envVarNumRowPerSubarray,
    m_envVarNumColPerSubarray,
    m_envVarMaxNumThreads,
    m_envVarAnalysisMode,
    m_envVarDebug,
    m_envVarLoadBalance,
    m_envVarBufferSize,
    m_envVarMatrixSize,
  };

  // Default values if not specified during init
  static constexpr int DEFAULT_NUM_RANKS = 2;
  static constexpr int DEFAULT_NUM_BANK_PER_RANK = 2;
  static constexpr int DEFAULT_NUM_SUBARRAY_PER_BANK = 128;
  static constexpr int DEFAULT_NUM_ROW_PER_SUBARRAY = 3;
  static constexpr int DEFAULT_NUM_COL_PER_SUBARRAY = 4096;
  static constexpr int DEFAULT_BUFFER_SIZE = 0;
  static constexpr int DEFAULT_MATRIX_SIZE = 4096;
  static constexpr PhotonicsDeviceEnum DEFAULT_SIM_TARGET = PHOTONICS_DEVICE_FULCRUM;

  //! @brief  Reset all member variables to default status
  inline void reset() {
    m_simConfigFile.clear();
    m_memConfigFile.clear();
    m_deviceType = PHOTONICS_DEVICE_NONE;
    m_simTarget = PHOTONICS_DEVICE_NONE;
    m_memoryProtocol = PHOTONICS_DEVICE_PROTOCOL_DDR;
    m_numRanks = 0;
    m_numBankPerRank = 0;
    m_numSubarrayPerBank = 0;
    m_numRowPerSubarray = 0;
    m_numColPerSubarray = 0;
    m_numThreads = 0;
    m_bufferSize = 0;
    m_matrixSize = 0;
    m_analysisMode = false;
    m_debug = 0;
    m_loadBalanced = false;
    m_envParams.clear();
    m_cfgParams.clear();
    m_isInit = false;
  }

  // PHOTONICS sim env variables
  std::string m_simConfigFile;
  std::string m_memConfigFile;
  PhotonicsDeviceEnum m_deviceType;
  PhotonicsDeviceEnum m_simTarget;
  PhotonicsDeviceProtocolEnum m_memoryProtocol;
  unsigned m_numRanks;
  unsigned m_numBankPerRank;
  unsigned m_numSubarrayPerBank;
  unsigned m_numRowPerSubarray;
  unsigned m_numColPerSubarray;
  unsigned m_numThreads;
  unsigned m_bufferSize;
  unsigned m_matrixSize;
  bool m_analysisMode;
  unsigned m_debug;
  bool m_loadBalanced;

  // Store original parameters for extension purpose
  std::unordered_map<std::string, std::string> m_envParams;
  std::unordered_map<std::string, std::string> m_cfgParams;
  bool m_isInit;
};

#endif

