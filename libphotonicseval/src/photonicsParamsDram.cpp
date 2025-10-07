// File: photonicsParamsDram.cc
// PHOTONICSeval Simulator - DRAM parameters
// Copyright (c) 2024 University of Virginia
// This file is licensed under the MIT License.
// See the LICENSE file in the root of this repository for more details.

#include "photonicsParamsDram.h"
#include "photonicsUtils.h"
#include "pimParamsDDRDram.h"
#include "pimParamsLPDDRDram.h"
#include "photonicsParamsHBMDram.h"
#include "pimParamsGDDRDram.h"
#include <string>
#include <algorithm>
#include <cctype>
#include <locale>
#include <stdexcept>

// Static factory method to create appropriate subclass based on protocol enum
std::unique_ptr<photonicsParamsDram> photonicsParamsDram::create(PhotonicsDeviceProtocolEnum deviceProtocol)
{
  if (deviceProtocol == PHOTONICS_DEVICE_PROTOCOL_DDR)
  {
    return std::make_unique<photonicsParamsDDRDram>();
  }
  else if (deviceProtocol == PHOTONICS_DEVICE_PROTOCOL_LPDDR)
  {
    return std::make_unique<photonicsParamsLPDDRDram>();
  }
  else if (deviceProtocol == PHOTONICS_DEVICE_PROTOCOL_HBM)
  {
    return std::make_unique<photonicsParamsHBMDram>();
  } 
  else if (deviceProtocol == PHOTONICS_DEVICE_PROTOCOL_GDDR)
  {
    return std::make_unique<photonicsParamsGDDRDram>();
  }
  else
  {
    std::string errorMessage("PHOTONICS-Error: Inavalid DRAM protocol parameter.\n");
    throw std::invalid_argument(errorMessage);
  }
}

// Static factory method to create appropriate subclass based on config file
std::unique_ptr<photonicsParamsDram> photonicsParamsDram::createFromConfig(const std::string& memConfigFilePath)
{
  std::unordered_map<std::string, std::string> params = photonicsUtils::readParamsFromConfigFile(memConfigFilePath);

  // Check if the "protocol" key exists
  if (params.find("protocol") == params.end())
  {
    std::string errorMessage("PHOTONICS-Error: Missing DRAM protocol parameter.\n");
    throw std::invalid_argument(errorMessage);
  }

  // Extract protocol from params
  std::string deviceProtocol = params["protocol"];

  // Instantiate the appropriate subclass based on the protocol
  if (deviceProtocol == "DDR3" || deviceProtocol == "DDR4" || deviceProtocol == "DDR5")
  {
    return std::make_unique<photonicsParamsDDRDram>(params);
  } 
  else if (deviceProtocol == "LPDDR3" || deviceProtocol == "LPDDR4") {
    return std::make_unique<photonicsParamsLPDDRDram>(params);
  } 
  else if (deviceProtocol == "HBM" || deviceProtocol == "HBM2") {
    return std::make_unique<photonicsParamsHBMDram>(params);
  }
  else if (deviceProtocol == "GDDR5" || deviceProtocol == "GDDR5X" || deviceProtocol == "GDDR6")
  {
    return std::make_unique<photonicsParamsGDDRDram>(params);
  }
  else
  {
    throw std::invalid_argument("Unknown protocol: " + deviceProtocol);
  }
}

