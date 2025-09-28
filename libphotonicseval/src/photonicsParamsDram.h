// File: photonicsParamsDram.h
// PHOTONICSeval Simulator - DRAM parameters
// Copyright (c) 2024 University of Virginia
// This file is licensed under the MIT License.
// See the LICENSE file in the root of this repository for more details.

#ifndef LAVA_PHOTONICS_PARAMS_DRAM_H
#define LAVA_PHOTONICS_PARAMS_DRAM_H

#include <string>
#include <memory>
#include "libphotonicseval.h"

//! @class  photonicsParamsDram
//! @brief  DRAM parameters (DRAMsim3 compatible)
class photonicsParamsDram
{
public:
  // Virtual destructor to ensure derived class destructors are called
  virtual ~photonicsParamsDram() = default;

  // Static factory method to create appropriate subclass based on protocol
  static std::unique_ptr<photonicsParamsDram> create(PhotonicsDeviceProtocolEnum deviceProtocol);

  // Static factory method to create appropriate subclass based on config file
  static std::unique_ptr<photonicsParamsDram> createFromConfig(const std::string& memConfigFilePath);

  // Virtual functions for protocol-specific implementation
  virtual int getDeviceWidth() const = 0;
  virtual int getBurstLength() const = 0;
  virtual int getNumChipsPerRank() const = 0;
  virtual double getNsRowRead() const = 0;
  virtual double getNsRowActivate() const = 0;
  virtual double getNsRowPrecharge() const = 0;
  virtual double getNsRowWrite() const = 0;
  virtual double getNsTCCD_S() const = 0;
  virtual double getNsTCCD_L() const = 0;
  virtual double getNsTCAS() const = 0;
  virtual double getNsAAP() const = 0;
  virtual double getTypicalRankBW() const = 0;
  virtual double getPjActPre() const = 0;
  virtual double getPjLogic() const = 0;
  virtual double getMwIDD2N() const = 0;
  virtual double getMwIDD3N() const = 0;
  virtual double getPjRead() const = 0;
  virtual double getPjWrite() const = 0;
  virtual double getPjPrecharge() const = 0;
  virtual double getPjActivate() const = 0;
  virtual double gettRCD() const = 0;
  virtual double gettRP() const = 0;
  virtual double gettCCD_L() const = 0;
  virtual double gettCCD_S() const = 0;
  virtual double gettRAS() const = 0; 
  virtual double gettCK() const = 0;
};

#endif

