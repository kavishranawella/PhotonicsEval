// File: photonicsPerfEnergyTables.h
// PHOTONICSeval Simulator - Performance and Energy Tables
// Copyright (c) 2024 University of Virginia
// This file is licensed under the MIT License.
// See the LICENSE file in the root of this repository for more details.

#ifndef LAVA_PHOTONICS_PERF_ENERGY_TABLES_H
#define LAVA_PHOTONICS_PERF_ENERGY_TABLES_H

#include "libphotonicseval.h"
#include "photonicsCmd.h"
#include <unordered_map>
#include <tuple>


namespace photonicsPerfEnergyTables
{
  // Perf-energy table of BitSIMD-V variants
  extern const std::unordered_map<PhotonicsDeviceEnum, std::unordered_map<PhotonicsDataType,
      std::unordered_map<PhotonicsCmdEnum, std::tuple<unsigned, unsigned, unsigned>>>> bitsimdPerfTable;
}

#endif

