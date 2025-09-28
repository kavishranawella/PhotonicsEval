// File: photonicsCmdFuse.cpp
// PHOTONICSeval Simulator - PHOTONICS API Fusion
// Copyright (c) 2024 University of Virginia
// This file is licensed under the MIT License.
// See the LICENSE file in the root of this repository for more details.

#include "photonicsCmdFuse.h"
#include <cstdio>


//! @brief  Photonics CMD: PHOTONICS API Fusion
bool
photonicsCmdFuse::execute()
{
  if (m_debugCmds) {
    std::printf("PHOTONICS-Cmd: API Fusion\n");
  }

  // Functional simulation
  // TODO: skip original updateStats
  bool success = true;
  for (auto& api : m_prog.m_apis) {
    PhotonicsStatus status = api();
    if (status != PHOTONICS_OK) {
      success = false;
      break;
    }
  }

  // Analyze API fusion opportunities
  success = success && updateStats();
  return success;
}

//! @brief  Photonics CMD: PHOTONICS API Fusion - update stats
bool
photonicsCmdFuse::updateStats() const
{
  // TODO: Parse m_prog and update stats
  return true;
}

