// File: photonicsCmdFuse.h
// PHOTONICSeval Simulator - PHOTONICS API Fusion
// Copyright (c) 2024 University of Virginia
// This file is licensed under the MIT License.
// See the LICENSE file in the root of this repository for more details.

#ifndef LAVA_PHOTONICS_CMD_FUSE_H
#define LAVA_PHOTONICS_CMD_FUSE_H

#include "libphotonicseval.h"
#include "photonicsCmd.h"


//! @class  photonicsCmdFuse
//! @brief  Photonics CMD: PHOTONICS API Fusion
class photonicsCmdFuse : public photonicsCmd
{
public:
  photonicsCmdFuse(PhotonicsProg prog) : photonicsCmd(PhotonicsCmdEnum::NOOP), m_prog(prog) {}
  virtual ~photonicsCmdFuse() {}
  virtual bool execute() override;
  virtual bool updateStats() const override;
private:
  PhotonicsProg m_prog;
};

#endif

