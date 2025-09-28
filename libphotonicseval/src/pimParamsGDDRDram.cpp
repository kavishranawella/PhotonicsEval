// File: photonicsParamsGDDRDram.cc
// PHOTONICSeval Simulator - GDDR DRAM parameters
// Copyright (c) 2024 University of Virginia
// This file is licensed under the MIT License.
// See the LICENSE file in the root of this repository for more details.

#include "pimParamsGDDRDram.h"
#include "photonicsUtils.h"
#include <sstream>
#include <string>
#include <algorithm>
#include <cctype>
#include <locale>
#include <stdexcept>

//! @brief  photonicsParamsGDDRDram ctor (based on GDDR6_8Gb_x16.ini from DRAMsim3)
photonicsParamsGDDRDram::photonicsParamsGDDRDram()
  : // [dram_structure]
    m_protocol("GDDR6"),
    m_bankgroups(4),
    m_banksPerGroup(4),
    m_rows(16384),
    m_columns(128),
    m_deviceWidth(16),
    m_BL(16),

    // [timing]
    m_tCK(0.66),
    m_CL(24),
    m_CWL(16),
    m_tRCDRD(24),
    m_tRCDWR(20),
    m_tRP(15),
    m_tRAS(32),
    m_tRFC(392),
    m_tREFI(8660),
    m_tRPRE(1),
    m_tWPRE(1),
    m_tRRD_S(9),
    m_tRRD_L(9),
    m_tWTR_S(7),
    m_tWTR_L(7),
    m_tFAW(32),
    m_tWR(16),
    m_tXS(132),
    m_tXP(12),
    m_tRTRS(1),
    m_tRTP_L(7),
    m_tRTP_S(7),
    m_tCCD_S(3),
    m_tCCD_L(4),
    m_tCKESR(7),
    m_tPPD(2),
    m_t32AW(420),
    m_RFCb(30),
    m_tREFIb(238),

    // [power]
    m_VDD(1.35),
    m_IDD0(71),
    m_IDD2P(45),
    m_IDD2N(60),
    m_IDD3P(50),
    m_IDD3N(61),
    m_IDD4W(231),
    m_IDD4R(248),
    m_IDD5AB(286),
    m_IDD5PB(45),
    m_IDD6x(35),
    
    // [system]
    m_channelSize(4096),
    m_channels(1),
    m_busWidth(128),
    m_addressMapping("rochrababgco"),
    m_queueStructure("PER_BANK"),
    m_rowBufPolicy("OPEN_PAGE"),
    m_cmdQueueSize(8),
    m_transQueueSize(32),
    // [other]
    m_epochPeriod(1499250),
    m_outputLevel(1)
{
}

//! @brief  photonicsParamsGDDRDram ctor with a config file
photonicsParamsGDDRDram::photonicsParamsGDDRDram(std::unordered_map<std::string, std::string> params)
{
  try {
    m_protocol = photonicsUtils::getParam(params, "protocol");
    m_bankgroups = std::stoi(photonicsUtils::getParam(params, "bankgroups"));
    m_banksPerGroup = std::stoi(photonicsUtils::getParam(params, "banks_per_group"));
    m_rows = std::stoi(photonicsUtils::getParam(params, "rows"));
    m_columns = std::stoi(photonicsUtils::getParam(params, "columns"));
    m_deviceWidth = std::stoi(photonicsUtils::getParam(params, "device_width"));
    m_BL = std::stoi(photonicsUtils::getParam(params, "BL"));

    m_tCK = std::stod(photonicsUtils::getParam(params, "tCK"));
    m_CL = std::stoi(photonicsUtils::getParam(params, "CL"));
    m_CWL = std::stoi(photonicsUtils::getParam(params, "CWL"));
    m_tRCDRD = std::stoi(photonicsUtils::getParam(params, "tRCDRD"));
    m_tRCDWR = std::stoi(photonicsUtils::getParam(params, "tRCDWR"));
    m_tRP = std::stoi(photonicsUtils::getParam(params, "tRP"));
    m_tRAS = std::stoi(photonicsUtils::getParam(params, "tRAS"));
    m_tRFC = std::stoi(photonicsUtils::getParam(params, "tRFC"));
    m_tREFI = std::stoi(photonicsUtils::getParam(params, "tREFI"));
    m_tRPRE = std::stoi(photonicsUtils::getParam(params, "tRPRE"));
    m_tWPRE = std::stoi(photonicsUtils::getParam(params, "tWPRE"));
    m_tRRD_S = std::stoi(photonicsUtils::getParam(params, "tRRD_S"));
    m_tRRD_L = std::stoi(photonicsUtils::getParam(params, "tRRD_L"));
    m_tWTR_S = std::stoi(photonicsUtils::getParam(params, "tWTR_S"));
    m_tWTR_L = std::stoi(photonicsUtils::getParam(params, "tWTR_L"));
    m_tFAW = std::stoi(photonicsUtils::getParam(params, "tFAW"));
    m_tWR = std::stoi(photonicsUtils::getParam(params, "tWR"));
    m_tCCD_S = std::stoi(photonicsUtils::getParam(params, "tCCD_S"));
    m_tCCD_L = std::stoi(photonicsUtils::getParam(params, "tCCD_L"));
    m_tCKESR = std::stoi(photonicsUtils::getParam(params, "tCKESR"));
    m_tXS = std::stoi(photonicsUtils::getParam(params, "tXS"));
    m_tXP = std::stoi(photonicsUtils::getParam(params, "tXP"));
    m_tRTRS = std::stoi(photonicsUtils::getParam(params, "tRTRS"));
    m_tRTP_L = std::stoi(photonicsUtils::getParam(params, "tRTP_L"));
    m_tRTP_S = std::stoi(photonicsUtils::getParam(params, "tRTP_S"));
    m_tPPD = std::stoi(photonicsUtils::getParam(params, "tPPD"));
    m_t32AW = std::stoi(photonicsUtils::getParam(params, "t32AW"));
    m_RFCb = std::stoi(photonicsUtils::getParam(params, "tRFCb"));
    m_tREFIb = std::stoi(photonicsUtils::getParam(params, "tREFIb"));

    m_VDD = std::stod(photonicsUtils::getParam(params, "VDD"));
    m_IDD0 = std::stoi(photonicsUtils::getParam(params, "IDD0"));
    m_IDD2P = std::stoi(photonicsUtils::getParam(params, "IDD2P"));
    m_IDD2N = std::stoi(photonicsUtils::getParam(params, "IDD2N"));
    m_IDD3P = std::stoi(photonicsUtils::getParam(params, "IDD3P"));
    m_IDD3N = std::stoi(photonicsUtils::getParam(params, "IDD3N"));
    m_IDD4W = std::stoi(photonicsUtils::getParam(params, "IDD4W"));
    m_IDD4R = std::stoi(photonicsUtils::getParam(params, "IDD4R"));
    m_IDD5AB = std::stoi(photonicsUtils::getParam(params, "IDD5AB"));
    m_IDD6x = std::stoi(photonicsUtils::getParam(params, "IDD6x"));

    m_channelSize = std::stoi(photonicsUtils::getParam(params, "channel_size"));
    m_channels = std::stoi(photonicsUtils::getParam(params, "channels"));
    m_busWidth = std::stoi(photonicsUtils::getParam(params, "bus_width"));
    m_addressMapping = photonicsUtils::getParam(params, "address_mapping");
    m_queueStructure = photonicsUtils::getParam(params, "queue_structure");
    m_rowBufPolicy = photonicsUtils::getParam(params, "row_buf_policy");
    m_cmdQueueSize = std::stoi(photonicsUtils::getParam(params, "cmd_queue_size"));
    m_transQueueSize = std::stoi(photonicsUtils::getParam(params, "trans_queue_size"));

    m_epochPeriod = std::stoi(photonicsUtils::getParam(params, "epoch_period"));
    m_outputLevel = std::stoi(photonicsUtils::getParam(params, "output_level"));
  } catch (const std::invalid_argument& e) {
    std::string errorMessage("PHOTONICS-Error: Missing or invalid parameter: ");
    errorMessage += e.what();
    errorMessage += "\n";
    throw std::invalid_argument(errorMessage);
  }
}

