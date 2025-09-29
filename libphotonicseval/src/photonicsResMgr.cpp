// File: photonicsResMgr.cpp
// PHOTONICSeval Simulator - PHOTONICS Resource Manager
// Copyright (c) 2024 University of Virginia
// This file is licensed under the MIT License.
// See the LICENSE file in the root of this repository for more details.

#include "photonicsResMgr.h"       // for photonicsResMgr
#include "photonicsDevice.h"       // for photonicsDevice
#include <cstdio>            // for printf
#include <algorithm>         // for sort, prev
#include <stdexcept>         // for throw, invalid_argument
#include <memory>            // for make_unique
#include <cassert>           // for assert
#include <string>            // for string


//! @brief  Print info of a PHOTONICS region
void
photonicsRegion::print(uint64_t regionId) const
{
  printf("{ PHOTONICS-Region %lu: CoreId = %d, Loc = (%u, %u), Size = (%u, %u) }\n",
         regionId, m_coreId, m_rowIdx, m_colIdx, m_numAllocRows, m_numAllocCols);
}

//! @brief  Print info of a PHOTONICS object
void
photonicsObjInfo::print() const
{
  std::printf("----------------------------------------\n");
  std::printf("PHOTONICS-Object: ObjId = %d, AllocType = %d, Regions =\n",
              m_objId, static_cast<int>(m_allocType));
  for (size_t i = 0; i < m_regions.size(); ++i) {
    m_regions[i].print(i);
  }
  std::printf("----------------------------------------\n");
}

//! @brief  Finalize obj info
void
photonicsObjInfo::finalize()
{
  std::unordered_map<PhotonicsCoreId, int> coreIdCnt;
  for (const auto& region : m_regions) {
    PhotonicsCoreId coreId = region.getCoreId();
    coreIdCnt[coreId]++;
    unsigned numRegionsPerCore = coreIdCnt[coreId];
    if (m_maxNumRegionsPerCore < numRegionsPerCore) {
      m_maxNumRegionsPerCore = numRegionsPerCore;
    }
  }
  m_numCoresUsed = coreIdCnt.size();
  m_numCoreAvailable = m_device->getNumCores();
  m_isLoadBalanced = m_device->getConfig().isLoadBalanced();

  const photonicsRegion& region = m_regions[0];
  m_maxElementsPerRegion = (uint64_t)region.getNumAllocRows() * region.getNumAllocCols() / m_bitsPerElementPadded;
  m_numColsPerElem = region.getNumColsPerElem();
}

//! @brief  Get number of bits per element
unsigned
photonicsObjInfo::getBitsPerElement(PhotonicsBitWidth bitWidthType) const
{
  switch (bitWidthType) {
    case PhotonicsBitWidth::ACTUAL:
    case PhotonicsBitWidth::SIM:
    case PhotonicsBitWidth::HOST:
      return photonicsUtils::getNumBitsOfDataType(m_dataType, bitWidthType);
    case PhotonicsBitWidth::PADDED:
      return m_bitsPerElementPadded;
    default:
      assert(0);
  }
  return 0;
}

//! @brief  Get all regions on a specific PHOTONICS core for current PHOTONICS object
std::vector<photonicsRegion>
photonicsObjInfo::getRegionsOfCore(PhotonicsCoreId coreId) const
{
  std::vector<photonicsRegion> regions;
  for (const auto& region : m_regions) {
    if (region.getCoreId() == coreId) {
      regions.push_back(region);
    }
  }
  return regions;
}

//! @brief  Copy data from host memory to PHOTONICS object data holder, with ref support
void
photonicsObjInfo::copyFromHost(void* src, uint64_t idxBegin, uint64_t idxEnd)
{
  // handle reference
  if (m_refObjId != -1) {
    photonicsObjInfo& refObj = m_device->getResMgr()->getObjInfo(m_refObjId);
    if (isDualContactRef()) {
      uint64_t numBytes = refObj.m_data.getNumBytes(idxBegin, idxEnd);
      std::vector<float> buffer(numBytes);
      // Cast src to uint8_t* and convert in one pass
      const uint8_t* srcBytes = static_cast<const uint8_t*>(src);
      for (size_t i = 0; i < numBytes; ++i) {
          buffer[i] = srcBytes[i] * 0.00390625f;  // normalize
      }
      refObj.m_data.copyFromHost(buffer.data(), idxBegin, idxEnd);
    } else {
      assert(0); // to be extended
    }
    return;
  }
  m_data.copyFromHost(src, idxBegin, idxEnd);
}

//! @brief  Copy data from PHOTONICS object data holder to host memory, with ref support
void
photonicsObjInfo::copyToHost(void* dest, uint64_t idxBegin, uint64_t idxEnd) const
{
  // handle reference
  if (m_refObjId != -1) {
    photonicsObjInfo &refObj = m_device->getResMgr()->getObjInfo(m_refObjId);
    if (isDualContactRef()) {
      uint64_t numBytes = refObj.m_data.getNumBytes(idxBegin, idxEnd);
      std::vector<float> buffer(numBytes);
      refObj.m_data.copyToHost(buffer.data(), idxBegin, idxEnd);
      // Cast destination pointer
      uint8_t* destBytes = static_cast<uint8_t*>(dest);
      for (size_t i = 0; i < numBytes; ++i) {
          destBytes[i] = static_cast<uint8_t>(buffer[i] * 256.0f);
      }
    } else {
      assert(0); // to be extended
    }
    return;
  }
  m_data.copyToHost(dest, idxBegin, idxEnd);
}

//! @brief  Copy data from a PHOTONICS object data holder to another, with ref support
void
photonicsObjInfo::copyToObj(photonicsObjInfo& destObj, uint64_t idxBegin, uint64_t idxEnd) const
{
  // handle reference
  if (m_refObjId != -1) {
    photonicsObjInfo &refObj = m_device->getResMgr()->getObjInfo(m_refObjId);
    if (isDualContactRef()) {
      uint64_t numBytes = m_data.getNumBytes(idxBegin, idxEnd);
      std::vector<uint8_t> buffer(numBytes);
      m_data.copyToHost(buffer.data(), idxBegin, idxEnd);
      for (auto& byte : buffer) { byte = ~byte; }
      refObj.m_data.copyFromHost(buffer.data(), idxBegin, idxEnd);
    } else {
      assert(0); // to be extended
    }
    return;
  }
  m_data.copyToObj(destObj.m_data, idxBegin, idxEnd);
}

//! @brief  Set an element at index with bit presentation, with ref support
void
photonicsObjInfo::setElementBits(uint64_t index, uint64_t bits)
{
  // handle reference
  if (m_refObjId != -1) {
    photonicsObjInfo& refObj = m_device->getResMgr()->getObjInfo(m_refObjId);
    if (isDualContactRef()) {
      bits = ~bits;
      refObj.m_data.setElementBits(index, bits);
    } else {
      assert(0); // to be extended
    }
    return;
  }
  m_data.setElementBits(index, bits);
}

//! @brief  Get bit representation of an element at index, with ref support
uint64_t
photonicsObjInfo::getElementBits(uint64_t index) const
{
  // handle reference
  if (m_refObjId != -1) {
    photonicsObjInfo& refObj = m_device->getResMgr()->getObjInfo(m_refObjId);
    if (isDualContactRef()) {
      uint64_t bits = 0;
      refObj.m_data.getElementBits(index, bits);
      bits = ~bits;
      return bits;
    } else {
      assert(0); // to be extended
    }
    return 0;
  }
  uint64_t bits = 0;
  m_data.getElementBits(index, bits);
  return bits;
}

//! @brief  Sync PHOTONICS object data from simulated memory
void
photonicsObjInfo::syncFromSimulatedMem()
{
  photonicsObjInfo &obj = (m_refObjId != -1 ? m_device->getResMgr()->getObjInfo(m_refObjId) : *this);
  unsigned numBits = getBitsPerElement(PhotonicsBitWidth::SIM);
  for (size_t i = 0; i < m_regions.size(); ++i) {
    photonicsRegion& region = m_regions[i];
    PhotonicsCoreId coreId = region.getCoreId();
    photonicsCore& core = m_device->getCore(coreId);
    uint64_t elemIdxBegin = region.getElemIdxBegin();
    uint64_t numElemInRegion = region.getNumElemInRegion();
    for (uint64_t j = 0; j < numElemInRegion; ++j) {
      auto [rowLoc, colLoc] = region.locateIthElemInRegion(j);
      uint64_t bits = isVLayout() ? core.getBitsV(rowLoc, colLoc, numBits)
                                  : core.getBitsH(rowLoc, colLoc, numBits);
      obj.m_data.setElementBits(elemIdxBegin + j, bits);
    }
  }
}

//! @brief  Sync PHOTONICS object data to simulated memory
void
photonicsObjInfo::syncToSimulatedMem() const
{
  const photonicsObjInfo &obj = (m_refObjId != -1 ? m_device->getResMgr()->getObjInfo(m_refObjId) : *this);
  unsigned numBits = getBitsPerElement(PhotonicsBitWidth::SIM);
  for (size_t i = 0; i < m_regions.size(); ++i) {
    const photonicsRegion& region = m_regions[i];
    PhotonicsCoreId coreId = region.getCoreId();
    photonicsCore& core = m_device->getCore(coreId);
    uint64_t elemIdxBegin = region.getElemIdxBegin();
    uint64_t numElemInRegion = region.getNumElemInRegion();
    for (uint64_t j = 0; j < numElemInRegion; ++j) {
      uint64_t bits = 0;
      obj.m_data.getElementBits(elemIdxBegin + j, bits);
      auto [rowLoc, colLoc] = region.locateIthElemInRegion(j);
      if (isVLayout()) {
        core.setBitsV(rowLoc, colLoc, bits, numBits);
      } else {
        core.setBitsH(rowLoc, colLoc, bits, numBits);
      }
    }
  }
}


//! @brief  photonicsResMgr ctor
photonicsResMgr::photonicsResMgr(photonicsDevice* device)
  : m_device(device),
    m_availObjId(0)
{
  unsigned numCores = m_device->getNumCores();
  unsigned numRowsPerCore = m_device->getNumRows();
  for (unsigned i = 0; i < numCores; ++i) {
    m_coreUsage[i] = std::make_unique<coreUsage>(numRowsPerCore);
  }
  m_debugAlloc = (m_device->getConfig().getDebug() & photonicsSimConfig::DEBUG_ALLOC);
}

//! @brief  photonicsResMgr dtor
photonicsResMgr::~photonicsResMgr()
{
}

//! @brief  Allocate a new PHOTONICS object
//!         For V layout, dataType determines the number of rows per region
//!         For H layout, dataType determines the number of bits per element
PhotonicsObjId
photonicsResMgr::photonicsAlloc(PhotonicsAllocEnum allocType, uint64_t numElements, PhotonicsDataType dataType)
{
  if (m_debugAlloc) {
    printf("PHOTONICS-Debug: photonicsAlloc: Request: %s %lu elements of type %s\n",
           photonicsUtils::photonicsAllocEnumToStr(allocType).c_str(), numElements,
           photonicsUtils::photonicsDataTypeEnumToStr(dataType).c_str());
  }

  if (numElements == 0) {
    printf("PHOTONICS-Error: photonicsAlloc: Invalid input parameter: 0 element\n");
    return -1;
  }

  unsigned bitsPerElement = photonicsUtils::getNumBitsOfDataType(dataType, PhotonicsBitWidth::SIM);

  std::vector<PhotonicsCoreId> sortedCoreId = getCoreIdsSortedByLeastUsage();
  photonicsObjInfo newObj(m_availObjId, dataType, allocType, numElements, bitsPerElement, m_device);
  m_availObjId++;

  unsigned numCores = m_device->getNumCores();
  unsigned numCols = m_device->getNumCols();
  unsigned numRowsToAlloc = 0;
  uint64_t numRegions = 0;
  unsigned numColsToAllocLast = 0;
  uint64_t numElemPerRegion = 0;
  uint64_t numElemPerRegionLast = 0;
  unsigned numColsPerElem = 0;
  if (allocType == PHOTONICS_ALLOC_V || allocType == PHOTONICS_ALLOC_V1) {
    // allocate one region per core, with vertical layout
    numRowsToAlloc = bitsPerElement;
    numRegions = (numElements - 1) / numCols + 1;
    numColsToAllocLast = numElements % numCols;
    if (numColsToAllocLast == 0) {
      numColsToAllocLast = numCols;
    }
    numElemPerRegion = numCols;
    numElemPerRegionLast = numColsToAllocLast;
    numColsPerElem = 1;
  } else if (allocType == PHOTONICS_ALLOC_H || allocType == PHOTONICS_ALLOC_H1) {
    // allocate one region per core, with horizontal layout
    numRowsToAlloc = 1;
    numRegions = (numElements * bitsPerElement - 1) / numCols + 1;
    numColsToAllocLast = (numElements * bitsPerElement) % numCols;
    if (numColsToAllocLast == 0) {
      numColsToAllocLast = numCols;
    }
    numElemPerRegion = numCols / bitsPerElement;
    numElemPerRegionLast = numColsToAllocLast / bitsPerElement;
    numColsPerElem = bitsPerElement;
  } else {
    printf("PHOTONICS-Error: photonicsAlloc: Unsupported allocation type %s\n",
           photonicsUtils::photonicsAllocEnumToStr(allocType).c_str());
    return -1;
  }

  if (m_debugAlloc) {
    printf("PHOTONICS-Debug: photonicsAlloc: Allocate %lu regions among %u cores\n",
           numRegions, numCores);
    printf("PHOTONICS-Debug: photonicsAlloc: Each region has %u rows x %u cols with %lu elements\n",
           numRowsToAlloc, numCols, numElemPerRegion);
    printf("PHOTONICS-Debug: photonicsAlloc: Last region has %u rows x %u cols with %lu elements\n",
           numRowsToAlloc, numColsToAllocLast, numElemPerRegionLast);
  }

  if (numRegions > numCores) {
    if (allocType == PHOTONICS_ALLOC_V1 || allocType == PHOTONICS_ALLOC_H1) {
      printf("PHOTONICS-Error: photonicsAlloc: Allocation type %s does not allow to allocate more regions (%lu) than number of cores (%u)\n",
             photonicsUtils::photonicsAllocEnumToStr(allocType).c_str(), numRegions, numCores);
      return -1;
    }
  }

  // create new regions
  bool success = true;
  for (unsigned i = 0; i < numCores; ++i) {
    m_coreUsage.at(i)->newAllocStart();
  }
  if (allocType == PHOTONICS_ALLOC_V || allocType == PHOTONICS_ALLOC_V1 || allocType == PHOTONICS_ALLOC_H || allocType == PHOTONICS_ALLOC_H1) {
    uint64_t elemIdx = 0;
    for (uint64_t i = 0; i < numRegions; ++i) {
      PhotonicsCoreId coreId = sortedCoreId[i % numCores];
      unsigned numColsToAlloc = (i == numRegions - 1 ? numColsToAllocLast : numCols);
      unsigned numElemInRegion = (i == numRegions - 1 ? numElemPerRegionLast : numElemPerRegion);
      photonicsRegion newRegion = findAvailRegionOnCore(coreId, numRowsToAlloc, numColsToAlloc);
      if (!newRegion.isValid()) {
        printf("PHOTONICS-Error: photonicsAlloc: Failed: Out of PHOTONICS memory\n");
        success = false;
        break;
      }
      newRegion.setElemIdxBegin(elemIdx);
      elemIdx += numElemInRegion;
      newRegion.setElemIdxEnd(elemIdx); // exclusive
      newRegion.setNumColsPerElem(numColsPerElem);
      newObj.addRegion(newRegion);

      // add to core usage map
      auto alloc = std::make_pair(newRegion.getRowIdx(), numRowsToAlloc);
      m_coreUsage.at(coreId)->addRange(alloc, newObj.getObjId());
    }
  }
  for (unsigned i = 0; i < numCores; ++i) {
    m_coreUsage.at(i)->newAllocEnd(success); // rollback if failed
  }

  if (!success) {
    return -1;
  }

  PhotonicsObjId objId = -1;
  if (newObj.isValid()) {
    objId = newObj.getObjId();
    newObj.finalize();
    // update new object to resource mgr
    m_objMap.insert(std::make_pair(newObj.getObjId(), newObj));
  }

  if (m_debugAlloc) {
    if (newObj.isValid()) {
      printf("PHOTONICS-Debug: photonicsAlloc: Allocated PHOTONICS object %d successfully\n", objId);
      newObj.print();
    } else {
      printf("PHOTONICS-Debug: photonicsAlloc: Failed\n");
    }
  }
  return objId;
}

//! @brief  Allocate a new PHOTONICS matrix object
PhotonicsObjId
photonicsResMgr::photonicsAllocMat(uint64_t numElements, PhotonicsDataType dataType)
{
  if (m_debugAlloc) {
    printf("PHOTONICS-Debug: photonicsAllocMat: Request: %lu elements of type %s\n",
           numElements, photonicsUtils::photonicsDataTypeEnumToStr(dataType).c_str());
  }

  if (numElements == 0) {
    printf("PHOTONICS-Error: photonicsAllocMat: Invalid input parameter: 0 element\n");
    return -1;
  }

  PhotonicsAllocEnum allocType = PHOTONICS_ALLOC_H1;

  unsigned bitsPerElement = photonicsUtils::getNumBitsOfDataType(dataType, PhotonicsBitWidth::SIM);

  std::vector<PhotonicsCoreId> sortedCoreId = getCoreIdsSortedByLeastUsage();
  photonicsObjInfo newObj(m_availObjId, dataType, allocType, numElements, bitsPerElement, m_device);
  m_availObjId++;

  unsigned numCores = m_device->getNumCores();
  unsigned numCols = m_device->getNumCols();
  unsigned numRowsToAlloc = 1;
  uint64_t numRegions = (numElements * bitsPerElement - 1) / numCols + 1;
  unsigned numColsToAllocLast = (numElements * bitsPerElement) % numCols;
  uint64_t numElemPerRegion = numCols / bitsPerElement;
  uint64_t numElemPerRegionLast = numColsToAllocLast / bitsPerElement;
  unsigned numColsPerElem = bitsPerElement;

  if (numColsToAllocLast == 0) {
    numColsToAllocLast = numCols;
  }

  if (m_debugAlloc) {
    printf("PHOTONICS-Debug: photonicsAllocMat: Allocate %lu regions among %u cores\n",
           numRegions, numCores);
    printf("PHOTONICS-Debug: photonicsAllocMat: Each region has %u rows x %u cols with %lu elements\n",
           numRowsToAlloc, numCols, numElemPerRegion);
    printf("PHOTONICS-Debug: photonicsAllocMat: Last region has %u rows x %u cols with %lu elements\n",
           numRowsToAlloc, numColsToAllocLast, numElemPerRegionLast);
  }

  if (numRegions > numCores) {
      printf("PHOTONICS-Error: photonicsAllocMat: Allocation type %s does not allow to allocate more regions (%lu) than number of cores (%u)\n",
             photonicsUtils::photonicsAllocEnumToStr(allocType).c_str(), numRegions, numCores);
      return -1;
  }

  // create new regions
  bool success = true;
  for (unsigned i = 0; i < numCores; ++i) {
    m_coreUsage.at(i)->newAllocStart();
  }

  uint64_t elemIdx = 0;
  for (uint64_t i = 0; i < numRegions; ++i) {
    PhotonicsCoreId coreId = sortedCoreId[i % numCores];
    unsigned numColsToAlloc = (i == numRegions - 1 ? numColsToAllocLast : numCols);
    unsigned numElemInRegion = (i == numRegions - 1 ? numElemPerRegionLast : numElemPerRegion);
    photonicsRegion newRegion = findAvailRegionOnCore(coreId, numRowsToAlloc, numColsToAlloc);
    if (!newRegion.isValid()) {
      printf("PHOTONICS-Error: photonicsAllocMat: Failed: Out of PHOTONICS memory\n");
      success = false;
      break;
    }
    newRegion.setElemIdxBegin(elemIdx);
    elemIdx += numElemInRegion;
    newRegion.setElemIdxEnd(elemIdx); // exclusive
    newRegion.setNumColsPerElem(numColsPerElem);
    newObj.addRegion(newRegion);

    // add to core usage map
    auto alloc = std::make_pair(newRegion.getRowIdx(), numRowsToAlloc);
    m_coreUsage.at(coreId)->addRange(alloc, newObj.getObjId());
  }

  for (unsigned i = 0; i < numCores; ++i) {
    m_coreUsage.at(i)->newAllocEnd(success); // rollback if failed
  }

  if (!success) {
    return -1;
  }

  PhotonicsObjId objId = -1;
  if (newObj.isValid()) {
    objId = newObj.getObjId();
    newObj.finalize();
    // update new object to resource mgr
    m_objMap.insert(std::make_pair(newObj.getObjId(), newObj));
  }

  if (m_debugAlloc) {
    if (newObj.isValid()) {
      printf("PHOTONICS-Debug: photonicsAllocMat: Allocated PHOTONICS object %d successfully\n", objId);
      newObj.print();
    } else {
      printf("PHOTONICS-Debug: photonicsAllocMat: Failed\n");
    }
  }
  return objId;
}


//! @brief  Allocate a new PHOTONICS object of type buffer
PhotonicsObjId
photonicsResMgr::photonicsAllocBuffer(uint32_t numElements, PhotonicsDataType dataType)
{
  if (m_debugAlloc) {
    printf("PHOTONICS-Debug: photonicsAlloc: Request: Global Buffer %u elements of type %s\n",
           numElements, photonicsUtils::photonicsDataTypeEnumToStr(dataType).c_str());
  }

  if (numElements == 0) {
    printf("PHOTONICS-Error: photonicsAlloc: Invalid input parameter: 0 element\n");
    return -1;
  }

  unsigned bitsPerElement = photonicsUtils::getNumBitsOfDataType(dataType, PhotonicsBitWidth::SIM);

  if (numElements * bitsPerElement > m_device->getBufferSize() * 8) {
    printf("PHOTONICS-Error: photonicsAlloc: Invalid input parameter: %u elements exceeds buffer size %u bytes\n",
           numElements, m_device->getBufferSize());
    return -1;
  }

  photonicsObjInfo newObj(m_availObjId, dataType, PHOTONICS_ALLOC_H, numElements, bitsPerElement, m_device, true);
  m_availObjId++;

  unsigned numCols = m_device->getNumCols();
  unsigned numRowsToAlloc = 1;
  uint64_t numRegions = 1; // AiM buffer size will be always size of one row; For UPMEM this will be different
  unsigned numColsToAllocLast = 0;
  uint64_t numElemPerRegion = 0;
  uint64_t numElemPerRegionLast = 0;
  unsigned numColsPerElem = 0;
  numColsToAllocLast = (numElements * bitsPerElement) % numCols;
  if (numColsToAllocLast == 0) {
    numColsToAllocLast = numCols;
  }
  numElemPerRegion = numCols / bitsPerElement;
  numElemPerRegionLast = numColsToAllocLast / bitsPerElement;
  numColsPerElem = bitsPerElement;

  if (m_debugAlloc) {
    printf("PHOTONICS-Debug: photonicsAlloc: Allocate %lu regions\n", numRegions);
    printf("PHOTONICS-Debug: photonicsAlloc: Each region has %u rows x %u cols with %lu elements\n",
           numRowsToAlloc, numCols, numElemPerRegion);
    printf("PHOTONICS-Debug: photonicsAlloc: Last region has %u rows x %u cols with %lu elements\n",
           numRowsToAlloc, numColsToAllocLast, numElemPerRegionLast);
  }

  // create new regions
  bool success = true;
  uint64_t elemIdx = 0;
  unsigned numColsToAlloc = numCols;
  unsigned numElemInRegion = numElemPerRegionLast;
  photonicsRegion newRegion;
  newRegion.setCoreId(0);  // Assign global buffer to core 0; this is fine for global buffers; for UPMEM, this will be different
  newRegion.setRowIdx(0);
  newRegion.setColIdx(0);
  newRegion.setNumAllocRows(numRowsToAlloc);
  newRegion.setNumAllocCols(numColsToAlloc);
  newRegion.setIsBuffer(true);
  newRegion.setElemIdxBegin(elemIdx);
  newRegion.setIsValid(true);
  elemIdx += numElemInRegion;
  newRegion.setElemIdxEnd(elemIdx); // exclusive
  newRegion.setNumColsPerElem(numColsPerElem);
  newObj.addRegion(newRegion);

  if (!success) {
    return -1;
  }

  PhotonicsObjId objId = -1;
  if (newObj.isValid()) {
    objId = newObj.getObjId();
    newObj.finalize();
    // update new object to resource mgr
    m_objMap.insert(std::make_pair(newObj.getObjId(), newObj));
  }

  if (m_debugAlloc) {
    if (newObj.isValid()) {
      printf("PHOTONICS-Debug: photonicsAlloc: Allocated PHOTONICS object of type Buffer %d successfully\n", objId);
      newObj.print();
    } else {
      printf("PHOTONICS-Debug: photonicsAlloc: Failed\n");
    }
  }
  return objId;
}

//! @brief  Allocate a PHOTONICS object associated with an existing object
//!         Number of elements must be identical between the two associated objects
//!         For V layout, no specific requirement on data type
//!         For H layout, data type of the new object must be equal to or narrower than the associated object
PhotonicsObjId
photonicsResMgr::photonicsAllocAssociated(PhotonicsObjId assocId, PhotonicsDataType dataType)
{
  if (m_debugAlloc) {
    printf("PHOTONICS-Debug: photonicsAllocAssociated: Request: Data type %s associated with PHOTONICS object ID %d\n",
           photonicsUtils::photonicsDataTypeEnumToStr(dataType).c_str(), assocId);
  }

  // check if assoc obj is valid
  if (m_objMap.find(assocId) == m_objMap.end()) {
    printf("PHOTONICS-Error: photonicsAllocAssociated: Invalid associated PHOTONICS object ID %d\n", assocId);
    return -1;
  }

  // associated object must not be a buffer
  const photonicsObjInfo& assocObj = m_objMap.at(assocId);
  if (assocObj.isBuffer()) {
    printf("PHOTONICS-Error: photonicsAllocAssociated: Associated PHOTONICS object ID %d is a buffer, which is not allowed.\n", assocId);
    return -1;
  }

  // get regions of the assoc obj
  unsigned numCores = m_device->getNumCores();

  // check if the request can be associated with ref
  PhotonicsAllocEnum allocType = assocObj.getAllocType();
  uint64_t numElements = assocObj.getNumElements();
  unsigned bitsPerElement = photonicsUtils::getNumBitsOfDataType(dataType, PhotonicsBitWidth::SIM);
  unsigned bitsPerElementAssoc = assocObj.getBitsPerElement(PhotonicsBitWidth::PADDED);
  if (allocType == PHOTONICS_ALLOC_V || allocType == PHOTONICS_ALLOC_V1) {
    if (m_debugAlloc) {
      printf("PHOTONICS-Debug: photonicsAllocAssociated: New object of data type %s (%u bits) is associated with object (%u bits) in V layout\n",
             photonicsUtils::photonicsDataTypeEnumToStr(dataType).c_str(), bitsPerElement, bitsPerElementAssoc);
    }
  } else if (allocType == PHOTONICS_ALLOC_H || allocType == PHOTONICS_ALLOC_H1) {
    if ((bitsPerElement > bitsPerElementAssoc) && (m_device->getSimTarget() != PHOTONICS_DEVICE_BANK_LEVEL && m_device->getSimTarget() != PHOTONICS_DEVICE_FULCRUM)) {
      printf("PHOTONICS-Error: photonicsAllocAssociated: New object data type %s (%u bits) is wider than associated object (%u bits), which is not supported in H layout\n",
            photonicsUtils::photonicsDataTypeEnumToStr(dataType).c_str(), bitsPerElement, bitsPerElementAssoc);
      return -1;
    } else if (bitsPerElement < bitsPerElementAssoc) {
      if (m_debugAlloc) {
        printf("PHOTONICS-Debug: photonicsAllocAssociated: New object of data type %s (%u bits) is padded to associated object (%u bits) in H layout\n",
                photonicsUtils::photonicsDataTypeEnumToStr(dataType).c_str(), bitsPerElement, bitsPerElementAssoc);
      }
      bitsPerElement = bitsPerElementAssoc;  // padding
    } else {
      // same bit width, no padding needed
      if (m_debugAlloc) {
        printf("PHOTONICS-Debug: photonicsAllocAssociated: New object of data type %s (%u bits) is associated with object (%u bits) in H layout\n",
                photonicsUtils::photonicsDataTypeEnumToStr(dataType).c_str(), bitsPerElement, bitsPerElementAssoc);
      }
    }
  } else {
    printf("PHOTONICS-Error: photonicsAllocAssociated: Unsupported allocation type %s\n",
           photonicsUtils::photonicsAllocEnumToStr(allocType).c_str());
    return -1;
  }

  // allocate associated regions
  photonicsObjInfo newObj(m_availObjId, dataType, allocType, numElements, bitsPerElement, m_device);
  m_availObjId++;

  unsigned numCols = m_device->getNumCols();
  uint64_t numRegions = 0;
  unsigned numColsToAllocLast = 0;
  uint64_t numElemPerRegion = 0;
  uint64_t numElemPerRegionLast = 0;
  unsigned numColsPerElem = 0;

  // The reason other horizontal bit-parallel (AiM, Aquabolt) PHOTONICS is not included in this condition is that
  // they support only 16-bit floats/ints.
  // If more bit-parallel PHOTONICSs are added, this condition should be extended.
  if ((allocType == PHOTONICS_ALLOC_H || allocType == PHOTONICS_ALLOC_H1) && (bitsPerElement > bitsPerElementAssoc) && (m_device->getSimTarget() == PHOTONICS_DEVICE_BANK_LEVEL || m_device->getSimTarget() == PHOTONICS_DEVICE_FULCRUM)) {
    // allocate one region per core, with horizontal layout
    numRegions = (numElements * bitsPerElement - 1) / numCols + 1;

    // This is a controversial design decision. I am not fully sold on this
    // TODO: discuss with professor before implementing the `non-controversial` design 
    if (numRegions > assocObj.getRegions().size()) {
      printf("PHOTONICS-Error: photonicsAllocAssociated: Allocation type %s does not allow to allocate more regions (%lu) than associated object (%lu)\n",
              photonicsUtils::photonicsAllocEnumToStr(allocType).c_str(), numRegions, assocObj.getRegions().size());
      return -1;
    }

    if (numRegions > numCores) {
      printf("PHOTONICS-Error: photonicsAllocAssociated: Allocation type %s does not allow to allocate more regions (%lu) than number of cores (%u)\n",
              photonicsUtils::photonicsAllocEnumToStr(allocType).c_str(), numRegions, numCores);
      return -1;
    }

    numColsToAllocLast = (numElements * bitsPerElement) % numCols;
    if (numColsToAllocLast == 0) {
      numColsToAllocLast = numCols;
    }
    numElemPerRegion = numCols / bitsPerElement;
    numElemPerRegionLast = numColsToAllocLast / bitsPerElement;
    numColsPerElem = bitsPerElement;  
  } 

  bool success = true;
  for (unsigned i = 0; i < numCores; ++i) {
    m_coreUsage.at(i)->newAllocStart();
  }

  unsigned regionIdx = 0;
  uint64_t elemIdx = 0;
  for (const photonicsRegion& region : assocObj.getRegions()) {
    if ((bitsPerElement > bitsPerElementAssoc) && (allocType == PHOTONICS_ALLOC_H || allocType == PHOTONICS_ALLOC_H1) && (m_device->getSimTarget() == PHOTONICS_DEVICE_BANK_LEVEL || m_device->getSimTarget() == PHOTONICS_DEVICE_FULCRUM)) {
      PhotonicsCoreId coreId = region.getCoreId();
      unsigned numAllocRows = region.getNumAllocRows() * bitsPerElement / bitsPerElementAssoc;
      unsigned numAllocCols = (regionIdx == numRegions - 1 ? numColsToAllocLast : numCols);
      photonicsRegion newRegion = findAvailRegionOnCore(coreId, numAllocRows, numAllocCols);
      if (!newRegion.isValid()) {
        printf("PHOTONICS-Error: photonicsAlloc: Failed: Out of PHOTONICS memory\n");
        success = false;
        break;
      }
      newRegion.setElemIdxBegin(elemIdx);
      elemIdx += (regionIdx == numRegions - 1 ? numElemPerRegionLast : numElemPerRegion);
      if (elemIdx != region.getElemIdxEnd()) {
        printf("PHOTONICS-Error: photonicsAllocAssociated: Mismatch in element index range: %lu vs %lu\n",
               elemIdx, region.getElemIdxEnd());
        success = false;
        break;
      }
      newRegion.setElemIdxEnd(region.getElemIdxEnd()); // exclusive
      newRegion.setNumColsPerElem(numColsPerElem);
      newObj.addRegion(newRegion);

      // add to core usage map
      auto alloc = std::make_pair(newRegion.getRowIdx(), numAllocRows);
      m_coreUsage.at(coreId)->addRange(alloc, newObj.getObjId());
    } else {
      PhotonicsCoreId coreId = region.getCoreId();
      unsigned numAllocRows = region.getNumAllocRows();
      unsigned numAllocCols = region.getNumAllocCols();
      if (allocType == PHOTONICS_ALLOC_V || allocType == PHOTONICS_ALLOC_V1) {
        numAllocRows = bitsPerElement;
      }
      photonicsRegion newRegion = findAvailRegionOnCore(coreId, numAllocRows, numAllocCols);
      if (!newRegion.isValid()) {
        printf("PHOTONICS-Error: photonicsAllocAssociated: Failed: Out of PHOTONICS memory\n");
        success = false;
        break;
      }
      newRegion.setElemIdxBegin(region.getElemIdxBegin());
      newRegion.setElemIdxEnd(region.getElemIdxEnd()); // exclusive
      newRegion.setNumColsPerElem(region.getNumColsPerElem());
      newObj.addRegion(newRegion);

      // add to core usage map
      auto alloc = std::make_pair(newRegion.getRowIdx(), numAllocRows);
      m_coreUsage.at(coreId)->addRange(alloc, newObj.getObjId());
    }
    regionIdx++;
  }
  for (unsigned i = 0; i < numCores; ++i) {
    m_coreUsage.at(i)->newAllocEnd(success); // rollback if failed
  }

  if (!success) {
    return -1;
  }

  PhotonicsObjId objId = -1;
  if (newObj.isValid()) {
    objId = newObj.getObjId();
    newObj.finalize();
    newObj.setAssocObjId(assocObj.getAssocObjId());
    // update new object to resource mgr
    m_objMap.insert(std::make_pair(newObj.getObjId(), newObj));
  }

  if (m_debugAlloc) {
    if (newObj.isValid()) {
      printf("PHOTONICS-Debug: photonicsAllocAssociated: Allocated PHOTONICS object %d successfully\n", objId);
      newObj.print();
    } else {
      printf("PHOTONICS-Debug: photonicsAllocAssociated: Failed\n");
    }
  }
  return objId;
}


//! @brief  Allocate a PHOTONICS object associated with an existing object
PhotonicsObjId
photonicsResMgr::photonicsAllocAssociatedSrcVec(PhotonicsObjId assocId, PhotonicsDataType dataType)
{
  if (m_debugAlloc) {
    printf("PHOTONICS-Debug: photonicsAllocAssociatedSrcVec: Request: Data type %s associated with PHOTONICS object ID %d\n",
           photonicsUtils::photonicsDataTypeEnumToStr(dataType).c_str(), assocId);
  }

  // check if assoc obj is valid
  if (m_objMap.find(assocId) == m_objMap.end()) {
    printf("PHOTONICS-Error: photonicsAllocAssociatedSrcVec: Invalid associated PHOTONICS object ID %d\n", assocId);
    return -1;
  }

  // associated object must not be a buffer
  const photonicsObjInfo& assocObj = m_objMap.at(assocId);
  if (assocObj.isBuffer()) {
    printf("PHOTONICS-Error: photonicsAllocAssociatedSrcVec: Associated PHOTONICS object ID %d is a buffer, which is not allowed.\n", assocId);
    return -1;
  }

  // get regions of the assoc obj
  unsigned numCores = m_device->getNumCores();
  unsigned numVCores = m_device->getNumRanks();
  unsigned numHCores = m_device->getNumBankPerRank();
  unsigned numVectorsPerCore = m_device->getNumSubarrayPerBank();

  // check if the request can be associated with ref
  PhotonicsAllocEnum allocType = assocObj.getAllocType();
  uint64_t numElements = assocObj.getNumElements()/(numVCores * numVectorsPerCore);
  unsigned bitsPerElement = photonicsUtils::getNumBitsOfDataType(dataType, PhotonicsBitWidth::SIM);
  unsigned bitsPerElementAssoc = assocObj.getBitsPerElement(PhotonicsBitWidth::PADDED);
  if ((bitsPerElement > bitsPerElementAssoc) && (m_device->getSimTarget() != PHOTONICS_DEVICE_BANK_LEVEL && m_device->getSimTarget() != PHOTONICS_DEVICE_FULCRUM)) {
    printf("PHOTONICS-Error: photonicsAllocAssociatedSrcVec: New object data type %s (%u bits) is wider than associated object (%u bits), which is not supported in H layout\n",
          photonicsUtils::photonicsDataTypeEnumToStr(dataType).c_str(), bitsPerElement, bitsPerElementAssoc);
    return -1;
  } else if (bitsPerElement < bitsPerElementAssoc) {
    if (m_debugAlloc) {
      printf("PHOTONICS-Debug: photonicsAllocAssociatedSrcVec: New object of data type %s (%u bits) is padded to associated object (%u bits) in H layout\n",
              photonicsUtils::photonicsDataTypeEnumToStr(dataType).c_str(), bitsPerElement, bitsPerElementAssoc);
    }
    bitsPerElement = bitsPerElementAssoc;  // padding
  } else {
    // same bit width, no padding needed
    if (m_debugAlloc) {
      printf("PHOTONICS-Debug: photonicsAllocAssociatedSrcVec: New object of data type %s (%u bits) is associated with object (%u bits) in H layout\n",
              photonicsUtils::photonicsDataTypeEnumToStr(dataType).c_str(), bitsPerElement, bitsPerElementAssoc);
    }
  }

  // allocate associated regions
  photonicsObjInfo newObj(m_availObjId, dataType, allocType, numElements, bitsPerElement, m_device);
  m_availObjId++;

  unsigned numCols = m_device->getNumCols();
  uint64_t numRegions = 0;
  unsigned numColsToAllocLast = 0;
  uint64_t numElemPerRegion = 0;
  uint64_t numElemPerRegionLast = 0;
  unsigned numColsPerElem = 0;


  // The reason other horizontal bit-parallel (AiM, Aquabolt) PHOTONICS is not included in this condition is that
  // they support only 16-bit floats/ints.
  // If more bit-parallel PHOTONICSs are added, this condition should be extended.
  if ((allocType == PHOTONICS_ALLOC_H || allocType == PHOTONICS_ALLOC_H1) && (bitsPerElement > bitsPerElementAssoc) && (m_device->getSimTarget() == PHOTONICS_DEVICE_BANK_LEVEL || m_device->getSimTarget() == PHOTONICS_DEVICE_FULCRUM)) {
    // allocate one region per core, with horizontal layout
    numRegions = (numElements * bitsPerElement - 1) / numCols + 1;

    // This is a controversial design decision. I am not fully sold on this
    // TODO: discuss with professor before implementing the `non-controversial` design 
    if (numRegions > assocObj.getRegions().size()) {
      printf("PHOTONICS-Error: photonicsAllocAssociatedSrcVec: Allocation type %s does not allow to allocate more regions (%lu) than associated object (%lu)\n",
              photonicsUtils::photonicsAllocEnumToStr(allocType).c_str(), numRegions, assocObj.getRegions().size());
      return -1;
    }

    if (numRegions > numCores) {
      printf("PHOTONICS-Error: photonicsAllocAssociatedSrcVec: Allocation type %s does not allow to allocate more regions (%lu) than number of cores (%u)\n",
              photonicsUtils::photonicsAllocEnumToStr(allocType).c_str(), numRegions, numCores);
      return -1;
    }

    numColsToAllocLast = (numElements * bitsPerElement) % numCols;
    if (numColsToAllocLast == 0) {
      numColsToAllocLast = numCols;
    }
    numElemPerRegion = numCols / bitsPerElement;
    numElemPerRegionLast = numColsToAllocLast / bitsPerElement;
    numColsPerElem = bitsPerElement;   
  }

  bool success = true;
  for (unsigned i = 0; i < numHCores; ++i) {
    m_coreUsage.at(i)->newAllocStart();
  }

  unsigned regionIdx = 0;
  uint64_t elemIdx = 0;
  for (unsigned i = 0; i < numHCores; ++i) {
    const photonicsRegion& region = assocObj.getRegions()[i];
    if ((bitsPerElement > bitsPerElementAssoc) && (allocType == PHOTONICS_ALLOC_H || allocType == PHOTONICS_ALLOC_H1) && (m_device->getSimTarget() == PHOTONICS_DEVICE_BANK_LEVEL || m_device->getSimTarget() == PHOTONICS_DEVICE_FULCRUM)) {
      PhotonicsCoreId coreId = region.getCoreId();
      unsigned numAllocRows = region.getNumAllocRows() * bitsPerElement / bitsPerElementAssoc;
      unsigned numAllocCols = (regionIdx == numRegions - 1 ? numColsToAllocLast : numCols);
      photonicsRegion newRegion = findAvailRegionOnCore(coreId, numAllocRows, numAllocCols);
      if (!newRegion.isValid()) {
        printf("PHOTONICS-Error: photonicsAlloc: Failed: Out of PHOTONICS memory\n");
        success = false;
        break;
      }
      newRegion.setElemIdxBegin(elemIdx);
      elemIdx += (regionIdx == numRegions - 1 ? numElemPerRegionLast : numElemPerRegion);
      if (elemIdx != region.getElemIdxEnd()) {
        printf("PHOTONICS-Error: photonicsAllocAssociatedSrcVec: Mismatch in element index range: %lu vs %lu\n",
               elemIdx, region.getElemIdxEnd());
        success = false;
        break;
      }
      newRegion.setElemIdxEnd(region.getElemIdxEnd()); // exclusive
      newRegion.setNumColsPerElem(numColsPerElem);
      newObj.addRegion(newRegion);

      // add to core usage map
      auto alloc = std::make_pair(newRegion.getRowIdx(), numAllocRows);
      m_coreUsage.at(coreId)->addRange(alloc, newObj.getObjId());
    } else {
      PhotonicsCoreId coreId = region.getCoreId();
      unsigned numAllocRows = region.getNumAllocRows();
      unsigned numAllocCols = region.getNumAllocCols();
      if (allocType == PHOTONICS_ALLOC_V || allocType == PHOTONICS_ALLOC_V1) {
        numAllocRows = bitsPerElement;
      }
      photonicsRegion newRegion = findAvailRegionOnCore(coreId, numAllocRows, numAllocCols);
      if (!newRegion.isValid()) {
        printf("PHOTONICS-Error: photonicsAllocAssociatedSrcVec: Failed: Out of PHOTONICS memory\n");
        success = false;
        break;
      }
      newRegion.setElemIdxBegin(region.getElemIdxBegin());
      newRegion.setElemIdxEnd(region.getElemIdxEnd()); // exclusive
      newRegion.setNumColsPerElem(region.getNumColsPerElem());
      newObj.addRegion(newRegion);

      // add to core usage map
      auto alloc = std::make_pair(newRegion.getRowIdx(), numAllocRows);
      m_coreUsage.at(coreId)->addRange(alloc, newObj.getObjId());
    }
    regionIdx++;
  }
  for (unsigned i = 0; i < numHCores; ++i) {
    m_coreUsage.at(i)->newAllocEnd(success); // rollback if failed
  }

  if (!success) {
    return -1;
  }

  PhotonicsObjId objId = -1;
  if (newObj.isValid()) {
    objId = newObj.getObjId();
    newObj.finalize();
    newObj.setAssocObjId(assocObj.getAssocObjId());
    // update new object to resource mgr
    m_objMap.insert(std::make_pair(newObj.getObjId(), newObj));
  }

  if (m_debugAlloc) {
    if (newObj.isValid()) {
      printf("PHOTONICS-Debug: photonicsAllocAssociatedSrcVec: Allocated PHOTONICS object %d successfully\n", objId);
      newObj.print();
    } else {
      printf("PHOTONICS-Debug: photonicsAllocAssociatedSrcVec: Failed\n");
    }
  }
  return objId;
}

//! @brief  Allocate a PHOTONICS object associated with an existing object
PhotonicsObjId
photonicsResMgr::photonicsAllocAssociatedDestVec(PhotonicsObjId assocId, PhotonicsDataType dataType)
{
  if (m_debugAlloc) {
    printf("PHOTONICS-Debug: photonicsAllocAssociatedDestVec: Request: Data type %s associated with PHOTONICS object ID %d\n",
           photonicsUtils::photonicsDataTypeEnumToStr(dataType).c_str(), assocId);
  }

  // check if assoc obj is valid
  if (m_objMap.find(assocId) == m_objMap.end()) {
    printf("PHOTONICS-Error: photonicsAllocAssociatedDestVec: Invalid associated PHOTONICS object ID %d\n", assocId);
    return -1;
  }

  // associated object must not be a buffer
  const photonicsObjInfo& assocObj = m_objMap.at(assocId);
  if (assocObj.isBuffer()) {
    printf("PHOTONICS-Error: photonicsAllocAssociatedDestVec: Associated PHOTONICS object ID %d is a buffer, which is not allowed.\n", assocId);
    return -1;
  }

  // check if the request can be associated with ref
  PhotonicsAllocEnum allocType = assocObj.getAllocType();
  unsigned bitsPerElement = photonicsUtils::getNumBitsOfDataType(dataType, PhotonicsBitWidth::SIM);
  unsigned bitsPerElementAssoc = assocObj.getBitsPerElement(PhotonicsBitWidth::PADDED);

  // get regions of the assoc obj
  unsigned numCores = m_device->getNumCores();
  unsigned numVCores = m_device->getNumRanks();
  unsigned numHCores = m_device->getNumBankPerRank();
  unsigned numVectorsPerCore = m_device->getNumSubarrayPerBank();
  unsigned numElementsPerVector = m_device->getNumCols()/bitsPerElement;

  uint64_t numElements = assocObj.getNumElements()/numElementsPerVector;
  
  if ((bitsPerElement > bitsPerElementAssoc) && (m_device->getSimTarget() != PHOTONICS_DEVICE_BANK_LEVEL && m_device->getSimTarget() != PHOTONICS_DEVICE_FULCRUM)) {
    printf("PHOTONICS-Error: photonicsAllocAssociatedDestVec: New object data type %s (%u bits) is wider than associated object (%u bits), which is not supported in H layout\n",
          photonicsUtils::photonicsDataTypeEnumToStr(dataType).c_str(), bitsPerElement, bitsPerElementAssoc);
    return -1;
  } else if (bitsPerElement < bitsPerElementAssoc) {
    if (m_debugAlloc) {
      printf("PHOTONICS-Debug: photonicsAllocAssociatedDestVec: New object of data type %s (%u bits) is padded to associated object (%u bits) in H layout\n",
              photonicsUtils::photonicsDataTypeEnumToStr(dataType).c_str(), bitsPerElement, bitsPerElementAssoc);
    }
    bitsPerElement = bitsPerElementAssoc;  // padding
  } else {
    // same bit width, no padding needed
    if (m_debugAlloc) {
      printf("PHOTONICS-Debug: photonicsAllocAssociatedDestVec: New object of data type %s (%u bits) is associated with object (%u bits) in H layout\n",
              photonicsUtils::photonicsDataTypeEnumToStr(dataType).c_str(), bitsPerElement, bitsPerElementAssoc);
    }
  }

  // allocate associated regions
  photonicsObjInfo newObj(m_availObjId, dataType, allocType, numElements, bitsPerElement, m_device);
  m_availObjId++;

  unsigned numCols = bitsPerElement;
  uint64_t numRegions = 0;
  unsigned numColsToAllocLast = 0;
  uint64_t numElemPerRegion = 0;
  uint64_t numElemPerRegionLast = 0;
  unsigned numColsPerElem = 0;


  // The reason other horizontal bit-parallel (AiM, Aquabolt) PHOTONICS is not included in this condition is that
  // they support only 16-bit floats/ints.
  // If more bit-parallel PHOTONICSs are added, this condition should be extended.
  if ((allocType == PHOTONICS_ALLOC_H || allocType == PHOTONICS_ALLOC_H1) && (bitsPerElement >= bitsPerElementAssoc) && (m_device->getSimTarget() == PHOTONICS_DEVICE_BANK_LEVEL || m_device->getSimTarget() == PHOTONICS_DEVICE_FULCRUM)) {
    // allocate one region per core, with horizontal layout
    numRegions = (numElements * bitsPerElement - 1) / numCols + 1;

    // This is a controversial design decision. I am not fully sold on this
    // TODO: discuss with professor before implementing the `non-controversial` design 
    if (numRegions > assocObj.getRegions().size()) {
      printf("PHOTONICS-Error: photonicsAllocAssociatedDestVec: Allocation type %s does not allow to allocate more regions (%lu) than associated object (%lu)\n",
              photonicsUtils::photonicsAllocEnumToStr(allocType).c_str(), numRegions, assocObj.getRegions().size());
      return -1;
    }

    if (numRegions > numCores) {
      printf("PHOTONICS-Error: photonicsAllocAssociatedDestVec: Allocation type %s does not allow to allocate more regions (%lu) than number of cores (%u)\n",
              photonicsUtils::photonicsAllocEnumToStr(allocType).c_str(), numRegions, numCores);
      return -1;
    }

    numColsToAllocLast = (numElements * bitsPerElement) % numCols;
    if (numColsToAllocLast == 0) {
      numColsToAllocLast = numCols;
    }
    numElemPerRegion = numCols / bitsPerElement;
    numElemPerRegionLast = numColsToAllocLast / bitsPerElement;
    numColsPerElem = bitsPerElement;   
  }

  bool success = true;
  for (unsigned i = 0; i < numCores; ++i) {
    m_coreUsage.at(i)->newAllocStart();
  }

  unsigned regionIdx = 0;
  uint64_t elemIdx = 0;
  for (const photonicsRegion& region : assocObj.getRegions()) {
    if ((bitsPerElement > bitsPerElementAssoc) && (allocType == PHOTONICS_ALLOC_H || allocType == PHOTONICS_ALLOC_H1) && (m_device->getSimTarget() == PHOTONICS_DEVICE_BANK_LEVEL || m_device->getSimTarget() == PHOTONICS_DEVICE_FULCRUM)) {
      PhotonicsCoreId coreId = region.getCoreId();
      unsigned numAllocRows = region.getNumAllocRows() * bitsPerElement / bitsPerElementAssoc;
      unsigned numAllocCols = (regionIdx == numRegions - 1 ? numColsToAllocLast : numCols);
      photonicsRegion newRegion = findAvailRegionOnCore(coreId, numAllocRows, numAllocCols);
      if (!newRegion.isValid()) {
        printf("PHOTONICS-Error: photonicsAlloc: Failed: Out of PHOTONICS memory\n");
        success = false;
        break;
      }
      newRegion.setElemIdxBegin(elemIdx);
      elemIdx += (regionIdx == numRegions - 1 ? numElemPerRegionLast : numElemPerRegion);
      if (elemIdx != region.getElemIdxEnd()) {
        printf("PHOTONICS-Error: photonicsAllocAssociatedDestVec: Mismatch in element index range: %lu vs %lu\n",
               elemIdx, region.getElemIdxEnd());
        success = false;
        break;
      }
      newRegion.setElemIdxEnd(region.getElemIdxEnd()); // exclusive
      newRegion.setNumColsPerElem(numColsPerElem);
      newObj.addRegion(newRegion);

      // add to core usage map
      auto alloc = std::make_pair(newRegion.getRowIdx(), numAllocRows);
      m_coreUsage.at(coreId)->addRange(alloc, newObj.getObjId());
    } else {
      PhotonicsCoreId coreId = region.getCoreId();
      unsigned numAllocRows = region.getNumAllocRows();
      unsigned numAllocCols = (regionIdx == numRegions - 1 ? numColsToAllocLast : numCols);
      if (allocType == PHOTONICS_ALLOC_V || allocType == PHOTONICS_ALLOC_V1) {
        numAllocRows = bitsPerElement;
      }
      photonicsRegion newRegion = findAvailRegionOnCore(coreId, numAllocRows, numAllocCols);
      if (!newRegion.isValid()) {
        printf("PHOTONICS-Error: photonicsAllocAssociatedDestVec: Failed: Out of PHOTONICS memory\n");
        success = false;
        break;
      }
      newRegion.setElemIdxBegin(elemIdx);
      elemIdx += (regionIdx == numRegions - 1 ? numElemPerRegionLast : numElemPerRegion);
      newRegion.setElemIdxEnd(elemIdx); // exclusive
      newRegion.setNumColsPerElem(region.getNumColsPerElem());
      newObj.addRegion(newRegion);

      // add to core usage map
      auto alloc = std::make_pair(newRegion.getRowIdx(), numAllocRows);
      m_coreUsage.at(coreId)->addRange(alloc, newObj.getObjId());
    }
    regionIdx++;
  }
  for (unsigned i = 0; i < numCores; ++i) {
    m_coreUsage.at(i)->newAllocEnd(success); // rollback if failed
  }

  if (!success) {
    return -1;
  }

  PhotonicsObjId objId = -1;
  if (newObj.isValid()) {
    objId = newObj.getObjId();
    newObj.finalize();
    newObj.setAssocObjId(assocObj.getAssocObjId());
    // update new object to resource mgr
    m_objMap.insert(std::make_pair(newObj.getObjId(), newObj));
  }

  if (m_debugAlloc) {
    if (newObj.isValid()) {
      printf("PHOTONICS-Debug: photonicsAllocAssociatedDestVec: Allocated PHOTONICS object %d successfully\n", objId);
      newObj.print();
    } else {
      printf("PHOTONICS-Debug: photonicsAllocAssociatedDestVec: Failed\n");
    }
  }
  return objId;
}

//! @brief  Free a PHOTONICS object
bool
photonicsResMgr::photonicsFree(PhotonicsObjId objId, unsigned numCores)
{
  if (m_objMap.find(objId) == m_objMap.end()) {
    printf("PHOTONICS-Error: photonicsFree: Invalid PHOTONICS object ID %d\n", objId);
    return false;
  }
  const photonicsObjInfo& obj = m_objMap.at(objId);

  if (!obj.isDualContactRef()) {
    for (unsigned i = 0; i < numCores; ++i) {
      m_coreUsage.at(i)->deleteObj(objId);
    }
  }
  m_objMap.erase(objId);

  // free all reference as well
  if (m_refMap.find(objId) != m_refMap.end()) {
    for (auto refId : m_refMap.at(objId)) {
      m_objMap.erase(refId);
    }
  }

  if (m_debugAlloc) {
    printf("PHOTONICS-Debug: photonicsFree: Deleted object %d\n", objId);
  }
  return true;
}

//! @brief  Create an obj referencing to a range of an existing obj
PhotonicsObjId
photonicsResMgr::photonicsCreateRangedRef(PhotonicsObjId refId, uint64_t idxBegin, uint64_t idxEnd)
{
  assert(0); // todo
  return -1;
}

//! @brief  Create an obj referencing to negation of an existing obj based on dual-contact memory cells
PhotonicsObjId
photonicsResMgr::photonicsCreateDualContactRef(PhotonicsObjId refId)
{
  // check if ref obj is valid
  if (m_objMap.find(refId) == m_objMap.end()) {
    std::printf("PHOTONICS-Error: Invalid ref object ID %d for PHOTONICS dual contact ref\n", refId);
    return -1;
  }

  const photonicsObjInfo& refObj = m_objMap.at(refId);
  if (refObj.isDualContactRef()) {
    std::printf("PHOTONICS-Error: Cannot create dual contact ref of dual contact ref %d\n", refId);
    return -1;
  }

  // The dual-contact ref has exactly same regions as the ref object.
  // The refObjId field points to the ref object.
  // The isDualContactRef field indicates that values need to be negated during read/write.
  photonicsObjInfo newObj = refObj;
  PhotonicsObjId objId = m_availObjId++;
  newObj.setObjId(objId);
  newObj.setRefObjId(refObj.getObjId());
  m_refMap[refObj.getObjId()].insert(objId);
  newObj.setIsDualContactRef(true);
  m_objMap.insert(std::make_pair(newObj.getObjId(), newObj));

  return objId;
}

//! @brief  Alloc resource on a specific core. Perform row allocation for now.
photonicsRegion
photonicsResMgr::findAvailRegionOnCore(PhotonicsCoreId coreId, unsigned numAllocRows, unsigned numAllocCols) const
{
  photonicsRegion region;
  region.setCoreId(coreId);
  region.setColIdx(0);
  region.setNumAllocRows(numAllocRows);
  region.setNumAllocCols(numAllocCols);

  // try to find an available slot
  unsigned prevAvail = m_coreUsage.at(coreId)->findAvailRange(numAllocRows);
  if (m_device->getNumRows() - prevAvail >= numAllocRows) {
    region.setRowIdx(prevAvail);
    region.setIsValid(true);
    return region;
  }

  return region;
}

//! @brief  Get a list of core IDs sorted by least usage
std::vector<PhotonicsCoreId>
photonicsResMgr::getCoreIdsSortedByLeastUsage() const
{
  std::vector<std::pair<unsigned, unsigned>> usages;
  for (unsigned coreId = 0; coreId < m_device->getNumCores(); ++coreId) {
    unsigned usage = m_coreUsage.at(coreId)->getTotRowsInUse();
    usages.emplace_back(usage, coreId);
  }
  std::sort(usages.begin(), usages.end());
  std::vector<PhotonicsCoreId> result;
  for (const auto& it : usages) {
    result.push_back(it.second);
  }
  return result;
}

//! @brief  Find next available range of rows with a given size
unsigned
photonicsResMgr::coreUsage::findAvailRange(unsigned numRowsToAlloc)
{
  unsigned prevAvail = 0;
  for (const auto& it : m_rangesInUse) {
    unsigned rowIdx = it.first.first;
    unsigned numRows = it.first.second;
    if (rowIdx - prevAvail >= numRowsToAlloc) {
      return prevAvail;
    }
    prevAvail = rowIdx + numRows;
  }
  return prevAvail;
}

//! @brief  Add a new range to core usage.
//! The new range will be aggregated with previous adjacent ragne if they are from same object
//! Returned range is after aggregation
void
photonicsResMgr::coreUsage::addRange(std::pair<unsigned, unsigned> range, PhotonicsObjId objId)
{
  // aggregate with the prev range
  if (!m_rangesInUse.empty()) {
    auto it = std::prev(m_rangesInUse.end());
    unsigned lastIdx = it->first.first;
    unsigned lastSize = it->first.second;
    PhotonicsObjId lastObjId = it->second;
    if (lastIdx + lastSize == range.first && lastObjId == objId) {
      m_newAlloc.erase(it->first);
      m_rangesInUse.erase(it);
      range = std::make_pair(lastIdx, lastSize + range.second);
    }
  }
  m_rangesInUse.insert(std::make_pair(range, objId));
  m_newAlloc.insert(range);
}

//! @brief  Delete an object from core usage
void
photonicsResMgr::coreUsage::deleteObj(PhotonicsObjId objId)
{
  for (auto it = m_rangesInUse.begin(); it != m_rangesInUse.end();) {
    if (it->second == objId) {
      it = m_rangesInUse.erase(it);
    } else {
      ++it;
    }
  }
}

//! @brief  Start a new allocation. This is preparing for rollback
void
photonicsResMgr::coreUsage::newAllocStart()
{
  m_newAlloc.clear();
}

//! @brief  End a new allocation. If failed, rollback all regions
void
photonicsResMgr::coreUsage::newAllocEnd(bool success)
{
  if (!success) {
    for (const auto &range : m_newAlloc) {
      m_rangesInUse.erase(range);
    }
  }
  m_newAlloc.clear();
}

//! @brief  If a PHOTONICS object uses vertical data layout
bool
photonicsResMgr::isVLayoutObj(PhotonicsObjId objId) const
{
  const photonicsObjInfo& obj = getObjInfo(objId);
  return obj.isVLayout();
}

//! @brief  If a PHOTONICS object uses horizontal data layout
bool
photonicsResMgr::isHLayoutObj(PhotonicsObjId objId) const
{
  const photonicsObjInfo& obj = getObjInfo(objId);
  return obj.isHLayout();
}

//! @brief  If a PHOTONICS object uses hybrid data layout
bool
photonicsResMgr::isHybridLayoutObj(PhotonicsObjId objId) const
{
  return false;
}

