// File: photonicsResMgr.h
// PHOTONICSeval Simulator - PHOTONICS Resource Manager
// Copyright (c) 2024 University of Virginia
// This file is licensed under the MIT License.
// See the LICENSE file in the root of this repository for more details.

#ifndef LAVA_PHOTONICS_RES_MGR_H
#define LAVA_PHOTONICS_RES_MGR_H

#include "libphotonicseval.h"      // for PhotonicsObjId, PhotonicsDataType
#include "photonicsUtils.h"        // for getNumBitsOfDataType, signExt, photonicsDataTypeEnumToStr, castTypeToBits
#include <vector>            // for vector
#include <unordered_map>     // for unordered_map
#include <set>               // for set
#include <map>               // for map
#include <string>            // for string
#include <memory>            // for unique_ptr
#include <cassert>           // for assert

class photonicsDevice;


//! @class  photonicsRegion
//! @brief  Represent a rectangle region in a PHOTONICS core
class photonicsRegion
{
public:
  photonicsRegion() {}
  ~photonicsRegion() {}

  void setCoreId(PhotonicsCoreId coreId) { m_coreId = coreId; }
  void setRowIdx(unsigned rowIdx) { m_rowIdx = rowIdx; }
  void setColIdx(unsigned colIdx) { m_colIdx = colIdx; }
  void setNumAllocRows(unsigned numAllocRows) { m_numAllocRows = numAllocRows; }
  void setNumAllocCols(unsigned numAllocCols) { m_numAllocCols = numAllocCols; }
  void setElemIdxBegin(uint64_t idx) { m_elemIdxBegin = idx; }
  void setElemIdxEnd(uint64_t idx) { m_elemIdxEnd = idx; }
  void setIsValid(bool val) { m_isValid = val; }
  void setNumColsPerElem(unsigned val) { m_numColsPerElem = val; }
  void setIsBuffer(bool val) { m_isBuffer = val; }

  PhotonicsCoreId getCoreId() const { return m_coreId; }
  unsigned getRowIdx() const { return m_rowIdx; }
  unsigned getColIdx() const { return m_colIdx; }
  unsigned getNumAllocRows() const { return m_numAllocRows; }
  unsigned getNumAllocCols() const { return m_numAllocCols; }
  uint64_t getElemIdxBegin() const { return m_elemIdxBegin; }
  uint64_t getElemIdxEnd() const { return m_elemIdxEnd; }
  uint64_t getNumElemInRegion() const { return m_elemIdxEnd - m_elemIdxBegin; }
  unsigned getNumColsPerElem() const { return m_numColsPerElem; }
  bool isBuffer() const { return m_isBuffer; }

  std::pair<unsigned, unsigned> locateIthElemInRegion(unsigned i) const {
    assert(i < getNumElemInRegion());
    unsigned rowIdx = m_rowIdx; // only one row of elements per region
    unsigned colIdx = m_colIdx + i * m_numColsPerElem;
    return std::make_pair(rowIdx, colIdx);
  }

  bool isValid() const { return m_isValid && m_coreId >= 0 && m_numAllocRows > 0 && m_numAllocCols > 0; }

  void print(uint64_t regionId) const;

private:
  PhotonicsCoreId m_coreId = -1;
  unsigned m_rowIdx = 0;        // starting row index
  unsigned m_colIdx = 0;        // starting col index
  unsigned m_numAllocRows = 0;  // number of rows of this region
  unsigned m_numAllocCols = 0;  // number of cols of this region
  uint64_t m_elemIdxBegin = 0;  // begin element index in this region
  uint64_t m_elemIdxEnd = 0;    // end element index in this region
  unsigned m_numColsPerElem = 0;  // number of cols per element
  bool m_isValid = false;
  bool m_isBuffer = false;  // true if this region is a buffer region
};

//! @class  photonicsDataHolder
//! @brief  A container holding raw data vector of a PHOTONICS object as a byte array
//! Assumption: Caller gurantees correct range and indices
class photonicsDataHolder
{
public:
  photonicsDataHolder(PhotonicsDataType dataType, uint64_t numElements)
    : m_dataType(dataType),
      m_numElements(numElements)
  {
    unsigned numBitsOfDataType = photonicsUtils::getNumBitsOfDataType(m_dataType, PhotonicsBitWidth::HOST);
    // Note: Each data element is stored as m_bytesPerElement bytes in this data holder.
    // This aligns with the number of bytes per element in the host void* ptr for memcpy.
    m_bytesPerElement = (numBitsOfDataType + 7) / 8;  // round up, e.g. 1 byte per bool
    m_data.resize(m_numElements * m_bytesPerElement);
  }
  ~photonicsDataHolder() {}

  // return the number of bytes within a given range
  uint64_t getNumBytes(uint64_t idxBegin, uint64_t idxEnd) const {
    uint64_t numElements = (idxEnd == 0 ? m_numElements : idxEnd - idxBegin);
    return numElements * m_bytesPerElement;
  }

  // copy data of range [idxBegin, idxEnd) from host ptr into holder
  // use full range if idxEnd is default 0
  bool copyFromHost(void* src, uint64_t idxBegin = 0, uint64_t idxEnd = 0) {
    uint64_t byteIndex = idxBegin * m_bytesPerElement;
    uint64_t numBytes = getNumBytes(idxBegin, idxEnd);
    std::memcpy(m_data.data() + byteIndex, src, numBytes);
    return true;
  }

  // copy data of range [idxBegin, idxEnd) from holder to host ptr
  // use full range if idxEnd is default 0
  bool copyToHost(void* dest, uint64_t idxBegin = 0, uint64_t idxEnd = 0) const {
    uint64_t byteIndex = idxBegin * m_bytesPerElement;
    uint64_t numBytes = getNumBytes(idxBegin, idxEnd);
    std::memcpy(dest, m_data.data() + byteIndex, numBytes);
    return true;
  }

  // copy data of range [idxBegin, idxEnd) from this holder to another holder
  // use full range if idxEnd is default 0
  bool copyToObj(photonicsDataHolder& dest, uint64_t idxBegin = 0, uint64_t idxEnd = 0) const {
    uint64_t byteIndex = idxBegin * m_bytesPerElement;
    uint64_t numBytes = getNumBytes(idxBegin, idxEnd);
    std::memcpy(dest.m_data.data() + byteIndex, m_data.data() + byteIndex, numBytes);
    return true;
  }

  // set an element at index from bit representation
  bool setElementBits(uint64_t index, uint64_t bits) {
    uint64_t byteIndex = index * m_bytesPerElement;
    std::memcpy(m_data.data() + byteIndex, &bits, m_bytesPerElement);
    return true;
  }

  // get bit representation of an element at index
  bool getElementBits(uint64_t index, uint64_t &bits) const {
    bits = 0;
    uint64_t byteIndex = index * m_bytesPerElement;
    std::memcpy(&bits, m_data.data() + byteIndex, m_bytesPerElement);
    bits = photonicsUtils::signExt(bits, m_dataType);
    return true;
  }

  // print all bytes for debugging
  void print() const {
    printf("PHOTONICS obj data holder: data-type = %s, num-elements = %lu, bytes-per-element = %u\n",
           photonicsUtils::photonicsDataTypeEnumToStr(m_dataType).c_str(), m_numElements, m_bytesPerElement);
    for (size_t i = 0; i < m_data.size(); ++i) {
      std::printf(" %02x", m_data[i]);
      if ((i + 1) % 64 == 0) { std::printf("\n"); }
    }
    std::printf("\n");
  }

private:
  std::vector<uint8_t> m_data;
  PhotonicsDataType m_dataType;
  uint64_t m_numElements;
  unsigned m_bytesPerElement;
};

//! @class  photonicsObjInfo
//! @brief  Meta data of a PHOTONICS object which includes
//!         - PHOTONICS object ID
//!         - One or more rectangle regions allocated in one or more PHOTONICS cores
//!         - Allocation type which specifies how data is stored in a region
class photonicsObjInfo
{
public:
  photonicsObjInfo(PhotonicsObjId objId, PhotonicsDataType dataType, PhotonicsAllocEnum allocType, uint64_t numElements, unsigned bitsPerElementPadded, photonicsDevice* device)
    : m_objId(objId),
      m_assocObjId(objId),
      m_dataType(dataType),
      m_allocType(allocType),
      m_data(dataType, numElements),
      m_numElements(numElements),
      m_bitsPerElementPadded(bitsPerElementPadded),
      m_device(device)
  {}
  photonicsObjInfo(PhotonicsObjId objId, PhotonicsDataType dataType, PhotonicsAllocEnum allocType, uint64_t numElements, unsigned bitsPerElementPadded, photonicsDevice* device, bool isBuffer)
    : m_objId(objId),
      m_assocObjId(objId),
      m_dataType(dataType),
      m_allocType(allocType),
      m_data(dataType, numElements),
      m_numElements(numElements),
      m_bitsPerElementPadded(bitsPerElementPadded),
      m_device(device),
      m_isBuffer(isBuffer)
  {}
  ~photonicsObjInfo() {}

  void addRegion(photonicsRegion region) { m_regions.push_back(region); }
  void setObjId(PhotonicsObjId objId) { m_objId = objId; }
  void setAssocObjId(PhotonicsObjId assocObjId) { m_assocObjId = assocObjId; }
  void setRefObjId(PhotonicsObjId refObjId) { m_refObjId = refObjId; }
  void setIsDualContactRef(bool val) { m_isDualContactRef = val; }
  void setNumColsPerElem(unsigned val) { m_numColsPerElem = val; }
  void finalize();

  PhotonicsObjId getObjId() const { return m_objId; }
  PhotonicsObjId getAssocObjId() const { return m_assocObjId; }
  PhotonicsObjId getRefObjId() const { return m_refObjId; }
  bool isDualContactRef() const { return m_isDualContactRef; }
  PhotonicsAllocEnum getAllocType() const { return m_allocType; }
  PhotonicsDataType getDataType() const { return m_dataType; }
  uint64_t getNumElements() const { return m_numElements; }
  unsigned getBitsPerElement(PhotonicsBitWidth bitWidthType) const;
  photonicsDevice* getDevice() { return m_device; }
  bool isValid() const { return m_numElements > 0 && m_bitsPerElementPadded > 0 && !m_regions.empty(); }
  bool isVLayout() const { return m_allocType == PHOTONICS_ALLOC_V || m_allocType == PHOTONICS_ALLOC_V1; }
  bool isHLayout() const { return m_allocType == PHOTONICS_ALLOC_H || m_allocType == PHOTONICS_ALLOC_H1; }
  bool isLoadBalanced() const { return m_isLoadBalanced; }
  bool isBuffer() const { return m_isBuffer; }

  const std::vector<photonicsRegion>& getRegions() const { return m_regions; }
  std::vector<photonicsRegion> getRegionsOfCore(PhotonicsCoreId coreId) const;
  unsigned getMaxNumRegionsPerCore() const { return m_maxNumRegionsPerCore; }
  unsigned getNumCoresUsed() const { return m_numCoresUsed; }
  unsigned getNumCoreAvailable() const { return m_numCoreAvailable; }
  unsigned getMaxElementsPerRegion() const { return m_maxElementsPerRegion; }
  unsigned getNumColsPerElem() const { return m_numColsPerElem; }

  void print() const;

  // Note: Below functions are wraper APIs to access PHOTONICS object data holder
  // For regular PHOTONICS objects:
  // - Support host-to-device, device-to-host, and device-to-device copying
  // - Use bit representation to set or get an element at specific element index
  // - Support ranges in [idxBegin, idxEnd). Use full range if idxEnd is 0
  // For reference PHOTONICS objects:
  // - A ref object directly access the data holder of the ref-to object
  // - Dual-contact ref negates all bits during operations
  void copyFromHost(void* src, uint64_t idxBegin = 0, uint64_t idxEnd = 0);
  void copyToHost(void* dest, uint64_t idxBegin = 0, uint64_t idxEnd = 0) const;
  void copyToObj(photonicsObjInfo& destObj, uint64_t idxBegin = 0, uint64_t idxEnd = 0) const;
  void setElementBits(uint64_t index, uint64_t bits);
  uint64_t getElementBits(uint64_t index) const;
  template <typename T> void setElement(uint64_t index, T val) {
    setElementBits(index, photonicsUtils::castTypeToBits(val));
  }

  // Note: Below two functions are for supporting mixed functional and micro-ops level simulation.
  // Functional simulation purely uses this PHOTONICS data holder for simulation speed,
  // while micro-ops level simulation uses simulated 2D memory arrays.
  // When a functional API is called during micro-ops level simulation, call below two functions
  // to sync the data between this PHOTONICS data holder and simulated memory arrays.
  void syncFromSimulatedMem();
  void syncToSimulatedMem() const;

private:
  PhotonicsObjId m_objId = -1;
  PhotonicsObjId m_assocObjId = -1;
  PhotonicsObjId m_refObjId = -1;
  PhotonicsDataType m_dataType;
  PhotonicsAllocEnum m_allocType;
  photonicsDataHolder m_data;
  uint64_t m_numElements = 0;
  unsigned m_bitsPerElementPadded = 0;
  unsigned m_numCoreAvailable = 0;
  std::vector<photonicsRegion> m_regions;  // a list of core ID and regions
  unsigned m_maxNumRegionsPerCore = 0;
  unsigned m_numCoresUsed = 0;
  unsigned m_maxElementsPerRegion = 0;
  unsigned m_numColsPerElem = 0; // number of cols per element
  bool m_isDualContactRef = false;
  photonicsDevice* m_device = nullptr; // for accessing simulated memory
  bool m_isLoadBalanced = true;
  bool m_isBuffer = false; // true if this is a global buffer
};


//! @class  photonicsResMgr
//! @brief  PHOTONICS resource manager
class photonicsResMgr
{
public:
  photonicsResMgr(photonicsDevice* device);
  ~photonicsResMgr();

  PhotonicsObjId photonicsAlloc(PhotonicsAllocEnum allocType, uint64_t numElements, PhotonicsDataType dataType);
  PhotonicsObjId photonicsAllocMat(uint64_t numElements, PhotonicsDataType dataType);
  PhotonicsObjId photonicsAllocAssociated(PhotonicsObjId assocId, PhotonicsDataType dataType);
  PhotonicsObjId photonicsAllocAssociatedSrcVec(PhotonicsObjId assocId, PhotonicsDataType dataType);
  PhotonicsObjId photonicsAllocAssociatedDestVec(PhotonicsObjId assocId, PhotonicsDataType dataType);
  PhotonicsObjId photonicsAllocBuffer(uint32_t numElements, PhotonicsDataType dataType);
  bool photonicsFree(PhotonicsObjId objId);
  PhotonicsObjId photonicsCreateRangedRef(PhotonicsObjId refId, uint64_t idxBegin, uint64_t idxEnd);
  PhotonicsObjId photonicsCreateDualContactRef(PhotonicsObjId refId);

  bool isValidObjId(PhotonicsObjId objId) const { return m_objMap.find(objId) != m_objMap.end(); }
  const photonicsObjInfo& getObjInfo(PhotonicsObjId objId) const { assert(objId != -1); return m_objMap.at(objId); }
  photonicsObjInfo& getObjInfo(PhotonicsObjId objId) { assert(objId != -1); return m_objMap.at(objId); }

  bool isVLayoutObj(PhotonicsObjId objId) const;
  bool isHLayoutObj(PhotonicsObjId objId) const;
  bool isHybridLayoutObj(PhotonicsObjId objId) const;

private:
  photonicsRegion findAvailRegionOnCore(PhotonicsCoreId coreId, unsigned numAllocRows, unsigned numAllocCols) const;
  std::vector<PhotonicsCoreId> getCoreIdsSortedByLeastUsage() const;
  
  //! @class  coreUsage
  //! @brief  Track row usage for allocation
  class coreUsage {
  public:
    coreUsage(unsigned numRowsPerCore) : m_numRowsPerCore(numRowsPerCore) {}
    ~coreUsage() {}
    unsigned getNumRowsPerCore() const { return m_numRowsPerCore; }
    unsigned getTotRowsInUse() const { return m_totRowsInUse; }
    unsigned findAvailRange(unsigned numRowsToAlloc);
    void addRange(std::pair<unsigned, unsigned> range, PhotonicsObjId objId);
    void deleteObj(PhotonicsObjId objId);
    void newAllocStart();
    void newAllocEnd(bool success);
  private:
    unsigned m_numRowsPerCore = 0;
    unsigned m_totRowsInUse = 0;
    std::map<std::pair<unsigned, unsigned>, PhotonicsObjId> m_rangesInUse;
    std::set<std::pair<unsigned, unsigned>> m_newAlloc;
  };

  photonicsDevice* m_device;
  PhotonicsObjId m_availObjId;
  std::unordered_map<PhotonicsObjId, photonicsObjInfo> m_objMap;
  std::unordered_map<PhotonicsCoreId, std::unique_ptr<photonicsResMgr::coreUsage>> m_coreUsage;
  std::unordered_map<PhotonicsObjId, std::set<PhotonicsObjId>> m_refMap;
  bool m_debugAlloc = 0;
};

#endif

