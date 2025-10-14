// Test: C++ version of lud
// Copyright (c) 2024 University of Virginia
// This file is licensed under the MIT License.
// See the LICENSE file in the root of this repository for more details.

#include "libphotonicseval.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <cstdlib>
#include <numeric>
#include <ctime>
#include <getopt.h>
#if defined(_OPENMP)
#include <omp.h>
#endif
#include "utilPhotonics.h"
#include <iomanip>
#include <chrono>
#include <cassert>
#include <Eigen/Dense>

#define VIENNACL_WITH_OPENCL
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <viennacl/matrix.hpp>
#include <viennacl/linalg/lu.hpp>
#include <viennacl/ocl/backend.hpp>
#include <viennacl/matrix_proxy.hpp>
#include <viennacl/linalg/prod.hpp>
using namespace viennacl::ocl;
using namespace viennacl::linalg;
using boost::numeric::ublas::matrix;
using boost::numeric::ublas::range;
using boost::numeric::ublas::project;

using namespace std;
using namespace Eigen;

#if USE_GPU_SOLVER
#include <cublas_v2.h>
#include <cuda_runtime.h>
#include <npp.h>

// Params ---------------------------------------------------------------------
typedef struct Params
{
  char *dramConfigFile;
  char *matAFile;
  char *vecUFile;
  bool shouldVerify;
  bool moreDebugPrints;
} Params;

void usage()
{
  fprintf(stderr,
          "\nUsage:  ./lud.out [options]"
          "\n"
          "\n    -v    should verify result with CPU"
          "\n    -i    input image file (default=generates matrix with random numbers)"
          "\n    -k    input csv file containing the kernel matrices (default=generates matrices with random numbers)"
          "\n    -c    input file containing dramsim config"
          "\n    -m    enable more debug prints (default = false)"
          "\n");
}

struct Params getInputParams(int argc, char **argv)
{
  struct Params p;
  p.dramConfigFile = nullptr;
  p.shouldVerify = false;
  p.moreDebugPrints = false;

  int opt;
  while ((opt = getopt(argc, argv, "h:c:v:m:")) >= 0)
  {
    switch (opt)
    {
    case 'h':
      usage();
      exit(0);
      break;
    case 'c':
      p.dramConfigFile = optarg;
      break; 
    case 'v':
      p.shouldVerify = (*optarg == 't') ? true : false;
      break;
    case 'm':
      p.moreDebugPrints = (*optarg == 't') ? true : false; 
      break;       
    default:
      fprintf(stderr, "\nUnrecognized option!\n");
      usage();
      exit(0);
    }
  }
  return p;
}

class LUFactorization {
public:
    LUFactorization(struct Params params)
        : params(params)
    {   
        if(initializeHardware()!=PHOTONICS_OK) {
          std::cout << "❌ Function: " << __func__ << "Abort: intializeHardware() failed" << std::endl;
          return;
        }

        std::cout << "Matrix size: " << N << " x " << N << "\n";

        A.resize(N, N);
        std::srand(static_cast<unsigned>(std::time(nullptr)));

        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j)
                A(i, j) = static_cast<float>(std::rand() % 10 + 1); // values 1–10

        // Optionally make diagonally dominant for stability:
        for (int i = 0; i < N; ++i)
            A(i, i) += static_cast<float>(50);

        std::cout << "\nGenerated matrix A (first 5×5 block):\n";
        for (int i = 0; i < std::min<int>(N, 5); ++i) {
            for (int j = 0; j < std::min<int>(N, 5); ++j)
                std::cout << A(i, j) << "\t";
            if (N > 5) std::cout << "...";
            std::cout << "\n";
        }

        tileCols = numHCores * numElementsPerVector;
        tileRows = numVCores * numVectorsPerCore;

        Akk_tile.resize(tileCols, tileCols);
        tempT.resize(tileCols, tileCols);

        matR.resize(tileSize);
        out.resize(numHCores * tileRows);
        vecU.resize(tileCols);

        // -------------------------------
        // Enumerate available platforms
        // -------------------------------
        std::vector<platform> platforms = get_platforms();
        std::cout << "Found " << platforms.size() << " OpenCL platform(s):\n";
        for (auto &p : platforms)
            std::cout << "  Platform: " << p.info() << std::endl;

        // -------------------------------
        // Pick first GPU device
        // -------------------------------
        std::vector<device> devices = platforms[0].devices(CL_DEVICE_TYPE_GPU);
        if (devices.empty()) {
            std::cerr << "No GPU devices found.\n";
            return;
        }
        std::cout << "Using device: " << devices[0].name() << "\n";

        setup_context(0, devices[0]);
        context &ctx = current_context();

        // Optional: safer kernel flags for NVIDIA OpenCL
        ctx.build_options("-cl-std=CL1.2 -cl-fast-relaxed-math");
    }

    ~LUFactorization(){
        photonicsFreeSrcVec(vecUObject);
        photonicsFreeDestVec(outObject);
        photonicsFreeMat(matAObject);
    }

    void solvePhotonics() {
        cout << "\nLU Factorization with Photonics+GPU..." << endl;

        A_pho.resize(N, N);
        A_pho_res.resize(N, N);
        A_pho_q.resize(N, N);
        A_pho_q_i.resize(N, N);
        A_pho_temp_i.resize(N, N);
        // cublasCreate(&handle);
        // cudaStreamCreate(&stream);
        // init_gpu_handles(stream);

        // -------------------------------
        // LU factorization on GPU + Photonics
        // -------------------------------
        quantize_matrix(minVal, maxVal);
        auto start = std::chrono::high_resolution_clock::now();
        project(A_pho_res, range(0, tileCols), range(0, N)) = project(A, range(0, tileCols), range(0, N));
        project(A_pho_res, range(tileCols, N), range(0, tileCols)) = project(A, range(tileCols, N), range(0, tileCols));

        // ---- Blocked LU: for k = 0.. in steps of tileCols ----
        for (int k = 0; k < N; k += tileCols)
        {
            start = std::chrono::high_resolution_clock::now();
            viennacl::range rK(k, k + tileCols);
            viennacl::range cK(k, k + tileCols);
            hostElapsedTimeGpu += (std::chrono::high_resolution_clock::now() - start);

            viennacl::range u_range(k, N);
            viennacl::matrix_range<viennacl::matrix<float>> AU_view(A_pho, rK, u_range);
            matrix<float> AU_block(tileCols, N - k);
            AU_block = project(A_pho_res, range(k, k+tileCols), range(k, N));
            viennacl::range l_range(k + tileCols, N);
            viennacl::matrix_range<viennacl::matrix<float>> AL_view(A_pho, l_range, cK);
            matrix<float> AL_block(tileCols, N - k);
            AL_block = project(A_pho_res, range(k+tileCols, N), range(k, k+tileCols));
            
            start = std::chrono::high_resolution_clock::now();
            viennacl::copy(AU_block, AU_view);
            viennacl::copy(AL_block, AL_view);
            transferElapsedTimeGpu += (std::chrono::high_resolution_clock::now() - start);

            start = std::chrono::high_resolution_clock::now();
            // Ranges for the panel (Akk), the row (Ak*), and the column (A*k)

            // --- 1) Panel LU on diagonal tile: Akk = Lkk * Ukk ---
            // We factor a temporary Akk_tile (b×b) on device, then write back.
            viennacl::matrix_range<viennacl::matrix<float>> Akk_view(A_pho, rK, cK);

            Akk_tile = Akk_view;                        // copy tile view -> dense tile
            viennacl::linalg::lu_factorize(Akk_tile);   // LU on the tile (GPU)
            Akk_view = Akk_tile;                        // write back LU factors into A

            // --- 2) Update U blocks to the right: U_kj = L_kk^{-1} A_kj  (left lower solve) ---
            for (int j = k + tileCols; j < N; j += tileCols)
            {
                viennacl::range cJ(j, j + tileCols);
                viennacl::matrix_range<viennacl::matrix<float>> Akj(A_pho, rK, cJ);

                // Solve L_kk * X = A_kj  (L is unit-lower from LU; use unit_lower_tag)
                viennacl::linalg::inplace_solve(Akk_tile, Akj, viennacl::linalg::unit_lower_tag());
                // After this: A_kj holds U_kj
            }

            // --- 3) Update L blocks below: L_ik = A_ik * U_kk^{-1} (right upper solve)
            // Implement right-side solve via transpose trick:
            //   X * U = B  <=>  U^T * X^T = B^T  (left lower solve on transposed system).
            for (int i = k + tileCols; i < N; i += tileCols)
            {
                viennacl::range rI(i, i + tileCols);
                viennacl::matrix_range<viennacl::matrix<float>> Aik(A_pho, rI, cK);

                // tempT := A_ik^T  (size bk x bi)
                tempT = trans(Aik);

                // Solve (U_kk^T) * tempT = A_ik^T.
                // U_kk is upper => U_kk^T is lower (non-unit).
                viennacl::linalg::inplace_solve(trans(Akk_tile), tempT, viennacl::linalg::lower_tag());

                // Write back: A_ik := (tempT)^T  (now holds L_ik)
                Aik = trans(tempT);
            }
            hostElapsedTimeGpu += (std::chrono::high_resolution_clock::now() - start);

            start = std::chrono::high_resolution_clock::now();
            viennacl::copy(AU_view, AU_block);
            transferElapsedTimeGpu += (std::chrono::high_resolution_clock::now() - start);
            project(A_pho_res, range(k, k+tileCols), range(k, N)) = AU_block;

            if(k >= N-tileCols) break;

            start = std::chrono::high_resolution_clock::now();
            viennacl::copy(AL_view, AL_block);
            transferElapsedTimeGpu += (std::chrono::high_resolution_clock::now() - start);
            project(A_pho_res, range(k+tileCols, N), range(k, k+tileCols)) = AL_block;

            quantize_portion(k, minVal, maxVal);

            calPhotonics(k);

            start = std::chrono::high_resolution_clock::now();
            A_pho_q_i.block(k+tileCols, k+tileCols, N-k-tileCols, N-k-tileCols) = A_pho_q_i.block(k+tileCols, k+tileCols, N-k-tileCols, N-k-tileCols) - A_pho_temp_i.block(k+tileCols, k+tileCols, N-k-tileCols, N-k-tileCols);
            hostElapsedTimeCpu += (std::chrono::high_resolution_clock::now() - start);
            A_pho_q = A_pho_q_i.cast<uint8_t>();

            dequantize_portion(k+tileCols, minVal, maxVal);
        }

        std::cout << "LU factorization completed in " << hostElapsedTimeGpu.count() << " ms on the GPU + Photonics\n";

        std::cout << "\nLU-factored matrix (first 5×5 block):\n";
        for (int i = 0; i < std::min<int>(N, 5); ++i) {
            for (int j = 0; j < std::min<int>(N, 5); ++j)
                std::cout << A_pho_res(i, j) << "\t";
            if (N > 5) std::cout << "...";
            std::cout << "\n";
        }
        // cublasDestroy(handle);
    }

    void solveGPU() {
        cout << "\nLU Factorizing with GPU..." << endl;

        A_gpu.resize(N, N);
        A_gpu_res.resize(N, N);

        // -------------------------------
        // LU factorization on GPU
        // -------------------------------
        auto start = std::chrono::high_resolution_clock::now();
        viennacl::copy(A, A_gpu);
        gpuTransferElapsedTime += (std::chrono::high_resolution_clock::now() - start);

        // ---- Blocked LU: for k = 0.. in steps of tileCols ----
        start = std::chrono::high_resolution_clock::now();
        for (int k = 0; k < N; k += tileCols)
        {

            // Ranges for the panel (Akk), the row (Ak*), and the column (A*k)
            viennacl::range rK(k, k + tileCols);
            viennacl::range cK(k, k + tileCols);

            // --- 1) Panel LU on diagonal tile: Akk = Lkk * Ukk ---
            // We factor a temporary Akk_tile (b×b) on device, then write back.
            viennacl::matrix_range<viennacl::matrix<float>> Akk_view(A_gpu, rK, cK);

            Akk_tile = Akk_view;                        // copy tile view -> dense tile
            viennacl::linalg::lu_factorize(Akk_tile);   // LU on the tile (GPU)
            Akk_view = Akk_tile;                        // write back LU factors into A

            // --- 2) Update U blocks to the right: U_kj = L_kk^{-1} A_kj  (left lower solve) ---
            for (int j = k + tileCols; j < N; j += tileCols)
            {
                viennacl::range cJ(j, j + tileCols);
                viennacl::matrix_range<viennacl::matrix<float>> Akj(A_gpu, rK, cJ);

                // Solve L_kk * X = A_kj  (L is unit-lower from LU; use unit_lower_tag)
                viennacl::linalg::inplace_solve(Akk_tile, Akj, viennacl::linalg::unit_lower_tag());
                // After this: A_kj holds U_kj
            }

            // --- 3) Update L blocks below: L_ik = A_ik * U_kk^{-1} (right upper solve)
            // Implement right-side solve via transpose trick:
            //   X * U = B  <=>  U^T * X^T = B^T  (left lower solve on transposed system).
            for (int i = k + tileCols; i < N; i += tileCols)
            {
                viennacl::range rI(i, i + tileCols);
                viennacl::matrix_range<viennacl::matrix<float>> Aik(A_gpu, rI, cK);

                // tempT := A_ik^T  (size bk x bi)
                tempT = trans(Aik);

                // Solve (U_kk^T) * tempT = A_ik^T.
                // U_kk is upper => U_kk^T is lower (non-unit).
                viennacl::linalg::inplace_solve(trans(Akk_tile), tempT, viennacl::linalg::lower_tag());

                // Write back: A_ik := (tempT)^T  (now holds L_ik)
                Aik = trans(tempT);
            }

            // --- 4) Schur update on trailing tiles: A_ij -= L_ik * U_kj  ---
            for (int i = k + tileCols; i < N; i += tileCols)
            {
                viennacl::range rI(i, i + tileCols);
                viennacl::matrix_range<viennacl::matrix<float>> Aik(A_gpu, rI, cK);  // L_ik

                for (int j = k + tileCols; j < N; j += tileCols)
                {
                    viennacl::range cJ(j, j + tileCols);
                    viennacl::matrix_range<viennacl::matrix<float>> Akj(A_gpu, rK, cJ);  // U_kj
                    viennacl::matrix_range<viennacl::matrix<float>> Aij(A_gpu, rI, cJ);  // update

                    // A_ij -= A_ik * A_kj
                    Aij = Aij - viennacl::linalg::prod(Aik, Akj);
                }
            }
        }
        gpuComputeElapsedTime = (std::chrono::high_resolution_clock::now() - start);
        std::cout << "LU factorization completed in " << gpuComputeElapsedTime.count() << " ms on the GPU\n";

        start = std::chrono::high_resolution_clock::now();
        viennacl::copy(A_gpu, A_gpu_res);
        gpuTransferElapsedTime += (std::chrono::high_resolution_clock::now() - start);


        std::cout << "\nLU-factored matrix (first 5×5 block):\n";
        for (int i = 0; i < std::min<int>(N, 5); ++i) {
            for (int j = 0; j < std::min<int>(N, 5); ++j)
                std::cout << A_gpu_res(i, j) << "\t";
            if (N > 5) std::cout << "...";
            std::cout << "\n";
        }
    }

    void printGpuPhotonicsBreakdown() { 
        cout << "\nPrinting GPU+Photonics Breakdown..." << endl;
        printCpuHostElapsedTime();
        printGpuHostElapsedTime();
        printGpuQuantElapsedTime();
        printGpuSyncElapsedTime();
        printGpuTransferElapsedTime();
    }

    void printGPUElapsedTime() { 
        cout << "GPU Compute elapsed time: " << std::fixed << std::setprecision(3) << gpuComputeElapsedTime.count() << " ms, Transfer elapsed time: " 
          << std::fixed << std::setprecision(3) << gpuTransferElapsedTime.count() << " ms."<< endl;
    }

    void compareWithPhotonics(){
        cout << "\nComparing Photonics and Software Solutions...\n";

        if (A_gpu_res.size1() != A_pho_res.size1() || A_gpu_res.size2() != A_pho_res.size2()) {
            std::cerr << "❌ Error: Matrix size mismatch: A_gpu_res(" << A_gpu_res.size1() << "x" << A_gpu_res.size2()
                      << ") vs A_pho_res(" << A_pho_res.size1() << "x" << A_pho_res.size2() << ")\n";
            return;
        }

        float max_abs_error = 0.0f;
        float total_abs_error = 0.0f;
        float normA = 0.0f;
        float norm_diff = 0.0f;
        int count_mismatch = 0;

        for (std::size_t i = 0; i < A_gpu_res.size1(); ++i) {
            for (std::size_t j = 0; j < A_gpu_res.size2(); ++j) {
                float diff = std::fabs(A_gpu_res(i, j) - A_pho_res(i, j));
                max_abs_error = std::max(max_abs_error, diff);
                total_abs_error += diff;
                normA += A_gpu_res(i, j) * A_gpu_res(i, j);
                norm_diff += diff * diff;

                if (diff > tolerance) {
                    ++count_mismatch;
                    if (count_mismatch <= 5) { // only show first few mismatches
                        std::cout << "  Mismatch at (" << i << "," << j << "): "
                                  << "A_gpu_res=" << A_gpu_res(i, j) << ", A_pho_res=" << A_pho_res(i, j)
                                  << ", diff=" << diff << "\n";
                    }
                }
            }
        }

        float mean_abs_error = total_abs_error / (A_gpu_res.size1() * A_gpu_res.size2());
        float relative_error = std::sqrt(norm_diff) / (std::sqrt(normA) + std::numeric_limits<float>::epsilon());

        std::cout << "\n=== Matrix Comparison Report ===\n";
        std::cout << "Max Absolute Error     : " << max_abs_error << "\n";
        std::cout << "Mean Absolute Error    : " << mean_abs_error << "\n";
        std::cout << "Relative Frobenius Err : " << relative_error << "\n";
        std::cout << "Values above tolerance : " << count_mismatch << "\n";
        std::cout << "================================\n";
        return;
    }

private:
    struct Params params;
    int N;          // total number of grid points
    int tileSize;
    int tileCols;
    int tileRows;
    const float tolerance = 1e-5f;
    const float minVal = -256.0f;
    const float maxVal = 256.0f;
    int num_bits = 8;
    float q_levels = (1 << num_bits) - 1;
    int numVCores, numHCores, numVectorsPerCore, numElementsPerVector;
    size_t bytes_R, bytes_u, bytes_u_q, bytes_out, bytes_out_q;

    matrix<float> A, A_gpu_res, A_pho_res;
    //matrix<uint8_t> A_pho_q;
    viennacl::matrix<float> A_gpu, A_pho;
    viennacl::matrix<float> Akk_tile, tempT;

    //MatrixXf A_pho_temp;
    //VectorXf u, u_new, u_temp, start_u, q_one_vec, u_gpu, u_temp_gpu, q_one_ext;
    ArrayXXi A_pho_q_i, A_pho_temp_i;
    Array<uint8_t, Dynamic, Dynamic> A_pho_q;
    //Array<uint8_t, Dynamic, 1> q_vec, q_vec_gpu;

    std::chrono::duration<float, std::milli> hostElapsedTimeCpu = std::chrono::duration<float, std::milli>::zero();
    std::chrono::duration<float, std::milli> hostElapsedTimeGpu = std::chrono::duration<float, std::milli>::zero();
    std::chrono::duration<float, std::milli> quantElapsedTimeGpu = std::chrono::duration<float, std::milli>::zero();
    std::chrono::duration<float, std::milli> syncElapsedTimeGpu = std::chrono::duration<float, std::milli>::zero();
    std::chrono::duration<float, std::milli> gpuComputeElapsedTime = std::chrono::duration<float, std::milli>::zero();
    std::chrono::duration<float, std::milli> gpuTransferElapsedTime = std::chrono::duration<float, std::milli>::zero();
    std::chrono::duration<float, std::milli> transferElapsedTimeGpu = std::chrono::duration<float, std::milli>::zero();

    PhotonicsObjId matAObject;
    PhotonicsObjId outObject;
    PhotonicsObjId vecUObject;
    PhotonicsDeviceProperties deviceParams;

    // cublasHandle_t handle;
    // cudaStream_t stream; 
    // NppStreamContext nppCtx;

    PhotonicsStatus status = PHOTONICS_OK;
    vector<uint8_t> vecU;
    vector<uint8_t> out;
    vector<uint8_t> matR;

    PhotonicsStatus initializeHardware(){
      if (!initPhotonicsAccel(params.dramConfigFile))
        return PHOTONICS_ERROR;

      if (!getHardwareSpecs(&deviceParams)){
        return PHOTONICS_ERROR;
      }

      numVCores = deviceParams.numRanks;
      numHCores = deviceParams.numBankPerRank;
      numVectorsPerCore = deviceParams.numSubarrayPerBank;
      numElementsPerVector = deviceParams.numColPerSubarray / 32;

      N = deviceParams.matrixSize;

      tileSize = numVCores * numHCores * numElementsPerVector * numVectorsPerCore;

      if ((numVCores > 1) || (numHCores >1))
      {
          cerr << "❌ Error: Chip synchronization not available in this benchmark! \n"
                  << "numVCores = " << numVCores << ", numHCores = " << numHCores << endl;
          return PHOTONICS_ERROR;
      }

      matAObject = photonicsAllocMat(tileSize, PHOTONICS_FP32);
      if (matAObject == -1)
      {
          std::cout << "❌ Function: " << __func__ << "Abort: photonicsAlloc failed for matrix A" << std::endl;
          return PHOTONICS_ERROR;
      }
      
      outObject = photonicsAllocAssociatedDestVec(matAObject, PHOTONICS_FP32);
      if (outObject == -1)
      {
          std::cout << "❌ Function: " << __func__ << "Abort: photonicsAllocAssociated failed for vector U = " << outObject << std::endl;
          return PHOTONICS_ERROR;
      }
      
      vecUObject = photonicsAllocAssociatedSrcVec(matAObject, PHOTONICS_FP32);
      if (vecUObject == -1)
      {
          std::cout << "❌ Function: " << __func__ << "Abort: photonicsAllocAssociated failed for vector U = " << vecUObject << std::endl;
          return PHOTONICS_ERROR;
      }

      return PHOTONICS_OK;

    }

    void calPhotonics(int k){

      for (int i = k + tileCols; i < N; i += tileRows)
      {
          for (int j = 0; j < tileRows; j++) {
              for (int l = 0; l < tileCols; l++) {
                matR[j * tileCols + l] = A_pho_q(i+j,l+k);
              }
          }
          status = photonicsCopyHostToDevice((void *)matR.data(), matAObject);
          if (status != PHOTONICS_OK)
          {
              std::cout << "❌ Function: " << __func__ << "Abort: photonicsCopyHostToDevice failed between matAObject and matA" << std::endl;
              return;
          }

          for (int j = k + tileCols; j < N; j += 1){
              for (int l = 0; l < tileCols; l++) {
                  vecU[l] = A_pho_q(l+k, j);
              }
              status = photonicsCopyHostToDevice((void *)vecU.data(), vecUObject);
              if (status != PHOTONICS_OK)
              {
                  std::cout << "❌ Function: " << __func__ << "Abort: photonicsCopyHostToDevice failed between vecUObject and vecU at i=" << i << std::endl;
                  return;
              }

              status = photonicsMvm(matAObject, vecUObject, outObject);
              if (status != PHOTONICS_OK)
              {
                  std::cout << "❌ Function: " << __func__ << "Abort: photonicsMvm Failed at i=" << i << std::endl;
                  return;
              }

              status = photonicsCopyDeviceToHost(outObject, out.data());
              if (status != PHOTONICS_OK)
              {
                  std::cout << "❌ Function: " << __func__ << "Abort: photonicsCopyDeviceToHost failed between outObject and out at i=" << i << std::endl;
                  return;
              }
              for (int l = 0; l < tileRows; l++) {
                  A_pho_temp_i(l+i, j) = static_cast<int>(out[l]);
              }
          }
      }
      return;
    }

    void quantize_matrix(float min_val, float max_val) {
        Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> 
              A_map(A.data().begin(), A.size1(), A.size2());
        auto start2 = std::chrono::high_resolution_clock::now();
        const float scale = q_levels / (max_val - min_val);
        A_pho_q.block(tileCols, tileCols, N-tileCols, N-tileCols) = 
                (((A_map.block(tileCols, tileCols, N-tileCols, N-tileCols).array() - min_val) * scale) + 0.5f)
                        .cwiseMin(static_cast<float>(q_levels))
                        .cwiseMax(0.0f)
                        .cast<uint8_t>();
        quantElapsedTimeGpu += (std::chrono::high_resolution_clock::now() - start2);
        A_pho_q_i = A_pho_q.cast<int>();
    }

    void quantize_portion(int k, float min_val, float max_val) {
        Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> 
              A_map(A_pho_res.data().begin(), A_pho_res.size1(), A_pho_res.size2());
        auto start2 = std::chrono::high_resolution_clock::now();
        const float scale = q_levels / (max_val - min_val);
        A_pho_q.block(k, k+tileCols, tileCols, N-k-tileCols) = 
                (((A_map.block(k, k+tileCols, tileCols, N-k-tileCols).array() - min_val) * scale) + 0.5f)
                        .cwiseMin(static_cast<float>(q_levels))
                        .cwiseMax(0.0f)
                        .cast<uint8_t>();
        A_pho_q.block(k+tileCols, k, N-k-tileCols, tileCols) = 
                (((A_map.block(k+tileCols, k, N-k-tileCols, tileCols).array() - min_val) * scale) + 0.5f)
                        .cwiseMin(static_cast<float>(q_levels))
                        .cwiseMax(0.0f)
                        .cast<uint8_t>();
        quantElapsedTimeGpu += (std::chrono::high_resolution_clock::now() - start2);
    }

    void dequantize_portion(int k, float min_val, float max_val) {
        Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> 
              A_map(A_pho_res.data().begin(), A_pho_res.size1(), A_pho_res.size2());
        auto start2 = std::chrono::high_resolution_clock::now();
        const float scale = (max_val - min_val) / q_levels;
        A_map.block(k, k, tileCols, N-k).array() = 
                (A_pho_q.block(k, k, tileCols, N-k).cast<float>() * scale) + min_val;
        if (k < N-tileCols){
            A_map.block(k+tileCols, k, N-k-tileCols, tileCols).array() = 
                    (A_pho_q.block(k+tileCols, k, N-k-tileCols, tileCols).cast<float>() * scale) + min_val;
        }
        quantElapsedTimeGpu += (std::chrono::high_resolution_clock::now() - start2);
    }

    // void init_gpu_handles(cudaStream_t stream) {
    //     cublasSetStream(handle, stream);

    //     // Fill NPP stream context (minimal set)
    //     memset(&nppCtx, 0, sizeof(nppCtx));
    //     nppCtx.hStream = stream;
    //     cudaDeviceProp prop{};
    //     int dev=0; cudaGetDevice(&dev); cudaGetDeviceProperties(&prop, dev);
    //     nppCtx.nCudaDeviceId                 = dev;
    //     nppCtx.nMultiProcessorCount          = prop.multiProcessorCount;
    //     nppCtx.nMaxThreadsPerMultiProcessor  = prop.maxThreadsPerMultiProcessor;
    //     nppCtx.nMaxThreadsPerBlock           = prop.maxThreadsPerBlock;
    //     nppCtx.nSharedMemPerBlock            = prop.sharedMemPerBlock;
    //     nppCtx.nCudaDevAttrComputeCapabilityMajor = prop.major;
    //     nppCtx.nCudaDevAttrComputeCapabilityMinor = prop.minor;
    // }

    void printGpuHostElapsedTime() { 
      cout << "GPU Host elapsed time: " << std::fixed << std::setprecision(3) << hostElapsedTimeGpu.count() << " ms." << endl;
    }
    
    void printCpuHostElapsedTime() { 
      cout << "CPU Host elapsed time: " << std::fixed << std::setprecision(3) << hostElapsedTimeCpu.count() << " ms." << endl;
    }

    void printGpuQuantElapsedTime() { 
      cout << "Quantization elapsed time: " << std::fixed << std::setprecision(3) << quantElapsedTimeGpu.count() << " ms." << endl;
    }

    void printGpuSyncElapsedTime() { 
      cout << "Photonics Core Synchronization elapsed time: " << std::fixed << std::setprecision(3) << syncElapsedTimeGpu.count() << " ms." << endl;
    }

    void printGpuTransferElapsedTime() { 
      cout << "GPU transfer elapsed time: " << std::fixed << std::setprecision(3) << transferElapsedTimeGpu.count() << " ms." << endl;
    }
};

int main(int argc, char *argv[])
{
  struct Params params = getInputParams(argc, argv);

  params.dramConfigFile = (char *) "./configs/photonics/lud.cfg";

  LUFactorization solver(params);
  
  solver.solveGPU();
  solver.solvePhotonics();

  solver.compareWithPhotonics();

  photonicsShowStats();
  solver.printGPUElapsedTime();

  solver.printGpuPhotonicsBreakdown();

  return 0;
}

#endif
