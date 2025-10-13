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
#include <Eigen/Dense>#

#define VIENNACL_WITH_OPENCL
#include <boost/numeric/ublas/matrix.hpp>
#include <viennacl/matrix.hpp>
#include <viennacl/linalg/lu.hpp>
#include <viennacl/ocl/backend.hpp>
#include <viennacl/matrix_proxy.hpp>
#include <viennacl/linalg/prod.hpp>
using namespace viennacl::ocl;
using namespace viennacl::linalg;
using boost::numeric::ublas::matrix;

using namespace std;
using namespace Eigen;

#if USE_GPU_SOLVER
#include <cublas_v2.h>
#include <cuda_runtime.h>
#include <npp.h>
#endif

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
            A(i, i) += static_cast<float>(N * 5);

        std::cout << "\nGenerated matrix A (first 5×5 block):\n";
        for (int i = 0; i < std::min<int>(N, 5); ++i) {
            for (int j = 0; j < std::min<int>(N, 5); ++j)
                std::cout << A(i, j) << "\t";
            if (N > 5) std::cout << "...";
            std::cout << "\n";
        }

        tileCols = numHCores * numElementsPerVector;
        tileRows = numVCores * numVectorsPerCore;

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
        viennacl::copy(A, A_pho);
        Akk_tile.resize(tileCols, tileCols);
        tempT.resize(tileCols, tileCols);

        // -------------------------------
        // LU factorization on GPU + Photonics
        // -------------------------------
        auto start = std::chrono::high_resolution_clock::now();

        // ---- Blocked LU: for k = 0.. in steps of tileCols ----
        for (int k = 0; k < N; k += tileCols)
        {

            // Ranges for the panel (Akk), the row (Ak*), and the column (A*k)
            viennacl::range rK(k, k + tileCols);
            viennacl::range cK(k, k + tileCols);

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

            // --- 4) Schur update on trailing tiles: A_ij -= L_ik * U_kj  ---
            for (int i = k + tileCols; i < N; i += tileCols)
            {
                viennacl::range rI(i, i + tileCols);
                viennacl::matrix_range<viennacl::matrix<float>> Aik(A_pho, rI, cK);  // L_ik

                for (int j = k + tileCols; j < N; j += tileCols)
                {
                    viennacl::range cJ(j, j + tileCols);
                    viennacl::matrix_range<viennacl::matrix<float>> Akj(A_pho, rK, cJ);  // U_kj
                    viennacl::matrix_range<viennacl::matrix<float>> Aij(A_pho, rI, cJ);  // update

                    // A_ij -= A_ik * A_kj
                    Aij = Aij - viennacl::linalg::prod(Aik, Akj);
                }
            }
        }

        hostElapsedTimeGpu = (std::chrono::high_resolution_clock::now() - start);
        std::cout << "LU factorization completed in " << hostElapsedTimeGpu.count() << " ms on the GPU + Photonics\n";

        // -------------------------------
        // Copy GPU → CPU and print
        // -------------------------------
        A_pho_res.resize(N, N);
        viennacl::copy(A_pho, A_pho_res);

        std::cout << "\nLU-factored matrix (first 5×5 block):\n";
        for (int i = 0; i < std::min<int>(N, 5); ++i) {
            for (int j = 0; j < std::min<int>(N, 5); ++j)
                std::cout << A_pho_res(i, j) << "\t";
            if (N > 5) std::cout << "...";
            std::cout << "\n";
        }

        //quantize_matrix(-1.0, 1.0);
        // status = photonicsCopyHostToDevice((void *)matR.data(), matAObject);
        // if (status != PHOTONICS_OK)
        // {
        //   std::cout << "❌ Function: " << __func__ << "Abort: photonicsCopyHostToDevice failed between matAObject and matA" << std::endl;
        //   return;
        // }


    }

    void solveGPU() {
        cout << "\nLU Factorizing with GPU..." << endl;

        A_gpu.resize(N, N);
        viennacl::copy(A, A_gpu);

        // -------------------------------
        // LU factorization on GPU
        // -------------------------------
        auto start = std::chrono::high_resolution_clock::now();

        viennacl::linalg::lu_factorize(A_gpu);

        gpuElapsedTime = (std::chrono::high_resolution_clock::now() - start);
        std::cout << "LU factorization completed in " << gpuElapsedTime.count() << " ms on the GPU\n";

        // -------------------------------
        // Copy GPU → CPU and print
        // -------------------------------
        A_gpu_res.resize(N, N);
        viennacl::copy(A_gpu, A_gpu_res);

        std::cout << "\nLU-factored matrix (first 5×5 block):\n";
        for (int i = 0; i < std::min<int>(N, 5); ++i) {
            for (int j = 0; j < std::min<int>(N, 5); ++j)
                std::cout << A_gpu_res(i, j) << "\t";
            if (N > 5) std::cout << "...";
            std::cout << "\n";
        }
    }

    const matrix<float>& getPhotonicsSolution() const { return A_pho_res; }

    const matrix<float>& getCPUSolutions() const { return A_cpu_res; }

    const matrix<float>& getGPUSolutions() const { return A_gpu_res; }

    void printGpuPhotonicsBreakdown() { 
        #if USE_GPU_SOLVER
            cout << "\nPrinting GPU+Photonics Breakdown..." << endl;
            printGpuHostElapsedTime();
            printGpuQuantElapsedTime();
            printGpuSyncElapsedTime();
            printGpuTransferElapsedTime();
        #endif
    }

    void printGPUElapsedTime() { 
        #if USE_GPU_SOLVER
            cout << "GPU elapsed time: " << std::fixed << std::setprecision(3) << gpuElapsedTime.count() << " ms." << endl;
        #endif
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
    int num_bits = 8;
    float q_levels = (1 << num_bits) - 1;
    int numVCores, numHCores, numVectorsPerCore, numElementsPerVector;
    size_t bytes_R, bytes_u, bytes_u_q, bytes_out, bytes_out_q;

    matrix<float> A, A_cpu, A_gpu_res, A_pho_res, A_cpu_res;
    matrix<uint8_t> A_pho_q;
    viennacl::matrix<float> A_gpu, A_pho;
    viennacl::matrix<float> Akk_tile, tempT;

    std::chrono::duration<float, std::milli> hostElapsedTimeGpu = std::chrono::duration<float, std::milli>::zero();
    std::chrono::duration<float, std::milli> quantElapsedTimeGpu = std::chrono::duration<float, std::milli>::zero();
    std::chrono::duration<float, std::milli> syncElapsedTimeGpu = std::chrono::duration<float, std::milli>::zero();
    std::chrono::duration<float, std::milli> gpuElapsedTime = std::chrono::duration<float, std::milli>::zero();
    std::chrono::duration<float, std::milli> transferElapsedTimeGpu = std::chrono::duration<float, std::milli>::zero();

    PhotonicsObjId matAObject;
    PhotonicsObjId outObject;
    PhotonicsObjId vecUObject;
    PhotonicsDeviceProperties deviceParams;

    #if USE_GPU_SOLVER
        float *d_R, *d_u, *d_one, *d_u_temp, *d_out;
        uint8_t *d_q, *d_out_q;
        cublasHandle_t handle;
        cudaStream_t stream; 
        NppStreamContext nppCtx;
    #endif

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

    // void calPhotonics(int i){

    //   quantize_vector(-1.0, 1.0);
    //   status = photonicsCopyHostToDevice((void *)vecU.data(), vecUObject);
    //   if (status != PHOTONICS_OK)
    //   {
    //     std::cout << "❌ Function: " << __func__ << "Abort: photonicsCopyHostToDevice failed between vecUObject and vecU at i=" << i << std::endl;
    //     return;
    //   }

    //   status = photonicsMvm(matAObject, vecUObject, outObject);
    //   if (status != PHOTONICS_OK)
    //   {
    //     std::cout << "❌ Function: " << __func__ << "Abort: photonicsMvm Failed at i=" << i << std::endl;
    //     return;
    //   }

    //   status = photonicsCopyDeviceToHost(outObject, out.data());
    //   if (status != PHOTONICS_OK)
    //   {
    //     std::cout << "❌ Function: " << __func__ << "Abort: photonicsCopyDeviceToHost failed between outObject and out at i=" << i << std::endl;
    //     return;
    //   }
    //   dequantize_vector(-1.0, 1.0);
    //   return;
    // }

    // // Vector quantization
    // void quantize_matrix(float min_val, float max_val) {
    //     const float scale = q_levels / (max_val - min_val);
    //     q_mat = (((R.array() - min_val) * scale) + 0.5f)
    //                 .cwiseMin(static_cast<float>(q_levels))
    //                 .cwiseMax(0.0f)
    //                 .cast<uint8_t>();
    //     for (int i = 0; i < N; i++) {
    //       for (int j = 0; j < N; j++) {
    //         matR[i * N + j] = q_mat(i,j);
    //       }
    //     }
    // }

    // void quantize_vector(float min_val, float max_val) {

    //     auto start2 = std::chrono::high_resolution_clock::now();
    //     const float scale = q_levels / (max_val - min_val);
    //     q_vec = (((u.array() - min_val) * scale) + 0.5f)
    //                 .cwiseMin(static_cast<float>(q_levels))
    //                 .cwiseMax(0.0f)
    //                 .cast<uint8_t>();
    //     quantElapsedTimeCpu += (std::chrono::high_resolution_clock::now() - start2);

    //     #if USE_GPU_SOLVER
    //         start2 = std::chrono::high_resolution_clock::now();
    //         Npp32f scale2 = q_levels / (max_val - min_val);
    //         Npp32f minVal = min_val;
    //         nppsSubC_32f_I(minVal, d_u, N);
    //         nppsMulC_32f_I(scale2, d_u, N);
    //         nppsConvert_32f8u_Sfs(d_u, d_q, N, NPP_RND_NEAR, 0);
    //         quantElapsedTimeGpu += (std::chrono::high_resolution_clock::now() - start2);
            
    //         start2 = std::chrono::high_resolution_clock::now();
    //         cudaMemcpy(q_vec_gpu.data(), d_q, bytes_u_q, cudaMemcpyDeviceToHost);
    //         transferElapsedTimeGpu += (std::chrono::high_resolution_clock::now() - start2);
    //     #endif

    //     std::memcpy(vecU.data(), q_vec.data(), bytes_u_q);
    // }

    // void dequantize_vector(float min_val, float max_val) {
    //     // Map raw 'out' data (uint8_t) into Eigen array of shape [N x numHCores]
    //     Map<const Array<uint8_t, Dynamic, Dynamic, RowMajor>>
    //         out_mat_temp(out.data(), N, numHCores);
    //     for (int i = 0; i < numHCores; i++) {
    //         out_mat_col[i] = out_mat_temp.col(i).cast<int>().matrix();
    //     }

    //     auto start2 = std::chrono::high_resolution_clock::now();
    //     for (int i = 1; i < numHCores; i++) {
    //         out_mat_col[0] = out_mat_col[0] + out_mat_col[i];
    //     }
    //     syncElapsedTimeCpu += (std::chrono::high_resolution_clock::now() - start2);

    //     start2 = std::chrono::high_resolution_clock::now();
    //     const float scale = (max_val - min_val) / q_levels;
    //     u_temp = (out_mat_col[0].array().cast<float>() * scale) + min_val;
    //     quantElapsedTimeCpu += (std::chrono::high_resolution_clock::now() - start2);

    //     #if USE_GPU_SOLVER
    //         // ----------------------
    //         // GPU version (Thrust)
    //         // ----------------------
    //         start2 = std::chrono::high_resolution_clock::now();
    //         cudaMemcpy(d_out_q, out.data(), bytes_out, cudaMemcpyHostToDevice);
    //         transferElapsedTimeGpu += (std::chrono::high_resolution_clock::now() - start2);

    //         int total = N * numHCores;
    //         start2 = std::chrono::high_resolution_clock::now();
    //         nppsConvert_8u32f(d_out_q, d_out, total);
    //         quantElapsedTimeGpu += (std::chrono::high_resolution_clock::now() - start2);
            
    //         start2 = std::chrono::high_resolution_clock::now();
    //         alpha_fused = alpha * scale;
    //         cublasSgemv(handle, CUBLAS_OP_T,
    //                     N, numHCores,
    //                     &alpha_fused,
    //                     d_out, N,
    //                     d_one, 1,
    //                     &beta,
    //                     d_u, 1);
    //         offset = alpha * min_val + gamma;
    //         cublasSaxpy(handle, N, &offset, d_one, 1, d_u, 1);
    //         hostElapsedTimeGpu += (std::chrono::high_resolution_clock::now() - start2);
    //     #endif
    // }

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
      cout << "Host elapsed time: " << std::fixed << std::setprecision(3) << hostElapsedTimeGpu.count() << " ms." << endl;
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

  // photonicsShowStats();
  solver.printGPUElapsedTime();

  solver.printGpuPhotonicsBreakdown();

  return 0;
}
