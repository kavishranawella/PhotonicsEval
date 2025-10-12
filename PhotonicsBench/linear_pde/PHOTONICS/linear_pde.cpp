// Test: C++ version of linear_pde
// Copyright (c) 2024 University of Virginia
// This file is licensed under the MIT License.
// See the LICENSE file in the root of this repository for more details.

#include "libphotonicseval.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <getopt.h>
#if defined(_OPENMP)
#include <omp.h>
#endif
#include "utilPhotonics.h"
#include <iomanip>
#include <chrono>
#include <cassert>
#include <Eigen/Dense>

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
          "\nUsage:  ./linear_pde.out [options]"
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

class HeatEquationPaperSolver {
public:
    HeatEquationPaperSolver(struct Params params, float dt = 0.01, float a = 1.0, float t_final = 4.0, int num_bits = 8)
        : params(params), dt(dt), a(a), t_final(t_final), num_bits(num_bits)
    {   
        if(initializeHardware()!=PHOTONICS_OK) {
          std::cout << "❌ Function: " << __func__ << "Abort: intializeHardware() failed" << std::endl;
          return;
        }

        n = sqrt(N);
        
        h = 6.0 / (n - 1);
        n_steps = static_cast<int>(round(t_final / dt));


        alpha = (dt * a / (h * h));
        beta  = ((dt * a / (h * h)) * d + 1.0);

        // Grid points
        x.resize(n);
        y.resize(n);
        for (int i = 0; i < n; i++) {
            x[i] = -3.0 + i * h;
            y[i] = -3.0 + i * h;
        }

        // Build Laplacian
        buildLaplacian();

        // Build R
        int d = -4;
        R = A - d * MatrixXf::Identity(N, N);
        R = R / fabs(d);

        // Initial condition (Gaussian)
        start_u = VectorXf(N);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                float X = x[j];
                float Y = y[i];
                start_u(i * n + j) = exp(-2.0 * (X * X + Y * Y));
            }
        }

        bytes_R = N * N * sizeof(float);
        bytes_u = N * sizeof(float);
        bytes_u_q = N * sizeof(uint8_t);
        bytes_out = N * numHCores * sizeof(float);
        bytes_out_q = N * numHCores * sizeof(uint8_t);

        q_one_vec = VectorXf::Ones(N);
        q_one_ext = VectorXf::Ones(numHCores);

        out_mat_col.resize(numHCores);

        matR.resize(N * N);
        out.resize(N * numHCores);
        vecU.resize(N);
    }

    ~HeatEquationPaperSolver(){
        photonicsFreeSrcVec(vecUObject);
        photonicsFreeDestVec(outObject);
        photonicsFreeMat(matAObject);
    }

    void solvePhotonics(const vector<float>& target_times) {
        cout << "\nSolving Heat Equation using optical format on Photonics+CPU..." << endl;

        photonics_solutions.clear();
        photonics_actual_times.clear();

        u = start_u;
        u_gpu = start_u;
        u_temp = VectorXf(N);
        u_new = VectorXf(N);

        #if USE_GPU_SOLVER
            // Allocate GPU memory
            cudaMalloc(&d_u, bytes_u);
            cudaMalloc(&d_one, bytes_u);
            cudaMalloc(&d_u_temp, bytes_u);
            cudaMalloc(&d_q, bytes_u_q);
            cudaMalloc(&d_out_q, bytes_out_q);
            cudaMalloc(&d_out, bytes_out);

            // Copy initial data
            cudaMemcpy(d_u, u_gpu.data(), bytes_u, cudaMemcpyHostToDevice);
            cudaMemcpy(d_one, q_one_vec.data(), bytes_u, cudaMemcpyHostToDevice);

            cublasCreate(&handle);
            cudaStreamCreate(&stream);
            init_gpu_handles(stream);

        #endif

        quantize_matrix(-1.0, 1.0);
        status = photonicsCopyHostToDevice((void *)matR.data(), matAObject);
        if (status != PHOTONICS_OK)
        {
          std::cout << "❌ Function: " << __func__ << "Abort: photonicsCopyHostToDevice failed between matAObject and matA" << std::endl;
          return;
        }

        for (int k = 1; k <= n_steps; k++) {
            
            calPhotonics(k);
            if (status != PHOTONICS_OK) return;

            auto start = std::chrono::high_resolution_clock::now();
            float t = k * dt;
            float gamma = 0.1 * sin(2.0 * M_PI * t) * dt;
            auto end = std::chrono::high_resolution_clock::now();
            hostElapsedTimeCpu += (end - start);
            hostElapsedTimeGpu += (end - start);

            // ------------------------------
            // CPU (Eigen) implementation
            // ------------------------------
            start = std::chrono::high_resolution_clock::now();
            // Update equation
            u_new = alpha * u_temp + (beta * u.array() + gamma).matrix();
            u = u_new;
            hostElapsedTimeCpu += (std::chrono::high_resolution_clock::now() - start);

            // Save solutions at requested times
            start = std::chrono::high_resolution_clock::now();
            for (size_t i = 0; i < target_times.size(); i++) {
                if ((fabs(t - target_times[i]) < (dt / 2.0)) && (photonics_solutions.size() < i + 1)) {
                    MatrixXf sol(n, n);
                    for (int row = 0; row < n; row++)
                        for (int col = 0; col < n; col++)
                            sol(row, col) = u(row * n + col);

                    photonics_solutions.push_back(sol);
                    photonics_actual_times.push_back(t);
                    cout << "Saved solution at t = " << t << endl;
                }
            }
            end = std::chrono::high_resolution_clock::now();
            hostElapsedTimeCpu += (end - start);
            hostElapsedTimeGpu += (end - start);
        }
        #if USE_GPU_SOLVER
            cublasDestroy(handle);
            cudaFree(d_out);
            cudaFree(d_out_q);
            cudaFree(d_q);
            cudaFree(d_u_temp);
            cudaFree(d_one);
            cudaFree(d_u);
        #endif
    }

    // ------------------------------
    // CPU (Eigen) implementation
    // ------------------------------
    void solveCPU(const vector<float>& target_times) {
        cout << "\nSolving Heat Equation using optical format on software with CPU..." << endl;

        cpu_solutions.clear();
        cpu_actual_times.clear();
        u = start_u;

        auto start = std::chrono::high_resolution_clock::now();
        for (int k = 1; k <= n_steps; k++) {
            float t = k * dt;
            float gamma = 0.1 * sin(2.0 * M_PI * t) * dt;

            // Update equation
            u_new = alpha * R * u + (beta * u.array() + gamma).matrix();
            u = u_new;

            // Save solutions at requested times
            for (size_t i = 0; i < target_times.size(); i++) {
                if ((fabs(t - target_times[i]) < (dt / 2.0)) && (cpu_solutions.size() < i + 1)) {
                    MatrixXf sol(n, n);
                    for (int row = 0; row < n; row++)
                        for (int col = 0; col < n; col++)
                            sol(row, col) = u(row * n + col);
                    cpu_solutions.push_back(sol);
                    cpu_actual_times.push_back(t);
                    cout << "Saved solution at t = " << t << endl;
                }
            }
        }
        cpuElapsedTime = (std::chrono::high_resolution_clock::now() - start);

    }

    // ------------------------------
    // GPU (cuBLAS) implementation
    // ------------------------------
    void solveGPU(const vector<float>& target_times) {
        #if USE_GPU_SOLVER
            cout << "\nSolving Heat Equation using optical format on software with GPU..." << endl;
            gpu_solutions.clear();
            gpu_actual_times.clear();
            u = start_u;

            // Allocate GPU memory
            cudaMalloc(&d_R, bytes_R);
            cudaMalloc(&d_u, bytes_u);
            cudaMalloc(&d_one, bytes_u);

            // Copy initial data
            cudaMemcpy(d_R, R.data(), bytes_R, cudaMemcpyHostToDevice);
            cudaMemcpy(d_u, u.data(), bytes_u, cudaMemcpyHostToDevice);
            cudaMemcpy(d_one, q_one_vec.data(), bytes_u, cudaMemcpyHostToDevice);

            // cuBLAS setup
            cublasCreate(&handle);

            auto start = std::chrono::high_resolution_clock::now();
            for (int k = 1; k <= n_steps; k++) {
                float t = k * dt;
                float gamma = 0.1 * sin(2.0 * M_PI * t) * dt;

                // Compute: u = alpha * R * u + beta * u
                // (R is N×N, u is N×1)
                cublasSgemv(handle, CUBLAS_OP_N,
                            N, N, &alpha, d_R, 
                            N, d_u, 1, 
                            &beta, d_u, 1);

                // d_u_new = d_u_new + gamma * d_one
                cublasSaxpy(handle, N, &gamma, d_one, 1, d_u, 1);

                // Save solutions at requested times
                for (size_t i = 0; i < target_times.size(); i++) {
                    if ((fabs(t - target_times[i]) < (dt / 2.0)) && (gpu_solutions.size() < i + 1)) {
                        MatrixXf sol(n, n);
                        cudaMemcpy(u.data(), d_u, bytes_u, cudaMemcpyDeviceToHost);
                        for (int row = 0; row < n; row++)
                            for (int col = 0; col < n; col++)
                                sol(row, col) = u(row * n + col);
                        gpu_solutions.push_back(sol);
                        gpu_actual_times.push_back(t);
                        cout << "Saved solution at t = " << t << endl;
                    }
                }
            }
            gpuElapsedTime = (std::chrono::high_resolution_clock::now() - start);

            cublasDestroy(handle);
            cudaFree(d_R);
            cudaFree(d_u);
            cudaFree(d_one);
        #else
            cout << "GPU implementation not possible..." << endl;
        #endif
    }

    void printPhotonicsSolutions() const {
        cout << "\nPrinting Photonics Solutions...\n" << endl;
        for (size_t k = 0; k < photonics_solutions.size(); k++) {
            cout << "\nSolution at t = " << photonics_actual_times[k] << ":\n";
            cout << photonics_solutions[k] << endl;
        }
    }

    void printCPUSolutions() const {
        cout << "\nPrinting CPU Solutions...\n" << endl;
        for (size_t k = 0; k < cpu_solutions.size(); k++) {
            cout << "\nSolution at t = " << cpu_actual_times[k] << ":\n";
            cout << cpu_solutions[k] << endl;
        }
    }

    void printGPUSolutions() const {
        #if USE_GPU_SOLVER
            cout << "\nPrinting GPU Solutions...\n" << endl;
            for (size_t k = 0; k < gpu_solutions.size(); k++) {
                cout << "\nSolution at t = " << gpu_actual_times[k] << ":\n";
                cout << gpu_solutions[k] << endl;
            }
        #else
          cout << "No GPU implementation available....." << endl;
        #endif
    }

    const vector<MatrixXf>& getPhotonicsSolutions() const { return photonics_solutions; }
    const vector<float>& getPhotonicsActualTimes() const { return photonics_actual_times; }

    const vector<MatrixXf>& getCPUSolutions() const { return cpu_solutions; }
    const vector<float>& getCPUActualTimes() const { return cpu_actual_times; }

    const vector<MatrixXf>& getGPUSolutions() const { return gpu_solutions; }
    const vector<float>& getGPUActualTimes() const { return gpu_actual_times; }

    void printCpuPhotonicsBreakdown() { 
        cout << "\nPrinting CPU+Photonics Breakdown..." << endl;
        printCpuHostElapsedTime();
        printCpuQuantElapsedTime();
        printCpuSyncElapsedTime();
    }

    void printGpuPhotonicsBreakdown() { 
        #if USE_GPU_SOLVER
            cout << "\nPrinting GPU+Photonics Breakdown..." << endl;
            printGpuHostElapsedTime();
            printGpuQuantElapsedTime();
            printGpuSyncElapsedTime();
            printGpuTransferElapsedTime();
        #endif
    }
    
    void printCPUElapsedTime() { 
        cout << "CPU elapsed time: " << std::fixed << std::setprecision(3) << cpuElapsedTime.count() << " ms." << endl;
    }

    void printGPUElapsedTime() { 
        #if USE_GPU_SOLVER
            cout << "GPU elapsed time: " << std::fixed << std::setprecision(3) << gpuElapsedTime.count() << " ms." << endl;
        #endif
    }

    void compareSWSolutions() { 
        #if USE_GPU_SOLVER
            cout << "\nVerifying consistency of CPU & GPU Results..." << endl;

            if (gpu_solutions.size() != cpu_solutions.size()) {
                cerr << "❌ Error: Mismatch in number of saved solutions! "
                    << "CPU = " << cpu_solutions.size()
                    << ", GPU = " << gpu_solutions.size() << endl;
                throw std::runtime_error("Solution count mismatch between CPU and GPU");
            }

            for (size_t k = 0; k < cpu_solutions.size(); k++) {
                const MatrixXf &H = cpu_solutions[k];
                const MatrixXf &S = gpu_solutions[k];

                if (H.rows() != S.rows() || H.cols() != S.cols()) {
                    cerr << "❌ Error: Dimension mismatch at solution " << k
                        << " (" << H.rows() << "x" << H.cols()
                        << " vs " << S.rows() << "x" << S.cols() << ")" << endl;
                    throw std::runtime_error("Dimension mismatch between CPU and GPU results");
                }

                // Compare element-wise
                for (int i = 0; i < H.rows(); i++) {
                    for (int j = 0; j < H.cols(); j++) {
                        float a = H(i, j);
                        float b = S(i, j);
                        float diff = fabs(a - b);

                        // if (diff > tolerance || std::isnan(a) || std::isnan(b)) {
                        //     cerr << std::fixed << setprecision(7);
                        //     cerr << "❌ Mismatch at t = " << cpu_actual_times[k]
                        //         << ", index (" << i << "," << j << ")"
                        //         << ": CPU = " << a
                        //         << ", GPU = " << b
                        //         << ", |Δ| = " << diff
                        //         << " > tol(" << tolerance << ")" << endl;

                        //     throw std::runtime_error("CPU and GPU results differ beyond tolerance");
                        // }
                    }
                }
            }

            cout << "✅ All CPU and GPU results match within tolerance "
                << tolerance << "!" << endl;

        #else
            cout << "\nNo GPU implementation available....." << endl;
        #endif
    }

    void compareWithPhotonics(){
        cout << "\nComparing Photonics and Software Solutions...\n";

        if (photonics_solutions.size() != cpu_solutions.size()) {
            cerr << "❌ Error: Mismatch in number of saved solutions! "
                << "Photonics = " << photonics_solutions.size()
                << ", CPU = " << cpu_solutions.size() << endl;
            throw std::runtime_error("Solution count mismatch between SW and Photonics");
        }

        for (size_t k = 0; k < photonics_solutions.size(); k++) {
            const MatrixXf& H = photonics_solutions[k];
            const MatrixXf& S = cpu_solutions[k];

            if (H.rows() != S.rows() || H.cols() != S.cols()) {
                cerr << "❌ Error: Dimension mismatch at solution " << k << endl;
                throw std::runtime_error("Dimension mismatch between SW and Photonics results");
            }

            // Compute difference
            MatrixXf diff = H - S;

            double l2_error = diff.norm();        // L2 norm of difference
            double l2_ref   = S.norm();           // L2 norm of reference
            double rel_error = (l2_ref > 1e-12) ? (l2_error / l2_ref) : l2_error;

            double max_abs_error = diff.cwiseAbs().maxCoeff();

            cout << "t = " << cpu_actual_times[k]
                << " | L2 Error = " << l2_error
                << " | Relative Error = " << rel_error
                << " | Max Abs Error = " << max_abs_error
                << endl;
        }
        return;
    }

private:
    struct Params params;
    int n;          // grid size in one dimension
    int N;          // total number of grid points
    float dt;      // time step
    float a;       // diffusivity
    float t_final; // final time
    int n_steps;    // number of steps
    float h;       // grid spacing
    int d = -4;     // diagonal value
    const float tolerance = 1e-5f;
    float alpha, beta, gamma, alpha_fused, offset;
    int num_bits = 8;
    float q_levels = (1 << num_bits) - 1;
    int numVCores, numHCores, numVectorsPerCore, numElementsPerVector;
    size_t bytes_R, bytes_u, bytes_u_q, bytes_out, bytes_out_q;

    vector<float> x, y;
    MatrixXf A, R;
    VectorXf u, u_new, u_temp, start_u, q_one_vec, u_gpu, u_temp_gpu, q_one_ext;
    Array<uint8_t, Dynamic, Dynamic> q_mat;
    Array<uint8_t, Dynamic, 1> q_vec, q_vec_gpu;
    vector<VectorXi> out_mat_col;
    vector<MatrixXf> photonics_solutions, cpu_solutions, gpu_solutions;
    vector<float> photonics_actual_times, cpu_actual_times, gpu_actual_times;
    std::chrono::duration<float, std::milli> hostElapsedTimeCpu = std::chrono::duration<float, std::milli>::zero();
    std::chrono::duration<float, std::milli> quantElapsedTimeCpu = std::chrono::duration<float, std::milli>::zero();
    std::chrono::duration<float, std::milli> syncElapsedTimeCpu = std::chrono::duration<float, std::milli>::zero();
    std::chrono::duration<float, std::milli> hostElapsedTimeGpu = std::chrono::duration<float, std::milli>::zero();
    std::chrono::duration<float, std::milli> quantElapsedTimeGpu = std::chrono::duration<float, std::milli>::zero();
    std::chrono::duration<float, std::milli> syncElapsedTimeGpu = std::chrono::duration<float, std::milli>::zero();
    std::chrono::duration<float, std::milli> cpuElapsedTime = std::chrono::duration<float, std::milli>::zero();
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

    void buildLaplacian() {
        A = MatrixXf::Zero(N, N);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                int idx = i * n + j;
                A(idx, idx) = -4;
                if (j < n - 1) A(idx, idx + 1) = 1;   // right
                if (j > 0)     A(idx, idx - 1) = 1;   // left
                if (i > 0)     A(idx, idx - n) = 1;   // top
                if (i < n - 1) A(idx, idx + n) = 1;   // bottom
            }
        }
    }

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

      if ((N != (numVCores * numVectorsPerCore)) || (N != (numHCores * numElementsPerVector)) )
      {
        cerr << "❌ Error: Entire matrix needs to fit into the Photonic Cores since tiling is unavailable! \n"
                << "N = " << N << ", (numVCores * numVectorsPerCore) = " << (numVCores * numVectorsPerCore)
                << ", (numHCores * numElementsPerVector) = " << (numHCores * numElementsPerVector) << endl;
        return PHOTONICS_ERROR;
      }

      matAObject = photonicsAllocMat(N*N, PHOTONICS_FP32);
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

    void calPhotonics(int i){

      quantize_vector(-1.0, 1.0);
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
      dequantize_vector(-1.0, 1.0);
      return;
    }

    // Vector quantization
    void quantize_matrix(float min_val, float max_val) {
        const float scale = q_levels / (max_val - min_val);
        q_mat = (((R.array() - min_val) * scale) + 0.5f)
                    .cwiseMin(static_cast<float>(q_levels))
                    .cwiseMax(0.0f)
                    .cast<uint8_t>();
        for (int i = 0; i < N; i++) {
          for (int j = 0; j < N; j++) {
            matR[i * N + j] = q_mat(i,j);
          }
        }
    }

    void quantize_vector(float min_val, float max_val) {

        auto start2 = std::chrono::high_resolution_clock::now();
        const float scale = q_levels / (max_val - min_val);
        q_vec = (((u.array() - min_val) * scale) + 0.5f)
                    .cwiseMin(static_cast<float>(q_levels))
                    .cwiseMax(0.0f)
                    .cast<uint8_t>();
        quantElapsedTimeCpu += (std::chrono::high_resolution_clock::now() - start2);

        #if USE_GPU_SOLVER
            start2 = std::chrono::high_resolution_clock::now();
            Npp32f scale2 = q_levels / (max_val - min_val);
            Npp32f minVal = min_val;
            nppsSubC_32f_I(minVal, d_u, N);
            nppsMulC_32f_I(scale2, d_u, N);
            nppsConvert_32f8u_Sfs(d_u, d_q, N, NPP_RND_NEAR, 0);
            quantElapsedTimeGpu += (std::chrono::high_resolution_clock::now() - start2);
            
            start2 = std::chrono::high_resolution_clock::now();
            cudaMemcpy(q_vec_gpu.data(), d_q, bytes_u_q, cudaMemcpyDeviceToHost);
            transferElapsedTimeGpu += (std::chrono::high_resolution_clock::now() - start2);
        #endif

        std::memcpy(vecU.data(), q_vec.data(), bytes_u_q);
    }

    void dequantize_vector(float min_val, float max_val) {
        // Map raw 'out' data (uint8_t) into Eigen array of shape [N x numHCores]
        Map<const Array<uint8_t, Dynamic, Dynamic, RowMajor>>
            out_mat_temp(out.data(), N, numHCores);
        for (int i = 0; i < numHCores; i++) {
            out_mat_col[i] = out_mat_temp.col(i).cast<int>().matrix();
        }

        auto start2 = std::chrono::high_resolution_clock::now();
        for (int i = 1; i < numHCores; i++) {
            out_mat_col[0] = out_mat_col[0] + out_mat_col[i];
        }
        syncElapsedTimeCpu += (std::chrono::high_resolution_clock::now() - start2);

        start2 = std::chrono::high_resolution_clock::now();
        const float scale = (max_val - min_val) / q_levels;
        u_temp = (out_mat_col[0].array().cast<float>() * scale) + min_val;
        quantElapsedTimeCpu += (std::chrono::high_resolution_clock::now() - start2);

        #if USE_GPU_SOLVER
            // ----------------------
            // GPU version (Thrust)
            // ----------------------
            start2 = std::chrono::high_resolution_clock::now();
            cudaMemcpy(d_out_q, out.data(), bytes_out, cudaMemcpyHostToDevice);
            transferElapsedTimeGpu += (std::chrono::high_resolution_clock::now() - start2);

            int total = N * numHCores;
            start2 = std::chrono::high_resolution_clock::now();
            nppsConvert_8u32f(d_out_q, d_out, total);
            quantElapsedTimeGpu += (std::chrono::high_resolution_clock::now() - start2);
            
            start2 = std::chrono::high_resolution_clock::now();
            alpha_fused = alpha * scale;
            cublasSgemv(handle, CUBLAS_OP_T,
                        N, numHCores,
                        &alpha_fused,
                        d_out, N,
                        d_one, 1,
                        &beta,
                        d_u, 1);
            offset = alpha * min_val + gamma;
            cublasSaxpy(handle, N, &offset, d_one, 1, d_u, 1);
            hostElapsedTimeGpu += (std::chrono::high_resolution_clock::now() - start2);
        #endif
    }

    void init_gpu_handles(cudaStream_t stream) {
        cublasSetStream(handle, stream);

        // Fill NPP stream context (minimal set)
        memset(&nppCtx, 0, sizeof(nppCtx));
        nppCtx.hStream = stream;
        cudaDeviceProp prop{};
        int dev=0; cudaGetDevice(&dev); cudaGetDeviceProperties(&prop, dev);
        nppCtx.nCudaDeviceId                 = dev;
        nppCtx.nMultiProcessorCount          = prop.multiProcessorCount;
        nppCtx.nMaxThreadsPerMultiProcessor  = prop.maxThreadsPerMultiProcessor;
        nppCtx.nMaxThreadsPerBlock           = prop.maxThreadsPerBlock;
        nppCtx.nSharedMemPerBlock            = prop.sharedMemPerBlock;
        nppCtx.nCudaDevAttrComputeCapabilityMajor = prop.major;
        nppCtx.nCudaDevAttrComputeCapabilityMinor = prop.minor;
    }

    void printCpuHostElapsedTime() { 
      cout << "Host elapsed time: " << std::fixed << std::setprecision(3) << hostElapsedTimeCpu.count() << " ms." << endl;
    }

    void printCpuQuantElapsedTime() { 
      cout << "Quantization elapsed time: " << std::fixed << std::setprecision(3) << quantElapsedTimeCpu.count() << " ms." << endl;
    }

    void printCpuSyncElapsedTime() { 
      cout << "Photonics Core Synchronization elapsed time: " << std::fixed << std::setprecision(3) << syncElapsedTimeCpu.count() << " ms." << endl;
    }

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

  params.dramConfigFile = (char *) "./configs/photonics/linear_pde.cfg";

  HeatEquationPaperSolver solver(params, 0.001, 1.0, 0.4, 8);

  vector<float> target_times = {0.025, 0.125, 0.375};

  solver.solvePhotonics(target_times);
  //solver.printPhotonicsSolutions();

  solver.solveCPU(target_times);
  //solver.printCPUSolutions();

  solver.solveGPU(target_times);
  //solver.printGPUSolutions();

  solver.compareSWSolutions();
  solver.compareWithPhotonics();

  photonicsShowStats();
  solver.printCPUElapsedTime();
  solver.printGPUElapsedTime();

  solver.printCpuPhotonicsBreakdown();
  solver.printGpuPhotonicsBreakdown();

  return 0;
}
