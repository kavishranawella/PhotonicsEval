// Test: C++ version of frozenk_nL_pde
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
          "\nUsage:  ./frozenk_nL_pde.out [options]"
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
    HeatEquationPaperSolver(struct Params params, int n = 16, float dt = 0.01, float a = 1.0, float t_final = 4.0, int num_bits = 8)
        : params(params), n(n), dt(dt), a(a), t_final(t_final), num_bits(num_bits)
    {
        N = n * n;
        h = 6.0 / (n - 1);
        n_steps = static_cast<int>(round(t_final / dt));

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

        matR.resize(N * N);
        vecU.resize(N);
    }

    void solveHardware(const vector<float>& target_times) {
        cout << "Solving Heat Equation using optical format on hardware..." << endl;

        hardware_solutions.clear();
        hardware_actual_times.clear();

        u=start_u;
        u_new = VectorXf(N);

        if(initializeHardware()!=PHOTONICS_OK) return;

        quantize_matrix(-1.0, 1.0);
        status = photonicsCopyHostToDevice((void *)matR.data(), matAObject);
        if (status != PHOTONICS_OK)
        {
          std::cout << "Function: " << __func__ << "Abort: photonicsCopyHostToDevice failed between matAObject and matA" << std::endl;
          return;
        }

        for (int k = 1; k <= n_steps; k++) {
            
            calPhotonics(k);
            if (status != PHOTONICS_OK) return;

            auto start = std::chrono::high_resolution_clock::now();

            float t = k * dt;

            // Update equation
            u = (dt * a / (h * h)) * u_new
                           + ((dt * a / (h * h)) * d + 1.0) * u
                           + qFunc(t) * dt;

            // Save solutions at requested times
            for (size_t i = 0; i < target_times.size(); i++) {
                if ((fabs(t - target_times[i]) < (dt / 2.0)) && (hardware_solutions.size() < i + 1)) {
                    MatrixXf sol(n, n);
                    for (int row = 0; row < n; row++)
                        for (int col = 0; col < n; col++)
                            sol(row, col) = u(row * n + col);

                    hardware_solutions.push_back(sol);
                    hardware_actual_times.push_back(t);
                    cout << "Saved solution at t = " << t << endl;
                }
            }
            auto end = std::chrono::high_resolution_clock::now();
            hostElapsedTime += (end - start);
        }
        photonicsFreeSrcVec(vecUObject);
        photonicsFreeDestVec(outObject);
        photonicsFreeMat(matAObject);
    }

    void solveSoftware(const vector<float>& target_times) {
        cout << "Solving Heat Equation using optical format on software..." << endl;

        software_solutions.clear();
        software_actual_times.clear();

        u=start_u;

        for (int k = 1; k <= n_steps; k++) {
            float t = k * dt;

            // Update equation
            VectorXf u_new = (dt * a / (h * h)) * R * u
                           + ((dt * a / (h * h)) * d + 1.0) * u
                           + qFunc(t) * dt;
            u = u_new;

            // Save solutions at requested times
            for (size_t i = 0; i < target_times.size(); i++) {
                if ((fabs(t - target_times[i]) < (dt / 2.0)) && (software_solutions.size() < i + 1)) {
                    MatrixXf sol(n, n);
                    for (int row = 0; row < n; row++)
                        for (int col = 0; col < n; col++)
                            sol(row, col) = u(row * n + col);

                    software_solutions.push_back(sol);
                    software_actual_times.push_back(t);
                    cout << "Saved solution at t = " << t << endl;
                }
            }
        }
    }

    void printHardwareSolutions() const {
        cout << "\nPrinting Hardware Solutions...\n" << endl;
        for (size_t k = 0; k < hardware_solutions.size(); k++) {
            cout << "\nSolution at t = " << hardware_actual_times[k] << ":\n";
            cout << hardware_solutions[k] << endl;
        }
    }

    void printSoftwareSolutions() const {
        cout << "\nPrinting Software Solutions...\n" << endl;
        for (size_t k = 0; k < software_solutions.size(); k++) {
            cout << "\nSolution at t = " << software_actual_times[k] << ":\n";
            cout << software_solutions[k] << endl;
        }
    }

    const vector<MatrixXf>& getHardwareSolutions() const { return hardware_solutions; }
    const vector<float>& getHardwareActualTimes() const { return hardware_actual_times; }

    const vector<MatrixXf>& getSoftwareSolutions() const { return software_solutions; }
    const vector<float>& getSoftwareActualTimes() const { return software_actual_times; }

    void printHostElapsedTime() { 
      cout << "Host elapsed time: " << std::fixed << std::setprecision(3) << hostElapsedTime.count() << " ms." << endl;
    }

    void compareSolutions(){
        cout << "\nComparing Hardware and Software Solutions...\n";

        if (hardware_solutions.size() != software_solutions.size()) {
            cerr << "Warning: Mismatch in number of saved solutions! "
                << "Hardware = " << hardware_solutions.size()
                << ", Software = " << software_solutions.size() << endl;
            return;
        }

        for (size_t k = 0; k < hardware_solutions.size(); k++) {
            const MatrixXf& H = hardware_solutions[k];
            const MatrixXf& S = software_solutions[k];

            if (H.rows() != S.rows() || H.cols() != S.cols()) {
                cerr << "Dimension mismatch at solution " << k << endl;
                continue;
            }

            // Compute difference
            MatrixXf diff = H - S;

            double l2_error = diff.norm();        // L2 norm of difference
            double l2_ref   = S.norm();           // L2 norm of reference
            double rel_error = (l2_ref > 1e-12) ? (l2_error / l2_ref) : l2_error;

            double max_abs_error = diff.cwiseAbs().maxCoeff();

            cout << "t = " << software_actual_times[k]
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
    int num_bits = 8;
    int q_levels = (1 << num_bits) - 1;

    vector<float> x, y;
    MatrixXf A, R;
    VectorXf u, u_new, start_u;
    vector<MatrixXf> hardware_solutions, software_solutions;
    vector<float> hardware_actual_times, software_actual_times;
    std::chrono::duration<float, std::milli> hostElapsedTime = std::chrono::duration<float, std::milli>::zero();

    PhotonicsObjId matAObject;
    PhotonicsObjId outObject;
    PhotonicsObjId vecUObject;

    PhotonicsStatus status = PHOTONICS_OK;
    std::vector<uint8_t> vecU;
    std::vector<uint8_t> matR;

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

    VectorXf qFunc(float t) {
        VectorXf q = VectorXf::Ones(N);
        return 0.1 * sin(2.0 * M_PI * t) * q;
    }

    PhotonicsStatus initializeHardware(){
      if (!initPhotonicsAccel(params.dramConfigFile))
        return PHOTONICS_ERROR;

      matAObject = photonicsAllocMat(N*N, PHOTONICS_FP32);
      if (matAObject == -1)
      {
        std::cout << "Function: " << __func__ << "Abort: photonicsAlloc failed for matrix A" << std::endl;
        return PHOTONICS_ERROR;
      }
      
      outObject = photonicsAllocAssociatedDestVec(matAObject, PHOTONICS_FP32);
      if (outObject == -1)
      {
        std::cout << "Function: " << __func__ << "Abort: photonicsAllocAssociated failed for vector U = " << outObject << std::endl;
        return PHOTONICS_ERROR;
      }
      
      vecUObject = photonicsAllocAssociatedSrcVec(matAObject, PHOTONICS_FP32);
      if (vecUObject == -1)
      {
        std::cout << "Function: " << __func__ << "Abort: photonicsAllocAssociated failed for vector U = " << vecUObject << std::endl;
        return PHOTONICS_ERROR;
      }

      return PHOTONICS_OK;

    }

    void calPhotonics(int i){

      quantize_vector(-1.0, 1.0);
      status = photonicsCopyHostToDevice((void *)vecU.data(), vecUObject);
      if (status != PHOTONICS_OK)
      {
        std::cout << "Function: " << __func__ << "Abort: photonicsCopyHostToDevice failed between vecUObject and vecU at i=" << i << std::endl;
        return;
      }

      status = photonicsIter(matAObject, vecUObject, outObject, 1);
      if (status != PHOTONICS_OK)
      {
        std::cout << "Function: " << __func__ << "Abort: photonicsIter Failed at i=" << i << std::endl;
        return;
      }

      status = photonicsCopyDeviceToHost(vecUObject, vecU.data());
      if (status != PHOTONICS_OK)
      {
        std::cout << "Function: " << __func__ << "Abort: photonicsCopyDeviceToHost failed between vecUObject and vecU at i=" << i << std::endl;
        return;
      }
      dequantize_vector(-1.0, 1.0);
      return;
    }

    // Vector quantization
    void quantize_matrix(float min_val, float max_val) {
        float scale = (max_val - min_val) / q_levels;
        for (int i = 0; i < N; i++) {
          for (int j = 0; j < N; j++) {
            int q = std::round((R(i,j) - min_val) / scale);
            q = std::max(0, std::min(q, q_levels));
            matR[i * N + j] = static_cast<uint8_t>(q);
          }
        }
    }

    void quantize_vector(float min_val, float max_val) {
        float scale = (max_val - min_val) / q_levels;
        for (int i = 0; i < N; i++) {
            int q = std::round((u(i) - min_val) / scale);
            q = std::max(0, std::min(q, q_levels));
            vecU[i] = static_cast<uint8_t>(q);
        }
    }

    void dequantize_vector(float min_val, float max_val) {
        float scale = (max_val - min_val) / q_levels;
        for (int i = 0; i < N; i++) {
            u_new(i) = vecU[i] * scale + min_val;
        }
    }
};

int main(int argc, char *argv[])
{
  struct Params params = getInputParams(argc, argv);

  params.dramConfigFile = (char *) "./configs/photonics/frozenk_nL_pde.cfg";

  HeatEquationPaperSolver solver(params, 16, 0.01, 1.0, 4.0, 8);

  vector<float> target_times = {0.25, 1.25, 3.75};

  solver.solveHardware(target_times);
  //solver.printHardwareSolutions();

  solver.solveSoftware(target_times);
  //solver.printSoftwareSolutions();

  solver.compareSolutions();

  photonicsShowStats();
  solver.printHostElapsedTime();

  return 0;
}
