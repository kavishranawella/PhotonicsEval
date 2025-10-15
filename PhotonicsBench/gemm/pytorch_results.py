import torch

# ----------------------------
# Config
# ----------------------------
alpha = 1
N, K1, K2 = 16384*alpha, 256, 2048
chunk1, chunk2 = 1024//(alpha**2), 1024//(alpha**2)
device = "cuda"
L1 = N // K1
L2 = N // K2
M1 = N//((K2 // K1)**2)

# ----------------------------
# CPU tensors
# ----------------------------
a_cpu = torch.randint(-1, 1, (N, N), dtype=torch.int8, pin_memory=True)
b_cpu = torch.randint(-1, 1, (N, N), dtype=torch.int8, pin_memory=True)
tensor_cpu_1 = torch.empty((N, N, L1), dtype=torch.int8, pin_memory=True)
tensor_cpu_2 = torch.empty((N, N, L2), dtype=torch.int8, pin_memory=True)

# ----------------------------
# First Approach: Direct GEMM on GPU
# ----------------------------
torch.cuda.synchronize()
start_transfer = torch.cuda.Event(enable_timing=True)
end_transfer = torch.cuda.Event(enable_timing=True)

start_transfer.record()
a_gpu = a_cpu.to(device=device, dtype=torch.float16, non_blocking=True)
b_gpu = b_cpu.to(device=device, dtype=torch.float16, non_blocking=True)
end_transfer.record()
torch.cuda.synchronize()
gpu_transfer_ms = start_transfer.elapsed_time(end_transfer)

# Warm-up
_ = a_gpu @ b_gpu
torch.cuda.synchronize()

start_compute = torch.cuda.Event(enable_timing=True)
end_compute = torch.cuda.Event(enable_timing=True)
start_compute.record()
c_gpu = a_gpu @ b_gpu
end_compute.record()
torch.cuda.synchronize()
gpu_compute_ms = start_compute.elapsed_time(end_compute)

result_cpu = c_gpu.to('cpu')

# Compute TFLOPs
flops = (2 * N - 1) * (N ** 2)
tflops = flops / (gpu_compute_ms / 1e3) / 1e12

print("\n=== Baseline GPU GEMM ===")
print(f"Matrix size: {N}x{N}")
print(f"Transfer time: {gpu_transfer_ms:.3f} ms")
print(f"Compute time: {gpu_compute_ms:.3f} ms")
print(f"Performance: {tflops:.2f} TFLOP/s")
print(f"Input shape: {a_cpu.shape}")
print(f"Result shape: {c_gpu.shape}")
print(f"Result dtype: {c_gpu.dtype}")

# ----------------------------
# Second Approach: Chunked “Photonics” style with Photonics MMM
# ----------------------------
print("\n=== Photonics-like Blocked GEMM + Reduction with Photonics MMM ===")

# Step 1: Block multiply and store each intermediate result to CPU
for j, i in enumerate(range(0, N, K1)):
    a_block_gpu = a_cpu[:, i:i + K1].to(device=device, dtype=torch.float16, non_blocking=True)
    b_block_gpu = b_cpu[i:i + K1, :].to(device=device, dtype=torch.float16, non_blocking=True)
    c_block_gpu = a_block_gpu @ b_block_gpu
    tensor_cpu_1[:, :, j] = c_block_gpu.to('cpu', dtype=torch.int8, non_blocking=False)
torch.cuda.synchronize()

# Step 2: Summation along 3rd dim in chunks
result_gpu_1 = torch.zeros((N, N), dtype=torch.float16, device=device)

# Warm up GPU context
_ = torch.sum(torch.zeros((1, 1, 1), device=device, dtype=torch.float16))
torch.cuda.synchronize()

pho_mmm_host_transfer_ms = 0.0
pho_mmm_host_compute_ms = 0.0

for i in range(0, N, chunk1):
    # --- Transfer timing ---
    start_t = torch.cuda.Event(enable_timing=True)
    end_t = torch.cuda.Event(enable_timing=True)

    start_t.record()
    chunk_gpu = tensor_cpu_1[i:i + chunk1, :, :].to(device, dtype=torch.float16, non_blocking=True)
    end_t.record()
    torch.cuda.synchronize()

    transfer_ms = start_t.elapsed_time(end_t)
    pho_mmm_host_transfer_ms += transfer_ms

    # --- Compute timing ---
    start_c = torch.cuda.Event(enable_timing=True)
    end_c = torch.cuda.Event(enable_timing=True)

    start_c.record()
    result_gpu_1[i:i + chunk1, :] += torch.sum(chunk_gpu, dim=2)
    end_c.record()
    torch.cuda.synchronize()

    compute_ms = start_c.elapsed_time(end_c)
    pho_mmm_host_compute_ms += compute_ms

result_cpu_1 = result_gpu_1.to('cpu')

pho_mmm_compute_ms = 50 * 1e-9 * L1 * L1 * M1 

print(f"Photonics tensor shape: {tensor_cpu_1.shape}")
print(f"Result shape: {result_gpu_1.shape}")
print(f"Host transfer time: {pho_mmm_host_transfer_ms:.3f} ms")
print(f"Host compute time: {pho_mmm_host_compute_ms:.3f} ms")
print(f"Photonics compute time: {pho_mmm_compute_ms:.3f} ms")

# ----------------------------
# Third Approach: Chunked “Photonics” style with Photonics MVM
# ----------------------------
print("\n=== Photonics-like Blocked GEMM + Reduction with Photonics MVM ===")

# Step 1: Block multiply and store each intermediate result to CPU
for j, i in enumerate(range(0, N, K2)):
    a_block_gpu = a_cpu[:, i:i + K2].to(device=device, dtype=torch.float16, non_blocking=True)
    b_block_gpu = b_cpu[i:i + K2, :].to(device=device, dtype=torch.float16, non_blocking=True)
    c_block_gpu = a_block_gpu @ b_block_gpu
    tensor_cpu_2[:, :, j] = c_block_gpu.to('cpu', dtype=torch.int8, non_blocking=False)
torch.cuda.synchronize()

# Step 2: Summation along 3rd dim in chunks
result_gpu_2 = torch.zeros((N, N), dtype=torch.float16, device=device)

# Warm up GPU context
_ = torch.sum(torch.zeros((1, 1, 1), device=device, dtype=torch.float16))
torch.cuda.synchronize()

pho_mvm_host_transfer_ms = 0.0
pho_mvm_host_compute_ms = 0.0

for i in range(0, N, chunk2):
    # --- Transfer timing ---
    start_t = torch.cuda.Event(enable_timing=True)
    end_t = torch.cuda.Event(enable_timing=True)

    start_t.record()
    chunk_gpu = tensor_cpu_2[i:i + chunk2, :, :].to(device, dtype=torch.float16, non_blocking=True)
    end_t.record()
    torch.cuda.synchronize()

    transfer_ms = start_t.elapsed_time(end_t)
    pho_mvm_host_transfer_ms += transfer_ms

    # --- Compute timing ---
    start_c = torch.cuda.Event(enable_timing=True)
    end_c = torch.cuda.Event(enable_timing=True)

    start_c.record()
    result_gpu_2[i:i + chunk2, :] += torch.sum(chunk_gpu, dim=2)
    end_c.record()
    torch.cuda.synchronize()

    compute_ms = start_c.elapsed_time(end_c)
    pho_mvm_host_compute_ms += compute_ms

result_cpu_2 = result_gpu_2.to('cpu')

pho_mvm_compute_ms = 50 * 1e-9 * N * L2 * L2

print(f"Photonics tensor shape: {tensor_cpu_2.shape}")
print(f"Result shape: {result_gpu_2.shape}")
print(f"Host transfer time: {pho_mvm_host_transfer_ms:.3f} ms")
print(f"Host compute time: {pho_mvm_host_compute_ms:.3f} ms")
print(f"Photonics compute time: {pho_mvm_compute_ms:.3f} ms")

# ----------------------------
# Result comparison
# ----------------------------

# 1. Compute absolute and relative differences
abs_diff = torch.abs(result_cpu - result_cpu_1)
max_diff = abs_diff.max().item()
mean_diff = abs_diff.mean().item()

# 2. L2 (Euclidean) norm of error
l2_error = torch.norm(result_cpu - result_cpu_1) / torch.norm(result_cpu)

# 3. Element-wise equality check within tolerance
is_close = torch.allclose(result_cpu, result_cpu_1, rtol=1e-3, atol=1e-2)

print("\n=== Comparison Results ===")
print(f"Shapes match? {result_cpu.shape == result_cpu_1.shape}")
print(f"Max absolute difference: {max_diff:.6f}")
print(f"Mean absolute difference: {mean_diff:.6f}")
print(f"Relative L2 error: {l2_error:.6e}")
print(f"Allclose within tolerance? {is_close}")