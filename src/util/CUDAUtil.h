// This file is in the public domain.

#ifndef LAB_CUDA_UTIL_H
#define LAB_CUDA_UTIL_H

#include "Exception.h"

#define exec(C) do { cudaError_t err = C; if (err != cudaSuccess) \
	THROW_EXCEPTION(CUDAException, "CUDA error: " << cudaGetErrorString(err)); } while (false)
#define checkKernelLaunchError() do { cudaError_t err = cudaGetLastError(); if (err != cudaSuccess) \
	THROW_EXCEPTION(CUDAException, "Error in CUDA kernel launch: " << cudaGetErrorString(err)); } while (false)

namespace Lab {

struct CUDAException : std::runtime_error {
	using std::runtime_error::runtime_error;
};

} // namespace Lab

#endif // LAB_CUDA_UTIL_H
