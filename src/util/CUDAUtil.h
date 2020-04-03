// This file is in the public domain.

#ifndef LAB_CUDA_UTIL_H
#define LAB_CUDA_UTIL_H

#include <cstddef> /* std::size_t */
#include <iostream>

#include "cuda.h"

#include "Exception.h"

#define exec(C) do { cudaError_t err = C; if (err != cudaSuccess) \
	THROW_EXCEPTION(CUDAException, "CUDA error: " << cudaGetErrorString(err)); } while (false)
#define checkKernelLaunchError() do { cudaError_t err = cudaGetLastError(); if (err != cudaSuccess) \
	THROW_EXCEPTION(CUDAException, "Error in CUDA kernel launch: " << cudaGetErrorString(err)); } while (false)

namespace Lab {

struct CUDAException : std::runtime_error {
	using std::runtime_error::runtime_error;
};

template<typename T>
struct CUDAHostDevMem {
	T* hostPtr;
	T* devPtr;
	std::size_t sizeInBytes;

	CUDAHostDevMem() : hostPtr(), devPtr(), sizeInBytes() {}
	CUDAHostDevMem(std::size_t size) {
		sizeInBytes = size * sizeof(T);
		exec(cudaMallocHost(&hostPtr, sizeInBytes));
		try {
			exec(cudaMalloc(&devPtr, sizeInBytes));
		} catch (...) {
			cudaError_t err = cudaFreeHost(hostPtr);
			if (err != cudaSuccess) {
				std::cerr << "[CUDAHostDevMem] Error in cudaFreeHost: " << cudaGetErrorString(err) << std::endl;
			}
			throw;
		}
	}
	CUDAHostDevMem(const CUDAHostDevMem&) = delete;
	CUDAHostDevMem& operator=(const CUDAHostDevMem&) = delete;
	CUDAHostDevMem(CUDAHostDevMem&&) = delete;
	CUDAHostDevMem& operator=(CUDAHostDevMem&& o) {
		if (hostPtr) THROW_EXCEPTION(InvalidStateException, "The CUDAHostDevMem object is already initialized.");
		hostPtr = o.hostPtr;
		devPtr = o.devPtr;
		sizeInBytes = o.sizeInBytes;
		o.hostPtr = nullptr;
		o.devPtr = nullptr;
		o.sizeInBytes = 0;
		return *this;
	}
	~CUDAHostDevMem() {
		cudaError_t err = cudaFree(devPtr);
		if (err != cudaSuccess) {
			std::cerr << "[~CUDAHostDevMem] Error in cudaFree: " << cudaGetErrorString(err) << std::endl;
		}
		err = cudaFreeHost(hostPtr);
		if (err != cudaSuccess) {
			std::cerr << "[~CUDAHostDevMem] Error in cudaFreeHost: " << cudaGetErrorString(err) << std::endl;
		}
	}

	cudaError_t copyHostToDevice() {
		return cudaMemcpy(devPtr, hostPtr, sizeInBytes, cudaMemcpyHostToDevice);
	}
	cudaError_t copyHostToDeviceAsync() {
		return cudaMemcpyAsync(devPtr, hostPtr, sizeInBytes, cudaMemcpyHostToDevice);
	}
	cudaError_t copyDeviceToHost() {
		return cudaMemcpy(hostPtr, devPtr, sizeInBytes, cudaMemcpyDeviceToHost);
	}
	cudaError_t copyDeviceToHostAsync() {
		return cudaMemcpyAsync(hostPtr, devPtr, sizeInBytes, cudaMemcpyDeviceToHost);
	}
};

template<typename T>
struct CUDADevMem {
	T* devPtr;
	std::size_t sizeInBytes;

	CUDADevMem() : devPtr(), sizeInBytes() {}
	CUDADevMem(std::size_t size) {
		sizeInBytes = size * sizeof(T);
		exec(cudaMalloc(&devPtr, sizeInBytes));
	}
	CUDADevMem(const CUDADevMem&) = delete;
	CUDADevMem& operator=(const CUDADevMem&) = delete;
	CUDADevMem(CUDADevMem&&) = delete;
	CUDADevMem& operator=(CUDADevMem&& o) {
		if (devPtr) THROW_EXCEPTION(InvalidStateException, "The CUDADevMem object is already initialized.");
		devPtr = o.devPtr;
		sizeInBytes = o.sizeInBytes;
		o.devPtr = nullptr;
		o.sizeInBytes = 0;
		return *this;
	}
	~CUDADevMem() {
		cudaError_t err = cudaFree(devPtr);
		if (err != cudaSuccess) {
			std::cerr << "[~CUDADevMem] Error in cudaFree: " << cudaGetErrorString(err) << std::endl;
		}
	}
};

namespace CUDAUtil {

inline
std::size_t
numberOfBlocks(std::size_t n, std::size_t blockSize) {
	return (n + (blockSize - 1U)) / blockSize;
}

inline
std::size_t
roundUpToMultipleOfBlockSize(std::size_t n, std::size_t blockSize) {
	return numberOfBlocks(n, blockSize) * blockSize;
}

} // namespace CUDAUtil
} // namespace Lab

#endif // LAB_CUDA_UTIL_H
