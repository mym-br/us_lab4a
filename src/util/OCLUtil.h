// This file is in the public domain.

#ifndef LAB_OCL_UTIL_H
#define LAB_OCL_UTIL_H

#include <cstddef> /* std::size_t */
#include <iostream>
#include <utility> /* move */

#include <CL/cl2.hpp>

#include "Exception.h"
#include "Log.h"



namespace Lab {

struct OCLException : std::runtime_error {
	using std::runtime_error::runtime_error;
};

template<typename T>
struct OCLPinnedHostMem {
	cl::CommandQueue& queue;
	std::size_t sizeInBytes;
	T* hostPtr;

	// mapFlag:
	//   CL_MAP_READ
	//   CL_MAP_WRITE
	//   CL_MAP_WRITE_INVALIDATE_REGION
	OCLPinnedHostMem(cl::Context& context, cl::CommandQueue& q, std::size_t size, int mapFlag)
			: queue(q)
			, sizeInBytes(size * sizeof(T)) {
		hostBuffer = cl::Buffer(context, CL_MEM_ALLOC_HOST_PTR, sizeInBytes);
		hostPtr = static_cast<T*>(queue.enqueueMapBuffer(
						hostBuffer, CL_TRUE /* blocking */, mapFlag,
						0 /* offset */, sizeInBytes));
	}
	OCLPinnedHostMem(const OCLPinnedHostMem&) = delete;
	OCLPinnedHostMem& operator=(const OCLPinnedHostMem&) = delete;
	OCLPinnedHostMem(OCLPinnedHostMem&& o)
			: queue(o.queue)
			, sizeInBytes(o.sizeInBytes)
			, hostPtr(o.hostPtr)
			, hostBuffer(std::move(o.hostBuffer)) {
		o.hostPtr = nullptr;
		o.sizeInBytes = 0;
	}
	OCLPinnedHostMem& operator=(OCLPinnedHostMem&& o) {
		clear();

		queue       = o.queue;
		sizeInBytes = o.sizeInBytes;
		hostPtr     = o.hostPtr;
		hostBuffer  = std::move(o.hostBuffer);

		o.hostPtr     = nullptr;
		o.sizeInBytes = 0;
		return *this;
	}
	~OCLPinnedHostMem() {
		try {
			clear();
		} catch (cl::Error& e) {
			std::cerr << "[~OCLPinnedHostMem] OpenCL error: " << e.what() << " (" << e.err() << ")." << std::endl;
		}
	}
	void clear() {
		if (!hostPtr) return;

		cl::Event event;
		queue.enqueueUnmapMemObject(hostBuffer, hostPtr, nullptr, &event);
		event.wait();

		hostBuffer = cl::Buffer();
		hostPtr = nullptr;
		sizeInBytes = 0;
	}
private:
	cl::Buffer hostBuffer;
};

namespace OCLUtil {

inline
cl::Context
initOpenCL()
{
	std::vector<cl::Platform> platforms;
	cl::Platform::get(&platforms);
	if (platforms.empty()) {
		THROW_EXCEPTION(UnavailableResourceException, "No OpenCL platforms available.");
	}

	if (Log::isDebugEnabled()) {
		for (auto& plat : platforms) {
			std::string name = plat.getInfo<CL_PLATFORM_NAME>();
			LOG_DEBUG << "OpenCL platform: " << name;

			std::string version = plat.getInfo<CL_PLATFORM_VERSION>();
			LOG_DEBUG << "OpenCL version: " << version;

			std::vector<cl::Device> devices;
			plat.getDevices(CL_DEVICE_TYPE_ALL, &devices);
			if (devices.empty()) {
				THROW_EXCEPTION(UnavailableResourceException, "No OpenCL devices available for platform " << name << '.');
			}

			for (auto& dev : devices) {
				cl_device_type deviceType = dev.getInfo<CL_DEVICE_TYPE>();
				std::string devName = dev.getInfo<CL_DEVICE_NAME>();
				LOG_DEBUG << "  device name: " << devName;
				switch (deviceType) {
				case CL_DEVICE_TYPE_CPU:
					LOG_DEBUG << "    type: CPU";
					break;
				case CL_DEVICE_TYPE_GPU:
					LOG_DEBUG << "    type: GPU";
					break;
				default:
					LOG_DEBUG << "    type: other (" << deviceType << ").";
				}
			}
		}
	}

	// Select platform.
	if (OPENCL_PLATFORM >= platforms.size()) {
		THROW_EXCEPTION(UnavailableResourceException, "Invalid OpenCL platform: " << OPENCL_PLATFORM << '.');
	}
	cl::Platform chosenPlatform = platforms[OPENCL_PLATFORM];

	// Select device.
	std::vector<cl::Device> devices;
	chosenPlatform.getDevices(CL_DEVICE_TYPE_ALL, &devices);
	if (devices.empty()) {
		THROW_EXCEPTION(UnavailableResourceException, "No OpenCL devices available.");
	}
	if (OPENCL_DEVICE >= devices.size()) {
		THROW_EXCEPTION(UnavailableResourceException, "Invalid OpenCL device: " << OPENCL_DEVICE << '.');
	}
	cl::Device chosenDevice = devices[OPENCL_DEVICE];

	// Create context.
	cl_context_properties contextProp[] = {CL_CONTEXT_PLATFORM, (cl_context_properties) chosenPlatform(), 0 /* end of list */};
	return cl::Context(chosenDevice, contextProp);
}

inline
std::size_t
numberOfGroups(std::size_t n, std::size_t groupSize) {
	return (n + (groupSize - 1U)) / groupSize;
}

inline
std::size_t
roundUpToMultipleOfGroupSize(std::size_t n, std::size_t groupSize) {
	return numberOfGroups(n, groupSize) * groupSize;
}

} // namespace OCLUtil
} // namespace Lab

#endif // LAB_OCL_UTIL_H
