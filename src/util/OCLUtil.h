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
						hostBuffer, CL_BLOCKING, mapFlag,
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

#define OCL_ERROR_CASE(A) case A: return #A;

const char*
getErrorString(cl_int error)
{
	switch (error) {
	OCL_ERROR_CASE(CL_SUCCESS)
	OCL_ERROR_CASE(CL_DEVICE_NOT_FOUND)
	OCL_ERROR_CASE(CL_DEVICE_NOT_AVAILABLE)
	OCL_ERROR_CASE(CL_COMPILER_NOT_AVAILABLE)
	OCL_ERROR_CASE(CL_MEM_OBJECT_ALLOCATION_FAILURE)
	OCL_ERROR_CASE(CL_OUT_OF_RESOURCES)
	OCL_ERROR_CASE(CL_OUT_OF_HOST_MEMORY)
	OCL_ERROR_CASE(CL_PROFILING_INFO_NOT_AVAILABLE)
	OCL_ERROR_CASE(CL_MEM_COPY_OVERLAP)
	OCL_ERROR_CASE(CL_IMAGE_FORMAT_MISMATCH)
	OCL_ERROR_CASE(CL_IMAGE_FORMAT_NOT_SUPPORTED)
	OCL_ERROR_CASE(CL_BUILD_PROGRAM_FAILURE)
	OCL_ERROR_CASE(CL_MAP_FAILURE)
	OCL_ERROR_CASE(CL_MISALIGNED_SUB_BUFFER_OFFSET)
	OCL_ERROR_CASE(CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST)
	OCL_ERROR_CASE(CL_COMPILE_PROGRAM_FAILURE)
	OCL_ERROR_CASE(CL_LINKER_NOT_AVAILABLE)
	OCL_ERROR_CASE(CL_LINK_PROGRAM_FAILURE)
	OCL_ERROR_CASE(CL_DEVICE_PARTITION_FAILED)
	OCL_ERROR_CASE(CL_KERNEL_ARG_INFO_NOT_AVAILABLE)
	OCL_ERROR_CASE(CL_INVALID_VALUE)
	OCL_ERROR_CASE(CL_INVALID_DEVICE_TYPE)
	OCL_ERROR_CASE(CL_INVALID_PLATFORM)
	OCL_ERROR_CASE(CL_INVALID_DEVICE)
	OCL_ERROR_CASE(CL_INVALID_CONTEXT)
	OCL_ERROR_CASE(CL_INVALID_QUEUE_PROPERTIES)
	OCL_ERROR_CASE(CL_INVALID_COMMAND_QUEUE)
	OCL_ERROR_CASE(CL_INVALID_HOST_PTR)
	OCL_ERROR_CASE(CL_INVALID_MEM_OBJECT)
	OCL_ERROR_CASE(CL_INVALID_IMAGE_FORMAT_DESCRIPTOR)
	OCL_ERROR_CASE(CL_INVALID_IMAGE_SIZE)
	OCL_ERROR_CASE(CL_INVALID_SAMPLER)
	OCL_ERROR_CASE(CL_INVALID_BINARY)
	OCL_ERROR_CASE(CL_INVALID_BUILD_OPTIONS)
	OCL_ERROR_CASE(CL_INVALID_PROGRAM)
	OCL_ERROR_CASE(CL_INVALID_PROGRAM_EXECUTABLE)
	OCL_ERROR_CASE(CL_INVALID_KERNEL_NAME)
	OCL_ERROR_CASE(CL_INVALID_KERNEL_DEFINITION)
	OCL_ERROR_CASE(CL_INVALID_KERNEL)
	OCL_ERROR_CASE(CL_INVALID_ARG_INDEX)
	OCL_ERROR_CASE(CL_INVALID_ARG_VALUE)
	OCL_ERROR_CASE(CL_INVALID_ARG_SIZE)
	OCL_ERROR_CASE(CL_INVALID_KERNEL_ARGS)
	OCL_ERROR_CASE(CL_INVALID_WORK_DIMENSION)
	OCL_ERROR_CASE(CL_INVALID_WORK_GROUP_SIZE)
	OCL_ERROR_CASE(CL_INVALID_WORK_ITEM_SIZE)
	OCL_ERROR_CASE(CL_INVALID_GLOBAL_OFFSET)
	OCL_ERROR_CASE(CL_INVALID_EVENT_WAIT_LIST)
	OCL_ERROR_CASE(CL_INVALID_EVENT)
	OCL_ERROR_CASE(CL_INVALID_OPERATION)
	OCL_ERROR_CASE(CL_INVALID_GL_OBJECT)
	OCL_ERROR_CASE(CL_INVALID_BUFFER_SIZE)
	OCL_ERROR_CASE(CL_INVALID_MIP_LEVEL)
	OCL_ERROR_CASE(CL_INVALID_GLOBAL_WORK_SIZE)
	OCL_ERROR_CASE(CL_INVALID_PROPERTY)
	OCL_ERROR_CASE(CL_INVALID_IMAGE_DESCRIPTOR)
	OCL_ERROR_CASE(CL_INVALID_COMPILER_OPTIONS)
	OCL_ERROR_CASE(CL_INVALID_LINKER_OPTIONS)
	OCL_ERROR_CASE(CL_INVALID_DEVICE_PARTITION_COUNT)
	OCL_ERROR_CASE(CL_PLATFORM_NOT_FOUND_KHR)
	default:
		return "Unknown OpenCL error.";
	}
}

} // namespace OCLUtil

std::ostream&
operator<<(std::ostream& out, const cl::Error& e)
{
	out << "OpenCL error " << e.err() << " in " << e.what() << ": " << OCLUtil::getErrorString(e.err()) << '.';
	return out;
}

} // namespace Lab

#endif // LAB_OCL_UTIL_H
