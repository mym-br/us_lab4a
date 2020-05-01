/***************************************************************************
 *  Copyright 2018, 2020 Marcelo Y. Matuda                                 *
 *                                                                         *
 *  This program is free software: you can redistribute it and/or modify   *
 *  it under the terms of the GNU General Public License as published by   *
 *  the Free Software Foundation, either version 3 of the License, or      *
 *  (at your option) any later version.                                    *
 *                                                                         *
 *  This program is distributed in the hope that it will be useful,        *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of         *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          *
 *  GNU General Public License for more details.                           *
 *                                                                         *
 *  You should have received a copy of the GNU General Public License      *
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.  *
 ***************************************************************************/
#ifndef SIMTRANSIENTACOUSTICFIELD_H
#define SIMTRANSIENTACOUSTICFIELD_H

#include <complex>
#include <condition_variable>
#include <mutex>
#include <thread>
#include <utility> /* pair, swap */
#include <vector>

#include <tbb/concurrent_queue.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/tbb.h>

#include "ArrayOfRectangularSourcesImpulseResponse.h"
#include "FFTWFilter2.h"
#include "IterationCounter.h"
#include "Matrix.h"
#include "Util.h"
#include "XYZValue.h"



namespace Lab {

template<typename TFloat>
class SimTransientAcousticField {
public:
	template<typename ImpulseResponse>
	static void getCircularSourceAcousticField(
			TFloat samplingFreq,
			TFloat propagationSpeed,
			TFloat sourceRadius,
			TFloat discretization,
			const std::vector<TFloat>& dvdt,
			Matrix<XYZValue<TFloat>>& gridData);

	template<typename ImpulseResponse>
	static void getCircularSourceAcousticFieldSingleThread(
			TFloat samplingFreq,
			TFloat propagationSpeed,
			TFloat sourceRadius,
			TFloat discretization,
			const std::vector<TFloat>& dvdt,
			Matrix<XYZValue<TFloat>>& gridData);

	template<typename ImpulseResponse>
	static void getRectangularSourceAcousticField(
			TFloat samplingFreq,
			TFloat propagationSpeed,
			TFloat sourceWidth,
			TFloat sourceHeight,
			TFloat discretization,
			const std::vector<TFloat>& dvdt,
			Matrix<XYZValue<TFloat>>& gridData);

	template<typename ImpulseResponse>
	static void getRectangularSourceAcousticFieldSingleThread(
			TFloat samplingFreq,
			TFloat propagationSpeed,
			TFloat sourceWidth,
			TFloat sourceHeight,
			TFloat discretization,
			const std::vector<TFloat>& dvdt,
			Matrix<XYZValue<TFloat>>& gridData);

	template<typename STImpulseResponse, typename MTImpulseResponse>
	static void getRectangularSourceAcousticFieldSTMT(
			TFloat samplingFreq,
			TFloat propagationSpeed,
			TFloat sourceWidth,
			TFloat sourceHeight,
			TFloat discretization,
			const std::vector<TFloat>& dvdt,
			Matrix<XYZValue<TFloat>>& gridData);

	template<typename ImpulseResponse>
	static void getArrayOfRectangularSourcesAcousticField(
			TFloat samplingFreq,
			TFloat propagationSpeed,
			TFloat sourceWidth,
			TFloat sourceHeight,
			TFloat discretization,
			const std::vector<TFloat>& dvdt,
			const std::vector<XY<TFloat>>& elemPos,
			const std::vector<TFloat>& focusDelay /* s */,
			Matrix<XYZValue<TFloat>>& gridData);

	template<typename ImpulseResponse>
	static void getArrayOfRectangularSourcesAcousticFieldDirect(
			TFloat samplingFreq,
			TFloat propagationSpeed,
			TFloat sourceWidth,
			TFloat sourceHeight,
			TFloat discretization,
			const std::vector<TFloat>& dvdt,
			const std::vector<XY<TFloat>>& elemPos,
			const std::vector<TFloat>& focusDelay /* s */,
			Matrix<XYZValue<TFloat>>& gridData);

	template<typename ImpulseResponse>
	static void getArrayOfRectangularSourcesAcousticFieldDirectSingleThread(
			TFloat samplingFreq,
			TFloat propagationSpeed,
			TFloat sourceWidth,
			TFloat sourceHeight,
			TFloat discretization,
			const std::vector<TFloat>& dvdt,
			const std::vector<XY<TFloat>>& elemPos,
			const std::vector<TFloat>& focusDelay /* s */,
			Matrix<XYZValue<TFloat>>& gridData);

	template<typename STImpulseResponse, typename MTImpulseResponse>
	static void getArrayOfRectangularSourcesAcousticFieldDirectSTMT(
			TFloat samplingFreq,
			TFloat propagationSpeed,
			TFloat sourceWidth,
			TFloat sourceHeight,
			TFloat discretization,
			const std::vector<TFloat>& dvdt,
			const std::vector<XY<TFloat>>& elemPos,
			const std::vector<TFloat>& focusDelay /* s */,
			Matrix<XYZValue<TFloat>>& gridData);

private:
	template<typename ImpulseResponse>
	struct CircularSourceThreadData {
		CircularSourceThreadData(
			TFloat samplingFreq,
			TFloat propagationSpeed,
			TFloat sourceRadius,
			TFloat discretization,
			const std::vector<TFloat>& dvdt)
				: ir(samplingFreq, propagationSpeed, sourceRadius, discretization)
		{
			filter.setCoefficients(dvdt, filterFreqCoeff);
		}
		ImpulseResponse ir;
		std::vector<std::complex<TFloat>> filterFreqCoeff;
		std::vector<TFloat> h;
		std::vector<TFloat> signal;
		FFTWFilter2<TFloat> filter;
	};

	template<typename ImpulseResponse>
	struct RectangularSourceThreadData {
		RectangularSourceThreadData(
			TFloat samplingFreq,
			TFloat propagationSpeed,
			TFloat sourceWidth,
			TFloat sourceHeight,
			TFloat discretization,
			const std::vector<TFloat>& dvdt)
				: ir(samplingFreq, propagationSpeed, sourceWidth, sourceHeight, discretization)
		{
			filter.setCoefficients(dvdt, filterFreqCoeff);
		}
		ImpulseResponse ir;
		std::vector<std::complex<TFloat>> filterFreqCoeff;
		std::vector<TFloat> h;
		std::vector<TFloat> signal;
		FFTWFilter2<TFloat> filter;
	};

	template<typename ImpulseResponse>
	struct ArrayOfRectangularSourcesThreadData {
		ArrayOfRectangularSourcesThreadData(
			TFloat samplingFreq,
			TFloat propagationSpeed,
			TFloat sourceWidth,
			TFloat sourceHeight,
			TFloat discretization,
			const std::vector<XY<TFloat>>& elemPos,
			const std::vector<TFloat>& focusDelay,
			const std::vector<TFloat>& dvdt)
				: ir(samplingFreq, propagationSpeed, sourceWidth, sourceHeight, discretization,
					elemPos, focusDelay)
		{
			filter.setCoefficients(dvdt, filterFreqCoeff);
		}
		ArrayOfRectangularSourcesImpulseResponse<TFloat, ImpulseResponse> ir;
		std::vector<std::complex<TFloat>> filterFreqCoeff;
		std::vector<TFloat> h;
		std::vector<TFloat> signal;
		FFTWFilter2<TFloat> filter;
	};

	template<typename ImpulseResponse>
	struct DirectArrayOfRectangularSourcesThreadData {
		DirectArrayOfRectangularSourcesThreadData(
			TFloat samplingFreq,
			TFloat propagationSpeed,
			TFloat sourceWidth,
			TFloat sourceHeight,
			TFloat discretization,
			const std::vector<XY<TFloat>>& elemPos,
			const std::vector<TFloat>& focusDelay,
			const std::vector<TFloat>& dvdt)
				: ir(samplingFreq, propagationSpeed, sourceWidth, sourceHeight, discretization,
					elemPos, focusDelay)
		{
			filter.setCoefficients(dvdt, filterFreqCoeff);
		}
		ImpulseResponse ir;
		std::vector<std::complex<TFloat>> filterFreqCoeff;
		std::vector<TFloat> h;
		std::vector<TFloat> signal;
		FFTWFilter2<TFloat> filter;
	};

	SimTransientAcousticField() = delete;

	template<typename T>
		static void exec(T& tls, Matrix<XYZValue<TFloat>>& gridData);
	template<typename T>
		static void execSingleThread(T& threadData, Matrix<XYZValue<TFloat>>& gridData);
	template<typename T, typename U>
		static void execSingleThreadMultiThread(T& stThreadData, U& tls, Matrix<XYZValue<TFloat>>& gridData);
};

template<typename TFloat>
template<typename T>
void
SimTransientAcousticField<TFloat>::exec(T& tls, Matrix<XYZValue<TFloat>>& gridData)
{
	IterationCounter::reset(gridData.n1());
	for (std::size_t i = 0, iEnd = gridData.n1(); i < iEnd; ++i) {
		tbb::parallel_for(tbb::blocked_range<std::size_t>(0, gridData.n2()),
			[&, i](const tbb::blocked_range<std::size_t>& r) {
				auto& local = tls.local();
				std::size_t hOffset;
				for (std::size_t j = r.begin(); j != r.end(); ++j) {
					XYZValue<TFloat>& point = gridData(i, j);
					local.ir.getImpulseResponse(point.x, point.y, point.z, hOffset, local.h);

					local.filter.filter(local.filterFreqCoeff, local.h, local.signal);

					TFloat minValue, maxValue;
					Util::minMax(local.signal, minValue, maxValue);
					point.value = maxValue - minValue;
				}
			}
		);
		IterationCounter::add(1);
	}
}

template<typename TFloat>
template<typename T>
void
SimTransientAcousticField<TFloat>::execSingleThread(T& threadData, Matrix<XYZValue<TFloat>>& gridData)
{
	IterationCounter::reset(gridData.n1());
	for (std::size_t i = 0, iEnd = gridData.n1(); i < iEnd; ++i) {
		std::size_t hOffset;
		for (std::size_t j = 0, jEnd = gridData.n2(); j < jEnd; ++j) {
			XYZValue<TFloat>& point = gridData(i, j);
			threadData.ir.getImpulseResponse(point.x, point.y, point.z, hOffset, threadData.h);

			threadData.filter.filter(threadData.filterFreqCoeff, threadData.h, threadData.signal);

			TFloat minValue, maxValue;
			Util::minMax(threadData.signal, minValue, maxValue);
			point.value = maxValue - minValue;
		}
		IterationCounter::add(1);
	}
}

template<typename TFloat>
template<typename T, typename U>
void
SimTransientAcousticField<TFloat>::execSingleThreadMultiThread(T& stThreadData, U& tls, Matrix<XYZValue<TFloat>>& gridData)
{
	const unsigned int totalThreads = std::thread::hardware_concurrency();
	const unsigned int numThreads = (totalThreads > 2) ? totalThreads - 2 : 0;

	IterationCounter::reset(gridData.n1() * gridData.n2());

	// Fill task queue.
	tbb::concurrent_queue<std::pair<unsigned int, unsigned int>> queue;
	for (std::size_t i = 0, iEnd = gridData.n1(); i < iEnd; ++i) {
		for (std::size_t j = 0, jEnd = gridData.n2(); j != jEnd; ++j) {
			queue.push(std::make_pair(j, i));
		}
	}

	// Create threads for the multithreaded impulse response.
	std::vector<std::thread> threadList(numThreads);
	for (unsigned int t = 0; t < numThreads; ++t) {
		threadList[t] = std::thread([&] () {
			auto& local = tls.local();
			std::pair<unsigned int, unsigned int> item;
			while (queue.try_pop(item)) {
				const auto [j, i] = item;

				XYZValue<TFloat>& point = gridData(i, j);
				std::size_t hOffset;
				local.ir.getImpulseResponse(point.x, point.y, point.z, hOffset, local.h);

				local.filter.filter(local.filterFreqCoeff, local.h, local.signal);

				TFloat minValue, maxValue;
				Util::minMax(local.signal, minValue, maxValue);
				point.value = maxValue - minValue;

				IterationCounter::add(1);
			}
		});
	}

	std::mutex m;
	std::condition_variable workAvailableOrExitingCond, workDoneCond;
	std::pair<unsigned int, unsigned int> procItem;
	std::vector<TFloat> hProc;
	enum class State {
		IDLE,
		PROCESSING,
		EXITING
	};
	State state = State::IDLE;

	// Create the thread that processes the result of the single-thread impulse response.
	std::thread procThread([&] () {
		while (true) {
			{
				std::unique_lock<std::mutex> locker(m);
				while (state == State::IDLE) {
					workAvailableOrExitingCond.wait(locker);
				}
				if (state == State::EXITING) return;
			}

			const auto [j, i] = procItem;
			XYZValue<TFloat>& point = gridData(i, j);

			stThreadData.filter.filter(stThreadData.filterFreqCoeff, hProc, stThreadData.signal);

			TFloat minValue, maxValue;
			Util::minMax(stThreadData.signal, minValue, maxValue);
			point.value = maxValue - minValue;

			{
				std::lock_guard<std::mutex> locker(m);
				state = State::IDLE;
			}
			workDoneCond.notify_all();

			IterationCounter::add(1);
		}
	});

	// Use the single-thread impulse response.
	std::pair<unsigned int, unsigned int> item;
	while (queue.try_pop(item)) {
		const auto [j, i] = item;

		const XYZValue<TFloat>& point = gridData(i, j);
		std::size_t hOffset;
		stThreadData.ir.getImpulseResponse(point.x, point.y, point.z, hOffset, stThreadData.h);

		{
			std::unique_lock<std::mutex> locker(m);
			while (state == State::PROCESSING) {
				workDoneCond.wait(locker);
			}
			std::swap(hProc, stThreadData.h);
			procItem = item;
			state = State::PROCESSING;
		}
		workAvailableOrExitingCond.notify_all();
	}

	{
		std::unique_lock<std::mutex> locker(m);
		while (state == State::PROCESSING) {
			workDoneCond.wait(locker);
		}
		state = State::EXITING;
	}
	workAvailableOrExitingCond.notify_all();

	for (unsigned int t = 0; t < numThreads; ++t) {
		threadList[t].join();
	}
	procThread.join();
}

template<typename TFloat>
template<typename ImpulseResponse>
void
SimTransientAcousticField<TFloat>::getCircularSourceAcousticField(
					TFloat samplingFreq,
					TFloat propagationSpeed,
					TFloat sourceRadius,
					TFloat discretization,
					const std::vector<TFloat>& dvdt,
					Matrix<XYZValue<TFloat>>& gridData)
{
	CircularSourceThreadData<ImpulseResponse> threadData{
		samplingFreq,
		propagationSpeed,
		sourceRadius,
		discretization,
		dvdt
	};
	tbb::enumerable_thread_specific<CircularSourceThreadData<ImpulseResponse>> tls(threadData);

	exec(tls, gridData);
}

template<typename TFloat>
template<typename ImpulseResponse>
void
SimTransientAcousticField<TFloat>::getCircularSourceAcousticFieldSingleThread(
					TFloat samplingFreq,
					TFloat propagationSpeed,
					TFloat sourceRadius,
					TFloat discretization,
					const std::vector<TFloat>& dvdt,
					Matrix<XYZValue<TFloat>>& gridData)
{
	CircularSourceThreadData<ImpulseResponse> threadData{
		samplingFreq,
		propagationSpeed,
		sourceRadius,
		discretization,
		dvdt
	};

	execSingleThread(threadData, gridData);
}

template<typename TFloat>
template<typename ImpulseResponse>
void
SimTransientAcousticField<TFloat>::getRectangularSourceAcousticField(
					TFloat samplingFreq,
					TFloat propagationSpeed,
					TFloat sourceWidth,
					TFloat sourceHeight,
					TFloat discretization,
					const std::vector<TFloat>& dvdt,
					Matrix<XYZValue<TFloat>>& gridData)
{
	RectangularSourceThreadData<ImpulseResponse> threadData{
		samplingFreq,
		propagationSpeed,
		sourceWidth,
		sourceHeight,
		discretization,
		dvdt
	};
	tbb::enumerable_thread_specific<RectangularSourceThreadData<ImpulseResponse>> tls(threadData);

	exec(tls, gridData);
}

template<typename TFloat>
template<typename ImpulseResponse>
void
SimTransientAcousticField<TFloat>::getRectangularSourceAcousticFieldSingleThread(
					TFloat samplingFreq,
					TFloat propagationSpeed,
					TFloat sourceWidth,
					TFloat sourceHeight,
					TFloat discretization,
					const std::vector<TFloat>& dvdt,
					Matrix<XYZValue<TFloat>>& gridData)
{
	RectangularSourceThreadData<ImpulseResponse> threadData{
		samplingFreq,
		propagationSpeed,
		sourceWidth,
		sourceHeight,
		discretization,
		dvdt
	};

	execSingleThread(threadData, gridData);
}

template<typename TFloat>
template<typename STImpulseResponse, typename MTImpulseResponse>
void
SimTransientAcousticField<TFloat>::getRectangularSourceAcousticFieldSTMT(
					TFloat samplingFreq,
					TFloat propagationSpeed,
					TFloat sourceWidth,
					TFloat sourceHeight,
					TFloat discretization,
					const std::vector<TFloat>& dvdt,
					Matrix<XYZValue<TFloat>>& gridData)
{
	RectangularSourceThreadData<STImpulseResponse> stThreadData{
		samplingFreq,
		propagationSpeed,
		sourceWidth,
		sourceHeight,
		discretization,
		dvdt
	};
	RectangularSourceThreadData<MTImpulseResponse> mtThreadData{
		samplingFreq,
		propagationSpeed,
		sourceWidth,
		sourceHeight,
		discretization,
		dvdt
	};
	tbb::enumerable_thread_specific<RectangularSourceThreadData<MTImpulseResponse>> tls(mtThreadData);

	execSingleThreadMultiThread(stThreadData, tls, gridData);
}

template<typename TFloat>
template<typename ImpulseResponse>
void
SimTransientAcousticField<TFloat>::getArrayOfRectangularSourcesAcousticField(
					TFloat samplingFreq,
					TFloat propagationSpeed,
					TFloat sourceWidth,
					TFloat sourceHeight,
					TFloat discretization,
					const std::vector<TFloat>& dvdt,
					const std::vector<XY<TFloat>>& elemPos,
					const std::vector<TFloat>& focusDelay,
					Matrix<XYZValue<TFloat>>& gridData)
{
	ArrayOfRectangularSourcesThreadData<ImpulseResponse> threadData{
		samplingFreq,
		propagationSpeed,
		sourceWidth,
		sourceHeight,
		discretization,
		elemPos,
		focusDelay,
		dvdt
	};
	tbb::enumerable_thread_specific<ArrayOfRectangularSourcesThreadData<ImpulseResponse>> tls(threadData);

	exec(tls, gridData);
}

template<typename TFloat>
template<typename ImpulseResponse>
void
SimTransientAcousticField<TFloat>::getArrayOfRectangularSourcesAcousticFieldDirect(
					TFloat samplingFreq,
					TFloat propagationSpeed,
					TFloat sourceWidth,
					TFloat sourceHeight,
					TFloat discretization,
					const std::vector<TFloat>& dvdt,
					const std::vector<XY<TFloat>>& elemPos,
					const std::vector<TFloat>& focusDelay,
					Matrix<XYZValue<TFloat>>& gridData)
{
	DirectArrayOfRectangularSourcesThreadData<ImpulseResponse> threadData{
		samplingFreq,
		propagationSpeed,
		sourceWidth,
		sourceHeight,
		discretization,
		elemPos,
		focusDelay,
		dvdt
	};
	tbb::enumerable_thread_specific<DirectArrayOfRectangularSourcesThreadData<ImpulseResponse>> tls(threadData);

	exec(tls, gridData);
}

template<typename TFloat>
template<typename ImpulseResponse>
void
SimTransientAcousticField<TFloat>::getArrayOfRectangularSourcesAcousticFieldDirectSingleThread(
					TFloat samplingFreq,
					TFloat propagationSpeed,
					TFloat sourceWidth,
					TFloat sourceHeight,
					TFloat discretization,
					const std::vector<TFloat>& dvdt,
					const std::vector<XY<TFloat>>& elemPos,
					const std::vector<TFloat>& focusDelay,
					Matrix<XYZValue<TFloat>>& gridData)
{
	DirectArrayOfRectangularSourcesThreadData<ImpulseResponse> threadData{
		samplingFreq,
		propagationSpeed,
		sourceWidth,
		sourceHeight,
		discretization,
		elemPos,
		focusDelay,
		dvdt
	};

	execSingleThread(threadData, gridData);
}

template<typename TFloat>
template<typename STImpulseResponse, typename MTImpulseResponse>
void
SimTransientAcousticField<TFloat>::getArrayOfRectangularSourcesAcousticFieldDirectSTMT(
					TFloat samplingFreq,
					TFloat propagationSpeed,
					TFloat sourceWidth,
					TFloat sourceHeight,
					TFloat discretization,
					const std::vector<TFloat>& dvdt,
					const std::vector<XY<TFloat>>& elemPos,
					const std::vector<TFloat>& focusDelay,
					Matrix<XYZValue<TFloat>>& gridData)
{
	DirectArrayOfRectangularSourcesThreadData<STImpulseResponse> stThreadData{
		samplingFreq,
		propagationSpeed,
		sourceWidth,
		sourceHeight,
		discretization,
		elemPos,
		focusDelay,
		dvdt
	};
	DirectArrayOfRectangularSourcesThreadData<MTImpulseResponse> mtThreadData{
		samplingFreq,
		propagationSpeed,
		sourceWidth,
		sourceHeight,
		discretization,
		elemPos,
		focusDelay,
		dvdt
	};
	tbb::enumerable_thread_specific<DirectArrayOfRectangularSourcesThreadData<MTImpulseResponse>> tls(mtThreadData);

	execSingleThreadMultiThread(stThreadData, tls, gridData);
}

} // namespace Lab

#endif // SIMTRANSIENTACOUSTICFIELD_H
