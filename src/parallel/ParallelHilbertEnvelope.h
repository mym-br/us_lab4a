/***************************************************************************
 *  Copyright 2014, 2017, 2018 Marcelo Y. Matuda                           *
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

#ifndef PARALLELHILBERTENVELOPE_H
#define PARALLELHILBERTENVELOPE_H

#include <cstddef> /* std::size_t */

#include <tbb/blocked_range.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>

#include "HilbertEnvelope.h"
#include "Matrix.h"
#include "Util.h"



namespace Lab {

template<typename TFloat>
class ParallelHilbertEnvelope {
public:
	template<typename T>
		class CalculateDim2;
	template<typename T, typename U>
		class GetAnalyticSignalDim2;

	template<typename T>
		static void calculateDim2(Matrix<T>& data);
	template<typename T>
		static void calculateDim2Value(Matrix<T>& data);
	template<typename T, typename U>
		static void getAnalyticSignalDim2(const Matrix<T>& origData, Matrix<U>& destData);
private:
	ParallelHilbertEnvelope() = delete;
};



template<typename TFloat>
template<typename T>
class ParallelHilbertEnvelope<TFloat>::CalculateDim2 {
public:
	CalculateDim2(tbb::enumerable_thread_specific<HilbertEnvelope<TFloat>>& envelopeTLS, Matrix<T>& data)
			: envelopeTLS_(envelopeTLS)
			, data_(data) { }
	void operator()(const tbb::blocked_range<std::size_t>& r) const {
		typename tbb::enumerable_thread_specific<HilbertEnvelope<TFloat>>::reference envelope = envelopeTLS_.local();
		for (std::size_t i = r.begin(); i != r.end(); ++i) {
			envelope.calculate(&data_(i, 0), data_.n2());
		}
	}
private:
	tbb::enumerable_thread_specific<HilbertEnvelope<TFloat>>& envelopeTLS_;
	Matrix<T>& data_;
};

template<typename TFloat>
template<typename T>
void
ParallelHilbertEnvelope<TFloat>::calculateDim2(Matrix<T>& data)
{
	tbb::enumerable_thread_specific<HilbertEnvelope<TFloat>> envelopeTLS;
	tbb::parallel_for(
		tbb::blocked_range<std::size_t>(0, data.n1()),
		CalculateDim2<typename Matrix<T>::ValueType>(envelopeTLS, data)
	);
}

template<typename TFloat>
template<typename T>
void
ParallelHilbertEnvelope<TFloat>::calculateDim2Value(Matrix<T>& data)
{
	Matrix<TFloat> aux;
	Util::copyValueToSimpleMatrix(data, aux);

	tbb::enumerable_thread_specific<HilbertEnvelope<TFloat>> envelopeTLS;
	tbb::parallel_for(
		tbb::blocked_range<std::size_t>(0, data.n1()),
		CalculateDim2<TFloat>(envelopeTLS, aux)
	);

	Util::copyValueFromSimpleMatrix(aux, data);
}

template<typename TFloat>
template<typename T, typename U>
class ParallelHilbertEnvelope<TFloat>::GetAnalyticSignalDim2 {
public:
	GetAnalyticSignalDim2(
			tbb::enumerable_thread_specific<HilbertEnvelope<TFloat>>& envelopeTLS,
			const Matrix<T>& origData,
			Matrix<U>& destData)
		: envelopeTLS_(envelopeTLS)
		, origData_(origData)
		, destData_(destData) { }
	void operator()(const tbb::blocked_range<std::size_t>& r) const {
		typename tbb::enumerable_thread_specific<HilbertEnvelope<TFloat>>::reference envelope = envelopeTLS_.local();
		for (std::size_t i = r.begin(); i != r.end(); ++i) {
			envelope.getAnalyticSignal(&origData_(i, 0), origData_.n2(), &destData_(i, 0));
		}
	}
private:
	tbb::enumerable_thread_specific<HilbertEnvelope<TFloat>>& envelopeTLS_;
	const Matrix<T>& origData_;
	Matrix<U>& destData_;
};

template<typename TFloat>
template<typename T, typename U>
void
ParallelHilbertEnvelope<TFloat>::getAnalyticSignalDim2(const Matrix<T>& origData, Matrix<U>& destData)
{
	destData.resize(origData.n1(), origData.n2());

	tbb::enumerable_thread_specific<HilbertEnvelope<TFloat>> envelopeTLS;
	tbb::parallel_for(
		tbb::blocked_range<std::size_t>(0, origData.n1()),
		GetAnalyticSignalDim2<T, U>(envelopeTLS, origData, destData)
	);
}

} // namespace Lab

#endif // PARALLELHILBERTENVELOPE_H
