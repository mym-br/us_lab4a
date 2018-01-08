#ifndef PARALLELHILBERTENVELOPE_H
#define PARALLELHILBERTENVELOPE_H

#include <cstddef> /* std::size_t */

#include <tbb/enumerable_thread_specific.h>
#include <tbb/tbb.h>

#include "HilbertEnvelope.h"
#include "Matrix2.h"



namespace Lab {

template<typename FloatType>
class ParallelHilbertEnvelope {
public:
	template<typename T>
		class CalculateDim2;
	template<typename T, typename U>
		class GetAnalyticSignalDim2;

	template<typename T>
		static void calculateDim2(Matrix2<T>& data);
	template<typename T, typename U>
		static void getAnalyticSignalDim2(const Matrix2<T>& origData, Matrix2<U>& destData);
private:
	ParallelHilbertEnvelope();
	~ParallelHilbertEnvelope();
};



template<typename FloatType>
template<typename T>
class ParallelHilbertEnvelope<FloatType>::CalculateDim2 {
public:
	CalculateDim2(tbb::enumerable_thread_specific<HilbertEnvelope<FloatType> >& envelopeTLS, Matrix2<T>& data)
			: envelopeTLS_(envelopeTLS)
			, data_(data) { }
	void operator()(const tbb::blocked_range<std::size_t>& r) const {
		typename tbb::enumerable_thread_specific<HilbertEnvelope<FloatType> >::reference envelope = envelopeTLS_.local();
		for (std::size_t i = r.begin(); i != r.end(); ++i) {
			envelope.calculate(&data_(i, 0), data_.n2());
		}
	}
private:
	tbb::enumerable_thread_specific<HilbertEnvelope<FloatType> >& envelopeTLS_;
	Matrix2<T>& data_;
};

template<typename FloatType>
template<typename T>
void
ParallelHilbertEnvelope<FloatType>::calculateDim2(Matrix2<T>& data)
{
	tbb::enumerable_thread_specific<HilbertEnvelope<FloatType> > envelopeTLS;
	tbb::parallel_for(
		tbb::blocked_range<std::size_t>(0, data.n1()),
		CalculateDim2<typename Matrix2<T>::ValueType>(envelopeTLS, data)
	);
}



template<typename FloatType>
template<typename T, typename U>
class ParallelHilbertEnvelope<FloatType>::GetAnalyticSignalDim2 {
public:
	GetAnalyticSignalDim2(
			tbb::enumerable_thread_specific<HilbertEnvelope<FloatType> >& envelopeTLS,
			const Matrix2<T>& origData,
			Matrix2<U>& destData)
		: envelopeTLS_(envelopeTLS)
		, origData_(origData)
		, destData_(destData) { }
	void operator()(const tbb::blocked_range<std::size_t>& r) const {
		typename tbb::enumerable_thread_specific<HilbertEnvelope<FloatType> >::reference envelope = envelopeTLS_.local();
		for (std::size_t i = r.begin(); i != r.end(); ++i) {
			envelope.getAnalyticSignal(&origData_(i, 0), origData_.n2(), &destData_(i, 0));
		}
	}
private:
	tbb::enumerable_thread_specific<HilbertEnvelope<FloatType> >& envelopeTLS_;
	const Matrix2<T>& origData_;
	Matrix2<U>& destData_;
};

template<typename FloatType>
template<typename T, typename U>
void
ParallelHilbertEnvelope<FloatType>::getAnalyticSignalDim2(const Matrix2<T>& origData, Matrix2<U>& destData)
{
	destData.resize(origData.n1(), origData.n2());

	tbb::enumerable_thread_specific<HilbertEnvelope<FloatType> > envelopeTLS;
	tbb::parallel_for(
		tbb::blocked_range<std::size_t>(0, origData.n1()),
		GetAnalyticSignalDim2<T, U>(envelopeTLS, origData, destData)
	);
}

} // namespace Lab

#endif // PARALLELHILBERTENVELOPE_H
