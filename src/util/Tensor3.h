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

#ifndef TENSOR3_H_
#define TENSOR3_H_

#include <limits>
#include <ostream>
#include <vector>

#include "Exception.h"



namespace Lab {

template<typename T, typename Alloc=std::allocator<T>>
class Tensor3 {
public:
	typedef typename std::vector<T, Alloc>::iterator Iterator;
	typedef typename std::vector<T, Alloc>::const_iterator ConstIterator;
	typedef typename std::vector<T, Alloc>::size_type SizeType;
	typedef T* Pointer;
	typedef const T* ConstPointer;
	typedef T& Reference;
	typedef const T& ConstReference;
	typedef T ValueType;

	class Dim1Iterator {
	public:
		Dim1Iterator(typename Tensor3<T, Alloc>::Pointer baseValuePtr, typename Tensor3<T, Alloc>::SizeType increment)
			: valuePtr_(baseValuePtr)
			, increment_(increment) {
		}
		Dim1Iterator(const Dim1Iterator& iter)
			: valuePtr_(iter.valuePtr_)
			, increment_(iter.increment_) {
		}
		Dim1Iterator& operator=(const Dim1Iterator& iter) {
			if (&iter != this) {
				valuePtr_ = iter.valuePtr_;
				increment_ = iter.increment_;
			}
			return *this;
		}
		bool operator==(const Dim1Iterator& iter) const { return !(valuePtr_ != iter.valuePtr_ || increment_ != iter.increment_); }
		bool operator!=(const Dim1Iterator& iter) const { return valuePtr_ != iter.valuePtr_ || increment_ != iter.increment_; }
		Dim1Iterator& operator++() { valuePtr_ += increment_; return *this; }
		Dim1Iterator operator++(int) { Dim1Iterator iter{*this}; valuePtr_ += increment_; return iter; }
		typename Tensor3<T, Alloc>::Reference operator*() { return *valuePtr_; }
	private:
		typename Tensor3<T, Alloc>::Pointer valuePtr_;
		const typename Tensor3<T, Alloc>::SizeType increment_;
	};

	class ConstDim1Iterator {
	public:
		ConstDim1Iterator(typename Tensor3<T, Alloc>::ConstPointer baseValuePtr, typename Tensor3<T, Alloc>::SizeType increment)
			: valuePtr_(baseValuePtr)
			, increment_(increment) {
		}
		ConstDim1Iterator(const ConstDim1Iterator& iter)
			: valuePtr_(iter.valuePtr_)
			, increment_(iter.increment_) {
		}
		ConstDim1Iterator& operator=(const ConstDim1Iterator& iter) {
			if (&iter != this) {
				valuePtr_ = iter.valuePtr_;
				increment_ = iter.increment_;
			}
			return *this;
		}
		bool operator==(const ConstDim1Iterator& iter) const { return !(valuePtr_ != iter.valuePtr_ || increment_ != iter.increment_); }
		bool operator!=(const ConstDim1Iterator& iter) const { return valuePtr_ != iter.valuePtr_ || increment_ != iter.increment_; }
		ConstDim1Iterator& operator++() { valuePtr_ += increment_; return *this; }
		ConstDim1Iterator operator++(int) { ConstDim1Iterator iter{*this}; valuePtr_ += increment_; return iter; }
		typename Tensor3<T, Alloc>::ConstReference operator*() { return *valuePtr_; }
	private:
		typename Tensor3<T, Alloc>::ConstPointer valuePtr_;
		const typename Tensor3<T, Alloc>::SizeType increment_;
	};

	typedef std::pair<Dim1Iterator, Dim1Iterator> Dim1Interval;
	typedef std::pair<ConstDim1Iterator, ConstDim1Iterator> ConstDim1Interval;

	typedef Dim1Iterator Dim2Iterator;
	typedef ConstDim1Iterator ConstDim2Iterator;
	typedef std::pair<Dim2Iterator, Dim2Iterator> Dim2Interval;
	typedef std::pair<ConstDim2Iterator, ConstDim2Iterator> ConstDim2Interval;

	typedef Iterator Dim3Iterator;
	typedef ConstIterator ConstDim3Iterator;
	typedef std::pair<Dim3Iterator, Dim3Iterator> Dim3Interval;
	typedef std::pair<ConstDim3Iterator, ConstDim3Iterator> ConstDim3Interval;

	Tensor3();
	Tensor3(SizeType n1, SizeType n2, SizeType n3);

	SizeType size() const { return data_.size(); }
	SizeType n1() const { return n1_; }
	SizeType n2() const { return n2_; }
	SizeType n3() const { return n3_; }

	Reference operator()(SizeType dim1, SizeType dim2, SizeType dim3);
	ConstReference operator()(SizeType dim1, SizeType dim2, SizeType dim3) const;

	Dim1Interval dim1Interval(SizeType dim2, SizeType dim3) {
		Pointer valuePtr = &data_[dim2 * n3_ + dim3];
		const SizeType n2n3 = n2_ * n3_;
		return std::make_pair(Dim1Iterator(valuePtr, n2n3), Dim1Iterator(valuePtr + n1_ * n2n3, n2n3));
	}

	ConstDim1Interval dim1Interval(SizeType dim2, SizeType dim3) const {
		ConstPointer valuePtr = &data_[dim2 * n3_ + dim3];
		const SizeType n2n3 = n2_ * n3_;
		return std::make_pair(ConstDim1Iterator(valuePtr, n2n3), ConstDim1Iterator(valuePtr + n1_ * n2n3, n2n3));
	}

	Dim2Interval dim2Interval(SizeType dim1, SizeType dim3) {
		const SizeType n2n3 = n2_ * n3_;
		Pointer valuePtr = &data_[dim1 * n2n3 + dim3];
		return std::make_pair(Dim2Iterator(valuePtr, n3_), Dim2Iterator(valuePtr + n2n3, n3_));
	}

	ConstDim2Interval dim2Interval(SizeType dim1, SizeType dim3) const {
		const SizeType n2n3 = n2_ * n3_;
		ConstPointer valuePtr = &data_[dim1 * n2n3 + dim3];
		return std::make_pair(ConstDim2Iterator(valuePtr, n3_), ConstDim2Iterator(valuePtr + n2n3, n3_));
	}

	Dim3Interval dim3Interval(SizeType dim1, SizeType dim2) {
		Iterator baseIter = data_.begin() + (dim1 * n2_ + dim2) * n3_;
		return std::make_pair(baseIter, baseIter + n3_);
	}

	ConstDim3Interval dim3Interval(SizeType dim1, SizeType dim2) const {
		ConstIterator baseIter = data_.begin() + (dim1 * n2_ + dim2) * n3_;
		return std::make_pair(baseIter, baseIter + n3_);
	}

	void resize(SizeType n1, SizeType n2, SizeType n3);
	void reset();
	void operator=(T value);
	bool empty() const { return data_.size() == 0; }

	Iterator begin() { return data_.begin(); }
	Iterator end() { return data_.end(); }
	ConstIterator begin() const { return data_.begin(); }
	ConstIterator end() const { return data_.end(); }
private:
	template<typename V>
	friend std::ostream& operator<<(std::ostream& out, const Tensor3<V>& m);

	void validateSize(SizeType n1, SizeType n2, SizeType n3);

	SizeType n1_;
	SizeType n2_;
	SizeType n3_;
	std::vector<T, Alloc> data_;
};

template<typename T, typename Alloc>
Tensor3<T, Alloc>::Tensor3() : n1_(), n2_(), n3_()
{
}

template<typename T, typename Alloc>
Tensor3<T, Alloc>::Tensor3(SizeType n1, SizeType n2, SizeType n3) : n1_(n1), n2_(n2), n3_(n3)
{
	validateSize(n1, n2, n3);

	data_.resize(n1 * n2 * n3);
}

template<typename T, typename Alloc>
typename Tensor3<T, Alloc>::Reference
Tensor3<T, Alloc>::operator()(SizeType dim1, SizeType dim2, SizeType dim3)
{
	return data_[(dim1 * n2_ + dim2) * n3_ + dim3];
}

template<typename T, typename Alloc>
typename Tensor3<T, Alloc>::ConstReference
Tensor3<T, Alloc>::operator()(SizeType dim1, SizeType dim2, SizeType dim3) const
{
	return data_[(dim1 * n2_ + dim2) * n3_ + dim3];
}

template<typename T, typename Alloc>
void
Tensor3<T, Alloc>::resize(SizeType n1, SizeType n2, SizeType n3)
{
	if (n1 == n1_ && n2 == n2_ && n3 == n3_) return;
	validateSize(n1, n2, n3);

	n1_ = n1;
	n2_ = n2;
	n3_ = n3;
	data_.resize(n1 * n2 * n3);
}

template<typename T, typename Alloc>
void
Tensor3<T, Alloc>::reset()
{
	n1_ = 0;
	n2_ = 0;
	n3_ = 0;
	data_.resize(0);
	// The memory is not deallocated.
}

template<typename T, typename Alloc>
void
Tensor3<T, Alloc>::validateSize(SizeType n1, SizeType n2, SizeType n3)
{
	if (n1 == 0) THROW_EXCEPTION(InvalidValueException, "n1 must be >= 1.");
	if (n2 == 0) THROW_EXCEPTION(InvalidValueException, "n2 must be >= 1.");
	if (n3 == 0) THROW_EXCEPTION(InvalidValueException, "n3 must be >= 1.");

	SizeType n[] = { n1, n2, n3 };
	// Sorts the array.
	if (n[1] > n[2]) std::swap(n[1], n[2]);
	if (n[0] > n[1]) std::swap(n[0], n[1]);
	if (n[1] > n[2]) std::swap(n[1], n[2]);

	const SizeType maxSize = std::numeric_limits<SizeType>::max() / sizeof(T);

	if (n[2] > maxSize / n[1]) THROW_EXCEPTION(InvalidValueException, "n1*n2 must be <= " << maxSize << '.');
	if (n[1] * n[2] > maxSize / n[0]) THROW_EXCEPTION(InvalidValueException, "n1*n2*n3 must be <= " << maxSize << '.');
}

template<typename T, typename Alloc>
void
Tensor3<T, Alloc>::operator=(T value)
{
	std::fill(data_.begin(), data_.end(), value);
}



template<typename T, typename Alloc>
std::ostream&
operator<<(std::ostream& out, const Tensor3<T, Alloc>& m)
{
	for (typename Tensor3<T, Alloc>::SizeType i = 0, n1 = m.n1(); i < n1; ++i) {
		for (typename Tensor3<T, Alloc>::SizeType j = 0, n2 = m.n2(); j < n2; ++j) {
			typename Tensor3<T, Alloc>::ConstDim3Interval dim3Interval = m.dim3Interval(i, j);
			for (typename Tensor3<T, Alloc>::ConstDim3Iterator iter = dim3Interval.first; iter != dim3Interval.second; ++iter) {
				out << ' ' << *iter;
			}
			out << '\n';
		}
		out << "\n\n";
	}

	return out;
}

} // namespace Lab

#endif /* TENSOR3_H_ */
