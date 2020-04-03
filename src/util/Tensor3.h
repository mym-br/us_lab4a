/***************************************************************************
 *  Copyright 2014, 2017, 2018, 2019 Marcelo Y. Matuda                     *
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

	template<typename U>
	class Range {
	public:
		Range(U begin, U end) : begin_(begin), end_(end) {}
		U& begin() { return begin_; }
		U& end() { return end_; }
		U& begin() const { return begin_; }
		U& end() const { return end_; }
	private:
		U begin_;
		U end_;
	};

	template<typename U>
	class Dim12Iterator {
	public:
		Dim12Iterator(U* baseValuePtr, typename Tensor3<T, Alloc>::SizeType increment)
			: valuePtr_(baseValuePtr)
			, increment_(increment) {
		}

		bool operator==(const Dim12Iterator<U>& iter) const {
			return valuePtr_ == iter.valuePtr_ && increment_ == iter.increment_;
		}
		bool operator!=(const Dim12Iterator<U>& iter) const {
			return !(*this == iter);
		}
		Dim12Iterator& operator++() {
			valuePtr_ += increment_;
			return *this;
		}
		Dim12Iterator operator++(int) {
			Dim12Iterator<U> iter{*this};
			valuePtr_ += increment_;
			return iter;
		}
		U& operator*() const {
			return *valuePtr_;
		}
	private:
		U* valuePtr_;
		const typename Tensor3<T, Alloc>::SizeType increment_;
	};

	typedef Iterator Dim3Iterator;
	typedef ConstIterator ConstDim3Iterator;

	Tensor3();
	Tensor3(SizeType n1, SizeType n2, SizeType n3);

	SizeType size() const { return data_.size(); }
	SizeType n1() const { return n1_; }
	SizeType n2() const { return n2_; }
	SizeType n3() const { return n3_; }
	const ValueType* data() const { return data_.data(); }
	ValueType* data() { return data_.data(); }

	Reference operator()(SizeType i1, SizeType i2, SizeType i3);
	ConstReference operator()(SizeType i1, SizeType i2, SizeType i3) const;

	Range<Dim12Iterator<T>> range1(SizeType i2, SizeType i3) {
		const SizeType n2n3 = n2_ * n3_;
		Pointer valuePtr = &data_[i2 * n3_ + i3];
		return Range<Dim12Iterator<T>>(
					Dim12Iterator<T>(valuePtr             , n2n3),
					Dim12Iterator<T>(valuePtr + n1_ * n2n3, n2n3));
	}
	Range<Dim12Iterator<const T>> range1(SizeType i2, SizeType i3) const {
		const SizeType n2n3 = n2_ * n3_;
		ConstPointer valuePtr = &data_[i2 * n3_ + i3];
		return Range<Dim12Iterator<const T>>(
					Dim12Iterator<const T>(valuePtr             , n2n3),
					Dim12Iterator<const T>(valuePtr + n1_ * n2n3, n2n3));
	}

	Range<Dim12Iterator<T>> range2(SizeType i1, SizeType i3) {
		const SizeType n2n3 = n2_ * n3_;
		Pointer valuePtr = &data_[i1 * n2n3 + i3];
		return Range<Dim12Iterator<T>>(
					Dim12Iterator<T>(valuePtr       , n3_),
					Dim12Iterator<T>(valuePtr + n2n3, n3_));
	}
	Range<Dim12Iterator<const T>> range2(SizeType i1, SizeType i3) const {
		const SizeType n2n3 = n2_ * n3_;
		ConstPointer valuePtr = &data_[i1 * n2n3 + i3];
		return Range<Dim12Iterator<const T>>(
					Dim12Iterator<const T>(valuePtr       , n3_),
					Dim12Iterator<const T>(valuePtr + n2n3, n3_));
	}

	Range<Dim3Iterator> range3(SizeType i1, SizeType i2) {
		auto baseIter = data_.begin() + (i1 * n2_ + i2) * n3_;
		return Range<Dim3Iterator>(baseIter, baseIter + n3_);
	}
	Range<ConstDim3Iterator> range3(SizeType i1, SizeType i2) const {
		auto baseIter = data_.cbegin() + (i1 * n2_ + i2) * n3_;
		return Range<ConstDim3Iterator>(baseIter, baseIter + n3_);
	}

	void resize(SizeType n1, SizeType n2, SizeType n3);
	void reset();
	void operator=(T value);
	bool empty() const { return data_.empty(); }

	Iterator begin() { return data_.begin(); }
	Iterator end() { return data_.end(); }
	ConstIterator begin() const { return data_.begin(); }
	ConstIterator end() const { return data_.end(); }
	ConstIterator cbegin() const { return data_.cbegin(); }
	ConstIterator cend() const { return data_.cend(); }

	bool operator==(const Tensor3<T, Alloc>& t) const {
		return n1_ == t.n1_ && n2_ == t.n2_ && n3_ == t.n3_ && data_ == t.data_;
	}
	bool operator!=(const Tensor3<T, Alloc>& t) const {
		return !(*this == t);
	}
private:
	template<typename V>
	friend std::ostream& operator<<(std::ostream& out, const Tensor3<V>& t);

	static void validateSize(SizeType n1, SizeType n2, SizeType n3);

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
Tensor3<T, Alloc>::operator()(SizeType i1, SizeType i2, SizeType i3)
{
	return data_[(i1 * n2_ + i2) * n3_ + i3];
}

template<typename T, typename Alloc>
typename Tensor3<T, Alloc>::ConstReference
Tensor3<T, Alloc>::operator()(SizeType i1, SizeType i2, SizeType i3) const
{
	return data_[(i1 * n2_ + i2) * n3_ + i3];
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
	// Sort the array.
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
operator<<(std::ostream& out, const Tensor3<T, Alloc>& t)
{
	for (typename Tensor3<T, Alloc>::SizeType i = 0, n1 = t.n1(); i < n1; ++i) {
		for (typename Tensor3<T, Alloc>::SizeType j = 0, n2 = t.n2(); j < n2; ++j) {
			auto range = t.range3(i, j);
			for (auto it = range.begin(); it != range.end(); ++it) {
				out << ' ' << *it;
			}
			out << '\n';
		}
		out << "\n\n";
	}

	return out;
}

} // namespace Lab

#endif /* TENSOR3_H_ */
