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

#ifndef MATRIX_H_
#define MATRIX_H_

#include <algorithm> /* swap */
#include <initializer_list>
#include <limits>
#include <ostream>
#include <vector>

#include "Exception.h"



namespace Lab {

template<typename T, typename Alloc=std::allocator<T>>
class Matrix {
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
	class Dim1Iterator {
	public:
		Dim1Iterator(U* baseValuePtr, typename Matrix<T, Alloc>::SizeType increment)
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
		~Dim1Iterator() = default;

		bool operator==(const Dim1Iterator<U>& iter) const {
			return valuePtr_ == iter.valuePtr_ && increment_ == iter.increment_;
		}
		bool operator!=(const Dim1Iterator<U>& iter) const {
			return !(*this == iter);
		}
		Dim1Iterator& operator++() {
			valuePtr_ += increment_;
			return *this;
		}
		Dim1Iterator operator++(int) {
			Dim1Iterator<U> iter{*this};
			valuePtr_ += increment_;
			return iter;
		}
		U& operator*() const {
			return *valuePtr_;
		}
	private:
		Dim1Iterator(Dim1Iterator&&) = delete;
		Dim1Iterator& operator=(Dim1Iterator&&) = delete;

		U* valuePtr_;
		const typename Matrix<T, Alloc>::SizeType increment_;
	};

	typedef Iterator Dim2Iterator;
	typedef ConstIterator ConstDim2Iterator;

	Matrix();
	Matrix(SizeType n1, SizeType n2);
	Matrix(std::initializer_list<std::vector<T>> il);

	SizeType size() const { return data_.size(); }
	SizeType n1() const { return n1_; }
	SizeType n2() const { return n2_; }

	Reference operator()(SizeType i1, SizeType i2);
	ConstReference operator()(SizeType i1, SizeType i2) const;

	Range<Dim1Iterator<T>> range1(SizeType i2) {
		Pointer valuePtr = &data_[i2];
		return Range<Dim1Iterator<T>>(
					Dim1Iterator<T>(valuePtr            , n2_),
					Dim1Iterator<T>(valuePtr + n1_ * n2_, n2_));
	}
	Range<Dim1Iterator<const T>> range1(SizeType i2) const {
		ConstPointer valuePtr = &data_[i2];
		return Range<Dim1Iterator<const T>>(
					Dim1Iterator<const T>(valuePtr            , n2_),
					Dim1Iterator<const T>(valuePtr + n1_ * n2_, n2_));
	}

	Range<Dim2Iterator> range2(SizeType i1) {
		auto baseIter = data_.begin() + i1 * n2_;
		return Range<Dim2Iterator>(baseIter, baseIter + n2_);
	}
	Range<ConstDim2Iterator> range2(SizeType i1) const {
		auto baseIter = data_.cbegin() + i1 * n2_;
		return Range<ConstDim2Iterator>(baseIter, baseIter + n2_);
	}

	void resize(SizeType n1, SizeType n2);
	void reset();
	void operator=(T value);
	bool empty() const { return data_.empty(); }
	void swap(Matrix& other);

	Iterator begin() { return data_.begin(); }
	Iterator end() { return data_.end(); }
	ConstIterator begin() const { return data_.begin(); }
	ConstIterator end() const { return data_.end(); }
	ConstIterator cbegin() const { return data_.cbegin(); }
	ConstIterator cend() const { return data_.cend(); }

	bool operator==(const Matrix<T, Alloc>& m) const {
		return n1_ == m.n1_ && n2_ == m.n2_ && data_ == m.data_;
	}
	bool operator!=(const Matrix<T, Alloc>& m) const {
		return !(*this == m);
	}
private:
	template<typename V>
	friend std::ostream& operator<<(std::ostream& out, const Matrix<V>& m);

	static void validateSize(SizeType n1, SizeType n2);

	SizeType n1_;
	SizeType n2_;
	std::vector<T, Alloc> data_;
};

template<typename T, typename Alloc>
Matrix<T, Alloc>::Matrix() : n1_(), n2_()
{
}

template<typename T, typename Alloc>
Matrix<T, Alloc>::Matrix(SizeType n1, SizeType n2) : n1_(n1), n2_(n2)
{
	validateSize(n1, n2);
	data_.resize(n1 * n2);
}

template<typename T, typename Alloc>
Matrix<T, Alloc>::Matrix(std::initializer_list<std::vector<T>> il)
{
	n1_ = il.size();
	if (n1_ == 0) {
		n2_ = 0;
		return;
	}
	n2_ = il.begin()->size();
	if (n2_ == 0) THROW_EXCEPTION(InvalidValueException, "Empty row.");
	validateSize(n1_, n2_);
	for (auto p = il.begin(); p != il.end(); ++p) {
		if (p->size() != n2_) THROW_EXCEPTION(InvalidValueException, "Inconsistent row size.");
		data_.insert(data_.end(), p->begin(), p->end());
	}
}

template<typename T, typename Alloc>
typename Matrix<T, Alloc>::Reference
Matrix<T, Alloc>::operator()(SizeType i1, SizeType i2)
{
	return data_[i1 * n2_ + i2];
}

template<typename T, typename Alloc>
typename Matrix<T, Alloc>::ConstReference
Matrix<T, Alloc>::operator()(SizeType i1, SizeType i2) const
{
	return data_[i1 * n2_ + i2];
}

template<typename T, typename Alloc>
void
Matrix<T, Alloc>::resize(SizeType n1, SizeType n2)
{
	if (n1 == n1_ && n2 == n2_) return;
	validateSize(n1, n2);

	n1_ = n1;
	n2_ = n2;
	data_.resize(n1 * n2);
}

template<typename T, typename Alloc>
void
Matrix<T, Alloc>::reset()
{
	n1_ = 0;
	n2_ = 0;
	data_.resize(0);
	// The memory is not deallocated.
}

template<typename T, typename Alloc>
void
Matrix<T, Alloc>::validateSize(SizeType n1, SizeType n2)
{
	if (n1 == 0) THROW_EXCEPTION(InvalidValueException, "n1 must be >= 1.");
	if (n2 == 0) THROW_EXCEPTION(InvalidValueException, "n2 must be >= 1.");

	const SizeType maxSize = std::numeric_limits<SizeType>::max() / sizeof(T);

	if (n1 < n2) {
		if (n2 > maxSize / n1) THROW_EXCEPTION(InvalidValueException, "n1*n2 must be <= " << maxSize << '.');
	} else {
		if (n1 > maxSize / n2) THROW_EXCEPTION(InvalidValueException, "n1*n2 must be <= " << maxSize << '.');
	}
}

template<typename T, typename Alloc>
void
Matrix<T, Alloc>::operator=(T value)
{
	std::fill(data_.begin(), data_.end(), value);
}

template<typename T, typename Alloc>
void
Matrix<T, Alloc>::swap(Matrix& other)
{
	std::swap(n1_, other.n1_);
	std::swap(n2_, other.n2_);
	data_.swap(other.data_);
}



template<typename T, typename Alloc>
std::ostream&
operator<<(std::ostream& out, const Matrix<T, Alloc>& m)
{
	for (typename Matrix<T, Alloc>::SizeType i = 0, n1 = m.n1(); i < n1; ++i) {
		auto range = m.range2(i);
		for (auto it = range.begin(); it != range.end(); ++it) {
			out << ' ' << *it;
		}
		out << '\n';
	}

	return out;
}

} // namespace Lab

#endif /* MATRIX_H_ */
