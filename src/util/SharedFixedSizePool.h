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

#ifndef SHAREDFIXEDSIZEPOOL_H_
#define SHAREDFIXEDSIZEPOOL_H_

#include <cstddef> /* std::size_t */
#include <vector>

#include "tbb/cache_aligned_allocator.h"
#include "tbb/concurrent_queue.h"



namespace Lab {

//TODO: remove this class. The allocator is useless.

//
// Note: the type T must be copy constructible.
template<typename T>
class SharedFixedSizePool {
public:
	SharedFixedSizePool(std::size_t size);
	virtual ~SharedFixedSizePool();

	std::size_t size() const { return pool_.size(); }

	T* get() {
		T* p;
		ptrQueue_.pop(p); // blocking
		return p;
	}
	void put(T* p) {
		ptrQueue_.push(p); // blocking
	}
private:
	SharedFixedSizePool(const SharedFixedSizePool&);
	SharedFixedSizePool& operator=(const SharedFixedSizePool&);

	std::vector<T, tbb::cache_aligned_allocator<T>> pool_;
	tbb::concurrent_bounded_queue<T*> ptrQueue_;
};

template<typename T>
SharedFixedSizePool<T>::SharedFixedSizePool(std::size_t size)
		: pool_(size)
{
	//Makes the container SLOW according to Intel.
	//ptrQueue_.set_capacity(size);

	// Fills the pointer queue.
	for (std::size_t i = 0; i < size; ++i) {
		ptrQueue_.push(&pool_[i]);
	}
}

template<typename T>
SharedFixedSizePool<T>::~SharedFixedSizePool()
{
}

} // namespace Lab

#endif /* SHAREDFIXEDSIZEPOOL_H_ */
