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

#ifndef RAWBUFFER_H_
#define RAWBUFFER_H_

#include <cstddef> /* std::size_t */
#include <cstring>
#include <limits>
#include <string>
#include <vector>

#include <boost/cstdint.hpp>

#include "Exception.h"



namespace Lab {

/*******************************************************************************
 *
 */
class RawBuffer {
public:
	typedef boost::uint8_t* iterator;
	typedef const boost::uint8_t* const_iterator;

	RawBuffer() : readIndex_() {
		buffer_.reserve(8192);
	}
	~RawBuffer() {}

	void reset();
	void reserve(std::size_t size);

	void putInt16(boost::int16_t value);
	boost::int16_t getInt16();

	void putUInt16(boost::uint16_t value);
	boost::uint16_t getUInt16();

	void putInt32(boost::int32_t value);
	boost::int32_t getInt32();

	void putUInt32(boost::uint32_t value);
	boost::uint32_t getUInt32();

	void putFloat(float value);
	float getFloat();

	void putString(const std::string& s);
	void getString(std::string& s);

	void putFloatArray(const std::vector<float>& a);
	void getFloatArray(std::vector<float>& a);

	template<typename T> void putInt16Array(const std::vector<T>& a);
	template<typename T> void putInt16Array(const T* a, std::size_t arraySize);
	template<typename T> void getInt16Array(std::vector<T>& a);
	template<typename T> void getInt16Array(T* a, std::size_t arraySize);

	boost::uint8_t& front() {
		return buffer_.front();
	}

	const boost::uint8_t& front() const {
		return buffer_.front();
	}

	std::size_t size() const {
		return buffer_.size();
	}
private:
	RawBuffer(const RawBuffer&);
	RawBuffer& operator=(const RawBuffer&);

	static void writeInt16(boost::int16_t value, std::vector<boost::uint8_t>& buffer, std::size_t& index);
	static boost::int16_t readInt16(const std::vector<boost::uint8_t>& buffer, std::size_t& index);
	static void writeUInt16(boost::uint16_t value, std::vector<boost::uint8_t>& buffer, std::size_t& index);
	static boost::uint16_t readUInt16(const std::vector<boost::uint8_t>& buffer, std::size_t& index);
	static void writeInt32(boost::int32_t value, std::vector<boost::uint8_t>& buffer, std::size_t& index);
	static boost::int32_t readInt32(const std::vector<boost::uint8_t>& buffer, std::size_t& index);
	static void writeUInt32(boost::uint32_t value, std::vector<boost::uint8_t>& buffer, std::size_t& index);
	static boost::uint32_t readUInt32(const std::vector<boost::uint8_t>& buffer, std::size_t& index);

	std::size_t readIndex_;
	std::vector<boost::uint8_t> buffer_; // big-endian
};


/*******************************************************************************
 *
 */
inline
void
RawBuffer::writeInt16(boost::int16_t value, std::vector<boost::uint8_t>& buffer, std::size_t& index)
{
	buffer[index++] = value >> 8U;
	buffer[index++] = value;
}

/*******************************************************************************
 *
 */
inline
boost::int16_t
RawBuffer::readInt16(const std::vector<boost::uint8_t>& buffer, std::size_t& index)
{
	union {
		boost::int16_t i;
		boost::uint16_t ui;
	};
	ui  = buffer[index++] << 8U;
	ui += buffer[index++];
	return i;
}

/*******************************************************************************
 *
 */
inline
void
RawBuffer::writeUInt16(boost::uint16_t value, std::vector<boost::uint8_t>& buffer, std::size_t& index)
{
	buffer[index++] = value >> 8U;
	buffer[index++] = value;
}

/*******************************************************************************
 *
 */
inline
boost::uint16_t
RawBuffer::readUInt16(const std::vector<boost::uint8_t>& buffer, std::size_t& index)
{
	boost::uint16_t value = buffer[index++] << 8U;
	value                += buffer[index++];
	return value;
}

/*******************************************************************************
 *
 */
inline
void
RawBuffer::writeInt32(boost::int32_t value, std::vector<boost::uint8_t>& buffer, std::size_t& index)
{
	buffer[index++] = value >> 24U;
	buffer[index++] = value >> 16U;
	buffer[index++] = value >> 8U;
	buffer[index++] = value;
}

/*******************************************************************************
 *
 */
inline
boost::int32_t
RawBuffer::readInt32(const std::vector<boost::uint8_t>& buffer, std::size_t& index)
{
	union {
		boost::int32_t i;
		boost::uint32_t ui;
	};
	ui  = buffer[index++] << 24U;
	ui += buffer[index++] << 16U;
	ui += buffer[index++] << 8U;
	ui += buffer[index++];
	return i;
}

/*******************************************************************************
 *
 */
inline
void
RawBuffer::writeUInt32(boost::uint32_t value, std::vector<boost::uint8_t>& buffer, std::size_t& index)
{
	buffer[index++] = value >> 24U;
	buffer[index++] = value >> 16U;
	buffer[index++] = value >> 8U;
	buffer[index++] = value;
}

/*******************************************************************************
 *
 */
inline
boost::uint32_t
RawBuffer::readUInt32(const std::vector<boost::uint8_t>& buffer, std::size_t& index)
{
	boost::uint32_t value = buffer[index++] << 24U;
	value                += buffer[index++] << 16U;
	value                += buffer[index++] << 8U;
	value                += buffer[index++];
	return value;
}

/*******************************************************************************
 *
 */
inline
void
RawBuffer::reset()
{
	readIndex_ = 0;
	buffer_.resize(0);
}

/*******************************************************************************
 *
 */
inline
void
RawBuffer::reserve(std::size_t size)
{
	readIndex_ = 0;
	buffer_.resize(size);
}

/*******************************************************************************
 *
 */
inline
void
RawBuffer::putInt16(boost::int16_t value)
{
	std::size_t endIndex = buffer_.size();
	buffer_.resize(endIndex + sizeof(value));
	writeInt16(value, buffer_, endIndex);
}

/*******************************************************************************
 *
 */
inline
boost::int16_t
RawBuffer::getInt16()
{
	if (readIndex_ > buffer_.size() - sizeof(boost::int16_t)) {
		THROW_EXCEPTION(EndOfBufferException, "Could not get an int16 from the buffer.");
	}

	return readInt16(buffer_, readIndex_);
}

/*******************************************************************************
 *
 */
inline
void
RawBuffer::putUInt16(boost::uint16_t value)
{
	std::size_t endIndex = buffer_.size();
	buffer_.resize(endIndex + sizeof(value));
	writeUInt16(value, buffer_, endIndex);
}

/*******************************************************************************
 *
 */
inline
boost::uint16_t
RawBuffer::getUInt16()
{
	if (readIndex_ > buffer_.size() - sizeof(boost::uint16_t)) {
		THROW_EXCEPTION(EndOfBufferException, "Could not get an uint16 from the buffer.");
	}

	return readUInt16(buffer_, readIndex_);
}

/*******************************************************************************
 *
 */
inline
void
RawBuffer::putInt32(boost::int32_t value)
{
	std::size_t endIndex = buffer_.size();
	buffer_.resize(endIndex + sizeof(value));
	writeInt32(value, buffer_, endIndex);
}

/*******************************************************************************
 *
 */
inline
boost::int32_t
RawBuffer::getInt32()
{
	if (readIndex_ > buffer_.size() - sizeof(boost::int32_t)) {
		THROW_EXCEPTION(EndOfBufferException, "Could not get a int32 from the buffer.");
	}

	return readInt32(buffer_, readIndex_);
}

/*******************************************************************************
 *
 */
inline
void
RawBuffer::putUInt32(boost::uint32_t value)
{
	std::size_t endIndex = buffer_.size();
	buffer_.resize(endIndex + sizeof(value));
	writeUInt32(value, buffer_, endIndex);
}

/*******************************************************************************
 *
 */
inline
boost::uint32_t
RawBuffer::getUInt32()
{
	if (readIndex_ > buffer_.size() - sizeof(boost::uint32_t)) {
		THROW_EXCEPTION(EndOfBufferException, "Could not get a uint32 from the buffer.");
	}

	return readUInt32(buffer_, readIndex_);
}

/*******************************************************************************
 *
 */
inline
void
RawBuffer::putFloat(float value)
{
	union {
		boost::uint32_t i;
		float f;
	};
	f = value;
	std::size_t endIndex = buffer_.size();
	buffer_.resize(endIndex + sizeof(value));
	writeUInt32(i, buffer_, endIndex);
}

/*******************************************************************************
 *
 */
inline
float
RawBuffer::getFloat()
{
	if (readIndex_ > buffer_.size() - sizeof(float)) {
		THROW_EXCEPTION(EndOfBufferException, "Could not get a float from the buffer.");
	}

	union {
		boost::uint32_t i;
		float f;
	};
	i = readUInt32(buffer_, readIndex_);
	return f;
}

/*******************************************************************************
 *
 */
inline
void
RawBuffer::putString(const std::string& s)
{
	std::size_t stringSize = s.size();
	if (stringSize > std::numeric_limits<boost::uint32_t>::max()) {
		THROW_EXCEPTION(InvalidValueException, "The string is too long.");
	}
	putUInt32(stringSize);

	std::size_t endIndex = buffer_.size();
	buffer_.resize(endIndex + stringSize);
	memcpy(&buffer_[endIndex], &s[0], stringSize);
}

/*******************************************************************************
 *
 */
inline
void
RawBuffer::getString(std::string& s)
{
	boost::uint32_t stringSize = getUInt32();

	if (readIndex_ > buffer_.size() - stringSize) {
		THROW_EXCEPTION(EndOfBufferException, "Could not get the string from the buffer.");
	}

	s.resize(stringSize);
	memcpy(&s[0], &buffer_[readIndex_], stringSize);
	readIndex_ += stringSize;
}

/*******************************************************************************
 *
 */
inline
void
RawBuffer::putFloatArray(const std::vector<float>& a)
{
	std::size_t arraySize = a.size();
	if (arraySize > std::numeric_limits<boost::uint32_t>::max()) {
		THROW_EXCEPTION(InvalidValueException, "The array is too long.");
	}
	putUInt32(arraySize);

	union {
		boost::uint32_t i;
		float f;
	};
	std::size_t endIndex = buffer_.size();
	buffer_.resize(endIndex + arraySize * sizeof(boost::uint32_t));
	for (boost::uint32_t j = 0; j < arraySize; ++j) {
		f = a[j];
		writeUInt32(i, buffer_, endIndex);
	}
}

/*******************************************************************************
 *
 */
inline
void
RawBuffer::getFloatArray(std::vector<float>& a)
{
	boost::uint32_t arraySize = getUInt32();

	if (readIndex_ > buffer_.size() - arraySize * sizeof(float)) {
		THROW_EXCEPTION(EndOfBufferException, "Could not get the float array from the buffer.");
	}

	a.resize(arraySize);
	union {
		boost::uint32_t i;
		float f;
	};
	for (boost::uint32_t j = 0; j < arraySize; ++j) {
		i = readUInt32(buffer_, readIndex_);
		a[j] = f;
	}
}

/*******************************************************************************
 *
 */
template<typename T>
void
RawBuffer::putInt16Array(const std::vector<T>& a)
{
	putInt16Array(&a[0], a.size());
}

/*******************************************************************************
 *
 */
template<typename T>
void
RawBuffer::putInt16Array(const T* a, std::size_t arraySize)
{
	if (arraySize > std::numeric_limits<boost::uint32_t>::max()) {
		THROW_EXCEPTION(InvalidValueException, "The array is too long.");
	}
	putUInt32(arraySize);

	std::size_t endIndex = buffer_.size();
	buffer_.resize(endIndex + arraySize * sizeof(boost::int16_t));
	for (boost::uint32_t j = 0; j < arraySize; ++j) {
		writeInt16(static_cast<boost::int16_t>(a[j]), buffer_, endIndex);
	}
}

/*******************************************************************************
 *
 */
template<typename T>
void
RawBuffer::getInt16Array(std::vector<T>& a)
{
	boost::uint32_t arraySize = getUInt32();

	if (readIndex_ > buffer_.size() - arraySize * sizeof(boost::int16_t)) {
		THROW_EXCEPTION(EndOfBufferException, "Could not get the int16 array from the buffer.");
	}

	a.resize(arraySize);
	for (boost::uint32_t j = 0; j < arraySize; ++j) {
		a[j] = static_cast<T>(readInt16(buffer_, readIndex_));
	}
}

/*******************************************************************************
 *
 */
template<typename T>
void
RawBuffer::getInt16Array(T* a, std::size_t arraySize)
{
	boost::uint32_t receivedArraySize = getUInt32();
	if (receivedArraySize != arraySize) {
		THROW_EXCEPTION(WrongBufferSizeException, "Wrong buffer size. Received=" << receivedArraySize << ", correct=" << arraySize << '.');
	}

	if (readIndex_ > buffer_.size() - arraySize * sizeof(boost::int16_t)) {
		THROW_EXCEPTION(EndOfBufferException, "Could not get the int16 array from the buffer.");
	}

	for (boost::uint32_t j = 0; j < arraySize; ++j) {
		a[j] = static_cast<T>(readInt16(buffer_, readIndex_));
	}
}

} // namespace Lab

#endif /* RAWBUFFER_H_ */
