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

#ifndef PARAMETERMAP_H
#define PARAMETERMAP_H

#include <memory>
#include <string>

#include <QHash>
#include <QString>

#include "Exception.h"
#include "KeyValueFileReader.h"



namespace Lab {

class ParameterMap {
public:
	explicit ParameterMap(const QString& filePath);
	~ParameterMap() = default;

	bool contains(const char* key) const;
	template<typename T> T value(const char* key) const;
	template<typename T> T value(const char* key, T minValue, T maxValue) const;

	void getValue(double& dest, const char* key, double minValue, double maxValue) const {
		dest = value<double>(key, minValue, maxValue);
	}
	void getValue(float& dest, const char* key, float minValue, float maxValue) const {
		dest = value<float>(key, minValue, maxValue);
	}
	void getValue(long& dest, const char* key, long minValue, long maxValue) const {
		dest = value<long>(key, minValue, maxValue);
	}
	void getValue(unsigned long& dest, const char* key, unsigned long minValue, unsigned long maxValue) const {
		dest = value<unsigned long>(key, minValue, maxValue);
	}
	void getValue(int& dest, const char* key, int minValue, int maxValue) const {
		dest = value<int>(key, minValue, maxValue);
	}
	void getValue(unsigned int& dest, const char* key, unsigned int minValue, unsigned int maxValue) const {
		dest = value<unsigned int>(key, minValue, maxValue);
	}
	void getValue(short& dest, const char* key, short minValue, short maxValue) const {
		dest = value<short>(key, minValue, maxValue);
	}
	void getValue(unsigned short& dest, const char* key, unsigned short minValue, unsigned short maxValue) const {
		dest = value<unsigned short>(key, minValue, maxValue);
	}
	void getValue(bool& dest, const char* key) const {
		dest = value<bool>(key);
	}
	void getValue(std::string& dest, const char* key) const {
		dest = value<std::string>(key);
	}
private:
	ParameterMap(const ParameterMap&) = delete;
	ParameterMap& operator=(const ParameterMap&) = delete;
	ParameterMap(ParameterMap&&) = delete;
	ParameterMap& operator=(ParameterMap&&) = delete;

	template<typename T> static T convertValue(const QString s);

	const QString filePath_;
	const KeyValueFileReader reader_;
};

typedef std::unique_ptr<const ParameterMap> ParamMapPtr;

template<typename T>
T
ParameterMap::value(const char* key) const
{
	const auto& map = reader_.map();

	if (!map.contains(key)) {
		THROW_EXCEPTION(InvalidValueException, "Key '" << key << "' not found in file " << filePath_.toStdString() << '.');
	}

	QString valueString = map.value(key);
	T value;
	try {
		value = convertValue<T>(valueString);
	} catch (const std::exception& e) {
		THROW_EXCEPTION(InvalidValueException, "Invalid value for key '" << key << "' in file " << filePath_.toStdString() << ": " << e.what());
	}
	return value;
}

template<typename T>
T
ParameterMap::value(const char* key, T minValue, T maxValue) const
{
	T v = value<T>(key);

	if (v > maxValue) {
		THROW_EXCEPTION(InvalidValueException, "The value for key '" << key << "' must be <= " << maxValue << " in file " << filePath_.toStdString() << '.');
	} else if (v < minValue) {
		THROW_EXCEPTION(InvalidValueException, "The value for key '" << key << "' must be >= " << minValue << " in file " << filePath_.toStdString() << '.');
	}

	return v;
}

} // namespace Lab

#endif // PARAMETERMAP_H
