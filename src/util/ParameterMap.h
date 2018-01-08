#ifndef PARAMETERMAP_H
#define PARAMETERMAP_H

#include <QHash>
#include <QString>

#include <boost/shared_ptr.hpp>
#include <boost/scoped_ptr.hpp>

#include "Exception.h"
#include "KeyValueFileReader.h"



namespace Lab {

class ParameterMap {
public:
	ParameterMap(const QString& filePath);

	bool contains(const char* key) const;
	template<typename T> T value(const char* key) const;
	template<typename T> T value(const char* key, T minValue, T maxValue) const;
private:
	ParameterMap(const ParameterMap&);
	ParameterMap& operator=(const ParameterMap&);

	template<typename T> T convertValue(const QString s) const;

	const QString filePath_;
	boost::scoped_ptr<const KeyValueFileReader> reader_;
};

typedef boost::shared_ptr<const ParameterMap> ConstParameterMapPtr;



template<typename T>
T
ParameterMap::value(const char* key) const
{
	const QHash<QString, QString>& map = reader_->map();

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
