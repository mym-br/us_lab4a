#include "ParameterMap.h"



namespace Lab {

ParameterMap::ParameterMap(const QString& filePath)
		: filePath_(filePath)
		, reader_(new KeyValueFileReader(filePath))
{
}

bool
ParameterMap::contains(const char* key) const
{
	const QHash<QString, QString>& map = reader_->map();
	return map.contains(key);
}

template<>
short
ParameterMap::convertValue<short>(const QString s) const
{
	bool ok;
	short v = s.toShort(&ok);
	if (!ok) {
		THROW_EXCEPTION(InvalidValueException, "Invalid short: " << s.toStdString() << '.');
	}
	return v;
}

template<>
unsigned short
ParameterMap::convertValue<unsigned short>(const QString s) const
{
	bool ok;
	unsigned short v = s.toUShort(&ok);
	if (!ok) {
		THROW_EXCEPTION(InvalidValueException, "Invalid unsigned short: " << s.toStdString() << '.');
	}
	return v;
}

template<>
int
ParameterMap::convertValue<int>(const QString s) const
{
	bool ok;
	int v = s.toInt(&ok);
	if (!ok) {
		THROW_EXCEPTION(InvalidValueException, "Invalid integer: " << s.toStdString() << '.');
	}
	return v;
}

template<>
unsigned int
ParameterMap::convertValue<unsigned int>(const QString s) const
{
	bool ok;
	unsigned int v = s.toUInt(&ok);
	if (!ok) {
		THROW_EXCEPTION(InvalidValueException, "Invalid unsigned integer: " << s.toStdString() << '.');
	}
	return v;
}

template<>
long
ParameterMap::convertValue<long>(const QString s) const
{
	bool ok;
	long v = s.toLong(&ok);
	if (!ok) {
		THROW_EXCEPTION(InvalidValueException, "Invalid long: " << s.toStdString() << '.');
	}
	return v;
}

template<>
unsigned long
ParameterMap::convertValue<unsigned long>(const QString s) const
{
	bool ok;
	unsigned long v = s.toULong(&ok);
	if (!ok) {
		THROW_EXCEPTION(InvalidValueException, "Invalid unsigned long: " << s.toStdString() << '.');
	}
	return v;
}

template<>
float
ParameterMap::convertValue<float>(const QString s) const
{
	bool ok;
	float v = s.toFloat(&ok);
	if (!ok) {
		THROW_EXCEPTION(InvalidValueException, "Invalid float: " << s.toStdString() << '.');
	}
	return v;
}

template<>
double
ParameterMap::convertValue<double>(const QString s) const
{
	bool ok;
	double v = s.toDouble(&ok);
	if (!ok) {
		THROW_EXCEPTION(InvalidValueException, "Invalid double: " << s.toStdString() << '.');
	}
	return v;
}

template<>
bool
ParameterMap::convertValue<bool>(const QString s) const
{
	if (s == "true") return true;
	if (s == "false") return false;
	THROW_EXCEPTION(InvalidValueException, "Invalid bool: " << s.toStdString() << '.');
}

template<>
std::string
ParameterMap::convertValue<std::string>(const QString s) const
{
	return s.toStdString();
}

} // namespace Lab
