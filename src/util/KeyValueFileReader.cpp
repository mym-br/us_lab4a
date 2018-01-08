#include "KeyValueFileReader.h"

#include <QFile>
#include <QStringList>
#include <QTextStream>

#include "Exception.h"



namespace Lab {

KeyValueFileReader::KeyValueFileReader(const QString& filePath)
{
	QFile file(filePath);
	if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
		THROW_EXCEPTION(IOException, "The file " << filePath.toStdString() << " could not be opened.");
	}

	QTextStream in(&file);
	in.setCodec("UTF-8");
	int lineNumber = 0;
	while (!in.atEnd()) {
		++lineNumber;

		QString line = in.readLine().trimmed();

		if (line.startsWith('#')) continue; // comment
		if (line.isEmpty()) continue;

		QStringList itemList = line.split('='/*, QString::SkipEmptyParts*/);
		if (itemList.size() == 1) {
			THROW_EXCEPTION(InvalidValueException, "Missing '=' separator at line " << lineNumber << " of file " << filePath.toStdString() << '.');
		}
		if (itemList.size() != 2) {
			THROW_EXCEPTION(InvalidValueException, "Invalid syntax at line " << lineNumber << " of file " << filePath.toStdString() << '.');
		}

		QString key = itemList.at(0).trimmed();
		if (key.isEmpty()) {
			THROW_EXCEPTION(InvalidValueException, "Empty key at line " << lineNumber << " of file " << filePath.toStdString() << '.');
		}
		QString value = itemList.at(1).trimmed();
		if (value.isEmpty()) {
			THROW_EXCEPTION(InvalidValueException, "Empty value at line " << lineNumber << " of file " << filePath.toStdString() << '.');
		}
		if (map_.contains(key)) {
			THROW_EXCEPTION(InvalidValueException, "Duplicated key '" << key.toStdString() << "' at line " << lineNumber << " of file " << filePath.toStdString() << '.');
		}
		map_.insert(key, value);
	}
}

KeyValueFileReader::~KeyValueFileReader()
{
}

} // namespace Lab
