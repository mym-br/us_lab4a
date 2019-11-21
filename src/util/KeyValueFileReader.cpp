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

} // namespace Lab
