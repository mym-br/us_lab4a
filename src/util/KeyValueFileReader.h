#ifndef KEYVALUEFILEREADER_H_
#define KEYVALUEFILEREADER_H_

#include <QHash>
#include <QString>



namespace Lab {

class KeyValueFileReader {
public:
	KeyValueFileReader(const QString& filePath);
	~KeyValueFileReader();

	const QHash<QString, QString>& map() const { return map_; }
private:
	KeyValueFileReader(const KeyValueFileReader&);
	KeyValueFileReader& operator=(const KeyValueFileReader&);

	QHash<QString, QString> map_;
};

} // namespace Lab

#endif /* KEYVALUEFILEREADER_H_ */
