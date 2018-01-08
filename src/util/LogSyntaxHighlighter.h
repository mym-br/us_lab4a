#ifndef LOGSYNTAXHIGHLIGHTER_H_
#define LOGSYNTAXHIGHLIGHTER_H_

#include <QSyntaxHighlighter>

#include <QRegExp>



QT_BEGIN_NAMESPACE
class QString;
class QTextDocument;
QT_END_NAMESPACE

namespace Lab {

class LogSyntaxHighlighter : public QSyntaxHighlighter {
	Q_OBJECT
public:
	LogSyntaxHighlighter(QTextDocument* parent);
protected:
	void highlightBlock(const QString& text);

	QRegExp errorBeginRegexp_;
	QRegExp errorEndRegexp_;
};

} // namespace Lab

#endif /* LOGSYNTAXHIGHLIGHTER_H_ */
