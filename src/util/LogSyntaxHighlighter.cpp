// Adapted from the QSyntaxHighlighter documentation in Qt Assistant.
// This file is in the public domain.

#include "LogSyntaxHighlighter.h"

#include <QString>

#include "Log.h"



namespace Lab {

LogSyntaxHighlighter::LogSyntaxHighlighter(QTextDocument* parent)
		: QSyntaxHighlighter(parent)
		, errorBeginRegexp_(ERROR_LOG_PREFIX)
		, errorEndRegexp_(ERROR_LOG_SUFFIX)
{
}

void
LogSyntaxHighlighter::highlightBlock(const QString& text)
{
	setCurrentBlockState(0);

	int beginIndex = (previousBlockState() == 1) ? 0 : text.indexOf(errorBeginRegexp_);
	while (beginIndex >= 0) {
		int endIndex = text.indexOf(errorEndRegexp_, beginIndex);

		int length;
		if (endIndex == -1) {
			setCurrentBlockState(1);
			length = text.length() - beginIndex;
		} else {
			length = endIndex - beginIndex + errorEndRegexp_.matchedLength();
		}
		setFormat(beginIndex, length, Qt::red);
		beginIndex = text.indexOf(errorBeginRegexp_, beginIndex + length);
	}
}

} // namespace Lab
