#ifndef FILEUTIL_H
#define FILEUTIL_H

#include <fstream>
#include <iterator> /* istreambuf_iterator */
#include <string>

#include "Exception.h"



namespace Lab {
namespace FileUtil {

std::string
loadASCIIFileToString(const char* fileName)
{
	std::ifstream in(fileName);
	if (!in) THROW_EXCEPTION(IOException, "File not found: " << fileName << '.');

	std::string s(
		(std::istreambuf_iterator<char>(in)),
		(std::istreambuf_iterator<char>()));
	if (!in) THROW_EXCEPTION(IOException, "An error occurred while reading the file " << fileName << '.');
	return s;
}

} // namespace FileUtil
} // namespace Lab

#endif // FILEUTIL_H
