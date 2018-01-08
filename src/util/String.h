#ifndef STRING_H
#define STRING_H

#include <sstream>
#include <string>



namespace Lab {
namespace String {

struct End {
};

class Begin {
public:
	Begin() { }

	template<typename T>
	Begin& operator<<(const T& v) {
		out_ << v;
		return *this;
	}

	std::string operator<<(const End& /*end*/) {
		return out_.str();
	}
private:
	std::ostringstream out_;
};

} // namespace String
} // namespace Lab

#endif // STRING_H
