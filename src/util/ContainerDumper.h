#ifndef CONTAINERDUMPER_H_
#define CONTAINERDUMPER_H_

#include <fstream>



namespace Lab {

class ContainerDumper {
public:
	template<typename InputIterator> static void save(const char* fileName, InputIterator first, InputIterator last);
private:
	ContainerDumper();
};

template<typename InputIterator>
void
ContainerDumper::save(const char* fileName, InputIterator first, InputIterator last)
{
	std::ofstream out(fileName);
	while (first != last) {
		out << *first++ << '\n';
	}
}

} // namespace Lab

#endif /* CONTAINERDUMPER_H_ */
