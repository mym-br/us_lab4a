#ifndef ARRAY_H_
#define ARRAY_H_

namespace Lab {

template<typename T, std::size_t N>
class Array {
public:
	typedef T value_type;
	typedef std::size_t size_type;

	Array() {}
	Array(const T& value)
	{
		fill(value);
	}
	~Array() {}

	T& operator[](size_type i)
	{
		return data_[i];
	}

	const T& operator[](size_type i) const
	{
		return data_[i];
	}

	static size_type size() { return N; }

	void fill(const T& value)
	{
		for (size_type i = 0; i < N; ++i) {
			data_[i] = value;
		}
	}
private:
	Array(const Array&);
	Array& operator=(const Array&);

	T data_[N];
};

} // namespace Lab

#endif /* ARRAY_H_ */
