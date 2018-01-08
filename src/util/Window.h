#ifndef WINDOW_H_
#define WINDOW_H_

#include <vector>

#include "Util.h"



namespace Lab {
namespace Window {

template<typename FloatType>
void
hamming(int n, std::vector<FloatType>& w)
{
	w.resize(n);
	for (int i = 0; i < n; ++i) {
		w[i] = 0.54 - 0.46 * std::cos((2.0 * PI * i) / (n - 1));
	}
}

} // namespace Window
} // namespace Lab

#endif /* WINDOW_H_ */
