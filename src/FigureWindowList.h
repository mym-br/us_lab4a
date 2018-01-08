#ifndef FIGUREWINDOWLIST_H_
#define FIGUREWINDOWLIST_H_

#include <vector>

#include <boost/shared_ptr.hpp>



namespace Lab {

template<typename W>
class FigureWindowList {
public:
	FigureWindowList() { }
	~FigureWindowList() { }

	W& get(int id);
	void clear();

private:
	struct FigureWindowData {
		int id;
		boost::shared_ptr<W> figure;
	};

	FigureWindowList(const FigureWindowList&);
	FigureWindowList& operator=(const FigureWindowList&);

	std::vector<FigureWindowData> list_;
};

template<typename W>
W&
FigureWindowList<W>::get(int id)
{
	for (typename std::vector<FigureWindowData>::iterator iter = list_.begin(), end = list_.end(); iter != end; ++iter) {
		if (iter->id == id) {
			return *(iter->figure);
		}
	}

	FigureWindowData data;
	data.id = id;
	data.figure.reset(new W);
	list_.push_back(data);
	return *(data.figure);
}

template<typename W>
void
FigureWindowList<W>::clear()
{
	list_.clear();
}

} // namespace Lab

#endif /* FIGUREWINDOWLIST_H_ */
