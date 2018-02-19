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
