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

#include <memory>
#include <vector>



namespace Lab {

template<typename W>
class FigureWindowList {
public:
	FigureWindowList() = default;
	~FigureWindowList() = default;

	W& get(int id);
	void clear();
private:
	struct FigureWindowData {
		int id;
		std::unique_ptr<W> figure;
	};

	FigureWindowList(const FigureWindowList&) = delete;
	FigureWindowList& operator=(const FigureWindowList&) = delete;
	FigureWindowList(FigureWindowList&&) = delete;
	FigureWindowList& operator=(FigureWindowList&&) = delete;

	std::vector<FigureWindowData> list_;
};

template<typename W>
W&
FigureWindowList<W>::get(int id)
{
	for (auto& item : list_) {
		if (item.id == id) {
			return *item.figure;
		}
	}

	list_.push_back(FigureWindowData{id, std::make_unique<W>()});
	return *list_.back().figure;
}

template<typename W>
void
FigureWindowList<W>::clear()
{
	list_.clear();
}

} // namespace Lab

#endif /* FIGUREWINDOWLIST_H_ */
