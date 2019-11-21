/***************************************************************************
 *  Copyright 2018 Marcelo Y. Matuda                                       *
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

#ifndef MULTILAYERIMAGEMETHOD_H
#define MULTILAYERIMAGEMETHOD_H

#include "Method.h"



namespace Lab {

class Project;

class MultiLayerImageMethod : public Method {
public:
	MultiLayerImageMethod(Project& project);
	virtual ~MultiLayerImageMethod() = default;

	virtual void execute();
private:
	MultiLayerImageMethod(const MultiLayerImageMethod&) = delete;
	MultiLayerImageMethod& operator=(const MultiLayerImageMethod&) = delete;
	MultiLayerImageMethod(MultiLayerImageMethod&&) = delete;
	MultiLayerImageMethod& operator=(MultiLayerImageMethod&&) = delete;

	Project& project_;
};

} // namespace Lab

#endif // MULTILAYERIMAGEMETHOD_H
