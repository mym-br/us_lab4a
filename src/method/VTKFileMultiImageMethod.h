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

#ifndef VTKFILEMULTIIMAGEMETHOD_H
#define VTKFILEMULTIIMAGEMETHOD_H

#include "Method.h"



namespace Lab {

class Project;

class VTKFileMultiImageMethod : public Method {
public:
	VTKFileMultiImageMethod(Project& project);
	virtual ~VTKFileMultiImageMethod();

	virtual void execute();

private:
	Project& project_;
};

} // namespace Lab

#endif // VTKFILEMULTIIMAGEMETHOD_H
