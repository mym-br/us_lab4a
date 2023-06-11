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

#ifndef PROCESSINGNODE_H_
#define PROCESSINGNODE_H_

#include <QObject>



namespace Lab {

class Controller;
class Project;

// This class does not throw exceptions.
class ProcessingNode : public QObject {
	Q_OBJECT
public:
	ProcessingNode(Controller& controller, Project& project);
Q_SIGNALS:
	void processingComplete();
	void error();
public Q_SLOTS:
	void process();
private:
	Project& project_;
};

} // namespace Lab

#endif /* PROCESSINGNODE_H_ */
