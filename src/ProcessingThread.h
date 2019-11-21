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

#ifndef PROCESSINGTHREAD_H_
#define PROCESSINGTHREAD_H_

#include <QThread>

namespace Lab {

class Controller;
class Project;

class ProcessingThread : public QThread {
	Q_OBJECT
public:
	ProcessingThread(Controller& controller, Project& sharedData);
protected:
	virtual void run();
private:
	Controller& controller_;
	Project& project_;
};

} // namespace Lab

#endif /* PROCESSINGTHREAD_H_ */
