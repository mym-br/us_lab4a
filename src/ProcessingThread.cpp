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

#include "ProcessingThread.h"

#include <memory>

#include "Log.h"
#include "ProcessingNode.h"



namespace Lab {

ProcessingThread::ProcessingThread(Controller& controller, Project& project)
		: controller_(controller)
		, project_(project)
{
}

ProcessingThread::~ProcessingThread()
{
}

void
ProcessingThread::run()
{
	//LOG_DEBUG << "Processing thread running...";

	auto processingNode = std::make_unique<ProcessingNode>(controller_, project_);

	exec();
}

} // namespace Lab
