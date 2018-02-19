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

#include "ProcessingNode.h"

#include <boost/scoped_ptr.hpp>

#include "Controller.h"
#include "Log.h"
#include "Method.h"
#include "Timer.h"



namespace Lab {

ProcessingNode::ProcessingNode(Controller& controller, Project& project)
		: project_(project)
{
	controller.connectProcessingNode(this);
}

ProcessingNode::~ProcessingNode()
{
}

void
ProcessingNode::process()
{
	LOG_INFO << "========== Processing started.";

	try {
		boost::scoped_ptr<Method> method(Method::get(project_));

		Timer procTimer;
		method->execute();
		LOG_INFO << "########## Processing time = " << procTimer.getTime() << " s";

	} catch (std::exception& e) {
		LOG_ERROR << "Exception caught during processing:\n" << e.what();
		emit error();
		return;
	} catch (...) {
		LOG_ERROR << "Unknown exception caught during processing.";
		emit error();
		return;
	}

	LOG_INFO << "========== Processing complete.";
	emit processingComplete();
}

} // namespace Lab
