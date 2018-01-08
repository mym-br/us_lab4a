#include "ProcessingThread.h"

#include <boost/scoped_ptr.hpp>

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
	LOG_DEBUG << "Processing thread starting...";

	boost::scoped_ptr<ProcessingNode> processingNode(new ProcessingNode(controller_, project_));

	exec();
}

} // namespace Lab
