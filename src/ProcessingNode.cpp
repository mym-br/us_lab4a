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
