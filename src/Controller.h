#ifndef CONTROLLER_H_
#define CONTROLLER_H_

#include <QObject>

#include "ProcessingNode.h"
#include "ProcessingThread.h"



namespace Lab {

class Project;

class Controller : public QObject {
	Q_OBJECT
public:
	Controller(Project& project);
	virtual ~Controller();

	void connectProcessingNode(const ProcessingNode* node);
	void exit();
	void enableProcessing();
	bool processingEnabled() const { return (state_ == STATE_PROCESSING_ENABLED); }
private:
	enum State {
		STATE_PROCESSING_STOPPED,
		STATE_PROCESSING_ENABLED
	};

	static void showTimestamp();

	State state_; // this state can only be read/changed by the main thread
	ProcessingThread processingThread_;
signals:
	void processingFinished();
	void processingRequested();

private slots:
	//void execAfterAcquisition(Acquisition::AcquisitionDataType* acqData, int step, int numberOfSteps);
	//void execAfterSleep();
	void execAfterProcessing();
	void execAfterError();
	//void execAfterAcquisitionTermination();
	//void handleAcquisitionError(AcquisitionNode::Function f);
	//void handleProcessingError(ProcessingNode::Function f);
	//void executeDelayedStateChanges();
};

} // namespace Lab

#endif /* CONTROLLER_H_ */
