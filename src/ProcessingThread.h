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
	virtual ~ProcessingThread();
protected:
	virtual void run();
private:
	Controller& controller_;
	Project& project_;
};

} // namespace Lab

#endif /* PROCESSINGTHREAD_H_ */
