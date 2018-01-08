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
	virtual ~ProcessingNode();
private:
	Project& project_;
signals:
	void processingComplete();
	void error();
public slots:
	void process();
};

} // namespace Lab

#endif /* PROCESSINGNODE_H_ */
