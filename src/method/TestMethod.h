#ifndef TESTMETHOD_H_
#define TESTMETHOD_H_

#include "Exception.h"
#include "Method.h"



namespace Lab {

struct TestException : public virtual Exception {};

class Project;

class TestMethod : public Method {
public:
	TestMethod(Project& project);
	virtual ~TestMethod();

	virtual void execute();

private:
	Project& project_;
};

} // namespace Lab

#endif /* TESTMETHOD_H_ */
