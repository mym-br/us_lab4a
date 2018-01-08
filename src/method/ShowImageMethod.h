#ifndef SHOWIMAGEMETHOD_H
#define SHOWIMAGEMETHOD_H

#include "Method.h"



namespace Lab {

class Project;

class ShowImageMethod : public Method {
public:
	ShowImageMethod(Project& project);
	virtual ~ShowImageMethod();

	virtual void execute();

private:
	Project& project_;
};

} // namespace Lab

#endif // SHOWIMAGEMETHOD_H
