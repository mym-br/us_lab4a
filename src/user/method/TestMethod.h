/***************************************************************************
 *  Copyright 2014, 2017, 2018, 2019 Marcelo Y. Matuda                     *
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

#ifndef TESTMETHOD_H_
#define TESTMETHOD_H_

#include "Exception.h"
#include "Method.h"



namespace Lab {

struct TestException : std::runtime_error {
	using std::runtime_error::runtime_error;
};

class Project;

class TestMethod : public Method {
public:
	explicit TestMethod(Project& project);
	virtual ~TestMethod() = default;

	virtual void execute();
private:
	TestMethod(const TestMethod&) = delete;
	TestMethod& operator=(const TestMethod&) = delete;
	TestMethod(TestMethod&&) = delete;
	TestMethod& operator=(TestMethod&&) = delete;

	template<typename T> void call(T pmf);

	void testAdd();
	void testAddElements();
	void testAddElements2();
	void testBessel();
	void testCentralDiff();
	void testDecimator();
	void testDirectFFTWFilter();
	void testFFT();
	void testFFTWFilter();
	void testFFTWFilter2();
	void testFillSequence();
	void testFillSequence2();
	void testFillSequence3();
	void testFillSequence4();
	void testHilbertTransform();
	void testInterpolator();
	void testInterpolator4X();
	void testKaiserWindow();
	void testLinearInterpolator();
	void testMatrix();
	void testMultiplyBy();
	void testSCF();
	void testStatistics();
	void testTensor3();

	Project& project_;
	unsigned int errorCount_;
	unsigned int figureNumber_;
};

} // namespace Lab

#endif /* TESTMETHOD_H_ */
