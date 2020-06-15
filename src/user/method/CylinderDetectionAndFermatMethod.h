/***************************************************************************
 *  Copyright 2014, 2017, 2018, 2020 Marcelo Y. Matuda                     *
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

#ifndef CYLINDERDETECTIONANDFERMATMETHOD_H_
#define CYLINDERDETECTIONANDFERMATMETHOD_H_

#include "Method.h"



namespace Lab {

class Project;

template<typename TFloat>
class CylinderDetectionAndFermatMethod : public Method {
public:
	CylinderDetectionAndFermatMethod(Project& project) : project_(project) { }
	virtual ~CylinderDetectionAndFermatMethod() = default;

	virtual void execute();
private:
	CylinderDetectionAndFermatMethod(const CylinderDetectionAndFermatMethod&) = delete;
	CylinderDetectionAndFermatMethod& operator=(const CylinderDetectionAndFermatMethod&) = delete;
	CylinderDetectionAndFermatMethod(CylinderDetectionAndFermatMethod&&) = delete;
	CylinderDetectionAndFermatMethod& operator=(CylinderDetectionAndFermatMethod&&) = delete;

	void acquireSignals();
	// Detects points using the peaks in STA images.
	void detectPoints();
	// Detects points using the peaks in STA images. The images are formed using sums of analytic signals.
	void detectPointsUsingVectorSum();
	// Detects points using the peaks in STA images. The images are formed using sums of the results of cross-correlations.
	void detectPointsUsingCrossCorrelation();
	// Detects points using geometry. The distances are calculated using one element at a time, as emitter and receiver.
	void detectPointsUsingCCBFPulseEcho();
	// Detects points using geometry. One element emits and others receive.
	void detectPointsUsingCCBFPitchCatch();
	// Detects points using arcs instead of images.
	void detectPointsInArcs();
	// Detects points using only geometry (envelope of arc and ellipses).
	void detectPointsUsingTangentCurveGeometry();
	void fitCircle();
	void execTwoMediumImaging();
	void execCombinedTwoMediumImaging();
	template<typename Proc> void execCombinedTwoMediumImagingCyl();

	void measureSpeed1();
	void measureSpeed1AndDistanceError();
	void measureSpeed1AndDistanceError2();
	void measureDistance();

	Project& project_;
};

} // namespace Lab

#endif /* CYLINDERDETECTIONANDFERMATMETHOD_H_ */
