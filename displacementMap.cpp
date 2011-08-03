/*
 * displacementMap.cpp
 *
 *  Created on: 3 aug. 2011
 *      Author: StaleyJW
 */

#include "vector3d.h"
#include "displacementMap.h"


DisplacementMap::DisplacementMap()
{

	initCommon();
}


DisplacementMap::DisplacementMap(std::string filename)
{

	initCommon();
}


// Loads a text file containing discrete displacement values at a given time step
// that were obtained from kWave simulation post-processed data.
void DisplacementMap::loadDisplacementMap(const int timeStep)
{

}


void DisplacementMap::initCommon()
{

}


// Returns a Vector3d object holding values for displacements in all axes.
Vector3d DisplacementMap::getDisplacements(const Vector3d &photonLocation)
{
	Vector3d result;

	return result;
}
