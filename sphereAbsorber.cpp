//
//  sphere.cpp
//  Xcode
//
//  Created by jacob on 6/20/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#include "sphere.h"
#include <cmath>
#include <iostream>
using std::cout;


Sphere::Sphere(const double r, const coords &center)
:Absorber(center)
{
    radius = r;
}


Sphere::~Sphere()
{
    // STUB
}


// FIXME: Need to check if the minimum distance from the line formed by the path
//        to the center of the absorber is within the radius bounds of the absorber.
//        If that is the case then the photon has moved through the absorber.
//        Based on the path length through the absorber a certain amount of 
//        energy needs to be dropped based on the local mu_a from the absorber.
bool Sphere::inAbsorber(const coords &photonCoords)
{
    // Convert the cartesian coordinates into a spherical coordinate structure.
    // sphereCoords scoords = cartesianToSpherical(center);
    
    // If the distance from the center of the absorber to the location of the photon is
    // larger than the radius of the absorber we know the photon hasn't crossed into the
    cout << "Checking if photon is in absorber\n";
    
}


// FIXME: Verify this is correct.
sphereCoords Sphere::cartesianToSpherical(const coords &center)
{
    sphereCoords temp;
    temp.r      = radius;
    temp.theta  = acos(center.z / temp.r);
    temp.phi    = atan2(center.y, center.x);
}