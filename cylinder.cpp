//
//  cylinder.cpp
//  Xcode
//
//  Created by jacob on 6/20/11.
//

#include "cylinder.h"



Cylinder::Cylinder(const double height, const double radius, const coords &center)
:Absorber(center)
{
    this->height = height;
    this->radius = radius;
}


Cylinder::~Cylinder()
{
    // STUB
}


bool Cylinder::inAbsorber(const coords &center)
{
    // STUB
}


void Cylinder::cartesianToCylindrical(void)
{
    // STUB
}