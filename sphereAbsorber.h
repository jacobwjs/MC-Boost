//
//  sphere.h
//  Xcode
//
//  Created by jacob on 6/20/11.
//

#ifndef SPHERE_H
#define SPHERE_H

#include "absorber.h"
#include "coordinates.h"

class Sphere : public Absorber
{
public:
    Sphere();
    Sphere(const double r, const coords &c);
    ~Sphere();
    
    virtual bool inAbsorber(const coords &c);
    sphereCoords    cartesianToSpherical(const coords &c);
    
    
private:
    double radius;
    
};

#endif // SPHERE_H