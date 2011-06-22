//
//  cylinder.h
//  Xcode
//
//  Created by jacob on 6/20/11.
//

#ifndef CYLINDER_H
#define CYLINDER_H

#include "absorber.h"
#include "coordinates.h"

class Cylinder : public Absorber 
{
public:
    Cylinder();
    Cylinder(const double height, const double radius, const coords &center); 
    ~Cylinder();
    
    virtual bool inAbsorber(const coords &c);
    void    cartesianToCylindrical(void);
    
private:
    // The height and radius of the cylindrical absorber.
    double height;
    double radius;
};

#endif