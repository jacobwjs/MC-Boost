//
//  absorber.h
//  Xcode
//
//  Created by jacob on 6/20/11.
//
#ifndef ABSORBER_H
#define ABSORBER_H

#include "coordinates.h"

class Absorber 
{
public:
    Absorber(const coords &center);
    ~Absorber();
    
    virtual bool inAbsorber(const coords &c) = 0;
    
    
    
protected:
    // The optical properties of the absorber.
    double mu_a, mu_s, anisotropy;
    
    // The coordinates of the center point of the absorber in the medium.
    coords center;
};

#endif