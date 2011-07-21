//
//  absorber.h
//  Xcode
//
//  Created by jacob on 6/20/11.
//
#ifndef ABSORBER_H
#define ABSORBER_H

#include "vector3D.h"
#include "vectorMath.h"
using namespace VectorMath;

class Absorber 
{
public:
    Absorber(const coords &location);
    Absorber(const boost::shared_ptr<Vector3d> location);
    virtual ~Absorber();
    
    virtual bool hitAbsorberBoundary(const boost::shared_ptr<Vector3d> photonVector) = 0;
    virtual bool inAbsorber(const boost::shared_ptr<Vector3d> photonVector) = 0;
    virtual bool crossedAbsorber(const boost::shared_ptr<Vector3d> currPoint, 
                                 const boost::shared_ptr<Vector3d> prevPoint) = 0;
    
    virtual double getAbsorberAbsorptionCoeff(void) = 0;
    virtual double getAbsorberScatteringCoeff(void) = 0;

    
    
protected:
    // The optical properties of the absorber.
    double mu_a, mu_s, anisotropy;
    
    // The coordinates of the center point of the absorber in the medium.
    boost::shared_ptr<Vector3d> center;
};

#endif