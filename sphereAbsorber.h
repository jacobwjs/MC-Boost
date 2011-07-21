//
//  sphere.h
//  Xcode
//
//  Created by jacob on 6/20/11.
//

#ifndef SPHEREABSORBER_H
#define SPHEREABSORBER_H

#include "absorber.h"


class SphereAbsorber : public Absorber
{
public:
    SphereAbsorber();
    SphereAbsorber(const double radius, const coords &c);
    SphereAbsorber(const double radius, const boost::shared_ptr<Vector3d> c);
    ~SphereAbsorber();
    
    virtual bool hitAbsorberBoundary(const boost::shared_ptr<Vector3d> photonVector);
    virtual bool inAbsorber(const boost::shared_ptr<Vector3d> photonVector);
    virtual bool crossedAbsorber(const boost::shared_ptr<Vector3d> A,
                                 const boost::shared_ptr<Vector3d> B);
    
    virtual double getAbsorberAbsorptionCoeff(void) {return this->mu_a;}
    virtual double getAbsorberScatteringCoeff(void) {return this->mu_s;}

    
    // Check if photon is within the radius of the absorber.
    bool inSphereVolume(const boost::shared_ptr<Vector3d> photonVector);
    sphereCoords    cartesianToSpherical(const coords &c);
    
    
private:
    double radius;
    
};

#endif // SPHERE_H