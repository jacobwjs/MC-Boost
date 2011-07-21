//
//  CylinderAbsorber.h
//  Xcode
//
//  Created by jacob on 6/20/11.
//

#ifndef CYLINDERABSORBER_H
#define CYLINDERABSORBER_H

#include "absorber.h"
#include <cmath>

class CylinderAbsorber : public Absorber 
{
public:
    CylinderAbsorber(const double radius, const coords &top, const coords &bottom); 
    ~CylinderAbsorber();
    
    virtual bool hitAbsorberBoundary(const boost::shared_ptr<Vector3d> photonVector);
    virtual bool inAbsorber(const boost::shared_ptr<Vector3d> photonVector);
    virtual bool crossedAbsorber(const boost::shared_ptr<Vector3d> photonVector);
    
    
    virtual double getAbsorberAbsorptionCoeff(void) {return this->mu_a;}
    virtual double getAbsorberScatteringCoeff(void) {return this->mu_s;}
    
    void    cartesianToCylindrical(void);
    
private:
    // The height and radius of the cylindrical absorber.
    double length;
    double radius;
    
    // Cartesian coordinates of the center location of one end
    // of the cyclinder.
    coords topCenter;       
    
    // Cartesian coordinates of center location of the other end
    // of the cyclinder.
    coords bottomCenter;
};

#endif