//
//  vectorMath.h
//  Xcode
//
//  Created by jacob on 7/11/11.
//  Copyright 2011 BMPI. All rights reserved.
//

#ifndef VECTOR3D_H
#define VECTOR3D_H

#include "coordinates.h"
#include "boost/smart_ptr.hpp"



class Vector3d
{
    
    
public:
    Vector3d();
    Vector3d(double x, double y, double z);
    Vector3d(const coords &location, const directionCos &dir);
    ~Vector3d();
    
    // Initializes the direction cosine private data.  This allows general use
    // of this vector class for objects that have position or position and direction.
    void withDirection(void) {direction = boost::shared_ptr<directionCos> (new directionCos);}
    
    
    // Overloaded operators for working with Vector3d.
    boost::shared_ptr<Vector3d> operator-(Vector3d &rhs);
    boost::shared_ptr<Vector3d> operator+(Vector3d &rhs);
    boost::shared_ptr<Vector3d> operator*(double num);
    boost::shared_ptr<Vector3d> operator*(Vector3d &rhs);
    
    // Get the coordinates.
    double getDirX(void);
    double getDirY(void);
    double getDirZ(void);
    
    // Set the coordinates.
    void setDirX(double x);
    void setDirY(double y);
    void setDirZ(double z);
    
    
    // Making these public to ease manipulation in vector math.
    // No longer need the getters and setters.
    coords location;            // Structure containing the cartesian coordinates of the photon.
    //coords previousLocation;    // Holds the coordinates of the previous location of the photon.
    //directionCos direction;      // Used in checking if the photon had passed through an absorber.
    
    
private:
    // This is private since not every object that needs location needs direction.
    boost::shared_ptr<directionCos> direction;     // The direction cosines of the photon.
    
};





#endif