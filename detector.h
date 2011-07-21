//
//  detector.h
//  Xcode
//
//  Created by jacob on 7/18/11.
//  Copyright 2011 BMPI. All rights reserved.
//

#ifndef DETECTOR_H
#define DETECTOR_H

#include "vector3D.h"
#include "vectorMath.h"
#include "logger.h"
using namespace VectorMath;

class Detector 
{
public:
    Detector(const Vector3d &centerPoint);
    Detector(const boost::shared_ptr<Vector3d> centerPoint);
    virtual ~Detector();
        
    virtual bool photonPassedThroughDetector(const boost::shared_ptr<Vector3d> p0,
                                             const boost::shared_ptr<Vector3d> p1) = 0;
    
    virtual void savePhotonExitCoordinates(const boost::shared_ptr<Vector3d> exitCoords) = 0;
    virtual void savePhotonExitWeight(void) = 0;
    
    virtual void setDetectorPlaneXY(void)
    {
        // Set which plane the detector resides.
        xz_plane = false;
        yz_plane = false;
        xy_plane = true;
        
        // Set the direction that the vector that is normal to the plane.
        normalVector.setDirX(0.0f);
        normalVector.setDirY(0.0f);
        normalVector.setDirZ(1.0f); normalVector.location.z = 1.0f;
        
    }
    
    virtual void setDetectorPlaneXZ(void)
    {
        // Set which plane the detector resides.
        yz_plane = false;
        xy_plane = false;
        xz_plane = true;
        
        // Set the direction that the vector that is normal to the plane.
        normalVector.setDirX(0.0f);
        normalVector.setDirY(1.0f); normalVector.location.y = 1.0f;
        normalVector.setDirZ(0.0f);
    }
    
    virtual void setDetectorPlaneYZ(void)
    {
        // Set which plane the detector resides.
        xz_plane = false;
        xy_plane = false;
        yz_plane = true;
        
        // Set the direction that the vector that is normal to the plane.
        normalVector.setDirX(1.0f); normalVector.location.x = 1.0f;
        normalVector.setDirY(0.0f);
        normalVector.setDirZ(0.0f);
    }
    
    
    
protected:
    // Center coordinates of the detector in the medium. [millimeters]
    Vector3d center;
    
    // Vector that is normal to the plane.
    Vector3d normalVector;
    
    // possible planes that the detector can be placed in 3D space.
    bool xy_plane;  
    bool xz_plane;
    bool yz_plane;
};


#endif