//
//  vectorMath.cpp
//  Xcode
//
//  Created by jacob on 7/11/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#include "vector3D.h"



Vector3d::Vector3d()
{
    location.x = location.y = location.z = 0.0f;
    //previousLocation.x = previousLocation.y = previousLocation. z = 0.0f;
    //direction.x = direction.y = direction.z = 0.0f;
}

Vector3d::Vector3d(double x, double y, double z)
{
    location.x = x;
    location.y = y;
    location.z = z;
    
    //direction.x = direction.y = direction.z = 0.0f;
}

Vector3d::Vector3d(const coords &c, const directionCos &dir)
{
    location.x = c.x;
    location.y = c.y;
    location.z = c.z;
    
    //previousLocation.x = previousLocation.y = previousLocation.z = 0.0f;
    
    // Initialize the direction cosine structure.
    withDirection();
    assert(direction.get() != NULL);
    direction->x = dir.x;
    direction->y = dir.y;
    direction->z = dir.z;
}


Vector3d::~Vector3d()
{
    
}

// Here we overload the - operator so we can subtract vectors 
boost::shared_ptr<Vector3d> Vector3d::operator-(Vector3d &rhs)
{
    boost::shared_ptr<Vector3d> result(new Vector3d);

    result->location.x = (this->location.x) - (rhs.location.x);
    result->location.y = (this->location.y) - (rhs.location.y);
    result->location.z = (this->location.z) - (rhs.location.z);
    
    return result;
}

// Here we overload the + operator so we can add vectors together 
boost::shared_ptr<Vector3d> Vector3d::operator+(Vector3d &rhs)
{
    boost::shared_ptr<Vector3d> result(new Vector3d);
    
    result->location.x = (this->location.x) + (rhs.location.x);
    result->location.y = (this->location.y) + (rhs.location.y);
    result->location.z = (this->location.z) + (rhs.location.z);
    
    return result;
}

boost::shared_ptr<Vector3d> Vector3d::operator*(double num)
{
    return boost::shared_ptr<Vector3d> (new Vector3d(location.x * num, 
                                                     location.y * num,
                                                     location.z * num));
}


// XXX: Is this correct?
boost::shared_ptr<Vector3d> Vector3d::operator*(Vector3d &rhs)
{
    boost::shared_ptr<Vector3d> result(new Vector3d);
    
    result->location.x = (this->location.x) * (rhs.location.x);
    result->location.y = (this->location.y) * (rhs.location.y);
    result->location.z = (this->location.z) * (rhs.location.z); 
    
    return result;
}


double Vector3d::getDirX(void)
{
    assert(direction.get() != NULL);
    return direction->x;
}


double Vector3d::getDirY(void)
{
    assert(direction.get() != NULL);
    return direction->y;
}


double Vector3d::getDirZ(void)
{
    assert(direction.get() != NULL);
    return direction->z;
}


void Vector3d::setDirX(double x)
{
    assert(direction.get() != NULL);
    direction->x = x;
}


void Vector3d::setDirY(double y)
{
    assert(direction.get() != NULL);
    direction->y = y;
}


void Vector3d::setDirZ(double z)
{
    assert(direction.get() != NULL);
    direction->z = z;
}

