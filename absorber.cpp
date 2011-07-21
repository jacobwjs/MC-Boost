//
//  absorber.cpp
//  Xcode
//
//  Created by jacob on 6/20/11.
//

#include "absorber.h"




Absorber::Absorber(const coords &c)
{
    center = boost::shared_ptr<Vector3d> (new Vector3d);
    center->location.x = c.x;
    center->location.y = c.y;
    center->location.z = c.z;
}

Absorber::Absorber(const boost::shared_ptr<Vector3d> c)
{
    center = boost::shared_ptr<Vector3d> (new Vector3d);
    center->location.x = c->location.x;
    center->location.y = c->location.y;
    center->location.z = c->location.z;
}

Absorber::~Absorber()
{
    // STUB
}

