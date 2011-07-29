//
//  logger.h
//  Xcode
//
//  Created by jacob on 7/18/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//


// Logger singleton.
// NOTE:  Construction is NOT thread-safe, must be initialized in main before any threads are spawned.
#ifndef LOGGER_H
#define LOGGER_H

#include "vector3D.h"
//#include "photon.h"
#include <fstream>
using std::ofstream;
#include <iostream>
using std::cout;
using std::endl;
#include <string>
#include <boost/thread/mutex.hpp>

class Photon;

class Logger 
{
public:
    static Logger * getInstance(void);
    
    void openExitFile(std::string filename);
    void openAbsorberFile(std::string filename);
    void openMetaData(std::string filename);
    
    
    void write(double val);
    void writeMetaData(const double absorberRadius, const double detectorRadius, 
                       const int Nphotons, const double detectorPlane, const Vector3d &absorberLocation);
    void writeExitData(const boost::shared_ptr<Vector3d> vectorCoords);
    void writeExitData(const boost::shared_ptr<Vector3d> vectorCoords,
                       const double weight,
                       bool tagged);
    void writeExitData(const boost::shared_ptr<Vector3d> photonVector,
                       const double weight);
    
    
    void writeAbsorberData(const double absorbedWeight);
    
    // XXX: Finish me
    void writeAbsorberData(const double absorbedWeight,
                           const double theta,
                           const double phi);

    
private:
    Logger();                            // default constructor is private
    Logger(Logger const&){};             // copy constructor is private
    ~Logger();
    
    Logger& operator=(Logger const&){};  // assignment operator is private
    
    static Logger * pInstance;
    
    // The output streams associated with data for the photon and data for
    // the absorbers.
    ofstream exit_data_stream;
    ofstream absorber_data_stream;
    ofstream meta_data_stream;
    
    boost::mutex m_mutex;
};

#endif