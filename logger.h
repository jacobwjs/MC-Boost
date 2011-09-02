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


#include <fstream>
using std::ofstream;
#include <iostream>
using std::cout;
using std::endl;
#include <string>
#include <boost/thread/mutex.hpp>


// Forward decleration of objects.
class Photon;
class Vector3d;




class Logger 
{
public:
    static Logger * getInstance(void);
    
    void openExitFile(std::string filename);
    void openAbsorberFile(std::string filename);
    
    void write(double val);
    void writeExitData(const boost::shared_ptr<Vector3d> vectorCoords);
    void writeExitData(const boost::shared_ptr<Vector3d> vectorCoords,
                       const double weight,
                       bool tagged);
    void writeExitData(const boost::shared_ptr<Vector3d> photonVector,
                       const double weight);
    void writeExitData(const boost::shared_ptr<Vector3d> photonVector,
                       const double weight,
                       const double transmissionAngle);
    
    
    void writeAbsorberData(const double absorbedWeight);
    
    // XXX: Finish me
    void writeAbsorberData(const double absorbedWeight,
                           const double theta,
                           const double phi);

    // Writes the photon's weight, transmission angle, modulated path length through the medium,
    // and its exit location on the exit detector window.
    void writeWeightAngleLengthCoords(const double exitWeight,
                                      const double transmissionAngle,
                                      const double modulatedPathLength,
                                      const boost::shared_ptr<Vector3d> photonVector);
    void writeWeightAngleLengthCoords(const double exitWeight,
    								  const double dirx,
    								  const double diry,
    								  const double dirz,
    								  const double modulatedPathLength,
    								  const boost::shared_ptr<Vector3d> photonVector);
    void writeWeightAngleLengthCoords(Photon &p);


    // XXX:
    // - Does this introduce race conditions by pointing to a threaded object that could
    //   potentially have data changing in obscure ways?  Unsure, but each object is given
    //   it's own CPU "core" to run on, which means any object's state between switches (which
    //   there should be none (i.e. context switching) in a perfect world since threads == cores)
    //   should be coherent.  
    void writePhoton(Photon *p);
    
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
    
    boost::mutex m_mutex;
};

#endif
