//
//  logger.cpp
//  Xcode
//
//  Created by jacob on 7/18/11.
//  Copyright 2011 BMPI. All rights reserved.
//

#include "vector3D.h"
#include "photon.h"
#include "logger.h"




Logger * Logger::pInstance = 0;




Logger::Logger()
{
        //cout << "logger dying";
}


Logger::~Logger()
{
    exit_data_stream.close();
    absorber_data_stream.close();
}

Logger * Logger::getInstance(void)
{
    if (!pInstance)
    {
        pInstance = new Logger();
    }
    return pInstance;
}


void Logger::openExitFile(std::string filename)
{
    // Ensure file stream is not already open.
    if (exit_data_stream.is_open())
        exit_data_stream.close();
    
    exit_data_stream.open(filename.c_str());
}

void Logger::openAbsorberFile(std::string filename)
{
    // Ensure file stream is not already open.
    if (absorber_data_stream.is_open())
        absorber_data_stream.close();
    
    absorber_data_stream.open(filename.c_str());
}


void Logger::write(double val)
{
    // Grab the lock to ensure that the logger doesn't get interrupted by a thread
    // in the middle of a write, causing the output to be corrupted.
    boost::mutex::scoped_lock lock(m_mutex);

    exit_data_stream << "val = " << val << endl;
}

void Logger::writeExitData(const boost::shared_ptr<Vector3d> photonVector)
{
    // Grab the lock to ensure that the logger doesn't get interrupted by a thread
    // in the middle of a write, causing the output to be corrupted.
    boost::mutex::scoped_lock lock(m_mutex);   
    
    exit_data_stream << photonVector;
    exit_data_stream.flush();
}


void Logger::writeExitData(const boost::shared_ptr<Vector3d> photonVector,
                   const double weight,
                   bool tagged)
{
    // Grab the lock to ensure that the logger doesn't get interrupted by a thread
    // in the middle of a write, causing the output to be corrupted.
    boost::mutex::scoped_lock lock(m_mutex);
    
    // Write out the location (x,y,z), exit angle (theta), weight of photon, and whether it was
    // tagged.
    exit_data_stream << tagged << "," 
                     << weight << ","
                     << photonVector->getDirZ() << ","
                     << photonVector << "\n";
    
    exit_data_stream.flush();
}



void Logger::writeExitData(const boost::shared_ptr<Vector3d> photonVector,
                           const double weight)
{
    // Grab the lock to ensure that the logger doesn't get interrupted by a thread
    // in the middle of a write, causing the output to be corrupted.
    boost::mutex::scoped_lock lock(m_mutex);
    
    // Write out the location (x,y,z), exit angle (theta), weight of photon, and whether it was
    // tagged.
    exit_data_stream << weight << "," 
                     << photonVector->getDirZ() << "," 
                     << photonVector << "\n";
    
}



void Logger::writeExitData(const boost::shared_ptr<Vector3d> photonVector,
                           const double exitWeight,
                           const double transmissionAngle
                           )
{
    // Grab the lock to ensure that the logger doesn't get interrupted by a thread
    // in the middle of a write, causing the output to be corrupted.
    boost::mutex::scoped_lock lock(m_mutex);
    
    // Write out the location (x,y,z), transmission angle (theta), weight of photon
    exit_data_stream << exitWeight << "," 
                     << transmissionAngle << "," 
                     << photonVector << "\n";
    
    exit_data_stream.flush();
}



void Logger::writeWeightAngleLengthCoords(const double exitWeight,
                                          const double transmissionAngle,
                                          const double modulatedPathLength,
                                          const boost::shared_ptr<Vector3d> photonVector)
{
    // Grab the lock to ensure that the logger doesn't get interrupted by a thread
    // in the middle of a write, causing the output to be corrupted.
    boost::mutex::scoped_lock lock(m_mutex);
    
    // Write out the location (x,y,z), transmission angle (theta), weight of photon
    exit_data_stream << exitWeight << "," 
                     << transmissionAngle << ","
                     << modulatedPathLength << ","
                     << photonVector << "\n";
    
    exit_data_stream.flush();
    
}

                                          
                                  
void Logger::writePhoton(Photon *p)
{
    // Grab the lock to ensure that the logger doesn't get interrupted by a thread
    // in the middle of a write, causing the output to be corrupted.
    boost::mutex::scoped_lock lock(m_mutex);
    
    cout << "Logger::writePhoton() stub\n";
}


void Logger::writeAbsorberData(const double absorbedWeight)
{
    absorber_data_stream << absorbedWeight << "\n";
    absorber_data_stream.flush();
}


