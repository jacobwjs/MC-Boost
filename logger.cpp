//
//  logger.cpp
//  Xcode
//
//  Created by jacob on 7/18/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

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
    exit_data_stream.open(filename.c_str());
}

void Logger::openAbsorberFile(std::string filename)
{
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


void Logger::writeAbsorberData(const double absorbedWeight)
{
    absorber_data_stream << absorbedWeight << "\n";
    absorber_data_stream.flush();
}


