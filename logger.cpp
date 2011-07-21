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
    m_stream.close();
}

Logger * Logger::getInstance(void)
{
    if (!pInstance)
    {
        pInstance = new Logger();
    }
    return pInstance;
}


void Logger::openFile(std::string filename)
{
    m_stream.open(filename.c_str());
}

void Logger::write(double val)
{
    // Grab the lock to ensure that the logger doesn't get interrupted by a thread
    // in the middle of a write, causing the output to be corrupted.
    boost::mutex::scoped_lock lock(m_mutex);

    m_stream << "val = " << val << endl;
}

void Logger::write(const boost::shared_ptr<Vector3d> vectorCoords)
{
    // Grab the lock to ensure that the logger doesn't get interrupted by a thread
    // in the middle of a write, causing the output to be corrupted.
    boost::mutex::scoped_lock lock(m_mutex);   
    
    m_stream << vectorCoords;
    m_stream.flush();
}