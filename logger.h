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
#include <fstream>
using std::ofstream;
#include <iostream>
using std::cout;
using std::endl;
#include <string>
#include <boost/thread/mutex.hpp>


class Logger 
{
public:
    static Logger * getInstance(void);
    
    void openFile(std::string filename);
    void write(double val);
    void write(const boost::shared_ptr<Vector3d> vectorCoords);
    

    
private:
    Logger();                            // default constructor is private
    Logger(Logger const&){};             // copy constructor is private
    ~Logger();
    
    Logger& operator=(Logger const&){};  // assignment operator is private
    
    static Logger * pInstance;
    ofstream m_stream;
    
    boost::mutex m_mutex;
};

#endif