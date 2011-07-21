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

void Logger::write(int val)
{
    m_stream << "val = " << val << endl;
    
}