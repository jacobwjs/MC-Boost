/*
 * loadPressureMap.cpp
 *
 *  Created on: 18 mei 2011
 *      Author: StaleyJW
 */

#include "pressureMap.h"

PressureMap::PressureMap()
{
	// Initialize the data structures and values for the pressure map.
	this->init();

	// FIXME: SHOULD USE BOOST FILESYSTEM TO GET initial_path AND DECIDE
	//        WHAT FILE TO LOAD BASED ON CURRENT WORKING DIRECTORY.
	// Assign the name to the pressure file.
	p_file.append("pressure-at-25us.txt");
}


PressureMap::PressureMap(const char *filename)
{
	// Initialize the data structures and values for the pressure map.
	this->init();

	// FIXME: SHOULD USE BOOST FILESYSTEM TO GET initial_path AND DECIDE
	//        WHAT FILE TO LOAD BASED ON CURRENT WORKING DIRECTORY.
	// Assign the name to the pressure file.
	p_file.append(filename);

}

void PressureMap::init(void)
{
	// FIXME:  THESE SHOULD NOT BE HARD CODED.
	Nx = Nz = Ny = 64;
	dx = dz = dy = 0.16; // [cm] (note: 100e-3/Nx in centimeters is 0.16;
	pgrid = new three_dim_array (boost::extents[Nx][Nz][Ny]);

	// Default location of where pressure files are located.  That is,
	// the current working directory of the executable.
	// FIXME: SHOULD USE BOOST FILESYSTEM TO REMAIN PLATFORM AGNOSTIC.
	p_file = "C:/Users/StaleyJW/Desktop/Software/MC-Boost/kWave-pressure/";

}

PressureMap::~PressureMap()
{
	delete pgrid;
}


void PressureMap::loadPressureFromKWave(void)
{


	// Open the file that contains the pressure values.
	pressure_file.open(p_file.c_str());

	if (!pressure_file) {
		cout << "!!! Error opening pressure map file !!!\n";
		exit(1);
	}
	else {
		cout << "Pressure map " << p_file.c_str() << "...  opened successfully\n";
	}


	double data = 0.0;

	for (array_index a = 0; a < Nx && pressure_file.good(); a++)
		for (array_index b = 0; b < Nz; b++)
			for (array_index c = 0; c < Ny; c++)
			{
				pressure_file >> data;
				(*pgrid)[a][b][c] = data;
			}

//	// Dump loaded pressure to stdout to see values.
//	for (array_index x = 0; x < Nx; x++)
//			for (array_index z = 0; z < Nz; z++)
//				for (array_index y = 0; y < Ny; y++)
//				{
//					cout << (*pgrid)[x][z][y] << endl;
//				}


	pressure_file.close();
}


// FIXME: ENSURE INDICES ARE WITHIN THE DIMENSIONS OF THE GRID.
double PressureMap::getPressureFromGrid(int a, int b, int c)
{
//	array_index x = x_location;
//	array_index z = z_location;
//	array_index y = y_location;
//	return (*pgrid)[x][z][y];
	cout << "PressureMap::getPressureFromGrid\n";
	cout << "a=" << a << ", b=" << b << ", c =" << c  << endl;
	return (*pgrid)[(array_index)a][(array_index)b][(array_index)c];
}

// FIXME: NEED TO DO ERROR CHECKING TO ENSURE BOUNDS OF THE GRID ARE RESPECTED.
// Returns the pressure from the grid based on supplied coordinates.
double PressureMap::getPressureCartCords(double a, double b, double c)
{

	int _x = floor(a/dx);
	int _z = floor(b/dz);
	int _y = floor(c/dy);
	cout << "PressureMap::getPressureCartCords\n";
	return getPressureFromGrid(_x, _z, _y);

}



