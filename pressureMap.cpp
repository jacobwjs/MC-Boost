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
}


PressureMap::PressureMap(const string &filename, const int Nx, const int Nz, const int Ny, const int dimension)
{
	// Assign the number of grid points (pixels in k-wave) used in the simulation.
	this->Nx = Nx;
	this->Ny = Ny;
	this->Nz = Nz;

	x_bound = y_bound = z_bound = dimension;

	// FIXME: SHOULD USE BOOST FILESYSTEM TO GET initial_path AND DECIDE
	//        WHAT FILE TO LOAD BASED ON CURRENT WORKING DIRECTORY.
	// Assign the name to the pressure file.
	p_file = filename;

	// Initialize the data structures and values for the pressure map.
	this->init();
}


void PressureMap::init(void)
{

	dx = x_bound/Nx; // [cm] (note: 20e-3/Nx in centimeters is 0.16;
	dy = y_bound/Ny;
	dz = z_bound/Nz;
	pgrid = new three_dim_array (boost::extents[Nx][Nz][Ny]);
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
//	cout << "PressureMap::getPressureFromGrid\n";
//	cout << "a=" << a << ", b=" << b << ", c =" << c  << endl;
	return (*pgrid)[(array_index)a][(array_index)b][(array_index)c];
}

// FIXME: NEED TO DO ERROR CHECKING TO ENSURE BOUNDS OF THE GRID ARE RESPECTED.
// Returns the pressure from the grid based on supplied coordinates.
double PressureMap::getPressureCartCords(double a, double b, double c)
{
	int _x = floor(a/dx);
	int _z = floor(b/dz);
	int _y = floor(c/dy);

	// Sanity check.
	assert((_x <= Nx && _x >= 0) &&
			(_y <= Ny && _y >= 0) &&
			(_z <= Nz && _z >= 0));

//	cout << "PressureMap::getPressureCartCords\n";
	//cout << "a=" << _x << ", b=" << _z << ", c=" << _y << endl;
	return getPressureFromGrid(_x, _z, _y);

}



