/*
 * loadPressureMap.cpp
 *
 *  Created on: 18 mei 2011
 *      Author: StaleyJW
 */

#include "pressureMap.h"

PressureMap::PressureMap()
{
	// FIXME:  THESE SHOULD NOT BE HARD CODED.
	Nx = Nz = Ny = 64;
	dx = dz = dy = 0.16; // [cm] (note: 100e-3/Nx in centimeters is 0.16;
	pgrid = new three_dim_array (boost::extents[Nx][Nz][Ny]);
}


PressureMap::~PressureMap()
{
	delete pgrid;
}


void PressureMap::loadPressureFromKWave(void)
{

	ifstream pressure_file;
	pressure_file.open("pressure-at-25us.txt");
	if (!pressure_file) {
		cout << "!!! Error opening pressure map file !!!\n";
		exit(1);
	}
	else {
		cout << "Pressure map opened successfully\n";
	}

	double data = 0.0;

	for (array_index x = 0; x < Nx && pressure_file.good(); x++)
		for (array_index z = 0; z < Nz; z++)
			for (array_index y = 0; y < Ny; y++)
			{
				pressure_file >> data;
				(*pgrid)[x][z][y] = data;
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
double PressureMap::getPressureFromGrid(int x, int z, int y)
{
//	array_index x = x_location;
//	array_index z = z_location;
//	array_index y = y_location;
//	return (*pgrid)[x][z][y];
	cout << "PressureMap::getPressureFromGrid\n";
	cout << "x=" << x << ", y =" << y << ", z=" << z << endl;
	return (*pgrid)[(array_index)x][(array_index)z][(array_index)y];
}

// FIXME: NEED TO DO ERROR CHECKING TO ENSURE BOUNDS OF THE GRID ARE RESPECTED.
// Returns the pressure from the grid based on supplied coordinates.
double PressureMap::getPressureCartCords(double x, double z, double y)
{

	int _x = floor(x/dx);
	int _z = floor(z/dz);
	int _y = floor(y/dy);
	cout << "PressureMap::getPressureCartCords\n";
	return getPressureFromGrid(_x, _z, _y);

}



