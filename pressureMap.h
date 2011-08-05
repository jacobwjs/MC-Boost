/*
 * loadPressureMap.h
 *
 *  Created on: 18 mei 2011
 *      Author: StaleyJW
 */

#ifndef LOADPRESSUREMAP_H_
#define LOADPRESSUREMAP_H_

#include <boost/multi_array.hpp>
#include <iostream>
using std::cout;
#include <fstream>
using std::fstream;
#include <string>



// Boost array that holds the pressure values after being
// loaded in from file.
typedef boost::multi_array<double, 3> three_dim_array;
typedef three_dim_array::index array_index;


// Forward decleration of Vector3d class.
class Vector3d;

class PressureMap
{
public:

	PressureMap(const std::string &filename, const int Nx, const int Nz, const int Ny, const int grid_size);
	PressureMap(const int Nx, const int Nz, const int Ny, const int grid_size);
	~PressureMap();

	void 	loadPressureMap(void);
	void	loadPressureMap(const std::string &filename, const int time_step);

	double 	getPressureFromGrid(int x, int z, int y);
	double 	getPressure(double x, double z, double y);
    double  getPressure(const Vector3d &location);
	int		getXBound(void) {return x_bound;}
	int		getYBound(void)	{return y_bound;}
	int		getZBound(void) {return z_bound;}
	double	getTransducerFreq(void) {return freq;}

	// Set the frequency of the transducer used to create the pressure map in k-Wave.
	void	setTransducerFreq(double freq) {this->freq = freq;}



private:
	// Ensure the default constructor is private.
	PressureMap();

	// Common init function for constructors of the class.
	void initCommon(void);

	// The bounds of the pressure grid [cm].
	int x_bound, y_bound, z_bound;

	// The number of voxels in the x, y, and z directions [cm].
	int Nx, Nz, Ny;

	// The voxel size [cm].
	double dx, dz, dy;

	// Frequency of the transducer that produced the pressure map.
	double freq;

	// Input stream
	std::ifstream pressure_file_stream;

	// Holds the name of the text file that contains the pressure values.
	std::string pressure_file;

	// Holds the pressure values obtained from k-Wave in a 3-dimensional grid
	// allowing us to index into the grid based on the coordinates of the phton
	// and retrieve the localized pressure.
	three_dim_array * pressure_grid;

};


#endif /* LOADPRESSUREMAP_H_ */
