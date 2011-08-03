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
#include <fstream>
#include <string>
using namespace std;


typedef boost::multi_array<double, 3> three_dim_array;
typedef three_dim_array::index array_index;


class PressureMap
{
public:
	PressureMap();
	PressureMap(const string &filename, const int x, const int z, const int y, const int dimension);
	~PressureMap();

	void 	loadPressureMap(void);
	double 	getPressureFromGrid(int x, int z, int y);
	double 	getPressureCartCoords(double x, double z, double y);
	int		getXBound(void) {return x_bound;}
	int		getYBound(void)	{return y_bound;}
	int		getZBound(void) {return z_bound;}
	double	getTransducerFreq(void) {return freq;}

	// Set the frequency of the transducer used to create the pressure map in k-Wave.
	void	setTransducerFreq(double freq) {this->freq = freq;}



private:
	void initCommon(void);	// Common init function for constructors of the class.

	// Holds the pressure values obtained from k-Wave in a 3-dimensional grid
	// allowing us to index into the grid based on the coordinates of the phton
	// and retrieve the localized pressure.
	three_dim_array * pgrid;

	// The bounds of the pressure grid [cm].
	int x_bound, y_bound, z_bound;

	// The number of voxels in the x, y, and z directions [cm].
	int Nx, Nz, Ny;

	// The voxel size [cm].
	double dx, dz, dy;

	// Frequency of the transducer that produced the pressure map.
	double freq;

	// Input stream
	ifstream pressure_file;

	// Holds the name of the text file that contains the pressure values.
	string p_file;

};


#endif /* LOADPRESSUREMAP_H_ */
