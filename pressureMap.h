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
	PressureMap(const char *filename);
	~PressureMap();

	void 	loadPressureFromKWave(void);
	double 	getPressureFromGrid(int x, int z, int y);
	double 	getPressureCartCords(double x, double z, double y);



private:
	void init(void);	// Common init function for constructors of the class.

	three_dim_array * pgrid;
	int Nx, Nz, Ny;
	double dx, dz, dy;
	ifstream pressure_file;
	string p_file;

};


#endif /* LOADPRESSUREMAP_H_ */
