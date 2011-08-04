/*
 * displacementMap.h
 *
 *  Created on: 3 aug. 2011
 *      Author: StaleyJW
 */

#ifndef DISPLACEMENTMAP_H_
#define DISPLACEMENTMAP_H_



#include <boost/multi_array.hpp>
#include <iostream>
using std::cout;
using std::endl;
#include <fstream>
using std::ifstream;
using std::ofstream;
#include <string>



// Boost array that holds the displacement values after being
// loaded in from file.
typedef boost::multi_array<double, 3> three_dim_array;
typedef three_dim_array::index array_index;


// Forward decleration of Vector3d class.
class Vector3d;


class DisplacementMap
{
public:
	DisplacementMap(const std::string &filename, const int Nx, const int Nz, const int Ny, const int grid_size);
    DisplacementMap(const int Nx, const int Nz, const int Ny, const int grid_size);

	~DisplacementMap();

	// Loads a text file containing discrete displacement values at a given time step
	// that were obtained from kWave simulation post-processed data.
	void	loadDisplacementMaps(const std::string &filename, const int timeStep);


	// Returns a Vector3d object holding values for displacements in all axes.
	// That is the returned Vector3d objects holds the values the coordinates of
	// the photon should be displaced accordingly.
	Vector3d	getDisplacements(const Vector3d &photonLocation);



private:
	// Ensure the default constructor can never be called.
	DisplacementMap();

	// Common initialization function.
	void	initCommon();

	// Input stream.
	std::ifstream disp_file_stream;

	// The bounds of the pressure grid [cm].
	int x_bound, y_bound, z_bound;

	// The number of voxels in the x, y, and z directions [cm].
	int Nx, Nz, Ny;

	// The voxel size [cm].
	double dx, dz, dy;

	// Holds the displacement values obtained from k-Wave in a 3-dimensional grid
	// allowing indexing into the grid based on the coordinates of the photon
	// and retrieve the localized displacement.
	// NOTE:
	// - Displacement happens on each dimension seperate from the other, based on the
	//   speed of sound in each direction.  Therefore post-processing of velocity data
	//   leaves displacement values in each direction, therefore we need 3 grids.
	three_dim_array * displacement_gridX;
	three_dim_array * displacement_gridY;
	three_dim_array * displacement_gridZ;



};

#endif /* DISPLACEMENTMAP_H_ */
