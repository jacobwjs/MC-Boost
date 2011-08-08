/*
 * displacementMap.cpp
 *
 *  Created on: 3 aug. 2011
 *      Author: StaleyJW
 */

#include "debug.h"
#include "vector3d.h"
#include "displacementMap.h"
#include <boost/lexical_cast.hpp>


// It's an error to create a DisplacementMap object without specifying attributes,
// therefore the default constructor should never be called.
DisplacementMap::DisplacementMap()
{
	cout << "Error: Default DisplacementMap() constructor called\n";
}


DisplacementMap::DisplacementMap(const std::string &filename, const int Nx, const int Nz, const int Ny, const int grid_size)
{
    // Assign the number of grid points (pixels in k-wave) used in the simulation.
	this->Nx = Nx;
	this->Ny = Ny;
	this->Nz = Nz;
    
	x_bound = y_bound = z_bound = grid_size;  // [cm]
	
    
    initCommon();
}


DisplacementMap::DisplacementMap(const int Nx, const int Nz, const int Ny, const int grid_size)
{
    // Assign the number of grid points (pixels in k-wave) used in the simulation.
	this->Nx = Nx;
	this->Ny = Ny;
	this->Nz = Nz;
    
	x_bound = y_bound = z_bound = grid_size;  // [cm]
	
    initCommon();
}



DisplacementMap::~DisplacementMap()
{
	if (displacement_gridX)
		delete displacement_gridX;

	if (displacement_gridY)
		delete displacement_gridY;

	if (displacement_gridZ)
		delete displacement_gridZ;
}



// Loads a text files containing discrete displacement values at a given time step, in all dimensions (i.e. x, y, z),
// that were obtained from kWave simulation post-processed data.
void DisplacementMap::loadDisplacementMaps(const std::string &filename, const int timeStep)
{



	// Assure memory has been allocated for the pressure values that
	// will be read in from file.  That is, initCommon() has already
	// been called.
	assert(displacement_gridX != NULL);
	assert(displacement_gridY != NULL);
	assert(displacement_gridZ != NULL);

	// A pointer to the string that is set for opening the displacement file.
	std::string file_to_open;

	// A pointer to one of the displacement grid arrays.  Set depending on
	// which grid should be filled below.
	three_dim_array *p_displacement_grid = NULL;

	// Load each data structure with their respective displacement data.
	for (int i = 0; i < 3; i++)
	{


		// Open the file that contains the pressure values for the specific
		// dimension based on the loop index.
		//
		if (i == 0)
		{  // X-displacement file.

			// Clear the string.
			file_to_open.clear();

			// Concatonate the values passed in to form a filename to read in.
			file_to_open = filename + "X-" + boost::lexical_cast<std::string>(timeStep);
			disp_file_stream.open(file_to_open.c_str());

			// The appropriate displacement grid is assigned to be filled below.
			p_displacement_grid = displacement_gridX;
		}
		else if (i == 1)
		{  // Y-displacement file.

			// Clear the string.
			file_to_open.clear();

			// Concatonate the values passed in to form a filename to read in.
			file_to_open = filename + "Y-" + boost::lexical_cast<std::string>(timeStep);
			disp_file_stream.open(file_to_open.c_str());

			// The appropriate displacement grid is assigned to be filled below.
			p_displacement_grid = displacement_gridY;
		}
		else
		{  // Z-displacement file.
			file_to_open.clear();

			// Concatonate the values passed in to form a filename to read in.
			file_to_open = filename + "Z-" + boost::lexical_cast<std::string>(timeStep);
			disp_file_stream.open(file_to_open.c_str());

			// The appropriate displacement grid is assigned to be filled below.
			p_displacement_grid = displacement_gridZ;
		}


        // Check for successful opening of the file.
		if (!disp_file_stream)
		{
			cout << "!!! Error opening displacement map file " << file_to_open.c_str() << "!!!\n";
			exit(1);
		}
		else
		{
			cout << "Displacement map " << file_to_open.c_str() << " opened successfully. ";
			cout << "Loading displacement values...\n";
		}


		double data = 0.0;
        // Read in data to the proper displacement array.
		for (array_index a = 0; a < Nx && disp_file_stream.good(); a++)
        {
			for (array_index b = 0; b < Nz; b++)
            {
				for (array_index c = 0; c < Ny; c++)
				{
					disp_file_stream >> data;
					(*p_displacement_grid)[a][b][c] = data;
                    //cout << (*p_displacement_grid)[a][b][c] << endl;
				}
            }
        }


		disp_file_stream.close();
	}
}


void DisplacementMap::initCommon()
{
    assert(Nx != 0 &&
           Ny != 0 &&
           Nz != 0);
    
	dx = (double)x_bound / (double)Nx; // [cm] (note: 20e-3/Nx in centimeters is 0.16;
	dy = (double)y_bound / (double)Ny;
	dz = (double)z_bound / (double)Nz;
	displacement_gridX = new three_dim_array (boost::extents[Nx][Nz][Ny]);
	displacement_gridY = new three_dim_array (boost::extents[Nx][Nz][Ny]);
	displacement_gridZ = new three_dim_array (boost::extents[Nx][Nz][Ny]);
}


// Returns a Vector3d object holding values for displacements in all axes.
boost::shared_ptr<Vector3d> DisplacementMap::getDisplacements(const Vector3d &photonLocation)
{

	const double THRESHOLD = 1.0e-12;
    boost::shared_ptr<Vector3d> result(new Vector3d);
    

    // Indices into the displacement grids.
    // NOTE:
    // - There is an effective transformation of the y (photon)
    //   and z (grid) transformation due to the grid orientation in kWave.
    int _x = floor(photonLocation.location.x/dx + THRESHOLD);
    int _z = floor(photonLocation.location.y/dz + THRESHOLD);
    int _y = floor(photonLocation.location.z/dy + THRESHOLD);
    

    // Sanity check.
	assert((_x < Nx && _x >= 0) &&
           (_y < Ny && _y >= 0) &&
           (_z < Nz && _z >= 0) ||
           assert_msg("_x=" << _x << " _y=" << _y << " _z=" << _z << "\n"
        		      << photonLocation.location.x << " "
        		      << photonLocation.location.y << " "
        		      << photonLocation.location.z));
    



    result->location.x = getDisplacementFromGridX(_x, _z, _y);
    result->location.y = getDisplacementFromGridY(_x, _z, _y);
    result->location.z = getDisplacementFromGridZ(_x, _z, _y);
    
	return result;

}



boost::shared_ptr<Vector3d> DisplacementMap::getDisplacements(const double x, const double y, const double z)
{
    
    
    // Indices into the displacement grids.
    // NOTE:
    // - There is an effective transformation of the y (photon)
    //   and z (grid) transformation due to the grid orientation in kWave.
    /*
    int _x = floor(x/dx);
	int _y = floor(z/dy);
    int _z = floor(y/dz);

    
    // Sanity check.
	assert((_x < Nx && _x >= 0) &&
           (_y < Ny && _y >= 0) &&
           (_z < Nz && _z >= 0));
    
    */
    
    boost::mutex::scoped_lock lock(m_displacement_mutex);

    boost::shared_ptr<Vector3d> result (new Vector3d);
    result->location.x = getDisplacementFromGridX(x, z, y);
    result->location.y = getDisplacementFromGridY(x, z, y);
    result->location.z = getDisplacementFromGridZ(x, z, y);
    
    return result;
}

// Returns the individual axis displacement value from their location in the grid.
double DisplacementMap::getDisplacementFromGridX(const int a, const int b, const int c)
{
	return (*displacement_gridX)[(array_index)a][(array_index)b][(array_index)c];
}

double DisplacementMap::getDisplacementFromGridY(const int a, const int b, const int c)
{
   	return (*displacement_gridY)[(array_index)a][(array_index)b][(array_index)c];
}

double DisplacementMap::getDisplacementFromGridZ(const int a, const int b, const int c)
{
	return (*displacement_gridZ)[(array_index)a][(array_index)b][(array_index)c];    
}
