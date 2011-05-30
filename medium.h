#ifndef MEDIUM_H
#define MEDIUM_H


#include "layer.h"
#include "pressureMap.h"
#include "photon.h" // Photon class is a friend of the Medium class.
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <boost/thread/mutex.hpp>
using namespace std;

// Maximum number of bins that hold absorption values.
const int MAX_BINS = 101;


// Medium is a container object that holds one or many layer objects that the
// photon is propagated through.  This allows easy simulation of heterogeneous 
// media with Monte Carlo simulations.
class Medium
{
	
public:

	friend class Photon;

	Medium();
    Medium(const int depth, const int x, const int y, const int z);
	~Medium();
	
	// Common values to be used with the constructors.
	void	initCommon(void);

	// Add some portion of the photon's energy that was lost at this interaction
	// point (i.e. due to absorption) to the medium's grid.
	void	absorbEnergy(const double z, const double energy);
	
	// Same as above, only the argument is an array of absorbed energy values
	// that is copied entirely to the Medium.
	void	absorbEnergy(const double *energy_array);

	// Print the grid for this medium.
	void	printGrid(const int num_photons);
	
	// Add a layer to the medium.
	void	addLayer(Layer *layer);
	
	// Add a pressure map object that holds pressure generated from K-Wave in order to simulate
	// acousto-optics.
	void 	addPressureMap(PressureMap *pmap);

	// Load the pressure data generated from K-Wave into the pressure map object.
	void 	loadPressure(void);

	// Return the pressure from the pressure grid based on cartesian coordinates
	// of the current location of the photon in the medium.
	double	getPressureFromCartCoords(double x, double z, double y);

	// Return the pressure from the pressure grid based on array index into
	// the matrix.
	double	getPressureFromGridCoords(int x, int z, int y);

	// Return the grid where absorption was accumulated.
	double * getPlanarGrid() {return Cplanar;}
	
	// Return the number of bins used in the grid.
	int		getBins() {return MAX_BINS;}
	
	// Return the radial size of the medium (cm).
	double	getRadialSize() {return radial_size;}

	// Return the bin size of the detector array (i.e. dr).
	double 	getRadialBinSize() {return radial_bin_size;}

	// Return the number of radial positions.
	double 	getNumRadialPos() {return num_radial_pos;}

	// Assign the array which will hold the planar absorbance values.
	void	setPlanarArray(double *planar);
	
	// Returns the absorption coefficient (mu_a) for a given depth (i.e. a layer).
	double	getLayerAbsorptionCoeff(double depth);
	
	// Returns the scattering coefficient (mu_s) for a given depth (i.e. a layer).
	double	getLayerScatterCoeff(double depth);
	
	// Return the anisotropy ('g') value for a given depth (i.e. a layer).
	double	getAnisotropyFromDepth(double depth);
	
	// Return layer from depth passed in.
	Layer * getLayerFromDepth(double depth);

	// Return the layer above the current layer.
	Layer * getLayerAboveCurrent(double depth);

	// Return the layer below the current layer.
	Layer * getLayerBelowCurrent(double depth);

    // Return the max depth of the medium.

    double 	getDepth() {return depth;}

    double 	getMediumZaxisBound(void) {return z_bound;}
    double 	getMediumXaxisBound(void) {return x_bound;}
    double 	getMediumYaxisBound(void) {return y_bound;}
	
    // Write photon coordinates to file.
    void 	writePhotonCoords(vector<double> &coords);

    // Write photon exit locations and phases to file.
    void	writeExitCoordsAndPhase(double x_disp, double y_disp, double displaced_path_length);

	
private:
	double	radial_size;			// Maximum radial size.
	int		num_radial_pos;			// Number of radial positions (i.e. NR).
	double	radial_bin_size;		// Radial bin size of the medium (i.e dr).
	
	// The arrays that hold the weights dropped during interaction points.
	//double	Cplanar[MAX_BINS];		// Planar photon concentration.
	double *Cplanar;
	//double	Ccylinder[MAX_BINS];	// Clindrical photon concentration.
	//double	Cspherical[MAX_BINS];	// Spherical photon concentration.
	
    // The total depth of the medium.
    double depth;
<<<<<<< Updated upstream
    double x_bound,
           y_bound,
           z_bound;
=======
    double x_bound, y_bound, z_bound;
>>>>>>> Stashed changes
	
	// Create a vector to hold the layers of the medium.
	vector<Layer *> p_layers;

	// Mutex to serialize access to the sensor array.
	boost::mutex m_sensor_mutex;

	// Mutex to serialize access to the data file that is written
	// by photons to store their coordinates during propagation.
	boost::mutex m_coords_mutex;

	// Mutex to serialize access to the data file that is written
	// by photons to store their coordinates during propagation.
	boost::mutex m_exit_phase_mutex;

	// Create a pointer to a PressureMap object.  For use with
	// modeling acousto-optics.
	PressureMap * pmap;


	// File for dumping photon paths to.  Used in the Photon class.
	ofstream coords_file;
	ofstream exit_location_and_phase_file;

};

#endif	// MEDIUM_H

