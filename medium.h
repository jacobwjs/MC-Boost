#ifndef MEDIUM_H
#define MEDIUM_H


#include "layer.h"
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
using namespace std;

// Maximum number of bins that hold absorption values.
const int MAX_BINS = 101;


// Medium is a container object that holds one or many layer objects that the
// photon is propagated through.  This allows easy simulation of heterogeneous 
// media with Monte Carlo simulations.
class Medium
{
	
public:

	Medium();
    Medium(const int depth);
	~Medium();
	
	// Add some portion of the photon's energy that was lost at this interaction
	// point (i.e. due to absorption) to the medium's grid.
	void	absorbEnergy(const double z, const double energy);
	
	// Print the grid for this medium.
	void	printGrid(const int num_photons);
	
	// Add a layer to the medium.
	void	addLayer(Layer *layer);
	
	// Return the grid where absorption was accumulated.
	double * getPlanarGrid() {return Cplanar;}
	
	// Return the number of bins used in the grid.
	int		getBins() {return MAX_BINS;}
	
	// Return the radial size of the medium (cm).
	double	getRadialSize() {return radial_size;}

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
	
    
    // Return the max depth of the medium.
    double getDepth() {return depth;}
	
	
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
	
	// Create a vector to hold the layers of the medium.
	vector<Layer *> p_layers;
};

#endif	// MEDIUM_H

