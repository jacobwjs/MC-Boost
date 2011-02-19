#ifndef MEDIUM_H
#define MEDIUM_H


#include "layer.h"
#include "photon.h"
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
using namespace std;

const int MAX_BINS = 101;

class Medium
{
	
public:

	Medium();
	~Medium();
	
	// Add some portion of the photon's energy that was lost at this interaction
	// point (i.e. due to absorption) to the medium's grid.
	void absorbEnergy(const double z, const double energy);
	
	// Print the grid for this medium.
	void printGrid(const int num_photons);
	
	// Add a layer to the medium.
	void addLayer(Layer *layer);
	
	// Inject a photon into the medium at the x and y coordinates.
	void injectPhoton(const int x_start, const int y_start, Photon *p);
	
	// Move the photon through the medium stochastically.
	void propogatePhoton(Photon *photon);
	
	// Return the grid where absorption was accumulated.
	double * getPlanarGrid() {return Cplanar;}
	
	// Return the number of bins used in the grid.
	int getBins() {return MAX_BINS;}
	
	// Return the radial size of the medium (cm).
	double getRadialSize() {return radial_size;}

	void setPlanarArray(double *planar);
	
	
	
private:
	double	radial_size;			// Maximum radial size.
	int		num_radial_pos;			// Number of radial positions (i.e. NR).
	double	radial_bin_size;		// Radial bin size of the medium (i.e dr).
	
	// The arrays that hold the weights dropped during interaction points.
	//double	Cplanar[MAX_BINS];		// Planar photon concentration.
	double *Cplanar;
	double	Ccylinder[MAX_BINS];	// Clindrical photon concentration.
	double	Cspherical[MAX_BINS];	// Spherical photon concentration.
	
	
	// Create a vector to hold the layers of the medium.
	vector<Layer *> p_layers;
};

#endif	// MEDIUM_H

