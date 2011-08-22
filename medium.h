#ifndef MEDIUM_H
#define MEDIUM_H

#include "photon.h" // Photon class is a friend of the Medium class.
#include <vector>
#include <string>
#include <iostream>
using std::cout;
using std::endl;
#include <iomanip>
#include <fstream>
#include <boost/thread/mutex.hpp>
#include <boost/shared_ptr.hpp>

// Maximum number of bins that hold absorption values.
const int MAX_BINS = 101;


// Forward declaration of PressureMap and DisplacementMap objects.
class Detector;
class Layer;
class Vector3d;






// Medium is a container object that holds one or many layer objects that the
// photon is propagated through.  This allows easy simulation of heterogeneous 
// media with Monte Carlo simulations.
class Medium
{
	
public:

	friend class Photon;


	// Constructor.
    Medium(const double x, const double y, const double z);
	~Medium();
    
    // Common initializations for the Medium object.  Called from constructors.
    void    initCommon(void);
	
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
    
    // Add a detector to the medium.
    void    addDetector(Detector *detector);
    
    // See if photon has crossed the detector plane.
    int    photonHitDetectorPlane(const boost::shared_ptr<Vector3d> p0);

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
	Layer * getLayerAboveCurrent(Layer *currentLayer);

	// Return the layer below the current layer.
	Layer * getLayerBelowCurrent(double depth);

    // Return the max depth of the medium.
    double 	getDepth() {return depth;}
    
    // Return the refractive index of the medium.
    //double  getRefractiveIndex(void) {return refractive_index;}
    

    // Return the bounds of the medium.
    double getXbound(void) {return x_bound;}
    double getYbound(void) {return y_bound;}
    double getZbound(void) {return z_bound;}
	
private:
    // Ensure the medium is defined with specific attributes, so we make the
    // default constructor private.
    Medium(void);

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
    double x_bound,
           y_bound,
           z_bound;
	
	// Create a STL vector to hold the layers of the medium.
    std::vector<Layer *> p_layers;
    
    // Create a STL vector to hold the detectors in the medium.
    std::vector<Detector *> p_detectors;
    
	// Mutex to serialize access to the sensor array.
	boost::mutex m_sensor_mutex;

	// Mutex to serialize access to the data file that is written
	// by photons.
	boost::mutex m_data_file_mutex;

    // The refrective index outside of the medium.  We assume air.
    double refractive_index;
    
};

#endif	// MEDIUM_H

