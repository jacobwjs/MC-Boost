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
class PressureMap;
class RefractiveMap;
class DisplacementMap;
class Detector;
class Layer;
class Vector3d;



typedef struct {
    // Create a pointer to a PressureMap object.  For use with
	// modeling acousto-optics.
	PressureMap * pmap;

	// Create a pointer to a RefractiveMap object.  For use with
	// modeling acousto-optics.
	RefractiveMap * nmap;
    
    // Create a pointer to a Displacement object.  For use with
    // modeling acousto-optics.
    DisplacementMap * dmap;
    
    // Frequency of the transducer used.
    double transducerFreq;
    
    // The wavenumber of the ultrasound.
    double waveNumber;

    // Number of time steps in the simulation.
    int totalTimeSteps;
    
} kWaveSim;



// Medium is a container object that holds one or many layer objects that the
// photon is propagated through.  This allows easy simulation of heterogeneous 
// media with Monte Carlo simulations.
class Medium
{
	
public:

	friend class Photon;

	// A structure that holds K-Wave simulation data and attributes
	// This is public since the members of the struct are objects with
	// their own private methods.  This is simply a convenient container.
	kWaveSim kwave;

	Medium(void);
    Medium(const double x, const double y, const double z);
	~Medium();
    
    // Common initializations for the Medium object.  Called from constructors.
    void    initCommon(void);

    // Set acoustic properties of medium.
    // NOTE: 'eta' is the pezio-optical coefficient
    void	setDensitySOSPezio(const double density, const double SOS, const double eta);

    // Set the density of the medium.
    void	setDensity(const double density) {this->density = density;}

    // Set the Speed-of-Sound of the medium.
    void	setSOS(const double sos) {this->speed_of_sound = sos;}

    // Set the pezio-optical coefficient of the medium.
    void	setPezioOpticCoeff(const double eta) {this->pezio_optical_coeff = eta;}

    // Set the background refractive index.
    void	setBackgroundRefractiveIndex(const double n_background) {this->background_refractive_index = n_background;}
	
	// Add some portion of the photon's energy that was lost at this interaction
	// point (i.e. due to absorption) to the medium's grid.
	void	absorbEnergy(const double z, const double energy);
	
	// Same as above, only the argument is an array of absorbed energy values
	// that is copied entirely to the Medium.
	void	absorbEnergy(const double *energy_array);

	// Print the grid for this medium.
	void	printGrid(const int num_photons);
    
    // Set the number of time steps that occurred in the K-Wave simulation.
    void    setKWaveTimeSteps(const int timeSteps) {kwave.totalTimeSteps = timeSteps;}
	
	// Add a layer to the medium.
	void	addLayer(Layer *layer);
    
    // Add a detector to the medium.
    void    addDetector(Detector *detector);
    
    // See if photon has crossed the detector plane.
    int    photonHitDetectorPlane(const boost::shared_ptr<Vector3d> p0);
	
	// Add a pressure map object that holds pressure generated from K-Wave
    void 	addPressureMap(PressureMap *p_map);
    
    // Add a refractive map object that holds refractive index values generated from k-Wave pressures.
    void	addRefractiveMap(RefractiveMap *n_map);

    // Add a displacement map object that holds pressure generated from k-Wave
    void    addDisplacementMap(DisplacementMap *d_map);
    
	// Load the pressure data generated from K-Wave (at simulation time step 'dt') into the pressure map object.
	void 	loadPressure(std::string &filename, const int dt);

	// Load the pressure data generated from K-Wave if only the file name is given.
	void	loadPressure(std::string &filename);
    
	// Load the displacement data generated from K-Wave (at simulation time step 'dt') into the displacement map object.
    void    loadDisplacements(std::string &filename, const int dt);
    
    // Loads a pressure file from k-Wave generated data and calculates the displacement, essentially converting the data to displacements
    // representative of simulated pressures.
    // TODO: Implement this method so that simulations can account for varying attributes over time due to heating from ultrasound.
//    void	loadDisplacementsFromPressure(std::string &filename, const int dt,
//    									  const double density,
//    									  const double speed_of_sound,
//    									  const double pezio_optical_coeff,
//    									  const double background_refractive_index);

    void	loadDisplacementsFromPressure(std::string &filename, const int dt);

    // Load the pressure data generated from k-Wave and convert it to refractive index values.
    // NOTE: eta is the pezio-optical coefficient of the medium.
    void 	loadRefractiveMap(std::string &filename, const double density, const double sos, const double n_background, const double eta);
    void    loadRefractiveMap(std::string &filename, const double density, const double sos, const double eta, const double n_background, const int dt);

	// Return the pressure from the pressure grid based on cartesian coordinates
	// of the current location of the photon in the medium.
	double	getPressureFromCartCoords(double x, double z, double y);
    
    // Return the pressure from the pressure grid based on the location of the photon.
    double  getPressureFromPhotonLocation(const boost::shared_ptr<Vector3d> photonCoords);
    

    // Return the displacement vector coordinates from the location of the photon in the medium.
    boost::shared_ptr<Vector3d> getDisplacementFromPhotonLocation(const boost::shared_ptr<Vector3d> photonCoords);

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
	Layer * getLayerAboveCurrent(Layer *currentLayer);

	// Return the layer below the current layer.
	Layer * getLayerBelowCurrent(double depth);

    // Return the max depth of the medium.
    double 	getDepth() {return depth;}
    
    // Return the refractive index of the medium.
    //double  getRefractiveIndex(void) {return refractive_index;}
    
    // Write photon coordinates to file.
    void 	writePhotonCoords(std::vector<double> &coords);

    // Write photon exit locations and phases to file.
    void	writeExitCoordsAndLength(std::vector<double> &coords_phase);

    // Write the photon exit locations, phase and weight to file.
    void	writeExitCoordsLengthWeight(std::vector<double> &coords_phase_weight);

    // Return the bounds of the medium.
    double getXbound(void) {return x_bound;}
    double getYbound(void) {return y_bound;}
    double getZbound(void) {return z_bound;}
	
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


	// File for dumping photon paths to.  Used in the Photon class.
    std::ofstream coords_file;

	// File for dumping data regarding exit location, path length, weight, etc.
	// to file for post processing in matlab.
    std::ofstream photon_data_file;
    
    // The refrective index outside of the medium.  We assume air.
    double refractive_index;
    
    // The unmodulated (i.e. background) refractive index of the medium.
    double background_refractive_index;


    // The voxel sizes of the medium.  This will match the size of the k-Wave simulation voxel size
    // and is only set here for convenience and later use in the 'Photon' class.
    double dx, dy, dz;
    double Nx, Ny, Nz;

    // The density of the medium.
    double density;

    // The acoustic speed-of-sound of the medium.
    double speed_of_sound;

    // The adiabatic piezo-optical coefficient of the medium.
    double pezio_optical_coeff;

};

#endif	// MEDIUM_H

