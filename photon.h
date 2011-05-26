// Class defines the properties of a photon.
#ifndef PHOTON_H
#define PHOTON_H

#include "medium.h"
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <fstream>
#include <iostream>
using namespace std;
//#include <boost/random/uniform_real.hpp>
//#include <boost/random/variate_generator.hpp>
//#include <boost/random/mersenne_twister.hpp>

#define ALIVE 1		// Value depicting Photon should continue propagation.
#define DEAD  0	    // Photon has lost all energy and failed roulette.
#define ONE_MINUS_COSZERO 1.0E-12
/* If 1-cos(theta) <= ONE_MINUS_COSZERO, fabs(theta) <= 1e-6 rad. */
/* If 1+cos(theta) <= ONE_MINUS_COSZERO, fabs(PI-theta) <= 1e-6 rad. */
#define THRESHOLD	0.01		// Threshold for determining if we should perform roulette
#define CHANCE      0.1  		// Used in roulette
#define PI			3.1415926
#define SIGN(x)           ((x)>=0 ? 1:-1)
//const int MAX_BINS = 101;



class Medium;


class Photon
{
public:
	// Constructors
	Photon(void);
	Photon(double x, double y, double z,
		   double dirx, double diry, double dirz);
	
	// Destructor
	~Photon(void);
	
	// Set the number of iterations this Photon (i.e. thread) will run.
	void	setIterations(const int n);

	// Move photon to new position
	void	hop(void);

	// Drop absorbed energy from photon weight due to absorption in medium.
	void	drop(void);

	// Change the trajectory of the photon due to scattering in the medium.
	void	spin(void);
	
	// Set the step size of the photon.
	void 	setStepSize(void);

	// Decide whether the photon should be transmitted to another layer
	// or internally reflected.
	void	transmitOrReflect(void);

	// Reset the Photon attributes so it can be propogated again.
	void	reset(void);
		
	// Give the photon a probabilistic chance of surviving or terminating
	// if its weight has dropped below a specified threshold.
	void	performRoulette(void);

	// Return the cartesian coordinates
	double	getX(void) {return x;}
	double	getY(void) {return y;}
	double	getZ(void) {return z;}

	// Return the direction cosines
	double	getDirX(void) {return dirx;}
	double	getDirY(void) {return diry;}
	double	getDirZ(void) {return dirz;}

	// Return the current weight of the photon
	double	getWeight(void) {return weight;}
	
	// Returns a random number 'n': 0 < n < 1
	double	getRandNum(void);
	
	// Return the status of the photon.
	bool	isAlive(void) {return status;}

	// Terminate the current photon.
	void	kill(void) {status = DEAD;}

	// Calculate the new location of the photon and
	// move it to those coordinates.
	void	updateLocation(void);

	// Update weight based on specular reflectance.
	void	specularReflectance(double n1, double n2);
	
	// Plot the photon's path.
	void	plotPath(void);
	
	// Inject the photon into the medium the given number of iterations.
	// 'state[1,2,3,4]' represent the random initial values for the state
	// of the random number generator.
	void	injectPhoton(Medium *m, const int num_iterations, unsigned int state1, unsigned int state2,
							unsigned int state3, unsigned int state4);
	
	// Sets initial trajectory values.
	void	initTrajectory(void);
	
	// Zero's out the local detection array.
	void	initDetectionArray(void);

	// Initialize the RNG.
	void	initRNG(unsigned int s1, unsigned int s2, unsigned int s3, unsigned int s4);

	// Routines related to the thread-safe RNG
	unsigned int TausStep(unsigned int &z, int s1, int s2, int s3, unsigned long M);
	unsigned int LCGStep(unsigned int &z, unsigned int A, unsigned long C);
	double	HybridTaus(void);


	// Check if photon is still within the medium.
	bool	isPhotonInMedium(void);

	// Check if photon has come into contact with a layer boundary.
	bool 	hitLayerBoundary(void);

	// Displace (i.e. update the location) the photon some distance
	// based on the pressure at that location.
	void 	displacePhotonFromPressure(void);

	// Write the coordinates of each scattering event to file for use
	// with plotting in matlab.
	void 	writeCoordsToFile(void);

	// Calculate path length of photon from point 'p1' to point 'p2'.
	double	getPathLength(double x, double y, double z);

	// Add the coordinates of the photon at it's current position.
	void	captureLocationCoords(void);

	// Check if exit location is through the aperture that will fall
	// on the detector.
	bool	didExitThroughDetectorAperture(void);

	
private:
	// Number of times this Photon (i.e., thread) will execute; where one execution
	// is the full cycle of photon propagation.
	int iterations;
	
	// Holds value of number of iterations thus far.
	int cnt;	
	
	// The number of photons that exit through the detector aperture.
	int cnt_through_aperture;

	// Location of the photon.
	double	x, y, z;
	
	// Location of the photon with displacement from ultrasound source taken into account.
	double x_disp, y_disp, z_disp;

	// Radial position.
	double r;
	
	// Direction of the photon.
	double	dirx, diry, dirz;
	
	// Weight of the photon.
	double	weight;
	
	// Step size for the photon.
	double	step;
	
	// Step size to boundary.  Used when calculating distance from layer
	// boundary to current position of the photon.  Specifically, it is
	// the remainder of the step size after calculating the distance to
	// the layer boundary.
	// i.e. (step_size - distance_to_boundary)/mu_t
	double step_remainder;

	// status for current photon - dead (false) or alive (true).
	bool	status;
	
	// cosine and sine theta.  Used for trajectory calculations.
	double	cos_theta, sin_theta;
	
	// The azimuthal angle
	double	psi;
	
	
	// The number of steps this photon has taken while propagating through
	// the medium.
	int num_steps;
	
	// Pointer to the medium which this photon will propagate through.
	Medium *m_medium;

	// The thread id associated with this photon object.  The value is passed
	// in from the loop that creates the threads in main.cpp.
	int thread_id;

	// Local absorption array that holds values during execution.  This array
	// is copied over to the global absorption array (i.e. in the medium) once
	// a photon has finished propagating in the medium.
	// FIXME: 101 SHOULD NOT BE HARD CODED VALUE.
	double local_Cplanar[101];
	double	radial_size;			// Maximum radial size.
	int		num_radial_pos;			// Number of radial positions (i.e. NR).
	double	radial_bin_size;		// Radial bin size of the medium (i.e dr).


	// Boost Random Number Library implementation of Mersenne-twister RNG.
	//boost::mt19937 gen;

	// Used with the state for the thread-safe RNG.
	unsigned int z1, z2, z3, z4;

	// Tracks the path length of the photon through the medium.
	double original_path_length;
	double displaced_path_length;

	// Holds the x,y,z coordinates of the photon for each scattering event.
	// Used to plot the photon's path.
	vector<double> coords;


}; 		

#endif // end PHOTON_H
