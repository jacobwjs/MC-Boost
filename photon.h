// Class defines the properties of a photon.
#ifndef PHOTON_H
#define PHOTON_H

#include "medium.h"
#include <cmath>
#include <ctime>
#include <cstdlib>

#define ALIVE 1		// Value depicting Photon should continue propagation.
#define DEAD  0	    // Photon has lost all energy and failed roulette.

#define ONE_MINUS_COSZERO 1.0E-12
/* If 1-cos(theta) <= ONE_MINUS_COSZERO, fabs(theta) <= 1e-6 rad. */
/* If 1+cos(theta) <= ONE_MINUS_COSZERO, fabs(PI-theta) <= 1e-6 rad. */

#define THRESHOLD	0.01		// Threshold for determining if we should perform roulette
#define CHANCE      0.1  		// Used in roulette
#define PI			3.141592

#define SIGN(x)           ((x)>=0 ? 1:-1)

class Photon
{
public:
	Photon(void);
	Photon(double x, double y, double z,
		   double dirx, double diry, double dirz);
	~Photon(void);
	
	 
	// Implement the virtual function 'run' from Thread.
	virtual void run();
	
	// Set the number of iterations this Photon (i.e. thread) will run.
	void	setIterations(const int n);

	// Move photon to new position
	void	hop();

	// Drop absorbed energy from photon weight due to absorption in medium.
	void	drop();

	// Change the trajectory of the photon due to scattering in the medium.
	void	spin();
	
	// Reset the Photon attributes so it can be propogated again.
	void	reset();
	
	// Set the layer for this photon.
	//void	setLayer(Medium *tissue) {this->medium = tissue;}
	
	// Set the medium which the photon will propogate through.
	// NOTE: can be made up of multiple layers.
	void	setMedium(Medium *tissue) {this->medium = tissue;}
	
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

	// Step the photon using stochastic methods.
	void	move(Layer &l);
	
	// Update weight based on specular reflectance.
	void	specularReflectance(double n1, double n2);
	
	
	
	
private:
	// Number of times this Photon (i.e., thread) will execute.
	int iterations;
	
	// Holds value of number of iterations thus far.
	int cnt;	
	
	// Location of the photon.
	double	x, y, z;
	
	// Radial position.
	double r;
	
	// Direction of the photon.
	double	dirx, diry, dirz;
	
	// Weight of the photon.
	double	weight;
	
	// Step size for the photon.
	double	step;
	
	// status for current photon - dead (false) or alive (true).
	bool	status;
	
	// cosine and sine theta.  Used for trajectory calculations.
	double	cos_theta, sin_theta;
	
	// The azimuthal angle
	double	psi;	
	
	// The Medium that the photon will be propogated through.
	Medium *medium;
	
}; 		

#endif // end PHOTON_H
