
#include "debug.h"
#include "photon.h"

#undef DEBUG

Photon::Photon(void)
{
#ifdef DEBUG
	cout << "Creating Photon...\n";
#endif
	// seed the random number generator.
	srand(time(0));

	// Photon just created, so it is alive.
	status = ALIVE;

	// Weight, cartesian coords, and step size default values before photon
	// is moved through the medium.
	weight = 1;
	x = y = z = 0;
	step = 0;
	
	// No interactions thus far.
	num_steps = 0;
	
	// Randomly set photon trajectory to yield isotropic source.
	cos_theta = (2.0 * getRandNum()) - 1;
	sin_theta = sqrt(1.0 - cos_theta*cos_theta);
	psi = 2.0 * PI * getRandNum();
	dirx = sin_theta * cos(psi);
	diry = sin_theta * sin(psi);
	dirz = cos_theta;
	
	// 'cnt' represents the number of times a photon has propogated
	// through the medium.
	cnt = 0;
}

Photon::Photon(double x, double y, double z,
		double dirx, double diry, double dirz)
{
	cout << "Constructor stub...\n";
}

Photon::~Photon(void)
{
#ifdef DEBUG	
	cout << "Destructing Photon...\n";
#endif
}


// Set the number of iterations this thread will run.
void Photon::setIterations(const int num)
{
	iterations = num;
}



// 1) Hop - move the photon
// 2) Drop - drop weight due to absorption
// 3) Spin - update trajectory accordingly
// 4) Roulette - test to see if photon should live or die.
void Photon::injectPhoton(Medium *medium, const int iterations)
{
	// Before propagation we set the medium which will be used by the photon,
	// which is passed in.
	this->m_medium = medium;
	
	int i;
	// Inject 'iterations' number of photons into the medium.
	for (i = 0; i < iterations; i++) 
	{
		// While the photon has not been terminated by absorption or leaving
		// the medium we propagate it through he medium.
		while (isAlive()) 
		{

			
			// Move the photon.
			hop();
			
			
			// FIXME: Should probably have all weight deposited at surface if
			//		  photon leaves medium by total internal reflection.  Also
            //        the case where it reaches the max depth.
            
			// Ensure the photon has not left the medium by either total internal
			// reflection or transmission (only looking at z-axis).
			if (z >= 0 && z <= m_medium->getDepth()) {
				
				// Drop weight of the photon due to an interaction with the medium.
				drop();
				
				// Calculate the new coordinates of photon propogation.
				spin();
				
				// Test whether the photon should continue propagation from the
				// Roulette rule.
				performRoulette();
				
			}
			else {
                // If we make it here the photon has hit a boundary.  We simply absorb
                // all energy at the boundary.
                // FIXME:  Take into account specular reflectance since photon might not
                //          leave medium.
                m_medium->absorbEnergy(z, weight);
                break;  // break from while loop and execute reset().
			}
			
		} // end while() loop
		
		// Reset the photon and start propogation over from the beginning.
		reset();
		
	} // end for() loop
}



void Photon::plotPath()
{
	// STUB
}


void Photon::reset()
{
#ifdef DEBUG
	cout << "Reseting Photon...\n";
#endif
	
	// Photon just created, so it is alive.
	status = ALIVE;
	
	// Set back to initial weight values.
	weight = 1;
	
	// FIXME: 
	// need to reset to the photons initial location.
	x = y = z = 0;
	
	r = 0;
	step = 0;
	
	// Reset the number of interactions back to zero.
	num_steps = 0;
	
	// Randomly set photon trajectory to yield isotropic source.
	cos_theta = (2.0 * getRandNum()) - 1;
	sin_theta = sqrt(1.0 - cos_theta*cos_theta);
	psi = 2.0 * PI * getRandNum();
	dirx = sin_theta * cos(psi);
	diry = sin_theta * sin(psi);
	dirz = cos_theta;
	
	
}

// Step photon to new position.
void Photon::hop()
{
#ifdef DEBUG
	cout << "Hopping...\n";
#endif	
	
	double rnd = getRandNum();
	
	// Update the current values of the absorption and scattering coefficients
	// based on the depth in the medium (i.e. which layer the photon is in).
	double mu_a = m_medium->getLayerAbsorptionCoeff(z);  // cm^-1
	double mu_s = m_medium->getLayerScatterCoeff(z);	  // cm^-1
	
	// Calculate the new step length of the photon.
	step = -log(rnd)/(mu_a	+ mu_s);
	
	// Update position of the photon.
	x += step*dirx;
	y += step*diry;
	z += step*dirz;
}


// Return absorbed energy from photon weight at this location.
double Photon::drop()
{
#ifdef DEBUG
	cout << "Dropping...\n";
#endif	
	
	// Update the current values of the absorption and scattering coefficients
	// based on the depth in the medium (i.e. which layer the photon is in).
	double mu_a = m_medium->getLayerAbsorptionCoeff(z);  // cm^-1
	double mu_s = m_medium->getLayerScatterCoeff(z);	  // cm^-1
	
	// Calculate the albedo and remove a portion of the photon's weight for this
	// interaction.
	double albedo = mu_s / (mu_a + mu_s);
	double absorbed = weight * (1 - albedo);
	
	// Remove the portion of energy lost due to absorption at this location.
	weight -= absorbed;
	
	// Deposit lost energy in the grid of the medium.
	m_medium->absorbEnergy(z, absorbed);
	
	// Return the energy that was absorbed.
	return absorbed;
}

// Calculate the new trajectory of the photon.
void Photon::spin()
{
#ifdef DEBUG
	cout << "Spinning...\n";
#endif	
	
	// Get the anisotropy factor from the layer that resides at depth 'z' in
	// the medium.
	double g = m_medium->getAnisotropyFromDepth(z);
	
	
	double rnd = getRandNum();
	
	if (g == 0.0) {
		cos_theta = (2.0 * rnd) - 1.0;
	} 
	else {
		double temp = (1.0 - g*g)/(1.0 - g + 2*g*rnd);
		cos_theta = (1.0 + g*g - temp*temp)/(2.0*g);
	}
	sin_theta = sqrt(1.0 - cos_theta*cos_theta); /* sqrt() is faster than sin(). */
		
	// Sample 'psi'.
	psi = 2.0 * PI * getRandNum();
	double cos_psi = cos(psi);
	double sin_psi = 0.0;
	if (psi < PI) {
		sin_psi = sqrt(1.0 - cos_psi*cos_psi);     /* sqrt() is faster than sin(). */
	} 
	else {
		sin_psi = -sqrt(1.0 - cos_psi*cos_psi);
	}
	
	double uxx, uyy, uzz;
	double temp;
	/* New trajectory. */
	if (1 - fabs(dirz) <= ONE_MINUS_COSZERO) {      /* close to perpendicular. */
		uxx = sin_theta * cos_psi;
		uyy = sin_theta * sin_psi;
		uzz = cos_theta * SIGN(dirz);   /* SIGN() is faster than division. */
	} 
	else {					/* usually use this option */
		temp = sqrt(1.0 - dirz * dirz);
		uxx = sin_theta * (dirx * dirz * cos_psi - diry * sin_psi) / temp + dirx * cos_theta;
		uyy = sin_theta * (diry * dirz * cos_psi + dirx * sin_psi) / temp + diry * cos_theta;
		uzz = -sin_theta * cos_psi * temp + dirz * cos_theta;
	}
	
	// Update trajectory.
	dirx = uxx;
	diry = uyy;
	dirz = uzz;
}



void Photon::performRoulette(void)
{
	if (weight < THRESHOLD) {
		if (getRandNum() <= CHANCE) {
			weight /= CHANCE;
		}
		else {
			status = DEAD;
		}
	}

}



double Photon::getRandNum(void)
{
	double rnd = (double)rand()/(double)RAND_MAX;
	while ((rnd == 0) || (rnd == 1)) { // produces 0 < rnd < 1
		rnd = (double)rand()/(double)RAND_MAX;
	}
	return rnd;
}


// FIXME: Currently not in use
void Photon::specularReflectance(double n1, double n2)
{
	// update the weight after specular reflectance has occurred.
	weight = weight - (pow((n1 - n2), 2) / pow((n1 + n2), 2)) * weight;
}

