
#include "debug.h"
#include "photon.h"

#undef DEBUG


double const Pi=4*atan(1);


Photon::Photon(void)
{
#ifdef DEBUG
	cout << "Creating Photon...\n";
#endif

	/*
	// Photon just created, so it is alive.
	status = ALIVE;

	// Weight, cartesian coords, and step size default values before photon
	// is moved through the medium.
	weight = 1;
	x = y = z = 0;
	step = 0;
	step_remainder = 0;
	
	// No interactions thus far.
	num_steps = 0;
	*/
	this->reset();
	


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


void Photon::initTrajectory()
{
	// Randomly set photon trajectory to yield isotropic source.
	cos_theta = (2.0 * getRandNum()) - 1;
	sin_theta = sqrt(1.0 - cos_theta*cos_theta);
	psi = 2.0 * PI * getRandNum();
	dirx = sin_theta * cos(psi);
	diry = sin_theta * sin(psi);
	//dirz = cos_theta;
	dirz = 1;

}

void Photon::initDetectionArray()
{
	// Zero out the local detection array since it will accumulate
		// values over many iterations.
		int i;
		for (i = 0; i < MAX_BINS; i++) {
			local_Cplanar[i] = 0;
		}
}


void Photon::initRNG(unsigned int state1, unsigned int state2,
						unsigned int state3, unsigned int state4)
{
	z1 = state1;
	z2 = state2;
	z3 = state3;
	z4 = state4;
}


// BOOST thread library starts execution here.
// 1) Hop - move the photon
// 2) Drop - drop weight due to absorption
// 3) Spin - update trajectory accordingly
// 4) Roulette - test to see if photon should live or die.
void Photon::injectPhoton(Medium *medium, const int iterations, unsigned int state1, unsigned int state2,
							unsigned int state3, unsigned int state4)
{
	// seed the random number generator.
	//srand(time(0) + thread_id);

	// Seed the Boost RNG (Random Number Generator).
	//gen.seed(time(0) + thread_id);

	// Initialize the photon's properties before propagation begins.
	initTrajectory();
	initDetectionArray();
	initRNG(state1, state2, state3, state4);

	// Before propagation we set the medium which will be used by the photon.
	this->m_medium = medium;
	

	// Assign local values of the detection grid from the Medium.
	radial_bin_size = m_medium->getRadialBinSize();
	num_radial_pos = m_medium->getNumRadialPos();



	// Inject 'iterations' number of photons into the medium.
	int i;
	for (i = 0; i < iterations; i++) 
	{
		// While the photon has not been terminated by absorption or leaving
		// the medium we propagate it through he medium.
		while (isAlive()) 
		{

			// Calculate and set the step size for the photon.
			setStepSize();
			
			// XXX: NOT CURRENTLY IN USE.
			// If the step size caused the photon to cross a boundary of the
			// layer, and therefore potentially needing to take into account
			// reflection, transmission and refraction, we account for these
			// phonomenon accordingly.
			/*
			if (hitLayerBoundary())
			{
				hop();	// Move the photon to the boundary.
				transmitOrReflect();	// Calculate whether to transmit the
										// photon or reflect it.

			}
			*/


			
			
			// FIXME: Should probably have all weight deposited at surface if
			//		  photon leaves medium by total internal reflection.  Also
            //        the case where it reaches the max depth.
            
			// Ensure the photon has not left the medium by either total internal
			// reflection or transmission (only looking at z-axis).
			//if (z >= 0 && z <= m_medium->getDepth()) {
				// Move the photon in the medium.
				hop();

				// Drop weight of the photon due to an interaction with the medium.
				drop();
				
				// Calculate the new coordinates of photon propagation.
				spin();
				
				// Test whether the photon should continue propagation from the
				// Roulette rule.
				performRoulette();
				
			//}
			//else {
                // If we make it here the photon has hit a boundary.  We simply absorb
                // all energy at the boundary.
                // FIXME:  Take into account specular reflectance since photon might not
                //          leave medium.
                //m_medium->absorbEnergy(z, weight);
                //break;  // break from while loop and execute reset().
			//}
			
		} // end while() loop
		
		cout << "Non-displaced path length: " << scientific << setprecision(12) <<  original_path_length << endl;
		cout << "Displaced path length: " << scientific << setprecision(12) << displaced_path_length << endl;
		cout << "Displaced - Original = " << scientific << setprecision(12)
					<< displaced_path_length - original_path_length << endl;

		// Write out the x,y,z coordinates of the photons path as it propagated through
		// the medium.
		writeCoordsToFile();

		// Reset the photon and start propogation over from the beginning.
		reset();
		
	} // end for() loop


	// This thread has executed all of it's photons, so now we update the global
	// absorption array in the medium.
	m_medium->absorbEnergy(local_Cplanar);
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
	
	// FIXME: ASSUMES GRID IS 10x10x10 cm.
	// need to reset to the photons initial location.
	x = 5;
	y = 5;
	z = 1;
	

	// Reset the displaced ('displaced path') coordinates for the photon.
	x_disp = 5;
	y_disp = 5;
	z_disp = 1;


	r = 0;
	step = 0;
	step_remainder = 0;
	
	// Reset the number of interactions back to zero.
	num_steps = 0;
	
	// Reset the path lengths of the photon.
	original_path_length = 0;
	displaced_path_length = 0;

	// 'cnt' represents the number of times a photon has propogated
	// through the medium.
	cnt = 0;

	// Randomly set photon trajectory to yield isotropic or anisotropic source.
	initTrajectory();
	
	
}


// Update the step size for the photon.
void Photon::setStepSize()
{

	// Update the current values of the absorption and scattering coefficients
	// based on the depth in the medium (i.e. which layer the photon is in).
	//double mu_a = m_medium->getLayerAbsorptionCoeff(z);  // cm^-1
	//double mu_s = m_medium->getLayerScatterCoeff(z);	  // cm^-1

	double mu_a = 1.0;		// cm^-1
	double mu_s = 100.0;		// cm^-1

	// If last step put the photon on the layer boundary
	// then we need a new step size.  Otherwise, the step
	// is set to remainder size and update step_remainder to zero.
	if (step_remainder == 0.0) {
		double rnd = getRandNum();
		// Calculate the new step length of the photon.
		step = -log(rnd)/(mu_a	+ mu_s);
	}
	else
	{
		step = step_remainder;
		step_remainder = 0.0;
	}
}


double Photon::getPathLength(double x_dist, double y_dist, double z_dist)
{
	return sqrt(pow(x_dist, 2) + pow(y_dist, 2) + pow(z_dist, 2));
}


void Photon::captureLocationCoords(void)
{
	// Add the coordinates to the STL vector for the displaced and
	// non-displaced case.
	//	coords.push_back(x);
	//	coords.push_back(y);
	//	coords.push_back(z);
	coords.push_back(x_disp);
	coords.push_back(y_disp);
	coords.push_back(z_disp);
}

// Step photon to new position.
void Photon::hop()
{
#ifdef DEBUG
	cout << "Hopping...\n";
#endif	
	

	// Record the location of the photon during this interaction.
	captureLocationCoords();

	// Locations before the photon is moved.
	double temp_x = x;
	double temp_y = y;
	double temp_z = z;

	
	// Update position of the photon.
	x += step*dirx;
	y += step*diry;
	z += step*dirz;

	// Update positon for the 'displaced' position of the photon.  Allows us to track
	// displaced and non-displaced separately.
	x_disp += step*dirx;
	y_disp += step*diry;
	z_disp += step*dirz;

	// Calculate the path length of the photon WITHOUT displacement.
	original_path_length += getPathLength((x - temp_x), (y - temp_y), (z - temp_z));

	// Move the photon to the new position based on the displacement of from
	// the ultrasound wave.
	this->displacePhotonFromPressure();

	// Calculate the path length of the photon WITH displacement.
	displaced_path_length += getPathLength((x_disp - temp_x), (y_disp - temp_y), (z_disp - temp_z));

}


// Return absorbed energy from photon weight at this location.
void Photon::drop()
{
#ifdef DEBUG
	cout << "Dropping...\n";
#endif	
	
	// Update the current values of the absorption and scattering coefficients
	// based on the depth in the medium (i.e. which layer the photon is in).
	//double mu_a = m_medium->getLayerAbsorptionCoeff(z);  // cm^-1
	//double mu_s = m_medium->getLayerScatterCoeff(z);	  // cm^-1
	
    double mu_a = 1.0;		// cm^-1
	double mu_s = 100.0;		// cm^-1
    
	// Calculate the albedo and remove a portion of the photon's weight for this
	// interaction.
	double albedo = mu_s / (mu_a + mu_s);
	double absorbed = weight * (1 - albedo);
	
	// Remove the portion of energy lost due to absorption at this location.
	weight -= absorbed;
	
	// Deposit lost energy in the grid of the medium.
	//m_medium->absorbEnergy(z, absorbed);


	double r = fabs(z);
	int ir = (r/radial_bin_size);
	if (ir >= num_radial_pos) {
		ir = num_radial_pos;
	}
	// Deposit lost energy of the photon into the local detection grid.
	local_Cplanar[ir] += absorbed;


}

// Calculate the new trajectory of the photon.
void Photon::spin()
{
#ifdef DEBUG
	cout << "Spinning...\n";
#endif	
	
	// Get the anisotropy factor from the layer that resides at depth 'z' in
	// the medium.
	//double g = m_medium->getAnisotropyFromDepth(z);
	double g = 0.90;
    
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


// Write the coordinates of each scattering event of the photon
// to file for postprocessing with matlab.
void Photon::writeCoordsToFile(void)
{
	ofstream outfile("photon-paths.txt");

	// Iterate over all the coordinates in the stl vector and write
	// to file.
	for (int i = 0; i < coords.size(); i++)
		outfile << coords[i] << " ";

	outfile.close();
}


// FIXME: CURRENTLY ONLY DISPLACING IN ONE DIRECTION.  SHOULD USE A TENSOR.
void Photon::displacePhotonFromPressure(void)
{
	// Get the local pressure from the grid based on the coordinate of the photon.
	double pressure = m_medium->getPressureFromCartCoords(x, z, y);

	// Impedance of the tissue.
	double impedance = 1.63e6;

	// FIXME: THIS SHOULD NOT BE HARD CODED.
	// Frequency of the ultrasound transducer.
	double freq = 5e6;

	// Displace in the z-axis based on the pressure.
	double displacement = (abs(pressure)*1e6)/(2*Pi*freq*impedance);
	z_disp += displacement;

}


// XXX: Currently not in use
void Photon::specularReflectance(double n1, double n2)
{
	// update the weight after specular reflectance has occurred.
	weight = weight - (pow((n1 - n2), 2) / pow((n1 + n2), 2)) * weight;
}


void Photon::transmitOrReflect(void)
{
	// XXX: FINISH ME
}


bool Photon::hitLayerBoundary(void)
{
	// 1)Find distance to layer boundary where there could potentially be
	//   refractive index mismatches.
	// 2) If distance to boundary is less than the current step_size of the
	//   photon, we move the photon to the boundary and keep track of how
	//   much distance is left over from step_size - distance_to_boundary.

	double distance_to_boundary = 0.0;
	Layer *l = m_medium->getLayerFromDepth(z);


	// If the direction the photon is traveling is towards the deeper boundary
	// of the layer, we execute the first clause.  Otherwise we are moving to
	// the more superficial boundary of the layer.
	if (dirz > 0.0)
	{
		distance_to_boundary = (l->getDepthEnd() - z) / dirz;
	}
	else if (dirz < 0.0)
	{
		distance_to_boundary = (l->getDepthStart() - z) / dirz;
	}

	// If the step size of the photon is larger than the distance
	// to the boundary and we are moving in some direction along
	// the z-axis (i.e. not parallel to the layer boundary) we calculate
	// the left over step size and then step the photon to the boundary.
	if (dirz != 0.0 && step > distance_to_boundary)
	{
		step_remainder = (step - distance_to_boundary)*l->getTotalAttenuationCoeff();
		step = distance_to_boundary;
		return true;
	}
	else
	{
		return false;
	}
}



// Everything below deals with the random number generator.
unsigned int Photon::TausStep(unsigned int &z, int s1, int s2, int s3, unsigned long M)
{
	unsigned int b=(((z << s1) ^ z) >> s2);
	z = (((z & M) << s3) ^ b);
	return z;
}


unsigned int Photon::LCGStep(unsigned int &z, unsigned int A, unsigned long C)
{
	return z=(A*z+C);
}

double Photon::HybridTaus(void)
{
	// Combined period is lcm(p1,p2,p3,p4)~ 2^121
	return 2.3283064365387e-10 * (              // Periods for the RNG.
			TausStep(z1, 13, 19, 12, 4294967294UL) 	^  // p1=2^31-1
			TausStep(z2, 2, 25, 4, 4294967288UL) 	^    // p2=2^30-1
			TausStep(z3, 3, 11, 17, 4294967280UL) 	^   // p3=2^28-1
			LCGStep(z4, 1664525, 1013904223UL)        // p4=2^32
	);
}


double Photon::getRandNum(void)
{
	// Thread safe RNG.
	return HybridTaus();




	// Non-thread safe RNG
//	double rnd = (double)rand()/(double)RAND_MAX;
//	while ((rnd == 0) || (rnd == 1)) { // produces 0 < rnd < 1
//		rnd = (double)rand()/(double)RAND_MAX;
//	}
//	return rnd;




	// FIXME:  Using the Boost Random Library is MUCH slower when generating
	//			random numbers.  Questions to answer,
	//			- Is it the algorithm used (i.e. Mersenne-twister)?
	//			- Does the creation and destruction of these objects
	//			  below slow things down?  That is, should there be
	//			  creation of them on the heap with pointers so they
	//			  stay in existence?  Even possible?
	//			Until these questions are answered and a speedup is found
	//			this version of the simulation will be using the built in
	//          RNG above.
	// Boost implementation of the Mersenne-twister RNG.
//	boost::uniform_real<> dist(0, 1);
//	boost::variate_generator<boost::mt19937&, boost::uniform_real<> > rand_num(gen, dist);
//	return rand_num();
}
