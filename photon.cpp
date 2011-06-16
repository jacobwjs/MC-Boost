
#include "debug.h"
#include "photon.h"

#undef DEBUG

Photon::Photon(void)
{
#ifdef DEBUG
	cout << "Creating Photon...\n";
#endif
    
	this->initCommon();
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


// Common initialization function.
void Photon::initCommon(void)
{
    // Photon just created, so it is alive.
	status = ALIVE;
    
	// Set back to initial weight values.
	weight = 1;
    
	r = 0;
	step = 0;
	step_remainder = 0;
    
	// Reset the number of interactions back to zero.
	num_steps = 0;
    
	// 'cnt' represents the number of times a photon has propogated
	// through the medium.
	cnt = 0;
    
    // Reset the flags for hitting a layer boundary.
	hit_x_bound = hit_y_bound = hit_z_bound = false;
    
    // Reset the transmission angle for a photon.
    transmission_angle = 0;
    
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
                          unsigned int state3, unsigned int state4, InjectionCoords &coords)
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
    
    // Set the current layer the photon starts propagating through.  This will
    // be updated as the photon moves through layers by checking 'hitLayerBoundary'.
    currLayer = m_medium->getLayerFromDepth(z);
    
    
	// Assign local values of the detection grid from the Medium.
	radial_bin_size = m_medium->getRadialBinSize();
	num_radial_pos = m_medium->getNumRadialPos();
    
    // Set the location of illumination source and the initial cartesian coordinates of the photon.
    this->x = this->illuminationCoords.x = coords.x;
    this->y = this->illuminationCoords.y = coords.y;
    this->z = this->illuminationCoords.z = coords.z;
    
    // Move the photon through the medium.
    propagatePhoton(iterations);
	
    
}


void Photon::propagatePhoton(const int iterations)
{
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
            
            
			// Make various checks on the photon to see if layer or medium boundaries
			// are hit and whether the photon should be transmitted or reflected.
			
			
            
            
            if (hitLayerBoundary())
			{
				cout << "Hit layer boundary\n";
                
				hop();	// Move the photon to the layer boundary.
				transmitOrReflect("layer");
			}
            else if (hitMediumBoundary())
			{
				
				cout << "Hit medium boundary\n";
				hop();  // Move the photon to the medium boundary.
				transmitOrReflect("medium");
			}
			else
			{
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
                
			}
            
            
		} // end while() loop
        
		// Reset the photon and start propagation over from the beginning.
		reset();
        
	} // end for() loop
    
    
	// This thread has executed all of it's photons, so now we update the global
	// absorption array in the medium.
	m_medium->absorbEnergy(local_Cplanar);
    
}



void Photon::plotPath(void)
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
	x = illuminationCoords.x;
    y = illuminationCoords.y;
    z = illuminationCoords.z;
    
	r = 0;
	step = 0;
	step_remainder = 0;
    
	// Reset the number of interactions back to zero.
	num_steps = 0;
    
	// 'cnt' represents the number of times a photon has propogated
	// through the medium.
	cnt = 0;
    
    // Reset the flags for hitting a layer boundary.
	hit_x_bound = hit_y_bound = hit_z_bound = false;
    
    // Reset the transmission angle for a photon.
    transmission_angle = 0;
    
	// Randomly set photon trajectory to yield isotropic source.
	initTrajectory();

    // Reset the current layer from the injection coordinates of the photon.
    currLayer = m_medium->getLayerFromDepth(z);

}


// Update the step size for the photon.
void Photon::setStepSize()
{
    
	// Update the current values of the absorption and scattering coefficients
	// based on the depth in the medium (i.e. which layer the photon is in).
	//double mu_a = m_medium->getLayerAbsorptionCoeff(z);  // cm^-1
	//double mu_s = m_medium->getLayerScatterCoeff(z);	  // cm^-1
    double mu_a = currLayer->getAbsorpCoeff();
    double mu_s = currLayer->getScatterCoeff();
    
	//double mu_a = 1.0;		// cm^-1
	//double mu_s = 100.0;		// cm^-1
    
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
		step = step_remainder/(mu_a + mu_s);
		step_remainder = 0.0;
	}
}

// Step photon to new position.
void Photon::hop()
{
#ifdef DEBUG
	cout << "Hopping...\n";
#endif	
    
	//setStepSize();
    
	// Update position of the photon.
	x += step*dirx;
	y += step*diry;
	z += step*dirz;
}


// Return absorbed energy from photon weight at this location.
void Photon::drop()
{
#ifdef DEBUG
	cout << "Dropping...\n";
#endif	
    
	// Update the current values of the absorption and scattering coefficients
	// based on the depth in the medium (i.e. which layer the photon is in).
	double mu_a = m_medium->getLayerAbsorptionCoeff(z);  // cm^-1
	double mu_s = m_medium->getLayerScatterCoeff(z);	  // cm^-1
    
	//double mu_a = 1.0;		// cm^-1
	//double mu_s = 100.0;		// cm^-1
    
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
            cout << "Photon died in Roulette\n";
			status = DEAD;
		}
	}
}



// XXX: Currently not in use
void Photon::specularReflectance(double n1, double n2)
{
	// update the weight after specular reflectance has occurred.
	weight = weight - (pow((n1 - n2), 2) / pow((n1 + n2), 2)) * weight;
}



// XXX: Finish me
void Photon::transmit(const char *type)
{
    if (strcmp("layer", type) == 0)
    {
        cout << "Transmitting through layer\n";
        dirz = cos(transmission_angle);
        // If we transmit through the layer to another we must update
        // the current layer pointer of the photon so it will correctly 
        // calculate the next step size.
        if (dirz > 0)
            currLayer = m_medium->getLayerBelowCurrent(z);
        else
            currLayer = m_medium->getLayerAboveCurrent(z);
        
    }
    else if (strcmp("medium", type) == 0)
    {
        cout << "Transmitting through medium boundary\n";
        this->status = DEAD;
    }
    
}


void Photon::transmitOrReflect(const char *type)
{
	// Test whether to transmit the photon or reflect it.  If the reflectance is
	// greater than the random number generated, we internally reflect, otherwise
	// transmit.
	// Case when interaction is with a layer boundary.
	if (strcmp("layer", type) == 0)
	{
		// Stochastically determine if the photon should be transmitted or reflected.
		if (getLayerReflectance() > getRandNum())
		{
            cout << "Internally reflecting\n";
			internallyReflectZ();
            
			// Since the photon has interacted with the tissue we deposit weight.
			drop();
            
			// Perform roulette to determine if the photon should continue propagation.
			// NOTE: spin() is not performed since the photon has been internally reflected
			//	   	and the new direction has been taken care of by internallyReflect().
			performRoulette();
            
		}
		else
		{
			transmit("layer");
		}
	}
	// Case when interaction is with a medium boundary.
	else if (strcmp("medium", type) == 0)
	{
		// Stochastically determine if the photon should be transmitted or reflected.
		if (getMediumReflectance() > getRandNum())
		{
            cout << "Reflecting photon on medium boundary\n";
			// Depending on which medium boundary was hit, we reflect on that axis,
            // change the direction of the direction cosign, and reset the boolean flag.
            if (hit_x_bound)
            {
				internallyReflectX();
                hit_x_bound = false;
            }
			else if (hit_y_bound)
            {
				internallyReflectY();
                hit_y_bound = false;
            }
			else if (hit_z_bound)
            {
				internallyReflectZ();
                hit_z_bound = false;
            }
            else
            {
                cout << "Error, no medium boundary hit\n";
            }
            
            

            
            
			// Since the photon has interacted with the tissue we deposit weight.
			drop();
            
			// Perform roulette to determine if the photon should continue propagation.
			// NOTE: spin() is not performed since the photon has been internally reflected
			//	   	and the new direction has been taken care of by internallyReflect().
			performRoulette();
		}
        else
        {
            transmit("medium");
        }
        
        
	}
	else
	{
		cout << "Error: transmitOrReflect()\n";
	}
}



// XXX: *** Need to verify the logic below is correct ***
double Photon::getMediumReflectance(void)
{
	// Sanity check.
	assert((hit_x_bound == true) ||
           (hit_y_bound == true) ||
           (hit_z_bound == true));
    
    
	//Layer *currLayer = m_medium->getLayerFromDepth(z);
	double refract_index_n1 = currLayer->getRefractiveIndex();	// Current layer's refractive index.
	double refract_index_n2 = 1.0;	// Outside of the medium is only air.
    
    
    
    
	double direction = 0.0;
	if (hit_x_bound)
		direction = dirx;
	else if (hit_y_bound)
		direction = diry;
	else
		direction = dirz;
    
	// Calculate the incident angle based on the axis in which the photon hit the medium.
	double incident_angle = acos(abs(direction));
    
	// Calculate the critical angle.
	double critical_angle = asin(refract_index_n2 / refract_index_n1);
    
	// If the incident angle is larger than the critical angle, the reflectance
	// is set to 1, otherwise the reflectance is calculated from Fresnel's equation.
	if (incident_angle > critical_angle)
	{
		reflectance = 1;
	}
	else
	{
        // Calculate the transmission angle through the layer boundary.
		this->transmission_angle = asin(refract_index_n1/refract_index_n2 * sin(incident_angle));
		
        // Calcualte the Fresnel reflection 
        reflectance = 0.5 * (pow(sin(incident_angle-transmission_angle), 2)/pow(sin(incident_angle+transmission_angle), 2)
                             + pow(tan(incident_angle-transmission_angle), 2)/pow(tan(incident_angle+transmission_angle), 2));
	}
    
    
	return reflectance;
    
}


double Photon::getLayerReflectance(void)
{
	double refract_index_n1 = 0.0;	// Current layer's refractive index.
	double refract_index_n2 = 0.0;	// Next layer's refractive index.
	//Layer *currLayer = m_medium->getLayerFromDepth(z);
	Layer *nextLayer;
    
	double incident_angle = acos(abs(dirz));
	refract_index_n1 = currLayer->getRefractiveIndex();
    
	// If the photon is moving towards a deeper layer.
	if (dirz > 0)
	{
		nextLayer = m_medium->getLayerBelowCurrent(z);
	}
	// If the photon is moving towards a more shallow layer.
	else if (dirz < 0)
	{
		nextLayer = m_medium->getLayerAboveCurrent(z);
	}
	// Perpendicular propagation.
	else
	{
		// FIXME:
		// This is where propagation is normal to the boundary/layer
		// and should be transmitted regardless.
		cout << "Photon should transmit...\n";
	}
    
    
    // If the layer above or below is outside of the medium
    // we assign the refractive index to be air, otherwise
    // use the value from the layer.
    if (nextLayer == NULL)
        refract_index_n2 = 1.0;
    else
        refract_index_n2 = nextLayer->getRefractiveIndex();

    
	// Calculate the critical angle.
    if (refract_index_n2 > refract_index_n1)
    {
        // For specular refection we always remove some portion of the weight and transmit the photon
        // to the next layer.
        transmission_angle = asin(refract_index_n1/refract_index_n2 * sin(incident_angle));
        specularReflectance(refract_index_n1, refract_index_n2);
        //transmit("layer");
        // Since we know we transmit, set reflectance to zero.
        reflectance = 0;
    }
    else
    {
        double critical_angle = asin(refract_index_n2 / refract_index_n1);
        // If the incident angle is larger than the critical angle, the reflectance
        // is set to 1, otherwise the reflectance is calculated from Fresnel's equation.
        if (incident_angle > critical_angle)
        {
            reflectance = 1;
        }
        else
        {
            transmission_angle = asin(refract_index_n1/refract_index_n2 * sin(incident_angle));
            
            reflectance = 0.5 * (pow(sin(incident_angle-transmission_angle), 2)/pow(sin(incident_angle+transmission_angle), 2)
                                 + pow(tan(incident_angle-transmission_angle), 2)/pow(tan(incident_angle+transmission_angle), 2));
        }
    }
    
	
    
    
	return reflectance;
}



// Note: We take the absolute value in the case where the direction 
//       cosine is negative, since the distance to boundary would be
//       negative otherwise, which is untrue.  This arises due to assuming
//       the lower axis bound in each dimension (x, y, z) begins at zero.
//       This could also be achieved by simply subtracting the current location
//       from zero (e.g. 0-y/diry), which would change the sign as well.
bool Photon::hitMediumBoundary(void)
{
	double distance_to_boundary = 0.0;
	//Layer *layer = m_medium->getLayerFromDepth(z);
	double mu_t = currLayer->getTotalAttenuationCoeff();
	double x_bound = m_medium->getXbound();
	double y_bound = m_medium->getYbound();
	double z_bound = m_medium->getZbound();
    
	// Check if the photon has been given a step size past the outer edges of the medium.
	// If so we calculate the distance to the boundary.
	
    
    if (step*dirx + x >= x_bound || step*dirx + x <= 0)
	{
		hit_x_bound = true;
		if (dirx > 0) // Moving towards positive x_bound
			distance_to_boundary = (x_bound - x) / dirx;
		else
			distance_to_boundary = abs(x / dirx);
	}
	else if	(step*diry + y >= y_bound || step*diry + y <= 0)
	{
		hit_y_bound = true;
		if (diry > 0) // Moving towards positive y_bound
			distance_to_boundary = (y_bound - y) / diry;
		else
			distance_to_boundary = abs(y / diry);
	}
	else if (step*dirz + z >= z_bound || step*dirz + z <= 0)
	{
		hit_z_bound = true;
		if (dirz > 0) // Moving towards positive y_bound
			distance_to_boundary = (z_bound - z) / dirz;
		else
			distance_to_boundary = abs(z / dirz);
	}
	// No boundaries have been crossed, so return false.
	else
	{
		return false;
	}

    
	// If the step size of the photon is larger than the distance
	// to the boundary and we are moving in some direction along
	// the axis we calculate the left over step size and then step
	// the photon to the boundary with the next call to 'hop()'.
	if (step > distance_to_boundary && distance_to_boundary != 0)
	{
		step_remainder = (step - distance_to_boundary)*mu_t;
		step = distance_to_boundary;
		return true;
	}
	else
	{
		return false;
	}
}


bool Photon::hitLayerBoundary(void)
{
	// 1)Find distance to layer boundary where there could potentially be
	//   refractive index mismatches.
	// 2) If distance to boundary is less than the current step_size of the
	//   photon, we move the photon to the boundary and keep track of how
	//   much distance is left over from step_size - distance_to_boundary.
    
	double distance_to_boundary = 0.0;
	//Layer *layer = m_medium->getLayerFromDepth(z);
	double mu_t = currLayer->getTotalAttenuationCoeff();
    
    
    
	// If the direction the photon is traveling is towards the deeper boundary
	// of the layer, we execute the first clause.  Otherwise we are moving to
	// the more superficial boundary of the layer.
	if (dirz > 0.0)
	{
		distance_to_boundary = (currLayer->getDepthEnd() - z) / dirz;
	}
	else if (dirz < 0.0)
	{
		distance_to_boundary = (currLayer->getDepthStart() - z) / dirz;
	}
    
    
	// If the step size of the photon is larger than the distance
	// to the boundary and we are moving in some direction along
	// the z-axis (i.e. not parallel to the layer boundary) we calculate
	// the left over step size and then step the photon to the boundary.
	if (dirz != 0.0 && step > distance_to_boundary)
	{
		step_remainder = (step - distance_to_boundary)*mu_t;
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
