
#include "debug.h"
#include "logger.h"
#include "absorber.h"
#include "vector3D.h"
#include "layer.h"
#include "medium.h"
#include "photon.h"
#include "displacementMap.h"
#include "pressureMap.h"



#undef DEBUG

static int cnt = 0;

Photon::Photon(void)
{
#ifdef DEBUG
	cout << "Creating Photon...\n";
#endif
    
    // The current location of the photon.
    currLocation = boost::shared_ptr<Vector3d> (new Vector3d);
    currLocation->withDirection();  // Enable direction for this vector.
    
    // The location before the photon hopped.
    prevLocation = boost::shared_ptr<Vector3d> (new Vector3d);  // Does not require direction.
    prevLocation->withDirection();
    
	this->initCommon();
}



Photon::Photon(double x, double y, double z,
               double dirx, double diry, double dirz)
{
#ifdef DEBUG
	cout << "Constructor stub...\n";
#endif
    
    coords location;
    directionCos direction;
    
    location.x = x;
    location.y = y;
    location.z = z;
    
    direction.x = dirx;
    direction.y = diry;
    direction.z = dirz;
    
    // Location and direction of photon.
    currLocation = boost::shared_ptr<Vector3d> (new Vector3d(location, direction));
    
    // The location before the photon hopped.
    prevLocation = boost::shared_ptr<Vector3d> (new Vector3d);  // Does not require direction.
    
    this->initCommon();

}


Photon::~Photon(void)
{
#ifdef DEBUG	
	cout << "Destructing Photon...\n";
#endif
    cout << "cnt = " << cnt << "\n";
}


// Common initialization function.
void Photon::initCommon(void)
{
    // Photon just created, so it is alive.
	status = ALIVE;
    
	// Set to initial weight values.
	weight = 1;
    
    // Default value of tagged to false.
    tagged = false;
    
	r = 0;
	step = 0;
	step_remainder = 0;
    
	// Set the number of interactions back to zero.
	num_steps = 0;
    
    
    // Set the flags for hitting a layer boundary.
	hit_x_bound = hit_y_bound = hit_z_bound = false;
    
    // Set the transmission angle for a photon.
    transmission_angle = 0.0f;
    
    // Set the path lengths during initialization.
    unmodulated_optical_path_length = 0.0f;
    displaced_optical_path_length = 0.0f;
    refractiveIndex_optical_path_length = 0.0f;
    
}


// Set the number of iterations this thread will run.
void Photon::setIterations(const int num)
{
	iterations = num;
}


void Photon::initTrajectory(void)
{
	// Randomly set photon trajectory to yield anisotropic source.
	cos_theta = (2.0 * getRandNum()) - 1;
	sin_theta = sqrt(1.0 - cos_theta*cos_theta);
	psi = 2.0 * PI * getRandNum();
    
    // Set the initial direction cosines for this photon.
    currLocation->setDirX(sin_theta * cos(psi));
    currLocation->setDirY(sin_theta * sin(psi));
    currLocation->setDirZ(1.0f);    
    
}

void Photon::initAbsorptionArray()
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
    
    // Initialize the seeds and also save them in case this photon
    // exits through the aperture and we want to write them out to
    // disk.
	seeds.s1 =  z1 = state1;
	seeds.s2 =  z2 = state2;
	seeds.s3 = 	z3 = state3;
	seeds.s4 = 	z4 = state4;
}


// BOOST thread library starts execution here.
// 1) Hop - move the photon
// 2) Drop - drop weight due to absorption
// 3) Spin - update trajectory accordingly
// 4) Roulette - test to see if photon should live or die.
void Photon::injectPhoton(Medium *medium, const int iterations, unsigned int state1, unsigned int state2,
                          unsigned int state3, unsigned int state4, coords &laser)
{
	// seed the random number generator.
	//srand(time(0) + thread_id);
    
	// Seed the Boost RNG (Random Number Generator).
	//gen.seed(time(0) + thread_id);
    
	// Initialize the photon's properties before propagation begins.
	initTrajectory();
	initAbsorptionArray();
	initRNG(state1, state2, state3, state4);
    
	// Before propagation we set the medium which will be used by the photon.
	this->m_medium = medium;
    
	// Assign local values of the detection grid from the Medium.
	radial_bin_size = m_medium->getRadialBinSize();
	num_radial_pos = m_medium->getNumRadialPos();
    
    // Set the location of illumination source and the initial cartesian coordinates of the photon
    // when it is first incident on the medium.
    currLocation->location.x = this->illuminationCoords.x = laser.x;
    currLocation->location.y = this->illuminationCoords.y = laser.y;
    currLocation->location.z = this->illuminationCoords.z = laser.z;
    
    // Set the current layer the photon starts propagating through.  This will
    // be updated as the photon moves through layers by checking 'hitLayerBoundary'.
    currLayer = m_medium->getLayerFromDepth(currLocation->location.z);
    
    // Move the photon through the medium. 'iterations' represents the number of photons this
    // object (which is a thread) will execute.
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
            
            // Flags for testing if a photon hit/passed through a layer
            // or medium boundary.
			//bool hitLayer = checkLayerBoundary();
            bool hitMedium = checkMediumBoundary();
            
            
            
			//if (!hitLayer && !hitMedium)
            if (!hitMedium)
			{
                // sanity check.
                assert(this->status == ALIVE);
                
				// Move the photon in the medium.
				hop();
                
                
                // Now displace the photon at its new location some distance depending on
                // how the pressure has moved scattering particles and/or due to the change
                // in the path of the photon due to refractive index gradient.
                //displacePhotonFromPressure();
                
				// Drop weight of the photon due to an interaction with the medium.
				drop();
                
				// Calculate the new coordinates of photon propagation.
				spin();
                
				// Test whether the photon should continue propagation from the
				// Roulette rule.
				performRoulette();
                
			}
            
            
		} // end while() loop
        
        
        
		// Write out the x,y,z coordinates of the photons path as it propagated through
		// the medium.
		//writeCoordsToFile();
        
        
		// Reset the photon and start propogation over from the beginning.
		reset();
        
	} // end for() loop
    
    
	// This thread has executed all of it's photons, so now we update the global
	// absorption array in the medium.
	//m_medium->absorbEnergy(local_Cplanar);
    
}



void Photon::plotPath(void)

{
	// STUB
}




void Photon::reset()
{
#ifdef DEBUG
	cout << "Resetting Photon...\n";
#endif
    
	// Photon just created, so it is alive.
	status = ALIVE;
    
	// Set back to initial weight values.
	weight = 1;
    
    // Set tagged boolean back to false.
    tagged = false;
    
    // Set the vector that contains the current location of the photon.
	currLocation->location.x = illuminationCoords.x;
    currLocation->location.y = illuminationCoords.y;
    currLocation->location.z = illuminationCoords.z;
    
	r = 0;
	step = 0;
	step_remainder = 0;
    
	// Reset the number of interactions back to zero.
	num_steps = 0;
    
	// Reset the path lengths of the photon.
	unmodulated_optical_path_length = 0;
	displaced_optical_path_length = 0;
	refractiveIndex_optical_path_length = 0;


    // Reset the flags for hitting a layer boundary.
	hit_x_bound = false;
    hit_y_bound = false;
    hit_z_bound = false;
    
    // Reset the transmission angle for a photon.
    transmission_angle = 0;
    
	// Randomly set photon trajectory to yield isotropic or anisotropic source.
	initTrajectory();
    
    // Reset the current layer from the injection coordinates of the photon.
    currLayer = m_medium->getLayerFromDepth(currLocation->location.z);
    
    // Since the photon is being restarted we save the current state of the RNG
    // as these are our 'seeds' for this run.
    //
    seeds.s1 = z1;
    seeds.s2 = z2;
    seeds.s3 = z3;
    seeds.s4 = z4;
}


// Update the step size for the photon.
void Photon::setStepSize()
{
	// Update the current values of the absorption and scattering coefficients
	// based on the depth in the medium (i.e. which layer the photon is in).
    double mu_a = currLayer->getAbsorpCoeff(currLocation);
    double mu_s = currLayer->getScatterCoeff(currLocation);
    
    
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



void Photon::captureLocationCoords(void)
{
	cout << "Photon::captureLocationCoords() stub\n";
	// Add the coordinates to the STL vector for the displaced scattering locations.
    //	coords.push_back(x_disp);
    //	coords.push_back(y_disp);
    //	coords.push_back(z_disp);
}


void Photon::captureExitCoordsAndLength(void)
{
	// Add the coordinates to the STL vector for the displaced scattering locations
	// and the displaced length.
    //	photon_exit_data.push_back(x_disp);
    //	photon_exit_data.push_back(y_disp);
    //	photon_exit_data.push_back(displaced_path_length);
    cout << "Photon::captureExitCoordsAndLength() stub\n";
}


void Photon::captureExitCoordsLengthWeight(void)
{
	// Add the coordinates to the STL vector for the displaced scattering locations
	// and the displaced length.
    //	photon_exit_data.push_back(x_disp);
    //	photon_exit_data.push_back(y_disp);
    //	photon_exit_data.push_back(displaced_path_length);
    //	photon_exit_data.push_back(weight);
    cout << "Photon::captureExitCoordsLengthWeight() stub\n";
}


// Tests if the photon will come into contact with a layer boundary
// after setting the new step size.  If so the process of transmitting or
// reflecting the photon begins.
bool Photon::checkLayerBoundary(void)
{
    if (hitLayerBoundary())
    {
#ifdef DEBUG
        cout << "Hit layer boundary\n";
#endif
        
        hop();	// Move the photon to the layer boundary.
        transmitOrReflect("layer");
        return true;
    }
    
    return false;
}


bool Photon::checkMediumBoundary(void)
{
    if (hitMediumBoundary())
    {
#ifdef DEBUG			
        cout << "Hit medium boundary\n";
#endif
        hop();  // Move the photon to the medium boundary.
        transmitOrReflect("medium");
        return true;
    }
    
    return false;
}


// Check if the photon passed through the detection window.
// Returns the number of detectors that were crossed in this
// hop in the case of multiple detectors present.
bool Photon::checkDetector(void)
{
    int cnt =  m_medium->photonHitDetectorPlane(currLocation);
    // If cnt > 0 the photon exited through the bounds of the detector.
    if (cnt > 0) 
    {
        return true;
    }
    else
        return false;
}


// Step photon to new position.
void Photon::hop()
{
#ifdef DEBUG
	cout << "Hopping...\n";
#endif	
    
	num_steps++;
    
    
	// Save the location before making the hop.
	prevLocation->location.x = currLocation->location.x;
    prevLocation->location.y = currLocation->location.y;
    prevLocation->location.z = currLocation->location.z;
    prevLocation->setDirX(currLocation->getDirX());
    prevLocation->setDirY(currLocation->getDirY());
    prevLocation->setDirZ(currLocation->getDirZ());

    
	// Update the location
	currLocation->location.x += step * currLocation->getDirX();
	currLocation->location.y += step * currLocation->getDirY();
	currLocation->location.z += step * currLocation->getDirZ();


	// Calculate the unmodulated path length, which ignores the
	// pressure and displacement maps.  This is simply the background
	// refractive index.
	//unmodulated_optical_path_length += VectorMath::Distance(currLocation, prevLocation)*currLayer->getRefractiveIndex();;

}


// Return absorbed energy from photon weight at this location.
void Photon::drop()
{
#ifdef DEBUG
	cout << "Dropping...\n";
#endif	
    
    if (this->status == DEAD) return;
    
    
    
    
    double mu_a = 0.0f;
    double mu_s = 0.0f;
    double albedo = 0.0f;
    double absorbed = 0.0f;
    
    
    Absorber * absorber = currLayer->getAbsorber(currLocation);
    // If an absorber was returned, then we get the absorption and
    // scattering coefficients from it.  Otherwise we use the values
    // from the background layer.
    if (absorber != NULL)
    {
        mu_a = absorber->getAbsorberAbsorptionCoeff();
        mu_s = absorber->getAbsorberScatteringCoeff();
        
        // Calculate the albedo and remove a portion of the photon's weight for this
        // interaction.
        albedo = mu_s / (mu_a + mu_s);
        absorbed = weight * (1 - albedo);
        
        // Update the absorbed weight in this absorber.
        absorber->updateAbsorbedWeight(absorbed);
        
        // If this photon hit an absorber we set tagged to true, which
        // assumes our tagging volume completely encompasses the absorber
        // and is the same shape.
        tagged = true;
    }
    else
    {
        // Update the current values of the absorption and scattering coefficients
        // based on the depth in the medium (i.e. which layer the photon is in).
    	// NOTE:
    	// - No need to index into the layer and see if absorption and scattering coefficients
    	//   should be pulled from absorber, because we verified above in the if() that this was
    	//   not the case.  Saves a small amount of time searching through the absorbers.
        mu_a = currLayer->getAbsorpCoeff();
        mu_s = currLayer->getScatterCoeff();
        
        // Calculate the albedo and remove a portion of the photon's weight for this
        // interaction.
        albedo = mu_s / (mu_a + mu_s);
        absorbed = weight * (1 - albedo);
    }
    
    
    
	// Remove the portion of energy lost due to absorption at this location.
	weight -= absorbed;
    
	// Deposit lost energy in the grid of the medium.
	//m_medium->absorbEnergy(z, absorbed);
    
    
	// Deposit lost energy.
    //updateLocalWeightArray(absorbed);
    
}


// Update the local absorbed energy array.
void Photon::updateLocalWeightArray(const double absorbed)
{
    double r = fabs(currLocation->location.z);
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
    
    if (this->status == DEAD) return;
    
	// Get the anisotropy factor from the layer that resides at depth 'z' in
	// the medium.
	// FIXME: Need to index into layer and check if absorber causes this to change.
    double g = currLayer->getAnisotropy();
    
	double rnd = getRandNum();
    
    double dirZ = currLocation->getDirZ();
    double dirY = currLocation->getDirY();
    double dirX = currLocation->getDirX();
    
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
	if (1 - fabs(dirZ) <= ONE_MINUS_COSZERO) {      /* close to perpendicular. */
		uxx = sin_theta * cos_psi;
		uyy = sin_theta * sin_psi;
		uzz = cos_theta * SIGN(dirZ);   /* SIGN() is faster than division. */
	} 
	else {					/* usually use this option */
        
        
		temp = sqrt(1.0 - dirZ * dirZ);
		uxx = sin_theta * (dirX * dirZ * cos_psi - dirY * sin_psi) / temp + dirX * cos_theta;
		uyy = sin_theta * (dirY * dirZ * cos_psi + dirX * sin_psi) / temp + dirY * cos_theta;
		uzz = -sin_theta * cos_psi * temp + dirZ * cos_theta;
	}
    
	// Update trajectory.
	currLocation->setDirX(uxx);
	currLocation->setDirY(uyy);
	currLocation->setDirZ(uzz);
}



void Photon::performRoulette(void)
{
    // Photon has already been killed, presumably by leaving the medium.
    if (this->status == DEAD) return;
    
    
	if (weight < THRESHOLD) {
		if (getRandNum() <= CHANCE) {
			weight /= CHANCE;
		}
		else {
#ifdef DEBUG
            cout << "Photon died in Roulette\n";
#endif
			status = DEAD;
		}
	}
}


// Write the coordinates of each scattering event of the photon
// to file for postprocessing with matlab.
void Photon::writeCoordsToFile(void)
{
	cout << "Photon::writeCoordsToFile() stub\n";
    
	//	cout << "Non-displaced path length: " << scientific << setprecision(12) <<  unmodulated_path_length << endl;
	//	cout << "Displaced path length: " << scientific << setprecision(12) << displaced_path_length << endl;
	//	cout << "Displaced - Original = " << scientific << setprecision(12)
	//		 << displaced_path_length - unmodulated_path_length << endl;
    
    
}



void Photon::displacePhotonFromPressure(void)
{    
    
    // Photon does not get displaced on boundaries of medium.
    if (hit_x_bound || hit_y_bound || hit_z_bound) return;
    
#ifdef DEBUG
    int Nx = m_medium->kwave.dmap->getnumVoxelsXaxis();
    int Ny = m_medium->kwave.dmap->getnumVoxelsYaxis();
    int Nz = m_medium->kwave.dmap->getnumVoxelsZaxis();
    
	assert((_x < Nx && _x >= 0) &&
           (_y < Ny && _y >= 0) &&
           (_z < Nz && _z >= 0) ||
           assert_msg("_x=" << _x << " _y=" << _y << " _z=" << _z << "\n"
        		      << currLocation->location.x << " "
        		      << currLocation->location.y << " "
        		      << currLocation->location.z));
#endif    
    
    

	// Transform the location of the photon in the medium to discrete locations in the grid.
	//
	double dx = m_medium->kwave.dmap->getDx();
	double Nx = m_medium->kwave.dmap->getNumVoxelsXaxis();

	double dy = m_medium->kwave.dmap->getDy();
	double Ny = m_medium->kwave.dmap->getNumVoxelsYaxis();

	double dz = m_medium->kwave.dmap->getDz();
	double Nz = m_medium->kwave.dmap->getNumVoxelsZaxis();

	int _x = currLocation->location.x/dx - (currLocation->location.x/dx)/Nx;
	int _y = currLocation->location.y/dy - (currLocation->location.y/dy)/Ny;
	int _z = currLocation->location.z/dz - (currLocation->location.z/dz)/Nz;

    


    // Subtle case where index into grid is negative because of rounding errors above.
#ifdef DEBUG
    if (_z < 0 || _x < 0 || _y < 0)
    {
        // FIXME:
        // - This should be removed when detection of line-plane intersection
        //   is implemented.  For now, it is a necessary evil.
    	cout << "Error in array index calculation: Photon::displacePhotonFromPressure()\n";
        this->status = DEAD;
        return;
    }
#endif
    
    
    // Calculate changes in the optical path length due to the refractive index gradient
    // produced due to pressure variations.
    // XXX: Does this need to happen before displacement of the photon due to pressure,
    //      or should this happen after?
    //      NOTE: I think before, because the arc of the path (due to the refractive gradient)
    //            would place the photon at a new location.  However, currently only the variation
    //            in the optical path length is calculated, not the change in the position of the photon.
    //alterPathLengthFromRefractiveChanges();
    
    
    
    // Calculate the displacement due to pressure.
    //
    // Index into the displacement grids to retrieve pre-calculated value of how much to
    // displace the photon from it's current location based upon the pressure in the voxel.
    currLocation->location.x += m_medium->kwave.dmap->getDisplacementFromGridZ(_x, _y, _z);
    currLocation->location.y += m_medium->kwave.dmap->getDisplacementFromGridY(_x, _y, _z);
    currLocation->location.z += m_medium->kwave.dmap->getDisplacementFromGridX(_x, _y, _z);
    
    // Update the optical path length of the photon through the medium by
    // calculating the distance between the two points and multiplying by the refractive index.
    displaced_optical_path_length += VectorMath::Distance(prevLocation, currLocation) * currLayer->getRefractiveIndex();
}



// Alter the optical path length of the photon through the medium as it encounters
// refractive index changes along it's step.
// Basic idea is that the line segment from the photon step encounters changes
// in the mediums refractive indeces.  From Fermat's theorem we know that the 
// arc (S) that the line segment will make is an extremum based on the values
// of refractive index it encounters on the physical path length (L).  The change
// in the optical length is given by L*refractive_index.
// By using the parametric equation of a line from the previous and current locations
// of the photon it is possible to move along this line, index into the K-Wave supplied
// grid of pre-calculated refractive index changes, and alter the pathlength accordingly.
void Photon::alterPathLengthFromRefractiveChanges(void)
{
    // Photon does not get its optical path length changed on boundaries of medium.
    if (hit_x_bound || hit_y_bound || hit_z_bound) return;
    
    // Current value of the change in refractive index (i.e. n(x,y,z) coordinates).
    double n_curr = 0.0f;  
    // Previous value of the change in refractive index.
    double n_prev = 0.0f; 
    // The value of the background index of refraction (i.e. no change due to pressure).
    double n_background = 0.0f;
    
    // Scales the position along the line segment from previous to current location.
    // i.e. P(t) = prevLocation + t*currLocation
    double t_curr = 0.0f;   
    double t_prev = 0.0f;
    
    // The adiabatic piezo-optical coefficient of the material.
    double eta = 0.3211;
    
    // The wave number of the acoustic wave.
    double k = 2*PI/m_medium->kwave.transducerFreq;
    
    // Create the new vector that is the line segment from the previous position of the photon (p0) to 
    // it's current position (p1).
    boost::shared_ptr<Vector3d> lineSegment = (*currLocation) - (*prevLocation);
    
    // These track the distance on the line segment the photon travels that contains it's respective
    // index of refraction.  That is, by tracking where a given point in space starts and ends with
    // a given value of refractive index, we know which portion (i.e. the distance from
    // 'previousPointOnSegment' to 'pointOnSegment') should be multiplied by it's corresponding index
    // of refraction, thus yeilding an optical path length dependent on the gradient changes of the
    // refraction of index due to variations in pressure.
    boost::shared_ptr<Vector3d> pointOnSegment;
    boost::shared_ptr<Vector3d> previousPointOnSegment = prevLocation;  // Set to last end point from hopping in medium.
    
    
    // Prime the check below by setting the previous and current value of the refractive index changes to the same values.
    //
    // The grid indices that are calculated form the photon's position.
    // Transform the location of the photon in the medium to discrete locations in the grid.
    //
    double dx = m_medium->kwave.pmap->getDx();
    double Nx = m_medium->kwave.pmap->getNumVoxelsXaxis();
    
    double dy = m_medium->kwave.pmap->getDy();
    double Ny = m_medium->kwave.pmap->getNumVoxelsYaxis();
    
    double dz = m_medium->kwave.pmap->getDz();
    double Nz = m_medium->kwave.pmap->getNumVoxelsZaxis();
    
    
    
    int _x = previousPointOnSegment->location.x/dx - (previousPointOnSegment->location.x/dx)/Nx;
    int _y = previousPointOnSegment->location.y/dy - (previousPointOnSegment->location.y/dy)/Ny;
    int _z = previousPointOnSegment->location.z/dz - (previousPointOnSegment->location.z/dz)/Nz;   
    
   
    // Prime the refractive index values so they can be compared as 'n_curr' is updated based on position.
    n_background = currLayer->getRefractiveIndex(); 
    n_curr = n_prev = n_background + n_background*eta*k*m_medium->kwave.pmap->getPressureFromGrid(_x, _y, _z);
    

    
    
    // Because we must make the integration of the refractive index changes discrete, we
    // have a step size (ds) along the arc (S).
    // More specifically, value of (steps) dictates how many portions (i.e. the resolution) we take along the arc (S).
    int steps = 8;
    double ds = 1/(double)steps;
    for (int i = 0; i < steps; i++, t_curr += ds)
    {
        // Update the (t_curr) value to move along the line segment.
        // 'prevLocation' is where the last hop moved the photon. 'lineSegment'
        // is the line formed form 'prevLocation' to the new location of the photon 'currLocation'.
        // 'pointOnSegment' is a point on that line ('lineSegment') based on the value of 't'.
        //t_curr += ds;
        pointOnSegment = (*prevLocation) + (*((*lineSegment)*t_curr));
        
        // Transform the location of the photon in the medium to discrete locations in the grid.
        //
        _x = pointOnSegment->location.x/dx - (pointOnSegment->location.x/dx)/Nx;
        _y = pointOnSegment->location.y/dy - (pointOnSegment->location.y/dy)/Ny;
        _z = pointOnSegment->location.z/dz - (pointOnSegment->location.z/dz)/Nz;  
        
        
#ifdef DEBUG
        int Nx = m_medium->kwave.pmap->getnumVoxelsXaxis();
        int Ny = m_medium->kwave.pmap->getnumVoxelsYaxis();
        int Nz = m_medium->kwave.pmap->getnumVoxelsZaxis();
        
        assert((_x < Nx && _x >= 0) &&
               (_y < Ny && _y >= 0) &&
               (_z < Nz && _z >= 0) ||
               assert_msg("_x=" << _x << " _y=" << _y << " _z=" << _z << "\n"
                          << currLocation->location.x << " "
                          << currLocation->location.y << " "
                          << currLocation->location.z));
#endif
        
        // Update the value of the refractive index based on coordinates that retrieve the 
        // current discrete location of the pressure from the pressure map grid.
        //
        // background index of refraction.  Might be different within some occlusion, so we pass
        // in coordinates to get the value at a specific location.  If simulating a homogenous medium
        // a constant can be placed here, which should speed up the simulation since no searching will occur.
        // n_background = currLayer->getRefractiveIndex(pointOnSegment);
        n_curr = n_background + n_background*eta*k*m_medium->kwave.pmap->getPressureFromGrid(_x, _y, _z);
        
        // If there was a change in refractive index value, due to a change in pressure found in the pressure
        // grid, we accumulate this portion of the step, which is multiplied by this portion's refractive index
        // value, and update the optical path length accordingly.
        //
        if (n_curr != n_prev)
        {
            
            // Displace the photon on its arced path due to the refractive index changes along the step.
            //displacePhotonFromRefractiveGradient(n_prev, n_curr);
            
            // The point on the line segment where the this refractive index change started.
            previousPointOnSegment = (*prevLocation) + (*((*lineSegment)*t_prev));
            
            // The current point on the segment has moved to another voxel in the K-Wave simulation, therefore
            // we need to step back a portion of 'ds' that was moved so we have the length over a given 
            // refractive index change, not just beyond.  Therefore we subtract off the portion that was added
            // to 'ds' when updating the current point.
            pointOnSegment = (*prevLocation) + (*((*lineSegment)*(t_curr)));
            
            refractiveIndex_optical_path_length += VectorMath::Distance(pointOnSegment, previousPointOnSegment) * n_curr;
            
            // Update the values for the next iterations.
            n_prev = n_curr;
            t_prev = t_curr;
            
        }

        
    } // end for() loop
    
                                
    // If the last portion of the line segment (pointOnSegment - previousPointOnSegment) didn't see an index of 
    // refraction change, its optical path length won't be updated.  Also, if the entire step fit into a voxel, this takes 
    // that into account as well.
     refractiveIndex_optical_path_length += VectorMath::Distance(pointOnSegment, previousPointOnSegment) * n_curr;
    
}



// FIXME:
// - This currently does NOT work correctly.
// - Need to calculate the change in direction as the photon moves into a refractive
//   index change.
// - Need to figure out the orientation of the plane that the photon passes through
//   which causes refraction.  Is it normal to the direction of photon propagation?
//   Is it normal to the direction of ultrasound.
void Photon::displacePhotonFromRefractiveGradient(const double n1, const double n2)
{

    boost::shared_ptr<Vector3d> n(new Vector3d(0.0f, 0.0f, 1.0f));
    boost::shared_ptr<Vector3d> I(new Vector3d);
    I->location.x = currLocation->location.x;
    I->location.y = currLocation->location.y;
    I->location.z = currLocation->location.z;
    

    VectorMath::Normalize(I);
    double cosT = -VectorMath::dotProduct(n, I);
    
    double cosT2 = sqrt(1 - pow((n1/n2),2)*(1-pow(cosT,2)));
    
    boost::shared_ptr<Vector3d> R1 = (*I)*(n1/n2);
    boost::shared_ptr<Vector3d> R2 = (*I)*(n1/n2*cosT - cosT2);
    boost::shared_ptr<Vector3d> R = (*R1) + (*R2);
    R = (*((*I)*(n1/n2))) + (*((*n)*((n1/n2)*cosT-cosT2)));
    
    boost::shared_ptr<Vector3d> reflected = (*I) + (*((*I)*(cosT*2)));
    
    
    /*
    double dirX = currLocation->getDirX();
    double dirY = currLocation->getDirY();
    double dirZ = currLocation->getDirZ();
    
    double incidentX = acos(abs(dirX));
    double incidentY = acos(abs(dirY));
    double incidentZ = acos(abs(dirZ));
    
    double t_x = asin(n1/n2*sin(incidentX));
    double t_y = asin(n1/n2*sin(incidentY));
    double t_z = asin(n1/n2*sin(incidentZ));

    t_x = cos(t_x);
    t_y = cos(t_y);
    t_z = cos(t_z);
    */
    
}



// FIXME:
// - The photon, if moving back up to the 'Air' layer, should be able to transmit out of the medium.
//   Currently this does not happen (maybe because of reflection on the layer boundary?).
// - If the photon hits the layer at the bottom of the medium, it gets stuck there until
//   it pobabilistically is transmitted through the medium.  It should hit the layer, check
//   the medium, and determine if it should transmit or reflect.  This bug only arises occasionally.
void Photon::transmit(const char *type)
{
    

    // 'tempLayer' is used 
    Layer *tempLayer = NULL;
    
    if (strcmp("layer", type) == 0)
    {
#ifdef DEBUG
        cout << "Transmitting through layer\n";
        cout << currLocation;
#endif
        
        
        // If we transmit through the layer to another we must update
        // the current layer pointer of the photon so it will correctly 
        // calculate the next step size.
        if (currLocation->getDirZ() > 0)
            tempLayer = m_medium->getLayerBelowCurrent(currLocation->location.z);
        else
            tempLayer = m_medium->getLayerAboveCurrent(currLayer);
        
        
        

        // If 'tempLayer' is NULL we are at the edge of the medium since
        // the layer above or below does not exist.  Therefore we decide if the
        // photon should reflect or transmit from the medium boundary.
        if (tempLayer == NULL) 
        {
            // Since 'tempLayer' is NULL we have tried to retrieve a layer outside of the
            // medium's bounds.  Therefore 'currLayer' still is valid (i.e. not NULL) and we check
            // if the photon should reflect or transmit through the medium.
            if (hitMediumBoundary())
            {
                hop(); // Move the photon to the medium boundary.
                transmitOrReflect("medium");
                return;
            }
        }
        else // Transmitted to new layer.
        {
        	// Store index of refraction from layer the photon is coming from.
        	double ni = currLayer->getRefractiveIndex();

        	// Get the new index of refraction from the layer the photon is moving to.
        	double nt = tempLayer->getRefractiveIndex();

        	// Since the photon is being transmitted, update to the pointer to the new layer that
        	// was found above.
            currLayer = tempLayer;

            // NOTE: Assumes layers are normal to initial direction, which is dirZ.
            // Set the direction cosines.
            currLocation->setDirZ(cos(this->transmission_angle));
            currLocation->setDirY(currLocation->getDirY()*ni/nt);
            currLocation->setDirX(currLocation->getDirX()*ni/nt);
        }
        


    }
    else if (strcmp("medium", type) == 0)
    {
#ifdef DEBUG
        cout << "Transmitting through medium boundary\n";
        cout << currLocation;
#endif
        // If we transmit through the medium, meaning the photon exits the medium,
        // we see if the exit location passed through the detector.  If so, the exit
        // data is written out to file.
        if (checkDetector())
        {
            // If we hit the detector when transmitting the photon, then we write the exit
            // data to file.
//            Logger::getInstance()->writeWeightAngleLengthCoords(this->weight,
//                                                                this->transmission_angle,
//                                                                this->displaced_optical_path_length,
//                                                                this->currLocation);
        	// Update the direction cosines upon leaving the medium so calculations can be mode
        	// to see if this photon makes it's way to the detector (i.e. CCD camera).
        	// NOTE:
        	// - We assume outside of the medium is nothing but air with an index of refraction 1.0
        	//   so no division for x & y direction cosines because nt = 1.0
        	// - It is also assumed that the photon is transmitted through the x-y plane.
        	double ni = currLayer->getRefractiveIndex();
        	currLocation->setDirZ(cos(this->transmission_angle));
        	currLocation->setDirY(currLocation->getDirY()*ni);
        	currLocation->setDirX(currLocation->getDirX()*ni);

        	// Write exit data from logger.
        	Logger::getInstance()->writeWeightAngleLengthCoords(*this);
            
            // Write out the seeds that caused this photon to hop, drop and spin its way out the
            // exit-aperture.
            //
            Logger::getInstance()->writeRNGSeeds(seeds.s1, seeds.s2, seeds.s3, seeds.s4);
            cnt++;
        }
        
        // The photon has left the medium, so kill it.
        this->status = DEAD;
    }
    else
    {
        cout << "Error:  Wrong paramater passed to Photon::transmit()\n";
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
#ifdef DEBUG
            cout << "Internally reflecting on layer boundary\n";
#endif
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
#ifdef DEBUG
            cout << "Transmitting through layer\n";
#endif
			transmit("layer");
		}
	}
	// Case when interaction is with a medium boundary.
	else if (strcmp("medium", type) == 0)
	{
		// Stochastically determine if the photon should be transmitted or reflected.
		if (getMediumReflectance() > getRandNum())
		{
#ifdef DEBUG
            cout << "Reflecting photon on medium boundary\n";
#endif
			// Depending on which medium boundary was hit, we reflect on that axis,
            // change the direction of the direction cosign, and reset the boolean flag.
            if (hit_x_bound)
            {
				internallyReflectX();
            }
			else if (hit_y_bound)
            {
				internallyReflectY();
            }
			else if (hit_z_bound)
            {
				internallyReflectZ();
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
    
    
	//Layer *currLayer = m_medium->getLayerFromDepth(currLocation->location.z);
	double refract_index_n1 = currLayer->getRefractiveIndex();	// Current layer's refractive index.
	double refract_index_n2 = 1.0;	// Outside of the medium is only air.
    
    
	double axis_direction = 0.0;
	if (hit_x_bound)
		axis_direction = currLocation->getDirX();
	else if (hit_y_bound)
		axis_direction = currLocation->getDirY();
	else
		axis_direction = currLocation->getDirZ();
    
	// Calculate the incident angle based on the axis in which the photon hit the medium.
	double incident_angle = acos(abs(axis_direction));
    
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
        // Calculate the transmission angle through the medium boundary.
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
	Layer *nextLayer = NULL;
    
	double incident_angle = acos(abs(currLocation->getDirZ()));
	refract_index_n1 = currLayer->getRefractiveIndex();
    
	// If the photon is moving towards a deeper layer.
	if (currLocation->getDirZ() > 0)
	{
		nextLayer = m_medium->getLayerBelowCurrent(currLocation->location.z);
	}
	// If the photon is moving towards a more shallow layer.
	else if (currLocation->getDirZ() < 0)
	{
		nextLayer = m_medium->getLayerAboveCurrent(currLayer);
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
        // For specular reflection we always remove some portion of the weight and transmit the photon
        // to the next layer.
        transmission_angle = asin(refract_index_n1/refract_index_n2 * sin(incident_angle));
        specularReflectance(refract_index_n1, refract_index_n2);
        
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




// FIXME:
// - Switch to planar bounds detection.  Below is so convoluted it's scary.
//   There is also an edge case of injection the photon, at its first step,
//   has directon cosignZ == 1.  If diffuse reflection occurs and it leaves
//   the medium boundary on its first run, this case is not accounted for.
//   THIS IS BAD WHEN INDEXING INTO A DISPLACEMENT ARRAY AND HAVING A NEGATIVE
//   VALUE CAUSES THINGS TO BLOW UP AND SPEND AN ENTIRE DAY FINDING A SUBTLE BUG.
//   THIS IS A REMINDER TO SWITCH TO PLANE INTERSECTION DETECTION....
// Note: We take the absolute value in the case where the direction 
//       cosine is negative, since the distance to boundary would be
//       negative otherwise, which is untrue.  This arises due to assuming
//       the lower axis bound in each dimension (x, y, z) begins at zero.
//       This could also be achieved by simply subtracting the current location
//       from zero (e.g. 0-y/diry), which would change the sign as well.
bool Photon::hitMediumBoundary(void)
{
	double distance_to_boundary = 0.0;
    double distance_to_boundary_X = 0.0;
    double distance_to_boundary_Y = 0.0;
    double distance_to_boundary_Z = 0.0;
    
	//Layer *layer = m_medium->getLayerFromDepth(currLocation->location.z);
	double mu_t = currLayer->getTotalAttenuationCoeff(currLocation);
	double x_bound = m_medium->getXbound();
	double y_bound = m_medium->getYbound();
	double z_bound = m_medium->getZbound();
    
	// Check if the photon has been given a step size past the outer edges of the medium.
	// If so we calculate the distance to the boundary.
	
    // The step that puts the photon in a new location based on its current location.
    double x_step = step*currLocation->getDirX() + currLocation->location.x;
    double y_step = step*currLocation->getDirY() + currLocation->location.y;
    double z_step = step*currLocation->getDirZ() + currLocation->location.z;    
    
    
    
    
    if (x_step >= x_bound || x_step <= 0.0f)
	{
		hit_x_bound = true;
		if (currLocation->getDirX() > 0.0f) // Moving towards positive x_bound
			distance_to_boundary_X = (x_bound - currLocation->location.x) / currLocation->getDirX();
		else
			distance_to_boundary_X = abs(currLocation->location.x / currLocation->getDirX());
	}
    
	if (y_step >= y_bound || y_step <= 0.0f)
	{
		hit_y_bound = true;
		if (currLocation->getDirY() > 0.0f) // Moving towards positive y_bound
			distance_to_boundary_Y = (y_bound - currLocation->location.y) / currLocation->getDirY();
		else
			distance_to_boundary_Y = abs(currLocation->location.y / currLocation->getDirY());
	}
    
	if (z_step >= z_bound || z_step <= 0.0f)
	{
		hit_z_bound = true;
		if (currLocation->getDirZ() > 0.0f) // Moving towards positive z_bound
			distance_to_boundary_Z = (z_bound - currLocation->location.z) / currLocation->getDirZ();
		else
			distance_to_boundary_Z = abs(currLocation->location.z / currLocation->getDirZ());
	}
    
    // If none were hit, we can return false and forego any further checking.
    if (!hit_x_bound && !hit_y_bound && !hit_z_bound)
        return false;
    
    
    // Need to take care of the rare case that photons can be at a corner edge and step
    // past the boundary of the medium in 2 or more dimensions.  If it is found that
    // the photon has passed through the bounds of the grid in more than 2 bounds we always
    // want to move it the smallest distance (i.e. to the closest bounds), therefore we should
    // never cross 2 bounds (or more!) simultaneously.
    if (hit_x_bound && hit_y_bound)
    {
        if (distance_to_boundary_X < distance_to_boundary_Y)
        {
            distance_to_boundary = distance_to_boundary_X;
            hit_x_bound = true; 
            hit_y_bound = false;
            hit_z_bound = false;
        }
        else
        {
            distance_to_boundary = distance_to_boundary_Y;
            hit_x_bound = false; 
            hit_y_bound = true;
            hit_z_bound = false;
        }
    }
    else if (hit_x_bound && hit_z_bound)
    {
        if (distance_to_boundary_X < distance_to_boundary_Z)
        {
            distance_to_boundary = distance_to_boundary_X;
            hit_x_bound = true;
            hit_y_bound = false; 
            hit_z_bound = false;
        }
        else
        {
            distance_to_boundary = distance_to_boundary_Z;
            hit_x_bound = false;
            hit_y_bound = false;
            hit_z_bound = true;
        } 
    }
    else if (hit_y_bound && hit_z_bound)
    {
        if (distance_to_boundary_Y < distance_to_boundary_Z)
        {
            distance_to_boundary = distance_to_boundary_Y;
            hit_x_bound = false;
            hit_y_bound = true;
            hit_z_bound = false;
        }
        else
        {
            distance_to_boundary = distance_to_boundary_Z;
            hit_x_bound = false;
            hit_y_bound = false;
            hit_z_bound = true;
        }
        
    }
    else if (hit_x_bound && hit_y_bound && hit_z_bound)
    {
    	cout << "ERROR: hit 3 boundaries\n";
    	cout.flush();
    	distance_to_boundary = distance_to_boundary_X < distance_to_boundary_Y ?
    							distance_to_boundary_X : distance_to_boundary_Y;
    	distance_to_boundary = distance_to_boundary < distance_to_boundary_Z ?
    	    							distance_to_boundary : distance_to_boundary_Z;
    	hit_x_bound = true;
    	hit_y_bound = hit_z_bound = false;
    }
    else
    {
        // If we make it here, we have hit only one boundary, which means the other 2 have values
        // of zero for their distances.  Therefore, without checking we can simply add them all together
        // to get the single distance to boundary.
        distance_to_boundary = distance_to_boundary_X + distance_to_boundary_Y + distance_to_boundary_Z;
    }
    
    
    
	// If the step size of the photon is larger than the distance
	// to the boundary and we are moving in some direction along
	// the axis we calculate the left over step size and then step
	// the photon to the boundary with the next call to 'hop()'.
	//if (step > distance_to_boundary && distance_to_boundary != 0)
	if (step > distance_to_boundary)
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
	//Layer *layer = m_medium->getLayerFromDepth(currLocation->location.z);
	double mu_t = currLayer->getTotalAttenuationCoeff(currLocation);
    
    
    
	// If the direction the photon is traveling is towards the deeper boundary
	// of the layer, we execute the first clause.  Otherwise we are moving to
	// the more superficial boundary of the layer.
	if (currLocation->getDirZ() > 0.0)
	{
		distance_to_boundary = (currLayer->getDepthEnd() - currLocation->location.z) / currLocation->getDirZ();
	}
	else if (currLocation->getDirZ() < 0.0)
	{
		distance_to_boundary = (currLayer->getDepthStart() - currLocation->location.z) / currLocation->getDirZ();
	}
    
    
	// If the step size of the photon is larger than the distance
	// to the boundary and we are moving in some direction along
	// the z-axis (i.e. not parallel to the layer boundary) we calculate
	// the left over step size and then step the photon to the boundary.
	if (currLocation->getDirZ() != 0.0 && step > distance_to_boundary)
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


void Photon::specularReflectance(double n1, double n2)
{
	// update the weight after specular reflectance has occurred.
	weight = weight - (pow((n1 - n2), 2) / pow((n1 + n2), 2)) * weight;
}


// Update the direction cosine when internal reflection occurs on z-axis.
void Photon::internallyReflectZ(void) 
{
    currLocation->setDirZ(-1*currLocation->getDirZ());
    
    // Reset the flag.
    hit_z_bound = false;
}

// Update the direction cosine when internal reflection occurs on y-axis.
void Photon::internallyReflectY(void) 
{
    currLocation->setDirY(-1*currLocation->getDirY());
    
    // Reset the flag.
    hit_y_bound = false;
}


// Update the direction cosine when internal reflection occurs on z-axis.
void Photon::internallyReflectX(void) 
{
    currLocation->setDirX(-1*currLocation->getDirX());
    
    // Reset the flag.
    hit_x_bound = false;
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
