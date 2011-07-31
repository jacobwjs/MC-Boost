
#include "debug.h"
#include "photon.h"

#undef DEBUG


double const Pi=4*atan(1);


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
    
    // Default value of tagged to false.
    tagged = false;
    
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
			bool hitLayer = checkLayerBoundary();
            bool hitMedium = checkMediumBoundary();

   
            
			if (!hitLayer && !hitMedium)
			{
                // sanity check.
                assert(this->status == ALIVE);

                
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
	original_path_length = 0;
	displaced_path_length = 0;


	// Randomly set photon trajectory to yield isotropic or anisotropic source.
	initTrajectory();

    // Reset the current layer from the injection coordinates of the photon.
    currLayer = m_medium->getLayerFromDepth(currLocation->location.z);

}


// Update the step size for the photon.
void Photon::setStepSize()
{
	// Update the current values of the absorption and scattering coefficients
	// based on the depth in the medium (i.e. which layer the photon is in).
    double mu_a = currLayer->getAbsorpCoeff(currLocation);
    double mu_s = currLayer->getScatterCoeff();
    
    
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


double Photon::getPathLength(double x_dist, double y_dist, double z_dist)
{
	return sqrt(pow(x_dist, 2) + pow(y_dist, 2) + pow(z_dist, 2));
}


void Photon::captureLocationCoords(void)
{
	// Add the coordinates to the STL vector for the displaced scattering locations.
    cout << "captureLocationCoords() stub\n";
}


void Photon::captureExitCoordsAndLength(void)
{
	// Add the coordinates to the STL vector for the displaced scattering locations
	// and the displaced length.
//	photon_exit_data.push_back(x_disp);
//	photon_exit_data.push_back(y_disp);
//	photon_exit_data.push_back(displaced_path_length);
    cout << "captureExitCoordsAndLength() stub\n";
}


void Photon::captureExitCoordsLengthWeight(void)
{
	// Add the coordinates to the STL vector for the displaced scattering locations
	// and the displaced length.
//	photon_exit_data.push_back(x_disp);
//	photon_exit_data.push_back(y_disp);
//	photon_exit_data.push_back(displaced_path_length);
//	photon_exit_data.push_back(weight);
    cout << "captureExitCoordsLengthWeight() stub\n";
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
    
    // We have to hop() once to move the currentLocation to what was calculated
    // from setStepSize(), since we made it here after verifying that the new
    // step size would have put the photon outside of the medium, but it was
    // never made to move to the new position.  Therefore prevLocation == currLocation
    // until hop is made.
    //hop();
    int cnt =  m_medium->photonHitDetectorPlane(this->currLocation);
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

	cnt++;

/*
	// Record the location of the photon during this interaction.
	//captureLocationCoords();

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

	// FIXME: Does this still need to be here???
	if (z_disp > 0 && z_disp < 10 && x_disp > 0 && x_disp < 10 && y_disp > 0 && y_disp < 10) {
		// Move the photon to the new position based on the displacement from
		// the ultrasound wave.
		this->displacePhotonFromPressure();
	}
	else {
		this->status = DEAD;
	}

	// Calculate the path length of the photon WITH displacement.
	displaced_path_length += getPathLength((x_disp - temp_x), (y_disp - temp_y), (z_disp - temp_z));
*/
    
	//setStepSize();
    
    // Save the coordinates of the photon before it is moved through
    // the medium.  This allows us to check if the photon passed through
    // an absorber.
    prevLocation->location.x = currLocation->location.x;
    prevLocation->location.y = currLocation->location.y;
    prevLocation->location.z = currLocation->location.z;
    
    // Move the photon to the new location in the medium.
    currLocation->location.x += step*currLocation->getDirX();
	currLocation->location.y += step*currLocation->getDirY();
	currLocation->location.z += step*currLocation->getDirZ();
}


// Return absorbed energy from photon weight at this location.
void Photon::drop()
{
#ifdef DEBUG
	cout << "Dropping...\n";
#endif	
    
    double mu_a = 0.0f;
    double mu_s = 0.0f;
    double albedo = 0.0f;
    double absorbed = 0.0f;
    
    
    Absorber *absorber;
    absorber = currLayer->getAbsorber(currLocation);
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

	//	cout << "Non-displaced path length: " << scientific << setprecision(12) <<  original_path_length << endl;
	//	cout << "Displaced path length: " << scientific << setprecision(12) << displaced_path_length << endl;
	//	cout << "Displaced - Original = " << scientific << setprecision(12)
	//		 << displaced_path_length - original_path_length << endl;


	//m_medium->writePhotonCoords(coords);
    cout << "writeCoordsToFile() stub\n";
}


// Write exit locations and path lenghts of photons to file.
void Photon::writeExitLocationsAndLength(void)
{
	//m_medium->writeExitCoordsAndLength(photon_exit_data);
    cout << "Photon::writeExitLocationsAndLength(void) stub\n";
}

void Photon::writeExitLocationsLengthWeight(void)
{
	//m_medium->writeExitCoordsLengthWeight(photon_exit_data);
    cout << "Photon::writeExitLocationsLengthWeight(void) stub\n";
}

// FIXME: CURRENTLY ONLY DISPLACING IN ONE DIRECTION.  SHOULD USE A TENSOR.
void Photon::displacePhotonFromPressure(void)
{
	// Get the local pressure from the grid based on the coordinate of the photon.
	double pressure = m_medium->getPressureFromCartCoords(x_disp, y_disp, z_disp);


	// Impedance of the tissue.
	double impedance = currLayer->getImpedance();

	// Frequency of the ultrasound transducer.
	static double transducer_freq = m_medium->pmap->getTransducerFreq();

	// Displace the photon in the x-axis based on the pressure.
	// Note: Since we are firing the photons into the medium normal
	//		 to the propagation of the ultrasound wave, we displace on the
	//		 x-axis.  Also, multiplication by 100 puts the displacement from meters to cm.
	double displacement = 100*(pressure)/(2*Pi*transducer_freq*impedance); // [cm]
	x_disp += displacement;

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
        
        // Set the direction cosine.
        currLocation->setDirZ(cos(transmission_angle));
        
        
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
            }
        }
        else
        {
            currLayer = tempLayer;
        }
                    
    }
    else if (strcmp("medium", type) == 0)
    {
#ifdef DEBUG
        cout << "Transmitting through medium boundary\n";
        cout << currLocation;
#endif
        // If we transmit through the medium, meaning the photon exits the medium
        // we see if the exit location passed through the detector.  If so, the exit
        // location and exit angle are written out to file, but only when this photon
        // has been tagged (i.e. interacted with an absorber).
        if (checkDetector() && tagged)
        {
            // If we hit the detector when transmitting the photon, then we write the exit
            // data to file.
            Logger::getInstance()->writeExitData(this->currLocation, 
                                                 this->weight);
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
	//Layer *currLayer = m_medium->getLayerFromDepth(currLocation->location.z);
	Layer *nextLayer;
    
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
        // For specular refection we always remove some portion of the weight and transmit the photon
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





// Note: We take the absolute value in the case where the direction 
//       cosine is negative, since the distance to boundary would be
//       negative otherwise, which is untrue.  This arises due to assuming
//       the lower axis bound in each dimension (x, y, z) begins at zero.
//       This could also be achieved by simply subtracting the current location
//       from zero (e.g. 0-y/diry), which would change the sign as well.
bool Photon::hitMediumBoundary(void)
{
	double distance_to_boundary = 0.0;
	//Layer *layer = m_medium->getLayerFromDepth(currLocation->location.z);
	double mu_t = currLayer->getTotalAttenuationCoeff();
	double x_bound = m_medium->getXbound();
	double y_bound = m_medium->getYbound();
	double z_bound = m_medium->getZbound();
    
	// Check if the photon has been given a step size past the outer edges of the medium.
	// If so we calculate the distance to the boundary.
	
    
    if (step*currLocation->getDirX() + currLocation->location.x >= x_bound || 
        step*currLocation->getDirX() + currLocation->location.x <= 0)
	{
		hit_x_bound = true;
		if (currLocation->getDirX() > 0) // Moving towards positive x_bound
			distance_to_boundary = (x_bound - currLocation->location.x) / currLocation->getDirX();
		else
			distance_to_boundary = abs(currLocation->location.x / currLocation->getDirX());
	}
	else if	(step*currLocation->getDirY() + currLocation->location.y >= y_bound || 
             step*currLocation->getDirY() + currLocation->location.y <= 0)
	{
		hit_y_bound = true;
		if (currLocation->getDirY() > 0) // Moving towards positive y_bound
			distance_to_boundary = (y_bound - currLocation->location.y) / currLocation->getDirY();
		else
			distance_to_boundary = abs(currLocation->location.y / currLocation->getDirY());
	}
	else if (step*currLocation->getDirZ() + currLocation->location.z >= z_bound || 
             step*currLocation->getDirZ() + currLocation->location.z <= 0)
	{
		hit_z_bound = true;
		if (currLocation->getDirZ() > 0) // Moving towards positive z_bound
			distance_to_boundary = (z_bound - currLocation->location.z) / currLocation->getDirZ();
		else
			distance_to_boundary = abs(currLocation->location.z / currLocation->getDirZ());
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
	double mu_t = currLayer->getTotalAttenuationCoeff();
    
    
    
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
