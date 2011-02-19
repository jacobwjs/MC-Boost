
#include "debug.h"
#include "photon.h"



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


// Here is where execution begins for the thread and we execute the  main steps,
// 1) Hop - move the photon
// 2) Drop - drop weight due to absorption
// 3) Spin - update trajectory accordingly
// 4) Roulette - test to see if photon should live or die.


void Photon::run()
{
	/*
	// Execute this photon (i.e. thread) the given number of iterations, thus
	// allowing the work to be split up.
	while (cnt < iterations) {
		
		while (isAlive()) {
#ifdef DEBUG
			cout << "Running thread( " << Thread::getThreadID() << ")...";
			cout << "..hop(), drop(), spin(), roulette().  Propogated " << cnt << " times.\n";
#endif
			
			// Move the photon.
			hop();
			
			// Drop weight of the photon due to an interaction with the medium.
			drop();
			
			// Calculate the new coordinates of photon propogation.
			spin();
			
			// Check if the photon has reached the threshold weight, and if so
			// perform roulette to see if it should continue (ALIVE) or end (DEAD).
			performRoulette();
		}
		
		cnt++;
		reset();
	}
	
	// The thread has done its portion of the work, time to exit.
	//Thread::exit();

*/	 
}


void Photon::plotPath()
{
	plot([


void Photon::reset()
{
#ifdef DEBUG
	cout << "Reseting Photon...\n";
#endif
	
	// Photon just created, so it is alive.
	status = ALIVE;
	
	// Set back to default values.
	weight = 1;
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
	
	// FIXME: Should check at which depth the photon resides in the medium
	//        (e.g. double mu_a = medium->getlayerAbsorption(z)), but hard coded values
	//        are used for now, which assumes homogeneous single layer throughout. 
	double mu_a = 1.0; // cm^-1
	double mu_s = 0.0; // cm^-1
	step = -log(rnd)/(mu_a	+ mu_s);
	
	
	
	// Update position of the photon.
	x += step*dirx;
	y += step*diry;
	z += step*dirz;
	
	

	// Store the x, y, and z location for this interaction.  Used with plotting.
	location_x[num_steps] = x;
	location_y[num_steps] = y;
	location_z[num_steps] = z;
	// Check to see if we have a valid index into the array to store the coords.
	num_steps += 1;
	if (num_steps > 10000) {
		cout << "Error, indexing past array bounds\n";
	}
	
}

// Return absorbed energy from photon weight at this location.
double Photon::drop()
{
#ifdef DEBUG
	cout << "Dropping...\n";
#endif	
	
	// FIXME:  Should return the albedo for the layer in which the photon is 
	//			currently being propogated through. (e.g. double val = medium->getLayerAlbedo(z))
	//		   For now assume single homogenous layer.
	double mu_s = 0.0;
	double mu_a = 1.0;
	double albedo = mu_s / (mu_a + mu_s);
	double absorbed = weight * (1 - albedo);
	
	// Remove the portion of energy lost due to absorption at this location.
	weight -= absorbed;
	
	// Return the energy that was absorbed.
	return absorbed;
}


void Photon::spin()
{
#ifdef DEBUG
	cout << "Spinning...\n";
#endif	
	
	// FIXME:  Should index into medium to find 'g' for given layer at this depth.
	//		   e.g. double g = medium->getAnisotropy(z);
	//		   Assuming static value now (ie single homogeneous layer).
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

