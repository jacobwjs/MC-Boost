/*
 * Copyright BMPI 2011
 * J.W. Staley - MIRA, 
 *               Biomedical Photonics Imaging Group (BMPI), 
 *				 University of Twente
 *
 */

//#define DEBUG 1



#include "stdio.h"
#include "photon.h"
#include "medium.h"
#include "layer.h"
#include "pressureMap.h"
#include <cmath>
#include <time.h>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <boost/thread/thread.hpp> 
using namespace std;




const int MAX_PHOTONS = 1000000;


int main()
{


	// Create the medium in which the photons will be fired.
	Medium *tissue = new Medium();
	
	// Add the layer to the medium.  NOTE:  destruction of the 'Layer' object is
	// handled in the 'tissue' object.
							//mu_a, mu_s,  n, start, end
	tissue->addLayer(new Layer(0.0, 30.0, 1.33, 0,  10));
	
	

	// Add the pressure map object to the medium and load the pressure data.
	//tissue->addPressureMap(new PressureMap("testing.txt"));
	tissue->addPressureMap(new PressureMap("pressure-at-25us.txt"));
	tissue->loadPressure();

										// x,  y,  z (photon)
	//cout << "pressure = " << tissue->getPressureFromGridCoords(31, 11, 31) << endl;
															// z , y , x (pressure)


	// Allocate the planar grid and set it in the tissue.
	double *Cplanar = (double*)malloc(sizeof(double) * 101);
	tissue->setPlanarArray(Cplanar);
		
	// Notify user of execution.
	cout << "Launching photons (" << MAX_PHOTONS << ")...\n";

	// Capture the time before launching photons into the medium.
	clock_t start, end;
	start = clock();
	

	// Let boost decide how many threads to run on this architecture.
	const int NUM_THREADS = boost::thread::hardware_concurrency();
	//const int NUM_THREADS = 2;
	
	// Each thread needs it's own photon object to run, so we need to create
	// an equal amount of photon objects as threads.
	const int NUM_PHOTON_OBJECTS = NUM_THREADS;

	Photon photons[NUM_PHOTON_OBJECTS];
	boost::thread threads[NUM_THREADS];

	// Init the random number generator.
	srand(time(0));

	// Used to seed the RNG.
	unsigned int s1, s2, s3, s4;

	// Create the threads and give them photon objects to run.
	// Each photon object is run MAX_PHOTONS/NUM_THREADS times, which essentially
	// splits up the work (i.e. photon propagation) amongst many workers.
	for (int i = 0; i < NUM_PHOTON_OBJECTS; i++)
	{
		s1 = rand() + 128;
		s2 = rand() + 128;
		s3 = rand() + 128;
		s4 = rand() + 128;

		cout << "Launching photon object" << i << " iterations: " << MAX_PHOTONS/NUM_THREADS << endl;
		threads[i] = boost::thread(&Photon::injectPhoton, &photons[i], tissue, MAX_PHOTONS/NUM_THREADS,
									s1, s2, s3, s4);

	}

	// Join all created threads once they have done their work.
	for (int i = 0; i < NUM_PHOTON_OBJECTS; i++)
	{
		threads[i].join();
	}

	
	// Print out the elapsed time it took from beginning to end.
	end = ((double)clock() - start) / CLOCKS_PER_SEC;
	cout << "Time elapsed: " << end << endl;
	
	// Print the matrix of the photon absorptions to file.
	tissue->printGrid(MAX_PHOTONS);
	
	// Clean up memory allocated memory on the heap.
	delete tissue;
	




	return 0;
}




