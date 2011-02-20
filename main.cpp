/*
 * Copyright BMPI 2010
 * J.W. Staley - MIRA, 
 *               Biomedical Photonics Imaging Group (BMPI), 
 *				 University of Twente
 *
 */




#include "stdio.h"
#include "photon.h"
#include "medium.h"
#include "layer.h"
#include <time.h>
#include <iostream>
#include <vector>


using namespace std;

const int MAX_PHOTONS = 1000000;

//#define DEBUG 1



	
int main()
{
	
	
	// Create the medium in which the photons will be fired.
	Medium *tissue = new Medium();
	
	// Add the layer to the medium.  NOTE:  destruction of the 'Layer' object is
	// handled in the 'tissue' object.
	tissue->addLayer(new Layer());
	
	// Create a photon object, which will be propogated through the medium
	// MAX_PHOTON times.
	Photon *photon = new Photon();
	
	
	
	// Allocate the planar grid and set it in the tissue.
	double *Cplanar = (double*)malloc(sizeof(double) * 101);
	tissue->setPlanarArray(Cplanar);
	
	// Initial injection location of a photon.
	int x = 0;
	int y = 0;
	
	// Capture the time before launching photons into the medium.
	clock_t start, end;
	start = clock();
	

	
	/*
	// Simulate photons being injected into the medium. 
	for (int i = 0; i < MAX_THREADS; i++) 
    {
		//tissue->injectPhoton(x, y, photon);
	}
	 */
	
	// Non-threaded case.
	photon->injectPhoton(tissue, MAX_PHOTONS);
	
	
	
	// Print out the elapsed time it took from beginning to end.
	end = ((double)clock() - start) / CLOCKS_PER_SEC;
	cout << "Time elapsed: " << end << endl;
	
	
	
	
	// Print the matrix of the photon absorptions to file.
	tissue->printGrid(MAX_PHOTONS);
	
	
	// Clean up memory allocated memory on the heap.
	delete tissue;
	delete photon;
	
	return 0;
}


