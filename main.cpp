/*
 * Copyright BMPI 2010
 * J.W. Staley - MIRA, 
 *               Biomedical Photonics Imaging Group (BMPI), 
 *				 University of Twente
 *
 */

/* TO DO:
 * ------
 * - Set default thread stack size
 * - Update printGrid to calculate fluences based on layers. 
 */



#include "stdio.h"
#include "photon.h"
#include "medium.h"
#include "layer.h"
#include <time.h>
#include <iostream>
#include <vector>

using namespace std;

const int MAX_PHOTONS = 5000;

//#define DEBUG 1

int main()
{
	
	// Create the medium in which the photons will be fired.
	Medium *tissue = new Medium();
	// Add the layer to the medium.  NOTE:  destruction of the 'Layer' object is
	// handled in the 'tissue' object.
	tissue->addLayer(new Layer());
	
	Photon *photon = new Photon();
	
	int i;
	for (i = 0; i < MAX_PHOTONS; i++) 
    {
        tissue->injectPhoton(0, 0, photon);

	
	// Capture the time before launching photons into the medium.
	clock_t start, end;
	start = clock();
	
	
	
	// Print out the elapsed time it took from beginning to end.
	end = ((double)clock() - start) / CLOCKS_PER_SEC;
	cout << "Time elapsed: " << end << endl;
	
	
	
	
	// Print the matrix of the photon absorptions to file.
	//Medium *ptrMedium = &tissue;
	tissue->printGrid(MAX_PHOTONS);
	delete tissue;
	
	
	return 0;
}
