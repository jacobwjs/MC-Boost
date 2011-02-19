/*
 * Copyright BMPI 2010
 * J.W. Staley - MIRA, 
 *               Biomedical Photonics Imaging Group (BMPI), 
 *				 University of Twente
 *
 */



// For matlab integration
#include "mex.h"

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



void createFluence(double *vals, const double, const double);




void mexFunction( int nlhs, mxArray *plhs[],
				 int nrhs, const mxArray *prhs[])
{
 
	
	
//int main()
//{
	
	
	// Create the medium in which the photons will be fired.
	Medium *tissue = new Medium();
	
	// Add the layer to the medium.  NOTE:  destruction of the 'Layer' object is
	// handled in the 'tissue' object.
	tissue->addLayer(new Layer());
	
	// Create a photon object, which will be propogated through the medium
	// MAX_PHOTON times.
	Photon *photon = new Photon();
	
	
	
	/* Matlab specific initialations */
	// ------------------------------------------------------
	// Number of rows and colums of the array given to Matlab.
	
	mwSize mrows,ncols;
	mrows = 1;
	ncols = 101;
	
	// set the output pointer to the output matrix 
	plhs[0] = mxCreateDoubleMatrix(mrows, ncols, mxREAL);
	double *Cplanar = mxGetPr(plhs[0]);
	
	tissue->setPlanarArray(Cplanar);
	
	
	
	//double *Cplanar = (double*)malloc(sizeof(double) * 101);
	//tissue->setPlanarArray(Cplanar);
	
	// Initial injection location of a photon.
	int x = 0;
	int y = 0;
	
	// Capture the time before launching photons into the medium.
	clock_t start, end;
	start = clock();
	
	// Simulate photons being injected into the medium. 
	for (int i = 0; i < MAX_PHOTONS; i++) 
    {
		tissue->injectPhoton(x, y, photon);
	}

	
	
	
	
	
	
	// Print out the elapsed time it took from beginning to end.
	end = ((double)clock() - start) / CLOCKS_PER_SEC;
	cout << "Time elapsed: " << end << endl;
	
	
	
	
	// Print the matrix of the photon absorptions to file.
	//Medium *ptrMedium = &tissue;
	//tissue->printGrid(MAX_PHOTONS);
	
	
	
	//double *Cplanar = tissue->getPlanarGrid();

	// Calculate the fluence.
	//createFluence(Cplanar, tissue->getBins(), tissue->getRadialSize());
	
	// Copy over values to the matlab "array".
	/*
	for (int i = 0; i < 101; i++)
	{
		to_matlab[i] = Cplanar[i];
	}
	*/
	
	
	// Clean up memory allocated memory on the heap.
	delete tissue;
	delete photon;
	
	//return 0;
}


void createFluence(double *Cplanar, const double BINS, const double rad_size)
{
	// FIXME:  Assuming homogenous medium.
	double mu_a = 1.0;
	
	double *temp = Cplanar;
	double radial_bin_size = BINS / rad_size;
	double fluencePlanar = 0;
	double r = 0;
	double shellVolume = 0;
	
	shellVolume = radial_bin_size;
	for (int ir = 0; ir <= BINS-1; ir++, Cplanar++) {
		r = (ir + 0.5)*radial_bin_size;
		
		*Cplanar = (*Cplanar)/MAX_PHOTONS/shellVolume/mu_a;
	}

}

