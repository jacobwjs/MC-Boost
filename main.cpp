/*
 * Copyright BMPI 2011
 * J.W. Staley - MIRA, 
 *               Biomedical Photonics Imaging Group (BMPI), 
 *				 University of Twente
 *
 */


//#define DEBUG 1

#include "photon.h"
#include "medium.h"
#include "layer.h"
#include "sphereAbsorber.h"
#include "cylinderAbsorber.h"
#include "coordinates.h"
#include "vector3D.h"
#include "vectorMath.h"
#include "logger.h"
#include "circularDetector.h"
#include <cmath>
#include <ctime>
#include <vector>
#include <boost/thread/thread.hpp> 
#include <boost/lexical_cast.hpp>
#include <string>
#include <iostream>
using std::cout;
using std::endl;




// Number of photons to simulate.
const int MAX_PHOTONS = 1000000;

// Used to append to saved data files.
time_t epoch;
struct tm *ptr_ts;
std::string getCurrTime(void);


// Testing routines.
void testVectorMath(void);



// Simulation routines.
void runMonteCarlo(void);



int main()
{
	
    //testVectorMath();
    
	runMonteCarlo();
    
    
	return 0;
}


std::string getCurrTime(void)
{
    
    // Set current time variable to be used with naming data files that are saved from the simulations.
    epoch = time(NULL);
    ptr_ts = localtime(&epoch);
    
    return (boost::lexical_cast<std::string>(ptr_ts->tm_hour) + "_" +
            boost::lexical_cast<std::string>(ptr_ts->tm_min) + "_" +
            boost::lexical_cast<std::string>(ptr_ts->tm_sec));
}






void runMonteCarlo(void)
{
    
    

    // The logger is a singleton.  To bypass any problems with using singletons in a multi-threaded applicaton
    // initialization occurs in main before any threads are spawned.
    std::string exit_data_file;
    //file = "Absorber-data.txt";
    //Logger::getInstance()->openAbsorberFile(file);
    
    
    // The dimensions of the medium.
    //
    double X_dim = 2.0f;      // [cm]
    double Y_dim = 2.0f;      // [cm]
    double Z_dim = 2.0f;      // [cm]
    
    
	// Create the medium in which the photons will be propagate.
    //
	Medium *tissue = new Medium(X_dim, Y_dim, Z_dim);
	
    
    // Define a layer in the tissue.
    //
    double mu_a = 1.0f;
    double mu_s = 10.0f;
    double refractive_index = 1.33f;
    double anisotropy = 0.9;
    double start_depth = 0.0f; // [cm]
    double end_depth = Z_dim; // [cm]
    Layer *tissueLayer0 = new Layer(mu_a, mu_s, refractive_index, anisotropy, start_depth, end_depth);
    
    
    
    // Define a spherical absorber.
    //
    //    SphereAbsorber *absorber0 = new SphereAbsorber(.5, X_dim/2, Y_dim/2, Z_dim/2);
    //    absorber0->setAbsorberAbsorptionCoeff(3.0f);
    //    absorber0->setAbsorberScatterCoeff(mu_s);
    //    tissueLayer0->addAbsorber(absorber0);
    
    // Create a spherical detector.
    Detector *detector;
    CircularDetector circularExitDetector(0.5f, Vector3d(X_dim/2, Y_dim/2, Z_dim));
    circularExitDetector.setDetectorPlaneXY();  // Set the plane the detector is orientated on.
    detector = &circularExitDetector;
    
    
    // Add the objects to the medium.
    //
    tissue->addLayer(tissueLayer0);
    tissue->addDetector(detector);
    
    // Define the initial location of injection of the photons.
    //
    coords injectionCoords;
    injectionCoords.x = X_dim/2; // Centered
    injectionCoords.y = Y_dim/2; // Centered
    injectionCoords.z = 1e-15f;   // Just below the surface of the top-most layer.
	
    
	// Allocate the planar fluence grid and set it in the tissue.
    //	double *Cplanar = (double*)malloc(sizeof(double) * 101);
    //	tissue->setPlanarArray(Cplanar);
    
	
	// Let boost decide how many threads to run on this architecture.
	const int NUM_THREADS = boost::thread::hardware_concurrency();
	//const int NUM_THREADS = 1;
    
	// Each thread needs it's own photon object to run, so we need to create
	// an equal amount of photon objects as threads.
	const int NUM_PHOTON_OBJECTS = NUM_THREADS;
    
    // Photon array.  Each object in the array will be assigned their own seperate CPU core to run on.
	Photon photons[NUM_PHOTON_OBJECTS];
	boost::thread threads[NUM_THREADS];
    
    
	// Used to seed the RNG.
    //
	unsigned int s1, s2, s3, s4;
    
	// Capture the time before launching photons into the medium.
    //
	clock_t start, end;

    
    

        
        // Init the random number generator.
        //
        srand(time(0));
        
        // Open a file for each time step which holds exit data of photons
        // when they leave the medium through the detector aperture.
        //
        exit_data_file = "./Log/Exit-data/exit-aperture-" + boost::lexical_cast<std::string>(time(0)) + ".txt";
        Logger::getInstance()->openExitFile(exit_data_file);
        
        // Grab the start time before the simulation runs.
        start = clock();
        
        
        // Create the threads and give them photon objects to run.
        // Each photon object is run MAX_PHOTONS/NUM_THREADS times, which essentially
        // splits up the work (i.e. photon propagation) amongst many workers.
        //
        for (int i = 0; i < NUM_PHOTON_OBJECTS; i++)
        {
            // The state variables need to be >= 128.
            s1 = rand() + 128;
            s2 = rand() + 128;
            s3 = rand() + 128;
            s4 = rand() + 128;
            
            cout << "Launching photon object" << i << " iterations: " << MAX_PHOTONS/NUM_THREADS << endl;
            threads[i] = boost::thread(&Photon::injectPhoton, &photons[i], tissue, MAX_PHOTONS/NUM_THREADS,
                                       s1, s2, s3, s4, injectionCoords);
            
        }
        
        // Join all created threads once they have done their work.
        for (int i = 0; i < NUM_PHOTON_OBJECTS; i++)
        {
            threads[i].join();
        }
        


	
	// Print out the elapsed time it took from beginning to end.
	end = ((double)clock() - start) / CLOCKS_PER_SEC;
	cout << "\n\nTotal time elapsed: " << end << endl;
    
    
	// Print the matrix of the photon absorptions to file.
	//tissue->printGrid(MAX_PHOTONS);
	
	// Clean up memory allocated memory on the heap.
	if (tissue)
		delete tissue;
	
}





// Simple routine to test the vectorMath library.
void testVectorMath(void)
{
    
    boost::shared_ptr<Vector3d> p0(new Vector3d(2.0f, 1.0f, 1.0f));
    boost::shared_ptr<Vector3d> p1(new Vector3d(3.5f, 1.5f, 11.0f));
    boost::shared_ptr<Vector3d> dir;
    boost::shared_ptr<Vector3d> c0(new Vector3d(0.0f, 0.0f, 11.0f));
    boost::shared_ptr<Vector3d> c1(new Vector3d(2.0f, 3.0f, 11.0f));
    boost::shared_ptr<Vector3d> c2(new Vector3d(11.0f, 13.5f, 11.0f));
    boost::shared_ptr<Vector3d> n;
    
    
    n = VectorMath::crossProduct((*c1 - *c0), (*c2 - *c0));
    //n.reset(new Vector3d(1.0f, 2.0f, 3.0f));
    VectorMath::Normalize(n);
    
    double u = VectorMath::dotProduct(n, (*c0 - *p0)) / VectorMath::dotProduct(n, (*p1 - *p0));
    double THRESH = 0.0000000000001;
    if (u < 0.0f || u > 1.0f + THRESH)
        cout << "FALSE\n";
    
    cout << "n = " << n;
    cout << "u = " << u << endl;
    
    
    double z0 = p0->location.z;
    double z1 = p1->location.z;
    
    double y0 = p0->location.y;
    double y1 = p1->location.y;
    
    double x0 = p0->location.x;
    double x1 = p1->location.x;
    
    double distToPlane = abs(VectorMath::dotProduct(n, (*c0-*p0)) / VectorMath::Length(n));
    //D = VectorMath::Distance(c0, p0);
    cout << "distance to plane = " << distToPlane << endl;
    cout << (*c0 - *p0);
    
    double z= z0 + (z1-z0)*u;
    double y = y0 + (y1-y0)*u;
    double x = x0 + (x1-x0)*u;
    
    boost::shared_ptr<Vector3d> intersectPoint(new Vector3d(x, y, z));
    cout << "intersection point = " << intersectPoint;
    
    
    CircularDetector detector(1.0f, Vector3d(1.0f, 1.0f, 11.0f));
    detector.setDetectorPlaneXY();  // Set the plane the detector is orientated on.
    bool hitDetector = detector.photonPassedThroughDetector(p0, p1);
    cout << "hitDetector = " << hitDetector << endl;
    
    
	//    z = (*p) - (*x);
	//    cout << "p - x = " << z->X() << endl;
	//
	//    using namespace VectorMath;
	//    z = VectorMath::crossProduct(p, x);
	//    double d = VectorMath::dotProduct(p, x);
	//    cout << "Dot product = " << d << endl;
	//    VectorMath::Length(z);
	//
	//    //VectorMath vmath;
	//    //z = vmath.crossProduct(p, x);
	//    cout << "cross x = " << z->location.x << "\ncross y = " << z->location.y << "\ncross z = " << z->location.z << endl;
	//
}



