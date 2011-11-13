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
#include "pressureMap.h"
#include "displacementMap.h"
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
const int MAX_PHOTONS = 5e5;

// Used to append to saved data files.
time_t epoch;
struct tm *ptr_ts;
std::string getCurrTime(void);


// Testing routines.
void testVectorMath(void);
void testDisplacements(void);
void testPressures(void);



// Simulation routines.
void runAcoustoOptics(void);



int main()
{
	
    //testVectorMath();
    //testDisplacements();
    //testPressures();
    
	runAcoustoOptics();
    
    
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






void runAcoustoOptics(void)
{
    
    

    // The logger is a singleton.  To bypass any problems with using singletons in a multi-threaded applicaton
    // initialization occurs in main before any threads are spawned.
    std::string exit_data_file;
    std::string rng_seed_file = "rng_exit_aperture_seeds.txt";
    Logger::getInstance()->createRNGSeedFile(rng_seed_file);
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
    double mu_s = 70.0f;
    double refractive_index = 1.0f;
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
    CircularDetector circularExitDetector(0.2f, Vector3d(X_dim/2, Y_dim/2, Z_dim));
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
    injectionCoords.z = 1e-15f;   // Just below the surface of the 'air' layer.
	
    
	
    
	// Create and add the pressure map object to the medium and load the pressure data.
	// tissue->addPressureMap(new PressureMap("testing.txt"));
	//
    
    const int pgrid_x = 64;  // Number of pixels in the kWave pressure grid.
	const int pgrid_y = 64;
	const int pgrid_z = 64;
/*
	string pressure_file = "./kWave-pressure/pressure";
    PressureMap *pmap = new PressureMap(pgrid_x, pgrid_z, pgrid_y, X_dim);
	tissue->addPressureMap(pmap);
	//cout << "pressure = " << tissue->getPressureFromGridCoords(31, 11, 31) << endl;
*/
    
	// Create and add the displacement map object to the medium and load the displacement data.
    //
	const int dgrid_x = pgrid_x; // Displacements are calculated from same simulation grid, therefore same size.
	const int dgrid_y = pgrid_z;
	const int dgrid_z = pgrid_y;
	string displacement_file = "./kWave-displacements/2cm/homo-medium/disp";
	//DisplacementMap *dmap = new DisplacementMap(dgrid_x, dgrid_z, dgrid_y, X_dim);
    //tissue->addDisplacementMap(dmap);
    
    
    // Set the value of the transducer frequency used in the K-Wave simulation.
    tissue->kwave.transducerFreq = 1.0e6;
    
    
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
	clock_t start, start_per_simulation, end;
	start = clock();
    
    
    // For each time step that K-Wave gave ultrasound data, propagate
    // photons through and track modulation due to the acoustic source.
	// Number of time steps that were executed in the K-Wave simulation
    // that produced displacement and pressure data.
	const int KWAVESIM_TIME_STEPS = 1;

	for (int dt = 1; dt <= KWAVESIM_TIME_STEPS; dt++)
    {
        // Capture the time at the beginning of this simulation step.
        start_per_simulation = clock();
        
        // Init the random number generator.
        //
        srand(13);
        
        // Open a file for each time step which holds exit data of photons
        // when they leave the medium through the detector aperture.
        //
        exit_data_file = "./Log/Exit-data/2cm/homo-medium/exit-aperture-" + boost::lexical_cast<std::string>(dt) + ".txt";
        Logger::getInstance()->openExitFile(exit_data_file);
        
        // Load a pressure map and displacement maps at time step number 'dt'.
        //tissue->loadPressure(pressure_file, dt);
        //tissue->loadDisplacements(displacement_file, dt);
        
        
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

//        	s1 = 292 + i*10;
//        	s2 = 338 + i*10;
//        	s3 = 2893 + i*10;
//        	s4 = 559232 + i*10;

            
            cout << "Launching photon object" << i << " iterations: " << MAX_PHOTONS/NUM_THREADS << endl;
            threads[i] = boost::thread(&Photon::injectPhoton, &photons[i], tissue, MAX_PHOTONS/NUM_THREADS,
                                       s1, s2, s3, s4, injectionCoords);
            
        }
        
        // Join all created threads once they have done their work.
        for (int i = 0; i < NUM_PHOTON_OBJECTS; i++)
        {
            threads[i].join();
        }
        
        // Print out the elapsed time it took for this simulation step.
        end = ((double)clock() - start_per_simulation) / CLOCKS_PER_SEC;
        cout << "Time elapsed for simulation (" << dt << "): " << end << endl;
        
    }
	
	// Print out the elapsed time it took from beginning to end.
	end = ((double)clock() - start) / CLOCKS_PER_SEC;
	cout << "\n\nTotal time elapsed: " << end << endl;
    
    
	// Print the matrix of the photon absorptions to file.
	//tissue->printGrid(MAX_PHOTONS);
	
	// Clean up memory allocated on the heap.
	if (tissue)
		delete tissue;
	
}


void testDisplacements(void)
{
    // Add the pressure map object to the medium and load the pressure data.
	//tissue->addPressureMap(new PressureMap("testing.txt"));
	const int dgrid_x = 64;  // Number of pixels in the kWave pressure grid.
	const int dgrid_y = 64;
	const int dgrid_z = 64;
    const int grid_size = 2;
	//string pressure_file = "C:/Users/StaleyJW/Desktop/Software/MC-Boost/kWave-pressure/pressure-at-25us.txt";
	string displacement_file = "d:/Displacement_Data/disp";
    DisplacementMap *dmap = new DisplacementMap(dgrid_x, dgrid_z, dgrid_y, grid_size);
    dmap->loadDisplacementMaps(displacement_file, 100);
    //tissue->addPressureMap(pmap);
	//tissue->loadPressure();
    
    
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



void testPressures(void)
{
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
	    double mu_s = 70.0f;
	    double refractive_index = 1.33f;
	    double anisotropy = 0.9;
	    double start_depth = 0.0f; // [cm]
	    double end_depth = Z_dim; // [cm]
	    Layer *tissueLayer0 = new Layer(mu_a, mu_s, refractive_index, anisotropy, start_depth, end_depth);



	    // Add the objects to the medium.
	    //
	    tissue->addLayer(tissueLayer0);


	    // Define the initial location of injection of the photons.
	    //
	    coords injectionCoords;
	    injectionCoords.x = X_dim/2; // Centered
	    injectionCoords.y = Y_dim/2; // Centered
	    injectionCoords.z = 1e-15f;   // Just below the surface of the 'air' layer.




		// Create and add the pressure map object to the medium and load the pressure data.
		// tissue->addPressureMap(new PressureMap("testing.txt"));
		//
	    const int pgrid_x = 64;  // Number of pixels in the kWave pressure grid.
		const int pgrid_y = 64;
		const int pgrid_z = 64;
		string pressure_file = "D:/kWave-pressure/pressure21.txt";
	    PressureMap *pmap = new PressureMap(pgrid_x, pgrid_z, pgrid_y, X_dim);
	    pmap->loadPressureMap(pressure_file);
		//tissue->addPressureMap(pmap);
		//tissue->loadPressure(pressure_file);


	    //cout << "Pressure(19,31,31) = " << pmap->getPressureFromGrid(19, 31, 31) << endl;
		//cout << "Pressure(54,31,31) = " << pmap->getPressureFromGrid(54, 31, 31) << endl;
		cout << "Pressure(64,11,31) = " << pmap->getPressureFromGrid(63, 10, 30) << endl;
		cout << "Pressure(64,32,64) = " << pmap->getPressureFromGrid(63, 31, 63) << endl;



		if (tissue)
			delete tissue;
}
