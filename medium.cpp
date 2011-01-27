#include "debug.h"
#include "medium.h"
#include <cmath>

Medium::Medium()
{
	radial_size = 3.0;	// Total range in which bins are extended (cm).
	num_radial_pos = MAX_BINS-1;	// Set the number of bins.
	radial_bin_size = radial_size / num_radial_pos;
	
	// Initialize all the bins to zero since they will serve as accumulators.
	for (int i = 0; i < MAX_BINS; i++) {
		Cplanar[i] = 0;
	}
	
	
}

Medium::~Medium()
{
	
	// Free the memory for layers that were added to the medium.
	for (vector<Layer *>::iterator i = p_layers.begin(); i != p_layers.end(); ++i)
		delete *i;
}

// Add the layer to the medium by pushing it onto the STL vector container.
void Medium::addLayer(Layer *layer)
{
	p_layers.push_back(layer);
}


void Medium::absorbEnergy(const double z, const double energy)
{
#ifdef DEBUG
	cout << "Updating bin...\n";
#endif
	
	double r = fabs(z);
	short ir = (short)(r/radial_bin_size);
	if (ir >= num_radial_pos) ir = num_radial_pos;


	Cplanar[ir] += energy;

}


void Medium::injectPhoton(const int x_start, const int y_start, Photon *photon)
{
     cout << "Stub... injectPhoton\n");
}




void Medium::printGrid(const int numPhotons)
{
	
	// Open the file we will write to.
	ofstream output;
	output.open("fluences.txt");
	
	// Print the header information to the file.
	//output << "r [cm] \t Fsph [1/cm2] \t Fcyl [1/cm2] \t Fpla [1/cm2]\n";
	//output << "r [cm] \t Fplanar[1/cm^2]\n";
	
	double mu_a = p_layers[0]->getAbsorpCoeff();
	double fluencePlanar = 0;
	double r = 0;
	double shellVolume = 0;
	
	for (int ir = 0; ir <= num_radial_pos; ir++) {
		r = (ir + 0.5)*radial_bin_size;
		shellVolume = radial_bin_size;
		fluencePlanar = Cplanar[ir]/numPhotons/shellVolume/mu_a;
		
		// Print to file with the value for 'r' in fixed notation and having a
		// precision of 5 decimal places, followed by the fluence in scientific
		// notation with a precision of 3 decimal places.
		output << fixed << setprecision(5) << r << "\t \t";
		output << scientific << setprecision(3) <<  fluencePlanar << "\n";
	}
	
	// close the file.
	output.close();
	
	
}
