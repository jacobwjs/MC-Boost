#include "debug.h"
#include "medium.h"
#include <cmath>
#include <assert.h>

Medium::Medium()
{
	radial_size = 3.0;	// Total range in which bins are extended (cm).
	num_radial_pos = MAX_BINS-1;	// Set the number of bins.
	radial_bin_size = radial_size / num_radial_pos;
	
	Cplanar = NULL;
}

Medium::~Medium()
{
	
	// Free the memory for layers that were added to the medium.
	for (vector<Layer *>::iterator i = p_layers.begin(); i != p_layers.end(); ++i)
		delete *i;
}



void Medium::setPlanarArray(double *array)
{
	Cplanar = array;
	// Initialize all the bins to zero since they will serve as accumulators.
	for (int i = 0; i < MAX_BINS; i++) {
		Cplanar[i] = 0;
	}
}

// Add the layer to the medium by pushing it onto the vector container.
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
	if (ir >= num_radial_pos) {
		ir = num_radial_pos;
	}

	Cplanar[ir] += energy;

}


// Return the layer in the medium at the passed in depth 'z'.
// We iterate through the vector which contains pointers to the layers.
// When the correct layer is found from the depth we return the layer object.
Layer * Medium::getLayerFromDepth(double z)
{
	vector<Layer *>::iterator it;
	for (it = p_layers.begin(); it != p_layers.end(); it++) {
		// Find the layer we are in within the medium based on the depth (i.e. z)
		// that was passed in.  Break from the loop when we find the correct layer.
		if ((*it)->getDepthStart() <= z && (*it)->getDepthEnd() >= z)
			break;
	}
	
	// Return layer based on the depth passed in.
	return *it;
}


double Medium::getLayerAbsorptionCoeff(double z)
{
	double absorp_coeff = 0;
	vector<Layer *>::iterator it;
	for (it = p_layers.begin(); it != p_layers.end(); it++) {
		// Find the layer we are it in the medium based on the depth (i.e. z)
		// that was passed in.  Break from the loop when we find the correct layer.
		if ((*it)->getDepthStart() <= z && (*it)->getDepthEnd() >= z) {
			absorp_coeff = (*it)->getAbsorpCoeff();
			break;
		}
	}
	
	// If not found, report error.
	assert(absorp_coeff != 0);
	
	// Return the absorption coefficient value.
	return absorp_coeff;
}


double Medium::getLayerScatterCoeff(double z)
{
	double scatter_coeff = 0;
	vector<Layer *>::iterator it;
	for (it = p_layers.begin(); it != p_layers.end(); it++) {
		// Find the layer we are it in the medium based on the depth (i.e. z)
		// that was passed in.  Break from the loop when we find the correct layer.
		if ((*it)->getDepthStart() <= z && (*it)->getDepthEnd() >= z) {
			scatter_coeff = (*it)->getScatterCoeff();
			break;
		}
	}
	
	// If not found, report error.
	assert(scatter_coeff != 0);
	
	// Return the scattering coefficient for the layer that resides at depth 'z'.
	return scatter_coeff;
}


double Medium::getAnisotropyFromDepth(double z)
{
	double anisotropy = 0;
	vector<Layer *>::iterator it;
	for (it = p_layers.begin(); it != p_layers.end(); it++) {
		// Find the layer we are it in the medium based on the depth (i.e. z)
		// that was passed in.  Break from the loop when we find the correct layer.
		if ((*it)->getDepthStart() <= z && (*it)->getDepthEnd() >= z) {
			anisotropy = (*it)->getAnisotropy();
			break;
		}
	}
	
	// If not found, report error.
	assert(anisotropy != 0);
	
	// Return the anisotropy value for the layer that resides at depth 'z'.
	return anisotropy;
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
