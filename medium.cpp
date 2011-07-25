#include "debug.h"
#include "medium.h"
#include <cmath>
#include <assert.h>

#undef DEBUG

Medium::Medium()
{
    x_bound = y_bound = z_bound = 10;   // defaults to 10x10x10 (cm)
   	this->initCommon();
}

Medium::Medium(const int x, const int y, const int z)
{
    
    this->x_bound = x;
    this->y_bound = y;
    this->z_bound = z;
	this->initCommon();	
}



void Medium::initCommon(void)
{
    // Assuming air surrounds tissue.
    refractive_index = 1.0;
    
    radial_size = 3.0;	// Total range in which bins are extended (cm).
	num_radial_pos = MAX_BINS-1;	// Set the number of bins.
	radial_bin_size = radial_size / num_radial_pos;
    Cplanar = NULL;
}

Medium::~Medium()
{

	// Free the memory for layers that were added to the medium.
	for (vector<Layer *>::iterator it = p_layers.begin(); it != p_layers.end(); it++)
    {
        (*it)->writeAbsorberData();
        delete *it;
    }
}



void Medium::setPlanarArray(double *array)
{
	Cplanar = array;
	// Initialize all the bins to zero since they will serve as accumulators.
	int i;
	for (i = 0; i < MAX_BINS; i++) {
		Cplanar[i] = 0;
	}
}

// Add the layer to the medium by pushing it onto the vector container.
void Medium::addLayer(Layer *layer)
{
	p_layers.push_back(layer);
}


void Medium::addDetector(Detector *detector)
{
    p_detectors.push_back(detector);
}


void Medium::absorbEnergy(const double z, const double energy)
{
#ifdef DEBUG
	cout << "Updating bin...\n";
#endif

	double r = fabs(z);
	int ir = (r/radial_bin_size);
	if (ir >= num_radial_pos) {
		ir = num_radial_pos;
	}

	Cplanar[ir] += energy;

}


void Medium::absorbEnergy(const double *energy_array)
{
	int i;
	boost::mutex::scoped_lock lock(m_mutex);
	for (i = 0; i < MAX_BINS; i++) {
		// Grab the lock to serialize threads when updating
		// the global planar detection array in the Medium.
		Cplanar[i] += energy_array[i];
	}
}


// See if photon has crossed the detector plane.
int Medium::photonHitDetectorPlane(const boost::shared_ptr<Vector3d> p0)
{
    bool hitDetectorNumTimes = 0;
    // Free the memory for layers that were added to the medium.
	for (vector<Detector *>::iterator it = p_detectors.begin(); it != p_detectors.end(); it++)
    {
		if ((*it)->photonHitDetector(p0))
            hitDetectorNumTimes++;
    }
    
    return hitDetectorNumTimes;
}

Layer * Medium::getLayerAboveCurrent(Layer *currentLayer)
{
	// Ensure that the photon's z-axis coordinate is sane.  That is,
	// it has not left the medium.
	assert(currentLayer != NULL);

	// If we have only one layer, no need to iterate through the vector.
	// And we should return NULL since there is no layer above us.
	if (p_layers.size() == 1)
		return NULL;
    


	// Otherwise we walk the vector and return 'trailer' since it is the
	// one before the current layer (i.e. 'it').
	vector<Layer *>::iterator it;
	vector<Layer *>::iterator trailer;
	it = p_layers.begin(); // Get the first layer from the array.
    
    // If we are at the top of the medium there is no layer above, so return NULL;
    if (currentLayer == (*it))
        return NULL;
    
	while(it != p_layers.end()) {
		trailer = it;  // Assign the trailer to the current layer.
		it++;         // Advance the iterator to the next layer.

		// Find the layer we are in within the medium based on the depth (i.e. z)
		// that was passed in.  Break from the loop when we find the correct layer
		// because trailer will be pointing to the previous layer in the medium.
		//if ((*it)->getDepthStart() <= z && (*it)->getDepthEnd() >= z)
		if ((*it) == currentLayer)
            break;
	}

	// Sanity check.  If the trailer has made it to the end, which means
	// the iterator made it past the end, then there
	// was no previous layer found, and something went wrong.
	if (trailer == p_layers.end())
		return NULL;

	// If we make it here, we have found the previous layer.
	return *trailer;
}


Layer * Medium::getLayerBelowCurrent(double z)
{
	// Ensure that the photon's z-axis coordinate is sane.  That is,
	// it has not left the medium.
	assert(z >= 0 && z <= z_bound);

	// If we have only one layer, no need to iterate through the vector.
	// And we should return NULL since there is no layer below us.
	if (p_layers.size() == 1)
		return NULL;
    
    // The case where there is no layer below is since we are at the bottom of the
    // medium.
    if (z == z_bound)
        return NULL;


	vector<Layer *>::iterator it;
	for (it = p_layers.begin(); it != p_layers.end(); it++) {
		// Find the layer we are in within the medium based on the depth (i.e. z)
		// that was passed in.  Break from the loop when we find the correct layer.
		if ((*it)->getDepthStart() <= z && (*it)->getDepthEnd() >= z) {
			return *(++it);
		}
	}

	// If the above loop never returned a layer it means we made it through the list
	// so there is no layer below us, therefore we return null.
	return NULL;


}


// Return the layer in the medium at the passed in depth 'z'.
// We iterate through the vector which contains pointers to the layers.
// When the correct layer is found from the depth we return the layer object.
Layer * Medium::getLayerFromDepth(double z)
{
	// Ensure that the photon's z-axis coordinate is sane.  That is,
	// it has not left the medium.
	assert(z >= 0 && z <= z_bound);

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
	// Ensure that the photon's z-axis coordinate is sane.  That is,
	// it has not left the medium.
	assert(z >= 0 && z <= z_bound);

	double absorp_coeff = -1;
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
	assert(absorp_coeff != -1);

	// Return the absorption coefficient value.
	return absorp_coeff;
}


double Medium::getLayerScatterCoeff(double z)
{
	// Ensure that the photon's z-axis coordinate is sane.  That is,
	// it has not left the medium.
	assert(z >= 0 && z <= z_bound);

	double scatter_coeff = -1;
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
	assert(scatter_coeff != -1);

	// Return the scattering coefficient for the layer that resides at depth 'z'.
	return scatter_coeff;
}


double Medium::getAnisotropyFromDepth(double z)
{
	// Ensure that the photon's z-axis coordinate is sane.  That is,
	// it has not left the medium.
	assert(z >= 0 && z <= z_bound);

	double anisotropy = -1;
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
	assert(anisotropy != -1);

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
