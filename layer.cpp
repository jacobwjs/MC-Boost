

#include "layer.h"

// Default constructor values if nothing is specified.
Layer::Layer(void)
{	
	// Set scattering and absorption properties of layer.
	mu_a = 1.0;		// cm^-1
	mu_s = 100.0;		// cm^-1

	
	albedo = mu_s/(mu_s + mu_a);
	g = 0.90;
	refractive_index = 1.33;

	// The depth at which this layer starts (cm).
	depth_start = 0;
	
	// The depth at which this layer ends (cm).
	depth_end = 10;
}

Layer::Layer(double mu_a, double mu_s, double refractive_index,
			 double depth_start, double depth_end)
{
	this->mu_a = mu_a;
	this->mu_s = mu_s;
	this->mu_t = mu_a + mu_s;
	albedo = mu_s/(mu_s + mu_a);
	g = 0.90;
	this->refractive_index = refractive_index;
	
	this->depth_start = depth_start;
	this->depth_end = depth_end;
}

Layer::~Layer(void)
{
		// Free any memory allocated on the heap by this object.
}

void Layer::setAbsorpCoeff(double mu_a)
{
	this->mu_a = mu_a;
	
	// If we ever update the absorption coefficient we need to update the
	// transmission coefficient and albedo similarly.
	this->mu_t = mu_a + mu_s;
	updateAlbedo();
}

void Layer::setScatterCoeff(double mu_s)
{
	this->mu_s = mu_s;

	// If we ever update the scattering coefficient we need to update the
	// transmission coefficient albedo similarly.
	this->mu_t = mu_a + mu_s;
	updateAlbedo();
}


void Layer::updateAlbedo()
{
	albedo = mu_s/(mu_s + mu_a);
}

		


