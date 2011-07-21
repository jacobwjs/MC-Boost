// Defines attributes of a layer.
#ifndef LAYER_H
#define LAYER_H


#include "absorber.h"
#include "coordinates.h"
#include <vector>



class Layer
{	

public:
	Layer(void);
	Layer(double mu_a, double mu_s, double ref_index, double anisotropy,
		  double depth_start, double depth_end);
	~Layer(void);


    // Returns the absorption coefficient of the layer.
	double	getAbsorpCoeff(void) const	{return mu_a;}
    // Returns the absorption coeffiecient of the layer based on the photon's coordinates
    // Checks are made to see if the photon has made it's way into an absorber as well.
    double  getAbsorpCoeff(const boost::shared_ptr<Vector3d> photonVector);
    
    // Returns the scattering coefficient of the layer.
	double	getScatterCoeff(void) const	{return mu_s;}
    double  getScatterCoeff(const boost::shared_ptr<Vector3d> photonVector) const;
    
    // Returns total interaction coefficient (mu_a + mu_s).
	double	getTotalAttenuationCoeff(void) const	{return mu_t;}
    double  getTotalAttenuationCoeff(const boost::shared_ptr<Vector3d> photonVector) const;
    
    // Return the albedo
	double	getAlbedo(void) const			{return albedo;}
    
	double	getAnisotropy(void) 		{return g;}
    double  getAnisotropy(const boost::shared_ptr<Vector3d> photonVector) const;
    
	double 	getDepthStart(void) const 		{return depth_start;}
	double  getDepthEnd(void)	const		{return depth_end;}
    
	double	getRefractiveIndex(void) const	{return refractive_index;}
    double  getRefractiveIndex(const boost::shared_ptr<Vector3d> photonVector) const;

	void	setAbsorpCoeff(const double mu_a);
	void	setScatterCoeff(const double mu_s);
	void	updateAlbedo();
    
    void    addAbsorber(Absorber * a);

	
private:
	// Anisotropy factor.
	double g;
	
	// Absorption coefficient
	double mu_a;
	
	// Scattering coefficient
	double mu_s;
	
	// Transmission coefficient
	double mu_t;
	
	// The refractive index of the layer
	double refractive_index;
	
	// The width of the layer.
	//double radial_size;
	
	// z-coordinate value at which the layer starts.
	double depth_start;
	
	// z-coordinate value at which the layer ends.
	double depth_end;
	
	// Albedo of the layer.
	double albedo;	
    
    // A vector that holds all the abosrbers in this layer.
    std::vector<Absorber *> p_absorbers;
	
};

#endif // end LAYER_H


