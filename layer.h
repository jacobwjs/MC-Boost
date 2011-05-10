// Defines attributes of a layer.
#ifndef LAYER_H
#define LAYER_H


class Layer
{	

public:
	Layer(void);
	Layer(double mu_a, double mu_s, double ref_index,
		  double depth_start, double depth_end);
	~Layer(void);


	double	getAbsorpCoeff(void) 	{return mu_a;}
	double	getScatterCoeff(void)	{return mu_s;}
	double	getTotalAttenuationCoeff(void)	{return mu_t;}
	double	getAlbedo(void) {return albedo;}
	double	getAnisotropy(void) {return g;}
	double 	getDepthStart(void) 	{return depth_start;}
	double  getDepthEnd(void)		{return depth_end;}

	void	setAbsorpCoeff(const double mu_a);
	void	setScatterCoeff(const double mu_s);
	void	updateAlbedo();

	
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
	
};

#endif // end LAYER_H


