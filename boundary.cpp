#include "boundary.h"

Boundary::Boundary(void)
{
	depth = 10;
	width = 20;
	boundary_energy = 0;
	
}

Boundary::Boundary(double depth, double width)
{
	this->depth = depth;
	this->width = width;
	this->boundary_energy = 0;
}


