/*
# This file is part of the Astrometry.net suite.
# Licensed under a 3-clause BSD style license - see LICENSE
*/

#ifndef PQUAD_H
#define PQUAD_H

#include <stdint.h>

/**
 This file is just required for testing purposes (of solver.c)
 */

struct potential_quad
{
	anbool scale_ok;
	int fieldA, fieldB;
	// distance-squared between A and B, in pixels^2.
	double scale;
	double costheta, sintheta;
	// (field pixel noise / quad scale in pixels)^2
	double rel_field_noise2;
	anbool* inbox;
	int ninbox;
	// Number of TRUE entries in the initialized prefix [0, ninbox).
	int eligible_count;
	// Optional immutable prefix counts for shared whole-field geometry.
	uint16_t* inbox_prefix;
	double* xy;
};
typedef struct potential_quad pquad;

#endif
