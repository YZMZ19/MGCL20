/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#include "StdAfx.h"
#include "mg/Box.h"
#include "mg/Position.h"
#include "mg/Geometry.h"
#include "mg/Curve.h"
#include "mg/Surface.h"
#include "mg/Point.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

// MGGeometry
// Implementation of MGGeometry.

//Compute direction unit vector of the geometry.
MGUnit_vector MGGeometry::direction(const MGPosition& param) const{
	return mgZ_UVEC;
}

//Error allowed in the parameter space of the geometry.
double MGGeometry::parameter_error() const{
	MGBox prange=parameter_range();
	int n=prange.sdim();
	double error=0., ework;
	for(int i=0; i<n; i++){
		ework=prange(i).relative_error();
		error+=ework*ework;
	}
	return sqrt(error);
}
