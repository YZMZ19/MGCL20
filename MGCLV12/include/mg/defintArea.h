/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#ifndef _MGDefintArea_HH_
#define _MGDefintArea_HH_
#include "mg/MGCL.h"

/// @cond

#define mgdefintlen 608

/// Generate points and weights of the DE formula for DEFINT
/// in double precision. These are static data and initialized once.
/// MGDefintArea is a proprietary class for mgDefint(integration).
class MG_DLL_DECLR MGDefintArea{

private:
	MGDefintArea();
	MGDefintArea(const MGDefintArea&);
	MGDefintArea& operator=(const MGDefintArea&);

public:

	///Get the static instance.
	static MGDefintArea& instance();

	///Initialize the data.
	void init();

	double m_am[mgdefintlen*2]	/* was [608][2] */,
		m_a0[2], m_ap[mgdefintlen*2]	/* was [608][2] */,
		m_b0, m_bb[mgdefintlen];

	int m_nend;
	int m_npow;
	double m_eps0;
};

/// @endcond

#endif
