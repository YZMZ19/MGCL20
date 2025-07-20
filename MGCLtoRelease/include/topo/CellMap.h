/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#ifndef _MGCellMap_HH_
#define _MGCellMap_HH_
#include <map>
#include "mg/MGCL.h"

///@cond

class MGCellNB;

class MGCellMap{
public:
	typedef std::map<const MGCellNB*, MGCellNB*> map;

	MGCellMap(){;};
	virtual ~MGCellMap(){;};

	std::map<const MGCellNB*, MGCellNB*> m_map;

};

///@endcond

#endif
