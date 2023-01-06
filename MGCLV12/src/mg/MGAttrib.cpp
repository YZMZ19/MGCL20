/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#include "StdAfx.h"
#include "mg/Attrib.h"
#include "mg/Ifstream.h"
#include "mg/Ofstream.h"
#include "mgGL/Context.h"
#include "mgGL/GLAttrib.h"
#include "mg/GelFactory.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

// MGAttrib
// Implementation of MGAttrib.

AUTO_GEL_REGISTER(MGContext, MGCONTEXT_TID);
AUTO_GEL_REGISTER(MGAppearance, MGAPPEARANCE_TID);

//Read all member data.
//void MGAttrib::ReadMembers(MGIfstream& buf){;}
//Write all member data
//void MGAttrib::WriteMembers(MGOfstream& buf)const{;}

