/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
//! @file
//!	@brief  Declaration for class MGIgesParamDELine.
//!	@author System fugen

#include "StdAfx.h"
#include "mgiges/IgesParamLine.h"

//!	@brief MGIgesParamLine describes a line of Parameter Data of an IGES file.
//This is used to output one line os Parameter Data Section.

//! Constructs an object of class MGIgesParamLine.
MGIgesParamLine::MGIgesParamLine()
:m_DE_back_pointer(0){;}

//Construct inputting only DEpoiter.
MGIgesParamLine::MGIgesParamLine(int DEpointer)
:m_DE_back_pointer(DEpointer){;}

//one_line.size() must be <=64.
MGIgesParamLine::MGIgesParamLine(std::string&& one_line, int DEpointer)
	:m_paramLine(std::move(one_line)),m_DE_back_pointer(DEpointer){;};
