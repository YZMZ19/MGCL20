/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/

#if !defined( __MGIGESPARAMLINE_H__)
#define __MGIGESPARAMLINE_H__

#include <memory>
#include <string>
#include "mg/MGCL.h"

///MGIgesParamLine describes a line of Parameter Data of an IGES file.

///This is used to output one line os Parameter Data Section.
class MGIgesParamLine{

// Constructors.

public:
	/// Constructs an object of class MGIgesParamLine.
	MGIgesParamLine();

	///Construct inputting only DEpoiter.
	MGIgesParamLine(int DEpointer);

	///one_line.size() must be <=64.
	MGIgesParamLine(std::string&& one_line, int DEpointer);

	///Get the reference of paramline string area m_paramLine.
	std::string& paramLine(){return m_paramLine;};

	///get de pointer
	int DEpointer()const{return m_DE_back_pointer;};

private:
	std::string m_paramLine;///<parameter data entry line data(from column 1 to 64)
					///<that does not contain data entry back pointer data.
	int m_DE_back_pointer;	///<Data entry back pointer. The line number is
					///<obtained by MGIges::DEpointer_to_lnumber().
};

#endif // __MGIGESPARAMLINE_H__
