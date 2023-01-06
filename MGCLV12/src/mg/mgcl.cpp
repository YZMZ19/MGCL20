/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#include "StdAfx.h"
#include "mg/tolerance.h"
#include "mgGL/color.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

// Implementation of non-categolized global functions.
//

//Version no. definition. This is used in MGIfstream/MGOfstream.cpp.
extern MG_DLL_DECLR const char* _MGCL_VER = "MGCL1203";
extern MG_DLL_DECLR const char* _MGCL_FILE = "File System fugen,Inc";
MGColors MGColor::m_colors=MGColors();

#ifndef _CONSOLE

#include "mgGL/OpenGLView.h"
#include "mgGL/GLSLProgram.h"
#include "mgGL/VBOElement.h"
#include "mgGL/VBOLeaf.h"
#include "mgGL/ConstructionPlane.h"

MGColor mgVBOElement::m_hilightColor;
GLfloat mgVBOElement::m_pointSize=6.f;//Default outer point size.
const GLfloat mgVBOLeaf::m_PolygonModePointSizeBase=1000.;
mgGLSLProgram* mgGLSLProgram::m_CurrrentGLSL=0;
MGDrawParam mgVBOElement::m_drawPara = MGDrawParam();
bool m_gdiplus_initialized = false;

#endif //_CONSOLE

// Compute radian angle from cosine and sine value.
//Function's return value is angle in radian, from zero to 2PAI.
double MGAngle(double ca	//Cosine value
			, double sa) { //Sine value
	double ang;
	if(ca>=0. && sa>=0.){		// 0<= angle <=HALFPAI.
		if(ca>=sa) ang=asin(sa);
		else       ang=acos(ca);
	}else if(ca<=0. && sa>=0.){ 	// HALFPAI<= angle <=PAI.
		if(sa>=-ca) ang=acos(ca);
		else        ang=mgPAI-asin(sa);
	}else if(ca<=0. && sa<=0.){	// PAI<= angle <=3*PAI/2..
		if(sa>=ca) ang=mgPAI-asin(sa);
		else       ang=mgDBLPAI-acos(ca);
	}else{						// 3*PAI/2.<= angle <=DBLPAI.
		if(ca<=-sa) ang=mgDBLPAI-acos(ca);
		else        ang=mgDBLPAI+asin(sa);
	}
	//normalize computing error.
	if(ang<0.)
		ang=0.;
	else if(ang>mgDBLPAI)
		ang=mgDBLPAI;
	return ang;
}

namespace MGCL{
	
//Get the MGCL_Version number.
const char* Version(){
	return _MGCL_VER;
};

//Get the MGCL File validity.
const char* File_validity(){
	return _MGCL_FILE;
}

//桁合わせ、デフォルトでは上から4桁を残し、後は切り捨てる
double MG_DLL_DECLR decimalAlign(double dValue, int nDigit){
	if(dValue == 0.)
		return 0.;

	int nValueExp = int(std::log10(std::abs(dValue)) + 1.);
	double pow = std::pow(10., nDigit - nValueExp);
	return int(dValue * pow) / pow;
}

//Round the input angel degree value degree to the multiples of step.
//The input degree and step are assumed to be angle in degree.
//step must be greater than 1.
//degree must be greater or equal to 0. and less than 360.
//Returned is a positibe value.
double MG_DLL_DECLR round_angle(double degree, double step){
	assert(step>1. && degree<=360.);
	const double BASE=360.;

	int m=int(degree/step);
	double rounded1=double(m)*step;
	double rounded2=double(++m)*step;
	if(rounded2-degree >= degree-rounded1)
		return rounded1;
	
	if(MGAEqual(rounded2, BASE))
		return 0.;
	return rounded2;
}

//Read in integer_string into intData.
//Function's return value is
//  true: when value specified.
//  false:when value not specified, intData be 0.
bool get_integer(
	char pDelimeter,	//parameter delimeter
	std::istringstream& istrm,	//Input string stream that contains integer data.
		//The stream pointer will be advanced to the start position of the next item.
	int& intData	//output integer data that is converted from the istrm data.
){
	intData = 0;
	bool specified = true;

	istrm>>intData;
	if(istrm.rdstate()){
		specified = false;
		istrm.clear();
	}
	std::string dummy;
	std::getline(istrm, dummy, pDelimeter);
	return specified;
}
bool get_integer(
	char pDelimeter,	//parameter delimeter
	std::istringstream& istrm,	//Input string stream that contains integer data.
		//The stream pointer will be advanced to the start position of the next item.
	short& shortData	//output integer data that is converted from the istrm data.
){
	int intData;
	bool specified = get_integer(pDelimeter, istrm, intData);
	shortData = intData;
	return specified;
}

//Read in real_string into realData
//Function's return value is
//  true: when value specified.
//  false:when value not specified, realData be 0.
bool get_real(
	char pDelimeter,	//parameter delimeter
	std::istringstream& istrm,	//Input string stream that contains real data.
		//The stream pointer will be advanced to the start position of the next item.
	double& realData	//converted real data from istrm will be output.
){
	bool specified = true;;
	realData = 0.;
	istrm>>realData;
	if(istrm.rdstate()){
		istrm.clear();
		specified = false;
	} else{
		char ch;
		istrm.get(ch);
		if(ch=='D' || ch=='d'){
			int expo = 0;
			istrm>>expo;
			if(expo!=0){
				realData *= pow(10., expo);
			}
		} else if(ch==pDelimeter)
			return specified;
	}
	std::string dummy;
	std::getline(istrm, dummy, pDelimeter);
	return specified;
}
bool get_real(
	char pDelimeter,	//parameter delimeter
	std::istringstream& istrm,	//Input string stream that contains real data.
		//The stream pointer will be advanced to the start position of the next item.
	float& floatData	//converted real data from istrm will be output.
){
	double dData;
	bool specified = get_real(pDelimeter, istrm, dData);
	floatData = float(dData);
	return specified;
}

#ifndef _CONSOLE

void start_up(bool need_to_GdiStartUp) {
	if (need_to_GdiStartUp) {
		//-> GDI+ startup
		Gdiplus::GdiplusStartupInput gdiplusStartupInput;
		Gdiplus::GdiplusStartup(&m_gdiplusToken, &gdiplusStartupInput, NULL);
		m_gdiplus_initialized = true;
		//<- GDI+ startup
	}
	else {
		m_gdiplus_initialized = false;
	}
}

void shut_down() {
	if (m_gdiplus_initialized)
		Gdiplus::GdiplusShutdown(m_gdiplusToken);
}

#endif //_CONSOLE

}

