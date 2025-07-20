/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#pragma once 

#ifndef _MGDIRECTIONAL_LIGHT_HH_
#define _MGDIRECTIONAL_LIGHT_HH_

#include "mgGL/Light.h"

class MGOfstream;
class MGIfstream;
class MGVector;

/** @addtogroup GLAttrib
 *  @{
 */

///MGDirectionalLight is a directional light source that approximates infinite light sources as the sun.

///MGDirectionalLight can improve rendering performance over other local light sources,
///such as MGPointLight and MGSpotLight. Use MGDirectionalLight to set the direction
///of general lighting for a scene.
class MG_DLL_DECLR MGDirectionalLight : public MGLight{

public:
	
///////////// Constructors & destructor /////////////

MGDirectionalLight();
MGDirectionalLight(
    float intensity,		///<applied to GL_DIFFUSE and GL_SPECULAR
    float ambientIntensity,	///<applied to GL_AMBIENT
    const float color[3],	///<applied to GL_AMBIENT, GL_DIFFUSE, GL_SPECULAR
	const MGVector& direction///<Light direction
);

///~MGDirectionalLight();

///Assignment
MGDirectionalLight& operator=(const MGGel& gel2);
MGDirectionalLight& operator=(const MGDirectionalLight& gel2);

///comparison
bool operator<(const MGDirectionalLight& gel2)const;
bool operator<(const MGGel& gel2)const;

///Generate a newed clone object.
MGDirectionalLight* clone()const;

///render GLAttribute process.
///Function's return value is the lightnumber of this light executed.
int exec()const;
	
///////// Set/Get  /////////////

///Set the direction data.
void setDirection(const MGVector& direction);
void setDirection(const float direction[3]){
	for(int i=0; i<3; i++) m_direction[i]=direction[i];
}
void setDirection(float v0, float v1, float v2){
	m_direction[0]=v0;
	m_direction[1]=v1;
	m_direction[2]=v2;
}

///Get the direction data.
void getDirection(MGVector& direction)const;
void getDirection(float direction[3])const{
	for(int i=0; i<3; i++) direction[i]=m_direction[i];
}
void getDirection(float& v0, float& v1, float& v2)const{
	v0=m_direction[0];
	v1=m_direction[1];
	v2=m_direction[2];
}

/// Return This object's typeID
long identify_type() const{return MGDIRECTIONAL_LIGHT_TID;};

///Get the name of the class.
std::string whoami()const{return "DirectionalLight";};

///Read all member data.
void ReadMembers(MGIfstream& buf);

///Write all member data
void WriteMembers(MGOfstream& buf)const;

/// Output function.
std::ostream& toString(std::ostream&) const;

///Transform the gel by the argument.

///translation
void transform(const MGVector& v);

///scaling.
void transform(double scale);

///matrix transformation.
void transform(const MGMatrix& mat);

///general transformation.
void transform(const MGTransf& tr);

private:

    float m_direction[3];	///<GL_SPOT_DIRECTION

};

/** @} */ // end of GLAttrib group
#endif // _MGDIRECTIONAL_LIGHT_HH_
