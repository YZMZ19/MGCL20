/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#pragma once 

#ifndef _MGLIGHT_HH_
#define _MGLIGHT_HH_

#include "mgGL/GLAttrib.h"

class MGOfstream;
class MGIfstream;
class MGMatrix;
class MGSpotLight;
class MGDirectionalLight;
class MGPointLight;

/** @file */
/** @addtogroup GLAttrib
 *  @{
 */

/// MGLight is an abstract base class for light sources.

///MGDirectionalLight,MGPointLight, or MGSpotLight are the inheritted classes.
class MG_DLL_DECLR MGLight:public MGGLAttrib{

public:
	
///Define light mode.
enum LIGHT_MODE{
	UNDEFINED= mgGLMode::UNDEFINED,
	OFF= mgGLMode::DISABLED,///Disabled is used as OFF.
	ON=1
};

///Write out the light data to ostream.
MG_DLL_DECLR friend std::ostream& operator<< (std::ostream&, const MGLight&);

//////////Constructor///////////////

MGLight();
MGLight(
    float intensity,		///<applied to GL_DIFFUSE and GL_SPECULAR
    float ambientIntensity,	///<applied to GL_AMBIENT
    const float color[3]	///<applied to GL_AMBIENT, GL_DIFFUSE, GL_SPECULAR
);

public:

///Assignment

///Assignment
virtual MGLight& operator=(const MGGel& gel2);
virtual MGLight& operator=(const MGLight& gel2);

///comparison
virtual bool operator<(const MGLight& gel2)const;
virtual bool operator<(const MGGel& gel2)const;

///Generate a newed clone object.
virtual MGLight* clone()const;

///render process.
///Function's return value is the lightnumber of this light executed.
virtual int exec()const;

///Reset the lightNo light to default valur(using setUniformLights).
static void ResetLight(GLint lightNo);

///draw GLAttribute process.
virtual void drawAttrib(
	mgVBO& vbo,///<The target graphic object.
	bool no_color=false	///<if true, color attribute will be neglected.
)const{exec();};

///Obtain the light number of this.
GLenum get_light_num()const;

///Test if this light is on, or off.
///True: light is on, false light is off or undefined.
bool light_is_on()const{return m_flag==ON;};

///render GLAttribute process.
virtual void render(mgVBO& vbo)const{exec();};

///Turn on the appropriate mask bit for this attribute. See glPushAttrib().
virtual void set_draw_attrib_mask(unsigned int& mask)const{
	set_Amask(mask,LIGHTING_BIT);
};

///Set light number.
void set_light_number(int lnum){m_lightNum=lnum;};

///Turn off the appropriate mask bit for this attribute. See glPushAttrib().
virtual void reset_draw_attrib_mask(unsigned int& mask)const{
	reset_Amask(mask,LIGHTING_BIT);
};

///Turn on the appropriate mask bit for this attribute. See glPushAttrib().
virtual void set_render_attrib_mask(unsigned int& mask)const{set_Amask(mask,LIGHTING_BIT);};

///Turn off the appropriate mask bit for this attribute. See glPushAttrib().
virtual void reset_render_attrib_mask(unsigned int& mask)const{reset_Amask(mask,LIGHTING_BIT);};

//////////// Set/Get  /////////////

///Set intensity.
void setIntensity(float intensity){m_intensity=intensity;};
///Get intensity.
float getIntensity()const{return m_intensity;};

///Set ambient intensity.
void setAmbientIntensity(float ambientIntensity){
	m_ambientIntensity=ambientIntensity;
};
///Get ambient intensity.
float getAmbientIntensity()const{return m_ambientIntensity;};

///Set color.
void setColor(const float color[3]){
	for(int i=0;i<3; i++) m_color[i]=color[i];
};
void getColor(float color[3]){
	for(int i=0;i<3; i++) color[i]=m_color[i];
};
void setColor(float v0, float v1, float v2){
	m_color[0]=v0;
	m_color[1]=v1;
	m_color[2]=v2;
}
void getColor(float& v0, float& v1, float& v2){
	v0=m_color[0];
	v1=m_color[1];
	v2=m_color[2];
}

///Get the name of the class.
virtual std::string whoami()const{return "Light";};

///Turn on this light, returned is the old state(UNDEFINED/OFF/ON).
int turn_off();

///Turn on this light, returned is the old state(UNDEFINED/OFF/ON).
int turn_on();

///Read all member data.
virtual void ReadMembers(MGIfstream& buf);
///Write all member data
virtual void WriteMembers(MGOfstream& buf)const;

/// Output virtual function.
virtual std::ostream& toString(std::ostream&) const;

protected:
	unsigned int m_lightNum;///Color number of this light.
    float m_intensity;	///<applied to GL_DIFFUSE and GL_SPECULAR
    float m_ambientIntensity;///<applied to GL_AMBIENT
    float m_color[3];	///<applied to GL_AMBIENT, GL_DIFFUSE, GL_SPECULAR

///assignment
MGLight& set_light(const MGLight& gel2);

};

/** @} */ // end of GLAttrib group
#endif // _MGLIGHT_HH_
