/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#pragma once 

#ifndef _MGLIGHTS_HH_
#define _MGLIGHTS_HH_

#include "mgGL/GLAttrib.h"

class MGOfstream;
class MGIfstream;
class MGLight;

/** @addtogroup GLAttrib
 *  @{
 */

///a container class for light sources(MGDirectionalLight, MGPointLight, or MGSpotLight).

///MGLights contains the light data as std::vector(a vector of newed MGLight's).
class MG_DLL_DECLR MGLights:public MGGLAttrib{
public:
	using MYVEC = std::vector<UniqueLight>;
	typedef MYVEC::iterator iterator ;
	typedef MYVEC::const_iterator const_iterator;
	typedef MYVEC::reverse_iterator reverse_iterator ;
	typedef MYVEC::const_reverse_iterator const_reverse_iterator;

private:
	MYVEC m_lights;///<vector of MGLight.

public:

////////Special member functions/////////
~MGLights()=default;//Destructor.
MGLights(const MGLights& rhs);//Copy constructor.
MGLights& operator=(const MGLights& gel2);//Copy assignment.
MGLights(MGLights&& rhs)=default;//Move constructor.
MGLights& operator=(MGLights&& rhs)=default;//Move assignment.

MGLights();

///Construct MGLights of one light.
MGLights(
	MGLight* light
):MGGLAttrib(1),m_lights(1){
	m_lights[0].reset(light);
};

///Assignment
MGLights& operator=(const MGGel& gel2);

///comparison
bool operator<(const MGLights& gel2)const;
bool operator<(const MGGel& gel2)const;

const UniqueLight& back()const{return m_lights.back();};
UniqueLight& back(){return m_lights.back();};
iterator begin(){return m_lights.begin();};
const_iterator begin()const{return m_lights.begin();};
void clear(){m_lights.clear();};
bool empty()const{return m_lights.empty();};
iterator end(){return m_lights.end();};
const_iterator end()const{return m_lights.end();};
const UniqueLight& front()const{return m_lights.front();};
UniqueLight& front(){return m_lights.front();};
const UniqueLight& operator[](int i)const{return m_lights[i];};
UniqueLight& operator[](int i){return m_lights[i];};
void pop_back(){m_lights.pop_back();};

///Generate a newed clone object.
virtual MGLights* clone()const;

///render GLAttribute process.
virtual void exec()const;

///draw GLAttribute process.
void drawAttrib(
	mgVBO& vbo,///<The target graphic object.
	bool no_color=false	///<<if true, color attribute will be neglected.
)const{exec();};

///add one light.
///Function's return value is the numbe of lights defined.
int push_back(MGLight* light);

///render GLAttribute process.
void render(mgVBO& vbo)const{exec();};

///Turn on the appropriate mask bit for this attribute. See glPushAttrib().
void set_render_attrib_mask(unsigned int& mask)const{set_Amask(mask,LIGHTING_BIT);};

///Turn off the appropriate mask bit for this attribute. See glPushAttrib().
void reset_render_attrib_mask(unsigned int& mask)const{reset_Amask(mask,LIGHTING_BIT);};

///Turn on the appropriate mask bit for this attribute. See glPushAttrib().
void set_draw_attrib_mask(unsigned int& mask)const{set_Amask(mask,LIGHTING_BIT);};

///Obtain the light number defined.
int size()const{return (int)m_lights.size();};

///Turn off the appropriate mask bit for this attribute. See glPushAttrib().
void reset_draw_attrib_mask(unsigned int& mask)const{;};

/// Return This object's typeID
virtual long identify_type() const{return MGLIGHTS_TID;};

///Get the name of the class.
std::string whoami()const{return "Lights";};

///Read all member data.
virtual void ReadMembers(MGIfstream& buf);

///Write all member data
virtual void WriteMembers(MGOfstream& buf)const;

/// Output virtual function.
virtual std::ostream& toString(std::ostream&) const;


};

/** @} */ // end of GLAttrib group
#endif // _MGLIGHTS_HH_
