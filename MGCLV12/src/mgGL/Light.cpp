/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#include "StdAfx.h"
#include "mg/Position.h"
#include "mg/Transf.h"
#include "mg/Ofstream.h"
#include "mg/Ifstream.h"
#include "mgGL/DirectionalLight.h"
#include "mgGL/SpotLight.h"
#include "mgGL/glslprogram.h"
#include "mgGL/VBO.h"
#include "mgGL/Appearance.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

////////////////////////////////////////////////////////////////////

MGLight::MGLight()
:MGGLAttrib(0), m_intensity(1.), m_ambientIntensity(0.){
	for(int i=0; i<3; i++) m_color[i]=1.;
}
MGLight::MGLight(
    float intensity,		//applied to GL_DIFFUSE and GL_SPECULAR
    float ambientIntensity,	//applied to GL_AMBIENT
    const float color[3]	//applied to GL_AMBIENT, GL_DIFFUSE, GL_SPECULAR
):MGGLAttrib(0), m_intensity(intensity),
 m_ambientIntensity(ambientIntensity){
	for(int i=0; i<3; i++) m_color[i]=color[i];
}

//assignment
MGLight& MGLight::operator=(const MGLight& gel2){
	set_light(gel2);
	return *this;
}
MGLight& MGLight::operator=(const MGGel& gel2){
	const MGLight* gel2_is_this=dynamic_cast<const MGLight*>(&gel2);
	if(gel2_is_this)
		operator=(*gel2_is_this);
	return *this;
}

MGLight* MGLight::clone()const{
	return new MGLight(*this);
}

//assignment
MGLight& MGLight::set_light(const MGLight& gel2){
	MGGLAttrib::operator=(gel2);

	m_intensity=gel2.m_intensity;
	m_ambientIntensity=gel2.m_ambientIntensity;
	for(int i=0; i<3;i++)
		m_color[i]=gel2.m_color[i];
	return *this;
}

bool MGLight::operator<(const MGLight& gel2)const{
	if(m_intensity==gel2.m_intensity)
		if(m_ambientIntensity==gel2.m_ambientIntensity)
			return m_color[0]<gel2.m_color[0];
		else
			return m_ambientIntensity<gel2.m_ambientIntensity;
	else
		return m_intensity<gel2.m_intensity;
}
bool MGLight::operator<(const MGGel& gel2)const{
	const MGLight* gel2_is_this=dynamic_cast<const MGLight*>(&gel2);
	if(gel2_is_this)
		return operator<(*gel2_is_this);
	return false;
}
//
int MGLight::exec()const{
	GLenum lightNo=get_light_num();
	// ライトの最大数は10こまで。mgclShader.fragより。
	ASSERT(lightNo>=0 && lightNo < mgGLSLProgram::LIGHT_NUM);

	mgGLSLProgram* glsl=mgGLSLProgram::getCurrentGLSLProgram();

	glsl->setUniformLights(lightNo, mgGLSLProgram::isEnabled , enabled());

	if(enabled()){
		int i;
		glm::vec4 color = glm::vec4(0.0f, 0.0f, 0.0f, 1.0f);

		//Ambient Color. GL_AMBIENT
		for(i=0; i<3; i++){
			color[i]=m_color[i]*m_ambientIntensity;
		}
		glsl->setUniformLights(lightNo, mgGLSLProgram::ambientColor, color);

		//Light Color. GL_DIFFUSE and GL_SPECULAR
		for(i=0; i<3; i++){
			color[i]=m_color[i]*m_intensity;
		}
		glsl->setUniformLights(lightNo, mgGLSLProgram::diffuseColor, color);
		glsl->setUniformLights(lightNo, mgGLSLProgram::specularColor, color);
	}
	return lightNo;
}

void MGLight::ResetLight(GLint lightNo){
	ASSERT(lightNo >= 0 && lightNo<10);

	mgGLSLProgram* glsl=mgGLSLProgram::getCurrentGLSLProgram();

	glsl->setUniformLights(lightNo, mgGLSLProgram::isEnabled, false);
	glsl->setUniformLights(lightNo, mgGLSLProgram::ambientColor, glm::vec4(0.0f, 0.0f, 0.0f, 1.0f));
	glsl->setUniformLights(lightNo, mgGLSLProgram::diffuseColor, glm::vec4(0.0f, 0.0f, 0.0f, 1.0f));
	glsl->setUniformLights(lightNo, mgGLSLProgram::specularColor, glm::vec4(0.0f, 0.0f, 0.0f, 1.0f));
	glsl->setUniformLights(lightNo, mgGLSLProgram::position, glm::vec4(0.0f, 0.0f, 1.0f, 0.0f));
	glsl->setUniformLights(lightNo, mgGLSLProgram::spotDirection, glm::vec3(0.0f, 0.0f, -1.0f));
	glsl->setUniformLights(lightNo, mgGLSLProgram::spotExponent, 0.0f);
	glsl->setUniformLights(lightNo, mgGLSLProgram::spotCutoff , 180.0f);
	glsl->setUniformLights(lightNo, mgGLSLProgram::constantAttenuation, 1.0f);
	glsl->setUniformLights(lightNo, mgGLSLProgram::linearAttenuation, 0.0f);
	glsl->setUniformLights(lightNo, mgGLSLProgram::quadraticAttenuation, 0.0f);
}

//Get light number available next.
GLenum MGLight::get_light_num()const{
//	return GL_LIGHT0+m_lightNum;
	return (GLenum)m_lightNum;
}

///Turn off this light, returned is the old state(UNDEFINED/OFF/ON).
int MGLight::turn_off(){
	int old=m_flag;
	m_flag=OFF;
	return old;
}

///Turn on this light, returned is the old state(UNDEFINED/OFF/ON).
int MGLight::turn_on(){
	int old=m_flag;
	m_flag=ON;
	return old;
}

// Output function.
std::ostream& MGLight::toString(std::ostream& ostrm) const{
	ostrm<<"Lgiht=";
	if(m_flag==ON)
		ostrm<<"ON";
	else if(m_flag==OFF)
		ostrm<<"OFF";
	else{
		MGGLAttrib::toString(ostrm);
		return ostrm;
	}

	ostrm<<",intensity="<<m_intensity<<",ambientIntensity="<<m_ambientIntensity;
	ostrm<<std::endl;
	ostrm<<"color=["<<m_color[0]<<","<<m_color[1]<<","<<m_color[2]<<"]";
	return ostrm;
}

std::ostream& operator<< (std::ostream& ostrm, const MGLight& light){
	return light.toString(ostrm);
}

/////////// MGDirectionalLight ///////////
/////////// Constructors & destructor ///////////

MGDirectionalLight::MGDirectionalLight()
:MGLight(){
	m_direction[0]=0.;
	m_direction[1]=0.;
	m_direction[2]=1.;
}
MGDirectionalLight::MGDirectionalLight(
    float intensity,		//applied to GL_DIFFUSE and GL_SPECULAR
    float ambientIntensity,	//applied to GL_AMBIENT
    const float color[3],	//applied to GL_AMBIENT, GL_DIFFUSE, GL_SPECULAR
	const MGVector& direction
):MGLight(intensity,ambientIntensity,color){
	for(int i=0;i<3;i++) m_direction[i]=float(direction[i]);
}

//MGDirectionalLight::~MGDirectionalLight(){;}

MGDirectionalLight* MGDirectionalLight::clone()const{
	return new MGDirectionalLight(*this);
}

//assignment
MGDirectionalLight& MGDirectionalLight::operator=(const MGDirectionalLight& gel2){
	if(this==&gel2)
		return *this;

	MGLight::operator=(gel2);
	for(int i=0; i<3;)
		m_direction[i]=gel2.m_direction[i];
	return *this;
}
MGDirectionalLight& MGDirectionalLight::operator=(const MGGel& gel2){
	const MGDirectionalLight* gel2_is_this=dynamic_cast<const MGDirectionalLight*>(&gel2);
	if(gel2_is_this)
		operator=(*gel2_is_this);
	return *this;
}

bool MGDirectionalLight::operator<(const MGDirectionalLight& gel2)const{
	if(MGLight::operator<(gel2))
		return true;
	return m_direction[0]<gel2.m_direction[0];
}
bool MGDirectionalLight::operator<(const MGGel& gel2)const{
	const MGDirectionalLight* gel2_is_this=dynamic_cast<const MGDirectionalLight*>(&gel2);
	if(gel2_is_this)
		return operator<(*gel2_is_this);
	return false;
}

//exec GLAttribute process.
int MGDirectionalLight::exec()const{
	GLenum lightNo=MGLight::exec();

	mgGLSLProgram* glsl=mgGLSLProgram::getCurrentGLSLProgram();

	////Light Positionは無限遠点(w=0なので), 位置はdirectionを負に
	// GL_POSITION 
	glm::vec4 pos = glm::vec4(-m_direction[0], -m_direction[1],-m_direction[2], 0.0f);
	glsl->setUniformLights(lightNo, mgGLSLProgram::position , pos);

	//GL_SPOT_DIRECTION
	glm::vec3 dir = glm::vec3(m_direction[0], m_direction[1], m_direction[2]);
	glsl->setUniformLights(lightNo, mgGLSLProgram::spotDirection , dir);

	// GL_SPOT_CUTOFF 
	glsl->setUniformLights(lightNo, mgGLSLProgram::spotCutoff, 180.0f);
	
	return lightNo;
}

//////// Set/Get  ///////////

void MGDirectionalLight::setDirection(const MGVector& direction){
	for(int i=0; i<3; i++) m_direction[i]=float(direction[i]);
}

void MGDirectionalLight::getDirection(MGVector& direction)const{
	direction.resize(3);
	for(int i=0; i<3; i++) direction(i)=m_direction[i];
}

// Output function.
std::ostream& MGDirectionalLight::toString(std::ostream& ostrm) const{
	ostrm<<std::endl<<"DirectionalLight="; MGLight::toString(ostrm);
	ostrm<<",m_direction=["<<m_direction[0];
	for(int i=1; i<3; i++) ostrm<<","<<m_direction[i];
	ostrm<<"]";
	return ostrm;
}

//Transform the gel by the argument.
void MGDirectionalLight::transform(const MGVector& v)//translation
{
	for(int i=0; i<3; i++) m_direction[i]+=float(v[i]);
}
void MGDirectionalLight::transform(double scale)//scaling.
{
	for(int i=0; i<3; i++) m_direction[i]+=float(scale);
}
void MGDirectionalLight::transform(const MGMatrix& mat)//matrix transformation.
{
	int i;
	MGPosition V(3);for(i=0; i<3; i++) V(i)=m_direction[i];
	V*=mat;
	for(i=0; i<4; i++) m_direction[i]=float(V[i]);
}
void MGDirectionalLight::transform(const MGTransf& tr)//general transformation.
{
	int i;
	MGPosition V(4);for(i=0; i<4; i++) V(i)=m_direction[i];
	V*=tr;
	for(i=0; i<4; i++) m_direction[i]=float(V[i]);
}

/////////// MGPointLight ///////////
//////// Constructors  /////////////
MGPointLight::MGPointLight()
:MGLight(), m_radius(1.){
	m_location[0]=0.;
	m_location[1]=0.;
	m_location[2]=1.;
	m_location[3]=1.;

	m_attenuation[0]=1.;
	m_attenuation[1]=0.;
	m_attenuation[2]=0.;
}

MGPointLight::MGPointLight(
    float intensity,		//applied to GL_DIFFUSE and GL_SPECULAR
    float ambientIntensity,	//applied to GL_AMBIENT
    const float color[3],	//applied to GL_AMBIENT, GL_DIFFUSE, GL_SPECULAR
	const MGPosition& location,//location of the light.
	float radius,
	const float attenuation[3]//[0]=GL_CONSTANT_ATTENUATION,
							//[1]=GL_LINEAR_ATTENUATION
							//[2]=GL_QUADRATIC_ATTENUATION
):MGLight(intensity,ambientIntensity,color),m_radius(radius){
	for(int i=0;i<3;i++){
		m_location[i]=float(location[i]);
		m_attenuation[i]=float(attenuation[i]);
	}
	m_location[3]=1.;
}

//MGPointLight::~MGPointLight(){;}

//assignment
MGPointLight& MGPointLight::operator=(const MGPointLight& gel2){
	if(this==&gel2)
		return *this;

	MGLight::operator=(gel2);
	for(int i=0; i<4;)
		m_location[i]=gel2.m_location[i];
	m_radius=gel2.m_radius;
	for(int i=0; i<3;)
		m_attenuation[i]=gel2.m_attenuation[i];
	return *this;
}
MGPointLight& MGPointLight::operator=(const MGGel& gel2){
	const MGPointLight* gel2_is_this=dynamic_cast<const MGPointLight*>(&gel2);
	if(gel2_is_this)
		operator=(*gel2_is_this);
	return *this;
}

bool MGPointLight::operator<(const MGPointLight& gel2)const{
	if(MGLight::operator<(gel2))
		return true;
	return m_radius<gel2.m_radius;
}
bool MGPointLight::operator<(const MGGel& gel2)const{
	const MGPointLight* gel2_is_this=dynamic_cast<const MGPointLight*>(&gel2);
	if(gel2_is_this)
		return operator<(*gel2_is_this);
	return false;
}

MGPointLight* MGPointLight::clone()const{
	return new MGPointLight(*this);
}

//exec GLAttribute process.
int MGPointLight::exec()const{
	GLenum lightNo = MGLight::exec();

	mgGLSLProgram* glsl=mgGLSLProgram::getCurrentGLSLProgram();

	// Lightの方向を設定
	// GL_POSITION
	glm::vec4 pos = glm::vec4(m_location[0], m_location[1], m_location[2], m_location[3]);
	glsl->setUniformLights(lightNo, mgGLSLProgram::position , pos);

	// GL_SPOT_CUTOFF cos(180) = -1.0f;
	glsl->setUniformLights(lightNo, mgGLSLProgram::spotCutoff, 180.0f);

	// GL_CONSTANT_ATTENUATION
	glsl->setUniformLights(lightNo, mgGLSLProgram::constantAttenuation, m_attenuation[0]);

	// GL_LINEAR_ATTENUATION
	glsl->setUniformLights(lightNo, mgGLSLProgram::linearAttenuation, m_attenuation[1]);

	// GL_QUADRATIC_ATTENUATION
	glsl->setUniformLights(lightNo, mgGLSLProgram::quadraticAttenuation, m_attenuation[2]);

	/// TODO by Tomoko.
	/// PointLightのm_radiusはOpenGL的にどういう意味なのかがわからない。

	return lightNo;
}

// Output function.
std::ostream& MGPointLight::toString(std::ostream& ostrm) const{
	ostrm<<"PointLight=";
	MGLight::toString(ostrm);
	ostrm<<std::endl<<
		"m_location=["<<m_location[0]<<","<<m_location[1]<<","<<m_location[2]<<"]";
	ostrm<<", m_attenuation=["<<
		m_attenuation[0]<<","<<m_attenuation[1]<<","<<m_attenuation[2]<<"]";
	ostrm<<", m_radius="<<m_radius;
	return ostrm;
}

//Transform the gel by the argument.
void MGPointLight::transform(const MGVector& v)//translation
{
	for(int i=0; i<4; i++) m_location[i]+=float(v[i]);
}
void MGPointLight::transform(double scale)//scaling.
{
	for(int i=0; i<4; i++) m_location[i]+=float(scale);
	m_radius*=float(scale);
}
void MGPointLight::transform(const MGMatrix& mat)//matrix transformation.
{
	int i;
	MGPosition P(4);for(i=0; i<4; i++) P(i)=m_location[i];
	P*=mat;
	for(i=0; i<4; i++) m_location[i]=float(P[i]);
	m_radius*=float(mat.scale());
}
void MGPointLight::transform(const MGTransf& tr)//general transformation.
{
	int i;
	MGPosition P(4);for(i=0; i<4; i++) P(i)=m_location[i];
	P*=tr;
	for(i=0; i<4; i++) m_location[i]=float(P[i]);
	m_radius*=float(tr.affine().scale());
}

////////// Set/Get ///////////

void MGPointLight::setLocation(const MGPosition& location){
	for(int i=0; i<3; i++) m_location[i]=float(location[i]);
}
void MGPointLight::getLocation(MGPosition& location)const{
	location.resize(3);
	for(int i=0; i<3; i++) location(i)=m_location[i];
}

/////////// MGSpotLight ///////////

////////// Constructors //////////////

MGSpotLight::MGSpotLight()
:m_cutOffAngle(180.), m_exponent(0.){
	m_direction[0]=0.;
	m_direction[1]=0.;
	m_direction[2]=1.;
}

MGSpotLight::MGSpotLight(
    float intensity,		//applied to GL_DIFFUSE and GL_SPECULAR
    float ambientIntensity,	//applied to GL_AMBIENT
    const float color[3],	//applied to GL_AMBIENT, GL_DIFFUSE, GL_SPECULAR
	const MGPosition& location,
	float radius,
	const float attenuation[3],//[0]=GL_CONSTANT_ATTENUATION,
							//[1]=GL_LINEAR_ATTENUATION
							//[2]=GL_QUADRATIC_ATTENUATION
    const MGVector& direction,//GL_SPOT_DIRECTION
    float exponent,		//GL_SPOT_EXPONENT
    float cutOffAngle	//GL_SPOT_CUTOFF
):MGPointLight(intensity,ambientIntensity,color,location,radius,attenuation)
,m_exponent(exponent),m_cutOffAngle(cutOffAngle){
	for(int i=0; i<3; i++) m_direction[i]=float(direction[i]);
}

//MGSpotLight::~MGSpotLight(){;}
	
MGSpotLight* MGSpotLight::clone()const{
	return new MGSpotLight(*this);
}

///////////// Set/Get  ///////////

void MGSpotLight::setDirection(const MGVector& direction){
	for(int i=0; i<3; i++) m_direction[i]=float(direction[i]);
}
void MGSpotLight::getDirection(MGVector& direction){
	direction.resize(3);
	for(int i=0; i<3; i++) direction(i)=m_direction[i];
}

//assignment
MGSpotLight& MGSpotLight::operator=(const MGSpotLight& gel2){
	if(this==&gel2)
		return *this;

	MGLight::operator=(gel2);
	for(int i=0; i<3;)
		m_direction[i]=gel2.m_direction[i];
	m_exponent=gel2.m_exponent;
	m_cutOffAngle=gel2.m_cutOffAngle;
	return *this;
}
MGSpotLight& MGSpotLight::operator=(const MGGel& gel2){
	const MGSpotLight* gel2_is_this=dynamic_cast<const MGSpotLight*>(&gel2);
	if(gel2_is_this)
		operator=(*gel2_is_this);
	return *this;
}

bool MGSpotLight::operator<(const MGSpotLight& gel2)const{
	if(MGLight::operator<(gel2))
		return true;
	return m_exponent<gel2.m_exponent;
}
bool MGSpotLight::operator<(const MGGel& gel2)const{
	const MGSpotLight* gel2_is_this=dynamic_cast<const MGSpotLight*>(&gel2);
	if(gel2_is_this)
		return operator<(*gel2_is_this);
	return false;
}

//exec GLAttribute process.
int MGSpotLight::exec()const{
	GLenum lightNo=MGPointLight::exec();

	mgGLSLProgram* glsl=mgGLSLProgram::getCurrentGLSLProgram();

	// GL_SPOT_DIRECTION
	glm::vec3 dir = glm::vec3(m_direction[0],m_direction[1],m_direction[2]);

	glsl->setUniformLights(lightNo, mgGLSLProgram::spotDirection, dir );
	glsl->setUniformLights(lightNo, mgGLSLProgram::spotCutoff, m_cutOffAngle);
	glsl->setUniformLights(lightNo, mgGLSLProgram::spotExponent, m_exponent);

	return lightNo;
}

// Output function.
std::ostream& MGSpotLight::toString(std::ostream& ostrm) const{
	ostrm<<"SpotLight=";
	MGPointLight::toString(ostrm);
	ostrm<<std::endl<<
		"m_direction=["<<m_direction[0]<<","<<m_direction[1]<<","<<m_direction[2]<<"]";
	ostrm<<", m_exponent="<<m_exponent<<", m_cutOffAngle="<<m_cutOffAngle;
	return ostrm;
}

//Transform the gel by the argument.
void MGSpotLight::transform(const MGVector& v)//translation
{
	MGPointLight::transform(v);
	for(int i=0; i<3; i++) m_direction[i]+=float(v[i]);
}
void MGSpotLight::transform(double scale)//scaling.
{
	MGPointLight::transform(scale);
	for(int i=0; i<3; i++) m_direction[i]+=float(scale);
}
void MGSpotLight::transform(const MGMatrix& mat)//matrix transformation.
{
	MGPointLight::transform(mat);
	int i;
	MGPosition V(3);for(i=0; i<3; i++) V(i)=m_direction[i];
	V*=mat;
	for(i=0; i<4; i++) m_direction[i]=float(V[i]);
}
void MGSpotLight::transform(const MGTransf& tr)//general transformation.
{
	MGPointLight::transform(tr);
	int i;
	MGPosition V(4);for(i=0; i<4; i++) V(i)=m_direction[i];
	V*=tr;
	for(i=0; i<4; i++) m_direction[i]=float(V[i]);
}

// Serialization fucntion.
void MGLight::WriteMembers(MGOfstream& buf)const{
	MGGLAttrib::WriteMembers(buf);
	buf<<m_lightNum;
	buf<<m_intensity;
	buf<<m_ambientIntensity;
	for(int i=0; i<3; i++) buf<<m_color[i];
}
void MGLight::ReadMembers(MGIfstream& buf){
	MGGLAttrib::ReadMembers(buf);
	buf>>m_lightNum;
	buf>>m_intensity;
	buf>>m_ambientIntensity;
	for(int i=0; i<3; i++) buf>>m_color[i];
}

void MGDirectionalLight::WriteMembers(MGOfstream& buf)const{
	MGLight::WriteMembers(buf);
	for(int i=0; i<3; i++) buf<<m_direction[i];
}
void MGDirectionalLight::ReadMembers(MGIfstream& buf){
	MGLight::ReadMembers(buf);
	for(int i=0; i<3; i++) buf>>m_direction[i];
}

void MGPointLight::WriteMembers(MGOfstream& buf)const{
	MGLight::WriteMembers(buf);
	int i;
	for(i=0; i<4; i++) buf<<m_location[i];
	buf<<m_radius;
	for(i=0; i<3; i++) buf<<m_attenuation[i];
}
void MGPointLight::ReadMembers(MGIfstream& buf){
	MGLight::ReadMembers(buf);
	int i;
	for(i=0; i<4; i++) buf>>m_location[i];
	buf>>m_radius;
	for(i=0; i<3; i++) buf>>m_attenuation[i];
}

void MGSpotLight::WriteMembers(MGOfstream& buf)const{
	MGPointLight::WriteMembers(buf);
	for(int i=0; i<3; i++) buf<<m_direction[i];
	buf<<m_exponent;
	buf<<m_cutOffAngle;
}
void MGSpotLight::ReadMembers(MGIfstream& buf){
	MGPointLight::ReadMembers(buf);
	for(int i=0; i<3; i++) buf>>m_direction[i];
	buf>>m_exponent;
	buf>>m_cutOffAngle;
}
