/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#include "StdAfx.h"
#include "mg/Tolerance.h"
#include "mg/EReal.h"
#include "mg/LBRep.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

// MGTolerance.cc
// Implementation of MGTolerance.
//

MGTolerance::MGTolerance(): 
	//  �����o�f�[�^�̏�����	
	m_mach_zero(1.0E-20),// �������Ƃ݂Ȃ��Q�_�Ԃ̋���-----Machine Zero Version	
	m_wc_zero(0.5E-3),// �������Ƃ݂Ȃ��Q�_�Ԃ̋���-----Absolute Version	
	m_wc_zero_sqr(m_wc_zero*m_wc_zero),// m_wc_zero �̂Q��-----Absolute Version	
	m_rc_zero(1.0E-6),// �������Ƃ݂Ȃ��Q�_�Ԃ̋���-----Relative Version	
	m_rc_zero_sqr(m_rc_zero*m_rc_zero),// m_rc_zero �̂Q��-----Relative Version	
	m_angle_zero(.0025),// �Q���������Ƃ݂Ȃ��p�x(in radina).	
	m_line_zero(m_wc_zero),// �Q�Ȑ����������Ƃ݂Ȃ��g�������X	
	m_max_knot_ratio(5.0E+2),// �ׂ荇��Knot�̔�̍ő�l	
	m_count(1)// �X�^�b�N�J�E���^
{
	std::fill_n(m_mach_zero_stack, MG_MAX_TOL_STACK_SIZE, 0.);//m_mach_zero
	std::fill_n(m_wc_zero_stack, MG_MAX_TOL_STACK_SIZE, 0.);// m_wc_zero
	std::fill_n(m_rc_zero_stack, MG_MAX_TOL_STACK_SIZE, 0.);// m_rc_zero
	std::fill_n(m_angle_zero_stack, MG_MAX_TOL_STACK_SIZE, 0.);// m_angle_zero
	std::fill_n(m_line_zero_stack, MG_MAX_TOL_STACK_SIZE, 0.);// m_line_zero
	std::fill_n(m_max_knot_ratio_stack, MG_MAX_TOL_STACK_SIZE, 0.);

	m_mach_zero_stack[0] = m_mach_zero;
	m_wc_zero_stack[0] = m_wc_zero;
	m_rc_zero_stack[0] = m_rc_zero;
	m_angle_zero_stack[0] = m_angle_zero;
	m_line_zero_stack[0] = m_line_zero;
	m_max_knot_ratio_stack[0] = m_max_knot_ratio;
}

MGTolerance::~MGTolerance(){
	pop();
}

//
// �����o�֐�
//
MGTolerance& MGTolerance::instance(){
	static MGTolerance theInst;
	return theInst;
}

//  �X�V
//  m_mach_zero��ύX����B
double MGTolerance::set_mach_zero(double mach_zero){
	MGTolerance& t = instance();
	double save=t.m_mach_zero;
    t.m_mach_zero = mach_zero;
	if(mach_zero<=0.0) t.m_mach_zero=1.0E-20;
	return save;
}

//  m_wc_zero��ύX����B
double MGTolerance::set_wc_zero(double wc_zero){
	MGTolerance& t = instance();
	double save=t.m_wc_zero;
    t.m_wc_zero = wc_zero;
	if(t.m_wc_zero<t.m_mach_zero) t.m_wc_zero=t.m_mach_zero;
    t.m_wc_zero_sqr = t.m_wc_zero * t.m_wc_zero;
	return save;
}

//  m_rc_zero��ύX����B
double MGTolerance::set_rc_zero(double rc_zero){
	MGTolerance& t = instance();
 	double save=t.m_rc_zero;
    t.m_rc_zero = rc_zero;
	if(t.m_rc_zero<t.m_mach_zero) t.m_rc_zero=t.m_mach_zero;
    t.m_rc_zero_sqr = t.m_rc_zero * t.m_rc_zero;
	return save;
}

//  m_angle_zero��ύX����B
double MGTolerance::set_angle_zero(double angle_zero){
	MGTolerance& t = instance();
 	double save=t.m_angle_zero;
    t.m_angle_zero = angle_zero;
	if(t.m_angle_zero<t.m_mach_zero) t.m_angle_zero=t.m_mach_zero;
	return save;
}

//  m_line_zero��ύX����B
double MGTolerance::set_line_zero(double line_zero){
	MGTolerance& t = instance();
 	double save=t.m_line_zero;
    t.m_line_zero = line_zero;
	if(t.m_line_zero<t.m_mach_zero) t.m_line_zero=t.m_mach_zero;
	return save;
}

//  m_max_knot_ratio��ύX����B
double MGTolerance::set_max_knot_ratio(double max_knot_ratio){
	MGTolerance& t = instance();
	double save=t.m_max_knot_ratio;
	t.m_max_knot_ratio = max_knot_ratio;
	if(t.m_max_knot_ratio<10.) t.m_max_knot_ratio=10.;
	return save;
}

//  �X�^�b�N�� push ����B
void MGTolerance::push(){
	MGTolerance& t = instance();
	assert(t.m_count < MG_MAX_TOL_STACK_SIZE );//*****Stack overflow******
	if(t.m_count < MG_MAX_TOL_STACK_SIZE ) {
		t.m_mach_zero_stack[ t.m_count ] = t.m_mach_zero;
		t.m_wc_zero_stack[ t.m_count ] = t.m_wc_zero;
		t.m_rc_zero_stack[ t.m_count ] = t.m_rc_zero;
		t.m_angle_zero_stack[ t.m_count ] = t.m_angle_zero;
		t.m_line_zero_stack[ t.m_count ] = t.m_line_zero;
		t.m_max_knot_ratio_stack[ t.m_count ] = t.m_max_knot_ratio;
		t.m_count += 1;
	}
}

//  �X�^�b�N�� pop ����B
void MGTolerance::pop(){
	MGTolerance& t = instance();
    assert(t.m_count >0 );//*****Stack underflow******
	if( t.m_count > 0 ){
		int id=t.m_count-1;
		t.m_mach_zero = t.m_mach_zero_stack[ id ];
		t.m_wc_zero = t.m_wc_zero_stack[ id ];
		t.m_wc_zero_sqr = t.m_wc_zero * t.m_wc_zero;
		t.m_rc_zero = t.m_rc_zero_stack[ id ];
		t.m_rc_zero_sqr = t.m_rc_zero * t.m_rc_zero;
		t.m_angle_zero = t.m_angle_zero_stack[ id ];
		t.m_line_zero = t.m_line_zero_stack[ id ];
		t.m_max_knot_ratio = t.m_max_knot_ratio_stack[ id ];
		t.m_count = id;
    }
}

//Set world coordinate zero tolerance.
void mgTolSetRCZero::restore(){
	if(m_RCzeroSave>0.)
		MGTolerance::set_rc_zero(m_RCzeroSave);//restore the saved error.
	m_RCzeroSave = -1.;//Value to indicate the content is invalid.
}

//Update world coordinate zero tolerance.
void mgTolSetRCZero::update(double errorNew){
	if (m_RCzeroSave <= 0.)//If the content is invalid.
		m_RCzeroSave = MGTolerance::rc_zero();//save the current.

	MGTolerance::set_rc_zero(errorNew);
}

//Set world coordinate zero tolerance.
void mgTolSetWCZero::restore() {
	if (m_WCzeroSave > 0.)
		MGTolerance::set_wc_zero(m_WCzeroSave);//restore the saved error.
	m_WCzeroSave = -1.;//Value to indicate the content is invalid.
}

//Update world coordinate zero tolerance.
void mgTolSetWCZero::update(double errorNew) {
	if (m_WCzeroSave <= 0.)//If the content is invalid.
		m_WCzeroSave = MGTolerance::wc_zero();//save the current.

	MGTolerance::set_wc_zero(errorNew);
}

//Set world coordinate zero tolerance.
void mgTolSetLineZero::restore() {
	if (m_LinezeroSave > 0.)
		MGTolerance::set_line_zero(m_LinezeroSave);//restore the saved error.
	m_LinezeroSave = -1.;//Value to indicate the content is invalid.
}

//Update world coordinate zero tolerance.
void mgTolSetLineZero::update(double errorNew) {
	if (m_LinezeroSave <= 0.)//If the content is invalid.
		m_LinezeroSave = MGTolerance::line_zero();//save the current.

	MGTolerance::set_line_zero(errorNew);
}

//Set world coordinate zero tolerance.
mgTolSetWCLineZero::mgTolSetWCLineZero(double wczero, double linezero){
	m_WCzeroSave = MGTolerance::set_wc_zero(wczero);
	m_LineZeroSave = MGTolerance::set_line_zero(linezero);
}

//Set world coordinate zero tolerance.
void mgTolSetWCLineZero::restore() {
	if (m_WCzeroSave > 0.){
		MGTolerance::set_line_zero(m_LineZeroSave);
		MGTolerance::set_wc_zero(m_WCzeroSave);
	}
	m_WCzeroSave = -1.;
}

//Set world coordinate zero tolerance.
void mgTolSetAngleZero::restore() {
	if (m_angleZeroSave > 0.)
		MGTolerance::set_angle_zero(m_angleZeroSave);//restore the saved error.
	m_angleZeroSave = -1.;//Value to indicate the content is invalid.
}

//Update world coordinate zero tolerance.
void mgTolSetAngleZero::update(double errorNew) {
	if (m_angleZeroSave <= 0.)//If the content is invalid.
		m_angleZeroSave = MGTolerance::angle_zero();//save the current.

	MGTolerance::set_angle_zero(errorNew);
}

// Global Functions
         
//  �g�������X���l�����ė^����ꂽ�l���O�����ׂ�-----Machine Zero Version.
bool MGMZero(double data) {
     return (fabs(data) <= MGTolerance::mach_zero());
}
//  �g�������X���l�����ė^����ꂽ�Q�� double ����v���邩���ׂ�B
//  �������� true(non zero) ��ԋp����-----Absolute Version.
bool MGAEqual (double data1, double data2) {
      return fabs(data1-data2) <= MGTolerance::wc_zero();
}
         
//  �g�������X���l�����ė^����ꂽ�l���O�����ׂ�-----Absolute Version.
bool MGAZero(double data) {
     return (fabs(data) <= MGTolerance::wc_zero());
}

//Test if difference of two data is less than MGTolerance::rc_zero()
//after changing data1 and data2 proportionally for data1 or 2 to be 1.
bool MGREqual2(double data1, double data2) {
      return MGREqual(data1,data2);
}

//Test if difference of two data is less than MGTolerance::rc_zero()
//after changing data1 and data2 proportionally for data1 or 2 to be 1.
bool MGREqual(double data1, double data2){
	double d2e=data2*MGTolerance::rc_zero();
	if(d2e<0.) d2e=-d2e;
	double d2md1=data2-data1;
	if(d2md1>d2e) return 0;
	if(d2md1<-d2e) return 0;
	return 1;
}

//Test if difference of two data is equal.
//Comparison is:
//test if abs(data1-data2)/base_length is less than MGTolerance::rc_zero().
bool MGREqual_base(double data1, double data2, double base_length){
	return MGRZero2(data1-data2,base_length);
}
bool MGREqual_baseEReal(MGEReal data1, MGEReal data2, const MGEReal& base_length){
	if(data1.finite() && data2.finite()){
		double dif=data1.value()-data2.value();
		if(base_length.finite())
			return MGRZero2(dif,base_length.value());
		if(dif>=0.)
			return dif<=MGTolerance::wc_zero();
		return (-dif)<= MGTolerance::wc_zero();
	}
	if(data1.minus_infinite() && data2.minus_infinite()) 
		return true;
	if(data1.plus_infinite() && data2.plus_infinite())
		return true;
	return false;
}

//  �g�������X���l�����ė^����ꂽ�l���O�����ׂ�1-----Relative Version.
bool MGRZero(double data) {
     return (fabs(data) <= MGTolerance::rc_zero());
}

//�g�������X���l�����ė^����ꂽ�l���O�����ׂ�2-----Relative Version
//Test if data is less or equal to rc_zero() compared to base_length.
//Comparison is done after data and base_length are so changed
//that base_length is 1.
//If base_length is zero, MGRZero2 returns always false.
bool MGRZero2(double data, double base_length){
	double e=base_length*MGTolerance::rc_zero();
	if(e<0.) e=-e;
	if(data<-e) return 0;
	if(data>e) return 0;
	return 1;
}
bool MGRZero2(double data, const MGEReal& base_length){
	if(base_length.finite())
		return MGRZero2(data,base_length.value());
	if(data>=0.)
		return data<=MGTolerance::wc_zero();
    return (-data)<= MGTolerance::wc_zero();
}

//  �g�������X���l�����ė^����ꂽ�p�x(Radian)�����p(�΁^�Q)�����ׂ�B
bool MGRight_angle(double cos_data) {
	return (fabs(cos_data) <= MGTolerance::angle_zero());
}

//  �g�������X���l�����ė^����ꂽ�p�x(Radian)���O�����ׂ�B
bool MGZero_angle (double data){
      return (fabs(data) <= MGTolerance::angle_zero());
}

void MGPrintBSpline(int k, int n, const double* t, const double* rcoef, int irc, int ncd){
	MGBPointSeq rc(n,ncd);
	MGKnotVector kntv(k,n,t);
	for(int i=0; i<n; i++){
		for(int j=0; j<ncd; j++){
			rc(i,j)=rcoef[i+irc*j];
		}
	}
	MGLBRep lb; lb.buildLBRepFromMemberData(std::move(kntv), std::move(rc));
	std::cout<<lb<<std::endl;
}

//Define C interface function
double bzrzro_(){return MGTolerance::rc_zero();}
double bzamin_(){return MGTolerance::angle_zero();}
double bzmzro_(){return MGTolerance::mach_zero();}
double bkmax_(){return MGTolerance::max_knot_ratio();}
void bzprintBspl(int k, int n,const double *t, const double *rcoef, int irc, int ncd){
	MGPrintBSpline(k,n,t,rcoef,irc,ncd);
}
