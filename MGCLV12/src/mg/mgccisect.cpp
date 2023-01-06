/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#include "StdAfx.h"
#include "mg/Tolerance.h"
#include "mg/Position.h"
#include "mg/CCisect.h"
#include "mg/CSisect.h"
#include "mg/SSisect.h"
//#include "topo/CFisect.h"
#include "topo/HHisect.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

// MGCCisect.cc
// MGCCisect の実装ファイル

//
// コンストラクタ 
//
// 初期化なしで交点を生成
MGCCisect::MGCCisect()
	:m_param1(0.), m_param2(0.), m_rel(MGCCRELATION::MGCCREL_UNKNOWN){;}

// 全てのコンポーネントを指定して交点を生成
MGCCisect::MGCCisect (
	const MGPosition & is, double t1, double t2,
	const MGCCRELATION rel) 
	: m_ipoint(is), m_param1(t1), m_param2(t2), m_rel(rel) {;}

//
// メンバ関数
//
bool MGCCisect::operator== (const MGCCisect& cci)const{
	return MGREqual(m_param1,cci.m_param1) && MGREqual(m_param2,cci.m_param2);
}

//Ordering functions.
bool MGCCisect::operator< (const MGisect& is)const{
	const MGCCisect* cci=dynamic_cast<const MGCCisect*>(&is);
	if(cci) return operator<(*cci);
	auto csi = dynamic_cast<const MGCSisect*>(&is);
	if(csi) return true;
	//auto cfi = dynamic_cast<const MGCFisect*>(&is);
	//if(cfi) return true;
	auto ssi = dynamic_cast<const MGSSisect*>(&is);
	if(ssi) return true;
	auto hhi = dynamic_cast<const MGHHisect*>(&is);
	if(hhi) return true;
	return true;
}

bool MGCCisect::operator== (const MGisect& is)const{
	const MGCCisect* cci =dynamic_cast<const MGCCisect*>(&is);
	if(!cci) return false;
	return operator==(*cci);
}

//Exchange 1st and 2nd order of the parameter line representation.
void MGCCisect::exchange12(){
	double param1=m_param1;
	m_param1=m_param2;
	m_param2=param1;
}

//Debug Function
// Output virtual function.
std::ostream& MGCCisect::toString(std::ostream& ostrm)const{
//	ostrm.setf ( ios::scientific, ios::floatfield );
//	ostrm.precision ( 10 );
	ostrm << "MGCCisect::m_ipoint="<<m_ipoint
		<<", m_param1="<<m_param1
		<<", m_param2="<<m_param2
		<<", m_rel="<< m_rel << std::endl;
	return ostrm;
}
