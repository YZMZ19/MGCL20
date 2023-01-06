/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#include "StdAfx.h"
#include "mg/FSurface.h"
#include "mg/isect.h"
#include "topo/CFisects.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//MGCFisects defines a vector of MGCFisect.
//The vector is implemeted using STL's vector.
//All the methods to handle the vector are available from the STL's vector class,
//and public member m_CFivector. Refer to STL vector class.
//MGCFisects is used to represent intersection lines of a shell with
//another shell, a face, or a surface.
//The behavior of MGCFisect is like an auto_ptr. Copy or assignment
//of MGCFisect means transfer of the ownership of all the included curves
//to copied or assigned MGCFisect and original MGCFisect does not have the
//ownership more. Users should be aware of this fact.
//

///Constructor of 1 MGCFisect.
MGCFisects::MGCFisects(
	const MGCurve* curve,MGCFisect&& cfi
):MGCSisects(curve){
	m_list.emplace_back(new MGCFisect(std::move(cfi)));
}

///Insert MGCFisect at the index position i.
void MGCFisects::insertAt(iterator i, MGCFisect&& isect){
	m_list.insert(i,std::unique_ptr<MGCFisect>(new MGCFisect(std::move(isect))));
}

/// Output virtual function.
std::ostream& MGCFisects::toString(std::ostream& ostrm)const{
	const MGCFisects& cfis=*this;
	size_t n=cfis.size();
	ostrm<<"MGCFisects::number of isect="<<n<<std::endl;
	MGCFisects::const_iterator i;
	int j=0;
	for(i=cfis.begin(); i!=cfis.end(); i++)
		ostrm<<j++<<":"<<(*i);
	return ostrm;
}

void MGCFisects::push_back(MGCFisect && isect){
	m_list.emplace_back(new MGCFisect(std::move(isect)));
}

void MGCFisects::emplace_back(const MGCSisect& is, const MGFSurface& f){
	m_list.emplace_back(new MGCFisect(is, f));
}
