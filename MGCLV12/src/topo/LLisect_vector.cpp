/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#include "StdAfx.h"
#include "topo/Cell.h"
#include "topo/Loop.h"
#include "topo/LLisect_vector.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

// MGLLisect_vector defines a vector of MGLLisect.
// Used to represent Intersection points of two loops.

// Constructor
MGLLisect_vector::MGLLisect_vector():m_error_square(-1.){;}

MGLLisect_vector::MGLLisect_vector(const MGLoop& loop)		//Loop
{
	double err=loop.error();
	m_error_square=err*err*9.;
}

// 交点の全てのコンポーネントを指定して，交点リストに追加
//Add one intersection point to the list.
void MGLLisect_vector::append(const MGLLisect& lli){
// Adds the MGLLisect to the end of the list.
	LLiterator itr=begin(), itrend=end();
	for(; itr!=itrend; itr++){
		if((*itr).distance_square(lli)<=m_error_square) return;
	}

	push_back(lli);
}

void MGLLisect_vector::append(
	const MGPosition& uv,	//Parameter (u,v) of the parent face.
	const MGLEPoint& lp1,	//First loop's point data.
	const MGLEPoint& lp2)	//Second loop's point data.
{
	append(MGLLisect(uv,lp1,lp2));
}

// Adds the MGLLisect_vector to the end of the list.
void MGLLisect_vector::append(const MGLLisect_vector& list){
// Adds the MGLLisect_vector to the end of the list.
	const_LLiterator i;
	for(i=list.begin(); i!=list.end(); i++) append(*i);
}

//Debug Function
std::ostream& operator<< (std::ostream& out, const MGLLisect_vector& vec){
	out<<"MGLLisect_vector::"<<" ,m_error_square="<<vec.m_error_square;
	int n=vec.entries();
	out<<", number of isect="<<n<<std::endl;
	MGLLisect_vector::const_LLiterator itr; int i=0;
	for(itr=vec.begin(); itr!=vec.end(); itr++)
		out<<i++<<":"<<(*itr);
	return out;
}
