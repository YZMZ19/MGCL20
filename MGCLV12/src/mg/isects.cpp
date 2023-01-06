/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#include "StdAfx.h"
#include "mg/isects.h"
#include "mg/CCisects.h"
#include "mg/CSisects.h"
#include "mg/SSisects.h"
#include "mg/FSurface.h"
#include "topo/HHisects.h"
#include "topo/Face.h"
#include "topo/Shell.h"
//#include "topo/FFisect.h"

//MGisects is used to represent an array of intersection lines of 
//two objects.
//The behavior of MGisects is like an auto_ptr. Copy or assignment
//of MGisects means transfer of the ownership of all the included MGisect
//to copied or assigned MGisects and original MGisects does not have the
//ownership any more. Users should be aware of this fact.
//

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

////////// Constructor //////////

//Void constructor(of size 0)
MGisects::MGisects(
	const MGObject* obj1,
	const MGObject* obj2
):m_object1(obj1), m_object2(obj2){
	if(obj1 && obj2)
		assert(obj1->manifold_dimension()<=obj2->manifold_dimension());
}

//Debug Function
std::ostream& operator << (std::ostream& ostrm, const MGisects& is){
	is.toString(ostrm);
	return ostrm;
}

std::ostream & MGisects::toString(std::ostream & ostrm) const{
	const MGisects& is=*this;
	ostrm<<"MGisects::"<<&is<<", m_object1="<<is.m_object1<<", m_object2="<<is.m_object2;
	size_t n=is.size();
	ostrm<<", number of isect="<<n<<std::endl;
	MGisects::const_iterator itr; int i=0;
	for(itr=is.begin(); itr!=is.end(); itr++) ostrm<<i++<<":"<<(*itr);
	ostrm<<std::endl;
	return ostrm;
}
////////// Member Function. //////////

//Replace first and second order of MGisect.
void MGisects::exchange12(){
	bool change=false;
	if(m_object2){
		int m1=m_object1->manifold_dimension(), m2=m_object2->manifold_dimension();
		if(m1==m2){
			const MGObject* objsave=m_object1;
			m_object1=m_object2;
			m_object2=objsave;
			change=true;
		}
	}
	if(size()&&change){
		iterator i=begin(), ie=end();
		for(; i!=ie; i++) (**i).exchange12();
	}
}

//append all the member of isects to the end of the list.
//Transfers the ownership of all the isect in isects to this list.
void MGisects::push_back(MGisects&& isects){
	m_list.splice(m_list.end(),std::move(isects.m_list));
}

///append all the member of isects to the fron of the list.
///Transfers the ownership of the isect in isects to this list.
void MGisects::push_front(MGisects&& isects){
	m_list.splice(m_list.begin(), std::move(isects.m_list));
}

