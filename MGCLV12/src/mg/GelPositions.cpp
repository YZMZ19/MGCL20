/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/

#include "StdAfx.h"
#include "mg/GelPositions.h"
#include "mg/PickObjects.h"
#include "mg/Gel.h"
#include "mg/isects.h"
#include "mg/Curve.h"
#include "mg/FSurface.h"
#include "topo/Face.h"
#include "topo/Shell.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//
//Implements MGGelPositions Class.

//////////Constructor//////////

///Construct MGGelPositions of a MGGelPosition.
///When gelp.is_null(), the gelp will not be set.
MGGelPositions::MGGelPositions(const MGGelPosition& gelp)
:m_vec(1,gelp){
	if(gelp.is_null())
		pop_back();
}

//Conversion constructor of MGPickObjects.
MGGelPositions::MGGelPositions(const MGPickObjects& gelp)
:m_vec(gelp.size()) {
	size_t n = gelp.size();
	for(size_t i=0; i<n; i++)
		m_vec[i] = gelp[i];
}

///////////////operator overloaded//////////////

//////////Member Function//////////

// Output virtual function.
std::ostream& operator<<(std::ostream& ostrm, const MGGelPositions& gelps){
	ostrm<<"MGGelPositions::number of gelps="<<gelps.size()<<std::endl;
	for(int j=0; auto& gelp:gelps)
		ostrm<<"gelp-"<<j++<<":"<<gelp<<std::endl;
	return ostrm;
}

// erase element MGGelPosition gelp. Function's return value is the following iterator
// of the erased element x.
void MGGelPositions::remove(const MGGelPosition& gelp){
	iterator i=find(gelp), ie=end();
	if(i!=ie)
		erase(i);
}

//Test if there is a MGPickObject that includes input objin
//in this MGPickObjects' member. If input objin is MGShell,
//and a member is_shell_face(), test is performed to the shell.
//Returns true if objin is included in this MGPickObjects.
MGGelPositions::iterator MGGelPositions::includes(const MGObject* objin){
	iterator i=begin(), ie=end();
	for(; i!=ie; i++){
		const MGGelPosition& gelpi=*i;
		const MGObject* obji=gelpi.leaf_object();
		if(obji == objin)
			return i;

		if(gelpi.is_shell_face()){
			obji=gelpi.top_object();
			if(obji==objin)
				return i;
		}
	}
	return ie;
}

//Remove objects of type from this.
void MGGelPositions::remove(const MGAbstractGels& types){
	MYVEC& gels=m_vec;
	int n=(int)size();
	for(int i=n-1; i>=0; i--){
		if(gels[i].is_type(types))
			gels.erase(begin()+i);
	}
}

//Remove gelps from this.
void MGGelPositions::remove(const MGGelPositions& gelps){
	size_t ngelp=gelps.size();
	for(size_t i=0; i<ngelp; i++){
		const MGGelPosition& po=gelps[i];
		iterator hit = find(po);
		if(hit != end()){//if po was already in current object, erase it.
			erase(hit);	
		}
	}
}

void MGGelPositions::reset(int i, const MGGelPosition& pobj){
	(*this)[i]=pobj;
}

//replace this with the common objects of this and pobjs2.
void MGGelPositions::reset_with_common(const MGGelPositions& pobjs2){
	const_iterator je=pobjs2.end();
	reverse_iterator ri, riend=rend();
	for(ri=rbegin(); ri!=riend; ri--){
		if(pobjs2.find(*ri)==je){
			auto i=ri;
			erase(++i.base());
		}
	}
}

//replace this with symmetric_differecne of this and pobj, that is;
//(1) remove the same MGGelPosition from this and pobjss.
//(2) append the result pobjs2 to this.
//On return, pobjs2 will have null sequence.
void MGGelPositions::reset_with_symmetric_difference(MGGelPositions& pobjs2){
	int n1=(int)size(), n2=(int)pobjs2.size();
	MGGelPositions* p1=this;
	MGGelPositions* p2=&pobjs2;
	int n=n1;
	if(n1>n2){
		n=n2;
		p1=&pobjs2; p2=this;
	}
	for(int i=n-1; i>=0; i--){
		iterator j=p2->find((*p1)[i]);
		if(j!=p2->end()){
			p2->erase(j);
			p1->erase(p1->begin()+i);
		}
	}
	append(pobjs2);
}

//Select objects of input type from this.
//Function's return value is MGGelPositions selected.
//This will be unchanged.
MGGelPositions MGGelPositions::select(const MGAbstractGels& types)const{
	MGGelPositions pobjs2;
	const_iterator i=begin(), ie=end();
	for(; i!=ie; i++){
		if((*i).is_type(types))
			pobjs2.push_back(*i);
	}
	return pobjs2;
}

//Select the 1st MGCurve from this.
//Function's return value is MGGelPosition of MGCurve 1st encountered in this
//MGGelPosition sequence. If this did not includes any MGCurve,
//null MGGelPosition will be returned.
//This will be unchanged.
MGGelPosition MGGelPositions::select_1st_curve()const{
	const_iterator i=begin(), ie=end();
	for(; i!=ie; i++){
		const MGObject* obj=(*i).leaf_object();
		const MGCurve* curve=dynamic_cast<const MGCurve*>(obj);
		if(curve)
			return *i;
	}
	return MGGelPosition();
}

//Select all the MGCurve from this.
//MGGelPosition of MGCurve encountered in this MGGelPosition sequence will be appended
//in curves.
//This will be unchanged.
void MGGelPositions::select_curves(MGGelPositions& curves)const{
	const_iterator i=begin(), ie=end();
	for(; i!=ie; i++){
		const MGObject* obj=(*i).leaf_object();
		const MGCurve* curve=dynamic_cast<const MGCurve*>(obj);
		if(curve){
			curves.push_back(*i);
		}
	}
}

//Select the 1st MGFSurface from this.
//Function's return value is MGGelPosition of MGFSurface 1st encountered in this
//MGGelPosition sequence. If this did not includes any MGFSurface,
//null MGGelPosition will be returned.
//This will be unchanged.
MGGelPosition MGGelPositions::select_1st_fsurface()const{
	const_iterator i=begin(), ie=end();
	for(; i!=ie; i++){
		const MGObject* obj=(*i).leaf_object();
		const MGFSurface* f=dynamic_cast<const MGFSurface*>(obj);
		if(f)
			return *i;
	}
	return MGGelPosition();
}

//Select all the MGFSurface from this.
//MGGelPositions of MGFSurface encountered in this MGGelPosition sequence will be appended
//in surfaces.
//This will be unchanged.
void MGGelPositions::select_fsurfaces(MGGelPositions& surfaces)const{
	const_iterator i=begin(), ie=end();
	for(; i!=ie; i++){
		const MGObject* obj=(*i).leaf_object();
		const MGFSurface* surface=dynamic_cast<const MGFSurface*>(obj);
		if(surface){
			surfaces.push_back(*i);
		}
	}
}

//Test if this is symmetric to gels2.
//Symmetric means:
//(1) number of gels included is the same.
//(2) all of the gels are MGObject and they have the same manifold dimension.
bool MGGelPositions::symmetric(
	 const MGGelPositions& gels2
)const{
	size_t n=size();
	if(gels2.size()!=n)
		return false;

	const MYVEC& gels=m_vec;
	for(size_t i=0; i<n; i++){
		if(gels[i].symmetric(gels2[i]))
			continue;
		return false;
	}
	return true;
}
