/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#include "StdAfx.h"

#include "mg/Curve.h"
#include "mg/FSurface.h"
#include "mg/PickObject.h"
#include "mg/PickObjects.h"
#include "mg/GelPositions.h"
#include "topo/Shell.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

////////////////////////////////////////////////////////////////////
//a container class for MGPickObject.

//Construct MGPickObjects of one pobj.
MGPickObjects::MGPickObjects(
	const MGPickObject& pobj
){
	m_PickObjects.emplace_back(pobj.clone());
}

//Copy constructor
MGPickObjects::MGPickObjects(const MGPickObjects& pobjs){
	size_t n=pobjs.size();
	for(size_t i=0; i<n; i++)
		push_back(pobjs[i]);
}

MGPickObjects& MGPickObjects::operator+=(const MGPickObject& gelp){
	m_PickObjects.emplace_back(gelp.clone());
	return *this;
}

//Operator overload.
MGPickObjects& MGPickObjects::operator=(const MGPickObjects& pobjs){
	m_PickObjects.clear();
	push_back(pobjs);
	return *this;
}

//append the current objects(MGGelPositions).
void MGPickObjects::append_object(const MGGelPositions& gelps){
	size_t n=gelps.size();
	for(size_t i=0; i<n; i++)
		m_PickObjects.emplace_back(new MGPickObject(gelps[i]));
}

//Replace this sequence with [first,last).
void MGPickObjects::assign(const_iterator first, const_iterator last){
	std::vector<UniquePickObject> tempobjs;
	for(const_iterator i=first; i!=last; i++)
		tempobjs.emplace_back((**i).clone());
	m_PickObjects=std::move(tempobjs);
}

// erase sequence [first, last).
void MGPickObjects::erase(iterator first, iterator last){
	m_PickObjects.erase(first,last);
}
MGPickObjects::iterator MGPickObjects::erase(iterator i){
	return m_PickObjects.erase(i);
}

// erase after the elments after the front().
//Resutl has length 1 sequence.
void MGPickObjects::erase_except_front(){
	if(size()<=1) return;
	iterator first=begin()+1, last=end();
	erase(first,last);
}

//find the same pobj in this objects.
MGPickObjects::iterator MGPickObjects::find(
	const MGPickObject& pobj
){
	iterator i=begin(), ie=end();
	for(; i!=ie; i++){
		if((**i)==pobj)
			return i;
	}
	return ie;
}

//find the same pobj in this objects.
MGPickObjects::const_iterator MGPickObjects::find(
	const MGPickObject& pobj
)const{
	const_iterator i=begin(), ie=end();
	for(; i!=ie; i++){
		if((**i)==pobj)
			return i;
	}
	return ie;
}

//Test if there is a MGPickObject that includes input objin
//in this MGPickObjects' member. If input objin is MGShell,
//and a member is_shell_face(), test is performed to the shell.
//Returns true if objin is included in this MGPickObjects.
MGPickObjects::iterator MGPickObjects::includes(const MGObject* objin){
	iterator i=begin(), ie=end();
	for(; i!=ie; i++){
		const MGPickObject& pobji=**i;
		const MGObject* obji=pobji.leaf_object();
		if(obji == objin)
			return i;

		if(pobji.is_shell_face()){
			obji=pobji.top_object();
			if(obji==objin)
				return i;
		}
	}
	return ie;
}

//add one pobj.
//Function's return value is the numbe of PickObjects defined.
void MGPickObjects::push_back(const MGPickObject& pobj){
	m_PickObjects.emplace_back(pobj.clone());
}
void MGPickObjects::push_back(const MGPickObjects& pobjs){
	size_t n=pobjs.size();
	for(size_t i=0; i<n; i++)
		push_back(pobjs[i]);
}
void MGPickObjects::push_back(MGPickObjects&& pobjs) {
	std::move(pobjs.begin(), pobjs.end(), std::back_inserter(m_PickObjects));
}

void MGPickObjects::remove(const MGPickObject& pobj){
	iterator i=find(pobj);
	if(i!=end())
		erase(i);
}

//Remove objects of type from this pickobjects.
void MGPickObjects::remove(const MGAbstractGels& types){
	int n=(int)size();
	for(int i=n-1; i>=0; i--){
		if(m_PickObjects[i]->is_type(types))
			erase(begin()+i);
	}

}
void MGPickObjects::remove(const MGPickObjects& pobjs){
	const_iterator jend=pobjs.end();
	int n= (int)size();
	for(int i=n-1; i>=0; i--){
		MGPickObject& pobji=*(m_PickObjects[i]);
		if(pobjs.find(pobji)!=jend)
			erase(begin()+i);
	}

}

//Remove gelps from this pickobjects.
void MGPickObjects::remove(const MGGelPositions& gelps){
	MGGelPositions::const_iterator jend=gelps.end();
	int n= (int)size();
	for(int i=n-1; i>=0; i--){
		MGPickObject& pobji=*(m_PickObjects[i]);
		if(gelps.find(pobji)!=jend)
			erase(begin()+i);
	}
}

//reserve the size n, which are all null.
void MGPickObjects::reserve(size_t n){
	m_PickObjects.reserve(n);
}

void MGPickObjects::reset(size_t i, const MGPickObject& pobj){
	m_PickObjects[i].reset(pobj.clone());
}

//Select objects of specified type from this and reset with them.
void MGPickObjects::reset_objects(const MGAbstractGels& types){
	int n=(int)size();
	for(int i=n-1; i>=0; i--){
		if(!(m_PickObjects[i]->is_type(types)))
			erase(begin()+i);
	}
}

//replace this with the common objects of this and pobjs2.
void MGPickObjects::reset_with_common(const MGPickObjects& pobjs2){
	const_iterator jend=pobjs2.end();
	for(int i= (int)size()-1; i>=0; i--){
		if(pobjs2.find(*(m_PickObjects[i]))==jend)
			erase(i);
	}
}

//replace this with symmetric_differecne of this and pobj, that is;
//(1) remove the same MGPickObject from this and pobjs2.
//(2) append the result pobjs2 to this.
void MGPickObjects::reset_with_symmetric_difference(const MGPickObjects& pobjs2){
	MGPickObjects p2dif;
	const_iterator j=pobjs2.begin(),jend=pobjs2.end();
	for(; j!=jend;j++){
		const MGPickObject& pobji=**j;
		iterator i=find(pobji);
		if(i==end())
			p2dif.push_back(pobji);//If pobji not found in this.
		else
			erase(i);//If pobji found in this, erase it.
	}
	push_back(p2dif);//push back not found objects.
}

//Select objects of input type from this.
//Function's return value is pickobjects selected.
//This will be unchanged.
MGPickObjects MGPickObjects::select(const MGAbstractGels& types)const{
	MGPickObjects pobjs2;
	const_iterator i=begin(), ie=end();
	for(auto i = begin(), ie = end(); i!=ie; i++){
		const MGPickObject& pobji=**i;
		if(pobji.is_type(types))
			pobjs2.push_back(pobji);
	}
	return pobjs2;
}

//Select the 1st MGCurve from this.
//Function's return value is MGPickObject of MGCurve 1st encountered in this
//MGPickObject sequence. If this did not includes any MGCurve,
//null MGPickOjbect will be returned.
//This will be unchanged.
MGPickObject MGPickObjects::select_1st_curve()const{
	const MGCurve* curve=0;
	const_iterator i=begin(), ie=end();
	for(; i!=ie; i++){
		const MGObject* obj=(**i).leaf_object();
		curve=dynamic_cast<const MGCurve*>(obj);
		if(curve)
			return **i;
	}
	return MGPickObject();
}

//Select all the MGCurve from this.
//MGPickObject of MGCurve encountered in this MGPickObject sequence will be appended
//in curves.
//This will be unchanged.
void MGPickObjects::select_curves(MGPickObjects& curves)const{
	const_iterator i=begin(), ie=end();
	for(; i!=ie; i++){
		const MGObject* obj=(**i).leaf_object();
		const MGCurve* curve=dynamic_cast<const MGCurve*>(obj);
		if(curve){
			curves.push_back(**i);
		}
	}
}

//Select the 1st MGFSurface from this.
//Function's return value is MGPickObject of MGFSurface 1st encountered in this
//MGPickObject sequence. If this did not includes any MGFSurface,
//null MGPickObject will be returned.
//This will be unchanged.
MGPickObject MGPickObjects::select_1st_fsurface()const{
	const_iterator i=begin(), ie=end();
	for(; i!=ie; i++){
		const MGObject* obj=(**i).leaf_object();
		const MGFSurface* f=dynamic_cast<const MGFSurface*>(obj);
		if(f)
			return **i;
	}
	return MGPickObject();
}

//Select all the MGFSurface from this.
//MGPickObjects of MGFSurface encountered in this MGPickObject sequence will be appended
//in surfaces.
//This will be unchanged.
void MGPickObjects::select_fsurfaces(MGPickObjects& surfaces)const{
	const_iterator i=begin(), ie=end();
	for(; i!=ie; i++){
		const MGObject* obj=(**i).leaf_object();
		const MGFSurface* surface=dynamic_cast<const MGFSurface*>(obj);
		if(surface){
			surfaces.push_back(**i);
		}
	}
}

//Get an 1st object to tessellate from this pick objects.
//Function's return value is the object pointer(MGFSurface).
const MGFSurface* MGPickObjects::get_object_to_tessellate()const{
	const MGPickObjects& picked=*this;

	//Get the tessellation object(MGFace, or MGSurface).
	size_t nobj=picked.size();
	for(size_t i=0; i<nobj; i++){
		const MGObject* obj=picked[i].leaf_object();
		const MGFSurface* fs=obj->fsurface();
		if(fs)
			return fs;
	}
	return 0;//If an object not found, return.
}
