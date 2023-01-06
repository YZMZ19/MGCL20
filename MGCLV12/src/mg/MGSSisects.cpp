/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#include "StdAfx.h"
#include "mg/FSurface.h"
#include "mg/SSisects.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

// MGSSisects defines singly linked list of MGSSisect.
// Used to represent intersection points of two surfaces.

MGSSisects::MGSSisects(const MGFSurface * s1, const MGFSurface * s2)
:MGisects(dynamic_cast<const MGObject*>(s1)
		  ,dynamic_cast<const MGObject*>(s2)){
}
void MGSSisects::append(MGSSisect&& isect){
// Adds the MGSSisect to the end of the list.
	emplace_back(new MGSSisect(std::move(isect)));
}
void MGSSisects::append(MGSSisects&& isectlist){
	assert(object1()==isectlist.object1());
	assert(object2()==isectlist.object2());
	m_list.splice(m_list.end(),std::move(isectlist.m_list));
}

// 全てのコンポーネントを指定して交線を追加
//Add one intersection line to the list.
//iline, param1, and param2 must be newed objects, and their ownership
//are transfered to MGSSisects.
void MGSSisects::append(
	MGCurve* iline,
	MGCurve* param1,
	MGCurve* param2,
	const MGSSRELATION r1){
// Adds the MGSSisect to the end of the list.
	emplace_back<MGSSisect>(iline,param1,param2,r1);
}

// 全てのコンポーネントを指定して交線を追加
//Add one intersection line to the list.
void MGSSisects::append(
	const MGCurve& iline,
	const MGCurve& param1,
	const MGCurve& param2,
	const MGSSRELATION r1){
// Adds the MGSSisect to the end of the list.
	emplace_back<MGSSisect>(iline,param1,param2,r1);
}

//Find where in this ssi2  have common parts (in line_zero()) in 
//their world representation.
//Fucntion's return value is the iterator of this that had the common.
//		!=end():have common part. 
//		==end():no common part(except a point) found.
MGSSisects::iterator MGSSisects::find_common(const MGSSisect& ssi2){
	iterator i=begin(), ie=end();
	for(; i!=ie; i++){
		if(ssi2.has_common(isectCast<MGSSisect>(i)))
			return i;
	}
	return ie;
}

std::unique_ptr<MGSSisect> MGSSisects::removeAt(iterator i){
//Remove the MGSSisect and return the MGSSisect. If i is no valid, 
// behavior is undefined.
	auto isect=release<MGSSisect>(i);
	return isect;
}

std::unique_ptr<MGSSisect> MGSSisects::removeFirst(){
//Remove the first MGSSisect int the list and return the MGSSisect.
//If i is not valid, behavior is undefined.
	auto isect=releaseFront<MGSSisect>();
	return isect;
}

std::unique_ptr<MGSSisect> MGSSisects::removeLast(){
//Remove the first MGSSisect int the list and return the MGSSisect.
//If i is not valid, behavior is undefined.
	auto isect=releaseBack<MGSSisect>();
	return isect;
}

///Return the pointer to surface1.
const MGFSurface* MGSSisects::surface1() const{
	return dynamic_cast<const MGFSurface*>(object1());
}

///Return the pointer to surface2.
const MGFSurface* MGSSisects::surface2() const{
	return dynamic_cast<const MGFSurface*>(object2());
}
