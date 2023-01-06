/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#include "StdAfx.h"
#include "topo/HHisects.h"
#include "topo/Shell.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//MGHHisects defines a vector of MGHHisect.
//The vector is implemeted using STL's vector.
//All the methods to handle the vector are available from the STL's vector class,
//and public member m_HHivector. Refer to STL vector class.
//MGHHisects is used to represent intersection lines of a shell with
//another shell, a face, or a surface.
//

////////// Constructor //////////
MGHHisects::MGHHisects(const MGShell * shel1, const MGShell * shel2)
:MGisects(shel1, shel2){
}
///Constructor of 1 MGHHisect.
MGHHisects::MGHHisects(MGHHisect&& hhi){
	m_list.emplace_back(new MGHHisect(std::move(hhi)));
}

void MGHHisects::append(
	const MGFSurface * face1,
	const MGFSurface * face2,
	std::unique_ptr<MGSSisect>&& ssi){
	emplace_back<MGHHisect>(face1,face2, std::move(ssi));
}

const MGShell* MGHHisects::shell1()const{
	return static_cast<const MGShell*>(object1()); 
}
const MGShell* MGHHisects::shell2()const{
	return static_cast<const MGShell*>(object2()); 
}
