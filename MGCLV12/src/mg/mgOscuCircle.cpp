/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#include "StdAfx.h"
#include "mg/OscuCircle.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//
// Implemetation of MGOscuCircle Class.

//Friend Function

//Constructor

//	MGOscuCircle(const MGOscuCircle&);	//Copy Constructor.

//Destructor
//	~MGOscuCircle();	

//Member Function
MGOscuCircle& MGOscuCircle::add(const MGOscuCircleData& OscuCircle)
	//Add to the end of list.
{
	int size=(int)m_circle.size();
	if(size <= m_n) m_circle.resize(size+10);
	m_circle[m_n]=OscuCircle;
	m_n +=1;
	return *this;
}

MGOscuCircle& MGOscuCircle::add(int index, double radious)
//Add to the end of list.
{
	MGOscuCircleData OscuCircle(index,radious);
	return add(OscuCircle);
}

MGOscuCircleData MGOscuCircle::remove(int i)
	//Remove i-th OscuCircle.
{
	assert(i<m_n);
	MGOscuCircleData OscuCircle=m_circle[i];
	for(int j=i; j<m_n-1; j++) m_circle[j]=m_circle[j+1];
	m_n -=1;
	return OscuCircle;
}

//Operator overload.
//	MGOscuCircle& operator =(MGOscuCircle&);//Assignment operator overload.

const MGOscuCircleData& MGOscuCircle::operator()(int i) const
	// Get i-th OscuCircleData
{
	assert(i<m_n);
	return m_circle[i];
}
