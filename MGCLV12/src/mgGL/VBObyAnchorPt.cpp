#include "StdAfx.h"

#include "mgGL/VBObyAnchorPt.h"


///�f�t�H���g�R���X�g���N�^.
VBObyAnchorPt::VBObyAnchorPt(mgGLSL::CoordinateType coordinateType)
:mgVBO(coordinateType){
}

///copy constructor.
VBObyAnchorPt::VBObyAnchorPt(const VBObyAnchorPt& vbo)
:mgVBO(vbo), m_anchorPointStatic(vbo.m_anchorPointStatic){
}

///Assignment.
VBObyAnchorPt& VBObyAnchorPt::operator=(const VBObyAnchorPt& vbo)
{
	mgVBO::operator=(vbo);
	m_anchorPointStatic = vbo.m_anchorPointStatic;
	return *this;
}

VBObyAnchorPt::~VBObyAnchorPt()
{
	;
}

void VBObyAnchorPt::setStaticAttribAnchorPoint(const MGPosition& pos)
{
	m_anchorPointStatic=pos;
}

const MGPosition& VBObyAnchorPt::anchorPointStatic()const
{
	return m_anchorPointStatic;
}

const MGPosition * VBObyAnchorPt::getAnchorPosition(){
	return &m_anchorPointStatic;
}