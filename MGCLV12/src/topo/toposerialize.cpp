/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#include "StdAfx.h"
#include "mg/Ifstream.h"
#include "mg/Ofstream.h"
#include "mg/SurfCurve.h"

#include "topo/Complex.h"
#include "topo/PVertex.h"
#include "topo/BVertex.h"
#include "topo/Edge.h"
#include "topo/Face.h"
#include "topo/Loop.h"
#include "topo/Shell.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//--------------------------------------------
// MGTopology�̃V���A���C�Y�֐��Q
// 

//identify_type implements.
long MGComplex::identify_type()const{return MGCOMPLEX_TID;}
long MGLoop::identify_type()const{return MGLOOP_TID;}
long MGShell::identify_type()const{return MGSHELL_TID;}
long MGPVertex::identify_type()const{return MGPVERTEX_TID;}
long MGBVertex::identify_type()const{return MGBVERTEX_TID;}
long MGEdge::identify_type()const{return MGEDGE_TID;}
long MGFace::identify_type()const{return MGFACE_TID;}

// �����o�f�[�^�������݊֐�
void MGCell::WriteMembers(MGOfstream& buf)const{
	MGObject::WriteMembers(buf);
	buf.WritePointer(m_parent_complex);
	buf<<m_perror;
	buf.WritePointer(m_extent.get());
}
// �����o�f�[�^��ǂݏo���֐�
void MGCell::ReadMembers(MGIfstream& buf){
	//�e�N���X�̃����o�f�[�^�̓ǂݏo���B
	MGObject::ReadMembers(buf);
	m_parent_complex=dynamic_cast<MGComplex*>(buf.ReadPointer());
	buf>>m_perror;
	m_extent=UniqueGeometry(dynamic_cast<MGGeometry*>(buf.ReadPointer()));
}

///Write Object's Member Data
void MGBCell::WriteMembers(MGOfstream& buf) const{
	int n =(int)number_of_partner_members();
	buf<<n;
	if(n)
		for(auto& pcell: *m_partners)
			buf.WritePointer(pcell);
}

///Read Object's member data.
void MGBCell::ReadMembers(MGIfstream& buf){
	int i, n;
	buf>>n;
	if (n) {
		m_partners=std::make_unique< PCellVec>(n,nullptr);
		for (i = 0; i < n; i++) {
			(*m_partners)[i] = dynamic_cast<const MGPCell*>(buf.ReadPointer());
		}
	}
}

///Write Object's Member Data
void MGPCell::WriteMembers(MGOfstream& buf) const{
	buf.WritePointer(m_binder.get());
}

///Read Object's member data.
void MGPCell::ReadMembers(MGIfstream& buf){
	MGBCell* bcel= dynamic_cast<MGBCell*>(buf.ReadPointer());
	if(bcel){
		SharedBCell* sb=buf.findSharedBCell(bcel);
		if(sb)
			m_binder = *sb;
		else{
			m_binder.reset(bcel);
			buf.insertSharedBCell(&m_binder);
		}
	} else
		m_binder.reset();
}

///Write Object's Member Data
void MGPVertex::WriteMembers(MGOfstream& buf) const{
	MGPCell::WriteMembers(buf);
	buf<<m_t;
	buf.WritePointer(m_edge);
}

///Read Object's member data.
void MGPVertex::ReadMembers(MGIfstream& buf){
	MGPCell::ReadMembers(buf);
	buf>>m_t;
	m_edge= dynamic_cast<MGEdge*>(buf.ReadPointer());
	assert(m_edge);
}

///Write Object's Member Data
void MGBVertex::WriteMembers(MGOfstream& ostrm) const{
	MGBCell::WriteMembers(ostrm);
	ostrm.WritePointer(m_point.get());
}

///Read Object's member data.
void MGBVertex::ReadMembers(MGIfstream& buf){
	MGBCell::ReadMembers(buf);
	m_point.reset(dynamic_cast<MGPoint*>(buf.ReadPointer()));
}

//Write Object's Member Data
void MGEdge::WriteMembers(MGOfstream& buf) const{
	//�e�N���X�̃����o�f�[�^�̏������݁B
	MGCell::WriteMembers(buf);
	MGPCell::WriteMembers(buf);
	MGBCell::WriteMembers(buf);
	buf<<m_equal_to_binder;	
	buf.WritePointer(m_vertex[0].get());
	buf.WritePointer(m_vertex[1].get());
}

//Read Object's member data.
void MGEdge::ReadMembers(MGIfstream& buf){
	//�e�N���X�̃����o�f�[�^�̓ǂݏo���B
	MGCell::ReadMembers(buf);
	MGPCell::ReadMembers(buf);
	MGBCell::ReadMembers(buf);
	buf>>m_equal_to_binder;
	MGGel* g0=buf.ReadPointer();
	MGGel* g1 = buf.ReadPointer();
	m_vertex[0].reset(dynamic_cast<MGPVertex*>(g0));
	m_vertex[1].reset(dynamic_cast<MGPVertex*>(g1));
}

//--------------------------------------------
// MGFace�̃V���A���C�Y�֐��Q

// �����o�f�[�^���������ފ֐�
void MGFace::WriteMembers(MGOfstream& buf)const{
	//�e�N���X�̃����o�f�[�^�̏������݁B
	MGCell::WriteMembers(buf);
	m_box_param.dump(buf);

	int n = (int)m_boundaries.size();
	buf<<n;
	for(int i = 0; i<n; i++)
		buf.WritePointer(m_boundaries[i].get());
}
// �����o�f�[�^��ǂݏo���֐�
void MGFace::ReadMembers(MGIfstream& buf){
	//�e�N���X�̃����o�f�[�^�̓ǂݏo���B
	MGCell::ReadMembers(buf);
	m_box_param.restore(buf);

	int n;
	buf>>n;
	for(int i = 0; i<n; i++)
		m_boundaries.emplace_back(dynamic_cast<MGLoop*>(buf.ReadPointer()));
}

//--------------------------------------------
// MGComplex�̃V���A���C�Y�֐��Q
// 

// �����o�f�[�^�������݊֐�
void MGComplex::WriteMembers(MGOfstream& buf)const{
	// �e�N���X�̃f�[�^�������݊֐����Ă�ł����B
	MGObject::WriteMembers(buf);
	buf.WritePointer(m_parent_cell);

	int n=(int)m_pcells.size();
	buf<<n;
	for(auto& pcell:m_pcells)
		buf.WritePointer(pcell.get());
}

// �����o�f�[�^�ǂݏo���֐�
void MGComplex::ReadMembers(MGIfstream& buf){
	//�e�N���X�̃����o�f�[�^�̓ǂݏo���B
	MGObject::ReadMembers(buf);
	m_parent_cell = dynamic_cast<MGCell*>(buf.ReadPointer());

	int n;
	buf>>n;
	for(int i=0; i<n; i++)
		m_pcells.emplace_back(dynamic_cast<MGCell*>(buf.ReadPointer()));
}

//--------------------------------------------
// MGLoop�̃V���A���C�Y�֐��Q

// �����o�f�[�^���������ފ֐�
void MGLoop::WriteMembers(MGOfstream& buf)const{
	//�e�N���X�̃����o�f�[�^�̏������݁B
	MGComplex::WriteMembers(buf);
	buf << m_area << int(m_kind)
		<< m_perim_num_s << m_perim_num_e;
}
// �����o�f�[�^��ǂݏo���֐�
void MGLoop::ReadMembers(MGIfstream& buf){
	//�e�N���X�̃����o�f�[�^�̓ǂݏo���B
	MGComplex::ReadMembers(buf);
	int kind;
	buf >> m_area >> kind
		>> m_perim_num_s >> m_perim_num_e;
	m_kind=LoopKind(kind);
}
