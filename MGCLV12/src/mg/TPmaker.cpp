/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/

/** @file TPmaker.cpp
 *  Implementation for class mgTPmaker.
 */
#include "stdafx.h"
#include "mg/LBRep.h"
#include "mg/SBRep.h"
#include "mg/SBRepVecTP.h"
#include "mg/TPOnFaceBoundary.h"
#include "mg/TPmaker.h"

mgTPmaker::mgTPmaker(const MGLBRep* perim[4], const std::vector<const MGFSurface*> faces[4]) {
	assert(is_valid_perim(perim));
	for (int i = 0; i < 4; i++) {
		m_perimeters[i]=perim[i];
		m_faces[i] = &faces[i];
		m_perimTP[i].build(this, i, m_perimeters[i]);
	}
	for (int i = 0; i < 4; i++) {
		const MGLBRep& curvei = *m_perimeters[i];
		const MGLBRep& curvePre = *m_perimeters[(i+3) % 4];
		MGVector V1, V2;

		switch (i) {
		case 0:{
			V1 = curvei.eval(curvei.param_s(), 1);
			V2 = curvePre.eval(curvePre.param_s(), 1);
			break;
		}
		case 1: {
			V1 = curvePre.eval(curvePre.param_e(), 1);
			V2 = curvei.eval(curvei.param_s(), 1);
			break;
		}
		case 2: {
			V1 = curvei.eval(curvei.param_e(), 1);
			V2 = curvePre.eval(curvePre.param_e(), 1);
			break;
		}
		default: {
			V1 = curvePre.eval(curvePre.param_s(), 1);
			V2 = curvei.eval(curvei.param_e(), 1);
			break;
		}
		}
		MGVector N = V1 * V2;
		m_cornerNormal[i] = N.normalize();
	}
}

// 4�ӂ�TP�Ȑ����\�z���AMGSBRepTP�ɃZ�b�g����֐�
void mgTPmaker::make_tp(
	MGSBRepTP& result // �������ɃZ�b�g�����
){
	for(int i = 0; i<4; ++i){		
		// �ӋȐ��ɑ΂���TP�Ȑ��̃Z�O�����g���ЂƂ�(mgPerimTP)�쐬����
		// (��ӓ��Ŕ�є�т�TP�Ȑ���(mgTPOnFaceBoundary)�����������)
		mgPerimTP& peri = m_perimTP[i];
		peri.buidTpVec(*m_faces[i]);
		result.set_TP(i,peri.create_tp());//Create the tp.
	}
}

// �ڑ��l�Ӗʗp�ɁA�א؂��TP�Ȑ����v�Z����֐�
void mgTPmaker::make_subtp(
	MGSBRepVecTP& vectp
){
	// �ӋȐ��ɑ΂���TP�Ȑ��̃Z�O�����g���ЂƂȏ�쐬����
	// (��ӓ��Ŕ�є�т�TP�Ȑ��Ђ����������)
	for(int i = 0; i<4; ++i)
		m_perimTP[i].buidTpVec(*m_faces[i]);

	// mgPerimTP�̎���mgTPOnFaceBoundary����MGLBRep���擾����subtps�Ɋi�[����
	for(int j = 0; j < 4; j++){// ��[j]�Ɋւ���mgTPOnFaceBoundary�̃R���N�V����
		std::vector<UniqueLBRep>        subtps;
		auto& cont = m_perimTP[j].subtp_container();
		auto first = cont.begin(), last = cont.end();

		// MGLBRep*���R�s�[����subtps�Ɋi�[����
		for(; first != last; ++first)
			subtps.emplace_back(static_cast<MGLBRep*>((*first)->tp().clone()));
	
		vectp.set_TP(j, std::move(subtps), m_perimeters[j]->param_range());
	}
}
