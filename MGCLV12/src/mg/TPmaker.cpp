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

// 4辺のTP曲線を構築し、MGSBRepTPにセットする関数
void mgTPmaker::make_tp(
	MGSBRepTP& result // ←ここにセットされる
){
	for(int i = 0; i<4; ++i){		
		// 辺曲線に対してTP曲線のセグメントをひとつ(mgPerimTP)作成する
		// (一辺内で飛び飛びのTP曲線片(mgTPOnFaceBoundary)が生成される)
		mgPerimTP& peri = m_perimTP[i];
		peri.buidTpVec(*m_faces[i]);
		result.set_TP(i,peri.create_tp());//Create the tp.
	}
}

// 接続四辺面用に、細切れのTP曲線を計算する関数
void mgTPmaker::make_subtp(
	MGSBRepVecTP& vectp
){
	// 辺曲線に対してTP曲線のセグメントをひとつ以上作成する
	// (一辺内で飛び飛びのTP曲線片が生成される)
	for(int i = 0; i<4; ++i)
		m_perimTP[i].buidTpVec(*m_faces[i]);

	// mgPerimTPの持つmgTPOnFaceBoundaryからMGLBRepを取得してsubtpsに格納する
	for(int j = 0; j < 4; j++){// 辺[j]に関するmgTPOnFaceBoundaryのコレクション
		std::vector<UniqueLBRep>        subtps;
		auto& cont = m_perimTP[j].subtp_container();
		auto first = cont.begin(), last = cont.end();

		// MGLBRep*をコピーしてsubtpsに格納する
		for(; first != last; ++first)
			subtps.emplace_back(static_cast<MGLBRep*>((*first)->tp().clone()));
	
		vectp.set_TP(j, std::move(subtps), m_perimeters[j]->param_range());
	}
}
