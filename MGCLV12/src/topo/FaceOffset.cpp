/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#include "StdAfx.h"
#include "topo/Edge.h"
#include "topo/Face.h"
#include "topo/Loop.h"
#include "mg/Curve.h"
#include "mg/Surface.h"
#include "mg/Straight.h"
#include "mg/CCisects.h"
#include "mg/CParam_list.h"
#include "mg/Tolerance.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

using namespace std;

//�{�b�N�X�g�Ɉ͂܂���_������UV�Ȑ��𐶐�����
void MGFace::getTrimCrv(
	double uerror, double verror,//u and v parameter error.
	const MGBox& box,	//parameter (u,v) box of the surface.
	std::vector<UniqueCurve>& vecCrv	//paremter curves will be output
) const{
	double u0=box[0].low_point(), u1=box[0].high_point();
	double v0=box[1].low_point(), v1=box[1].high_point();

//	if(u0>159. && v0>-1. && u1<223. && v1<102.)
//		u0=u0;
	const int n = number_of_boundaries();
	for(int i = 0; i < n; i++){

	const MGLoop& lpi = *(loop(i));
	if(!lpi.active())
		continue;
	const int m = lpi.number_of_edges();
	for(int j = 0; j < m; j++){
		const MGEdge& ei = *(lpi.edge(j));
		//if(ei.on_surface_perimeter())continue;
		//�y�����[�^�ȊO�̃o�E���_���̏���
		MGBox boxEdge(ei.box());
		if((boxEdge&box).empty())
			continue;//���ʃ{�b�N�X�Ȃ�
		std::unique_ptr<MGCurve> atpCrv(ei.curve_limitted());

		//�G�b�W���^����ꂽ�{�b�N�X�Ɋ܂܂��Ƃ��̏���
		if(box >> boxEdge){
			//When boxEdge is included in box.
			vecCrv.emplace_back(atpCrv.release());
			continue;
		}

		MGCParam_list isectParamList(atpCrv.get());
		isectParamList.append(atpCrv->param_s());
		isectParamList.append(atpCrv->param_e());

		double dPreWcTol = MGTolerance::set_wc_zero(uerror);
		isectParamList.append(atpCrv->isect_1D(u0,0));
		MGTolerance::set_wc_zero(verror);
		isectParamList.append(atpCrv->isect_1D(v0,1));
		isectParamList.append(atpCrv->isect_1D(v1,1));
		MGTolerance::set_wc_zero(dPreWcTol);

		isectParamList.sort();
		//�p�����[�^�̒��ԓ_�̂t�u�l���{�b�N�X���̂Ƃ��Ȑ��𕪒f���ăx�N�g���ɑ}������
		MGCParam_list::const_iterator  param_iter = isectParamList.begin();
		for(; param_iter != isectParamList.end(); param_iter++){
			MGCParam_list::const_iterator next_param_iter = param_iter;
			next_param_iter++;  //���̃p�����[�^
			if(next_param_iter == isectParamList.end())
				break;				//�I���`�F�b�N
			double t0=*param_iter, t1=*next_param_iter;
			double tmid=(t0+t1)*.5;//���ԃp�����[�^
			MGPosition UVPos = atpCrv->eval(tmid);//�ʏ�p�����[�^
			if(box >> UVPos)
				vecCrv.emplace_back(atpCrv->part(t0,t1));
		}
	}
	
	}
}

//Offset.
//distance is plus value if the direction is toward normal vector of the
//face. Minus if against the normal vector.
//�G���[�R�[�h 0:���� -1:�ȗ����a�ȏ�̃I�t�Z�b�g�s�� -2:�ʐ����R���X�g���N�^�G��>>
int MGFace::offset(double distance, std::vector<UniqueFace>& vecOfsFace)const{
	vecOfsFace.clear();
	const MGSurface *psrf = surface();
	int err;
	std::vector<UniqueSurface> vecOfsSrf = psrf->offset(distance, err);
	if(err < 0)
		return err;

	if(vecOfsSrf.size()==1){
		//���E�������o��
		const std::vector<UniqueLoop>& bounds = boundaries();
		//Face�𐶐�����
		MGFace* f=new MGFace(vecOfsSrf.front().release(), bounds);
		vecOfsFace.emplace_back(f);
		return 0;
	}

	//trim surface
	const MGBox& uvbox=box_param();
	double uerror=uvbox[0].relative_error()*.2;
	double verror=uvbox[1].relative_error()*.2;
	for(auto i = vecOfsSrf.begin(), ie = vecOfsSrf.end(); i!=ie; i++){
		UniqueSurface& srfi=*i;
		std::vector<UniqueCurve> uvcurves;
		getTrimCrv(uerror, verror, srfi->box_param(), uvcurves);
		int n=(int)uvcurves.size();
		if(!n){//If no trim curves.
			if(in_range(srfi->center_param()))
				vecOfsFace.emplace_back(new MGFace(srfi.release()));
			continue;
		}
		//If there were trim curves.
		int pnum,l;
		vector<UniqueCurve>::const_iterator j=uvcurves.begin();
		int oldpnum=-1;
		for(l=0;l<n; j++, l++){
			if(srfi->on_perimeter(**j,pnum)){
				if(oldpnum==-1)
					oldpnum=pnum;
				else if(oldpnum!=pnum)
					break;
			}else
				break;
		}
		if(l==n){//If all the uvcurves were on the one perimeter
			if(in_range(srfi->center_param()))
				vecOfsFace.emplace_back(new MGFace(srfi.release()));
		}else{
			MGFace *pOfsFace = new MGFace(srfi.release());
			j=uvcurves.begin();
			for(l=0;l<n; j++, l++)
				pOfsFace->trim(**j);
			vecOfsFace.emplace_back(pOfsFace);
		}
	}
	return 0;
}

//Offset.
//distance is plus value if the direction is toward normal vector of the
//face. Minus if against the normal vector.
//�G���[�R�[�h 0:���� -1:�ʂɂ��ꂪ���� -2:�ȗ����a�ȏ�̃I�t�Z�b�g�s�� -3:�ʐ����R���X�g���N�^�G���[
MGFace MGFace::offset(double distance, int& error)const
{
	//Face�̖ʂ��I�t�Z�b�g����
	const MGSurface *srf = surface();
	std::unique_ptr<MGSurface> pofsSrf = srf->offset_c1(distance, error);
	if(error < 0){
		if(error==-1){	//Recover the surface of multiplicity order()-1
						//by remove_knot().
			MGSurface* srf2=srf->copy_surface();
			srf2->remove_knot();
			pofsSrf = std::unique_ptr<MGSurface>(srf2->offset_c1(distance, error));
			delete srf2;
			if(error)
				return MGFace();
		}else
			return MGFace();
	}
	//���E�������o��
	const std::vector<UniqueLoop>& bounds = boundaries();
	//Face�𐶐�����
	MGFace f2(pofsSrf.release(), bounds);
	return f2;
}

//Offset.
//distance is plus value if the direction is toward normal vector of the
//FSurface. Minus if against the normal vector.
//�G���[�R�[�h 0:���� -1:�ȗ����a�ȏ�̃I�t�Z�b�g�s�� -3:�ʐ����R���X�g���N�^�G���[
int MGFace::offset_fs(double distance, std::vector<UniqueFSurface>& vecOfsFSurface)const{
	std::vector<UniqueFace> faces;
	int error=offset(distance,faces);
	if(error)
		return error;
	std::move(faces.begin(), faces.end(), std::back_inserter(vecOfsFSurface));
	return 0;
}
