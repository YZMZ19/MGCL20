/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#include "StdAfx.h"
#include "mg/Bisection.h"
#include "mg/Box.h"
#include "mg/Knot.h"
#include "mg/Unit_vector.h"
#include "mg/CParam_list.h"
#include "mg/Position_list.h"
#include "mg/Straight.h"
#include "mg/LBRep.h"
#include "mg/TrimmedCurve.h"
#include "mg/CompositeCurve.h"
#include "mg/CSisects.h"
#include "mg/SSisects.h"
#include "mg/SurfCurve.h"
#include "mg/FSurface.h"
#include "mg/Plane.h"
#include "mg/SBRep.h"
#include "mg/RSBRep.h"
#include "mg/Tolerance.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

using namespace std;

//����������2�{��B�\���Ȑ���ڑ�����(������ނ̂Ƃ�)
MGLBRep* join2LBRep(const MGLBRep& crv1, const MGLBRep& crv2) {
	int which, cont;
	double ratio;
	cont = crv1.continuity(crv2, which, ratio);
	if (cont < 0 || which == 0 || which == 3)
		return NULL;

	auto lb = new MGLBRep(crv1);
	lb->connect(cont, which, crv2);
	return lb;
}

//����������B�\���Ȑ����X�g��ڑ�����(LBRep�̂�)�Bjoin_crvl�ɐڑ������Ȑ����X�g������B
//�߂�l�́A�����̋Ȑ����X�g�̌������Ⴄ�Ƃ��A����B�\�����m�łȂ������Ƃ�false���Ԃ�B
int join(
	std::vector<UniqueCurve>& crvl,
	std::vector<UniqueCurve>& join_crvl
) {
	int num = (int)crvl.size(), rc = 0;
	MGCurve* pre_pcrv = 0, * cur_pcrv, * next_pcrv;
	if (!num)
		return false;

	if (num == 1) {	//�Ȑ����P�{�̂Ƃ��̏���
		cur_pcrv = crvl[0]->clone();
		cur_pcrv->remove_knot();
		join_crvl.emplace_back(cur_pcrv);
		return 1;
	}

	cur_pcrv = crvl[0]->clone();
	for (int i = 0; i < num - 1; i++) {
		next_pcrv = crvl[i + 1].get();
		if (!cur_pcrv || !next_pcrv)
			return false;

		MGLBRep* lb1 = dynamic_cast<MGLBRep*>(cur_pcrv);
		MGLBRep* lb2 = dynamic_cast<MGLBRep*>(next_pcrv);
		if (!lb1 || !lb2)
			return false;

		pre_pcrv = join2LBRep(*lb1, *lb2);
		if (!pre_pcrv) {
			rc++;
			cur_pcrv->remove_knot();
			join_crvl.emplace_back(cur_pcrv);
			delete pre_pcrv;
			cur_pcrv = next_pcrv->clone();
			continue;
		}
		delete cur_pcrv;
		cur_pcrv = pre_pcrv;
	}
	cur_pcrv->remove_knot();
	join_crvl.emplace_back(cur_pcrv);
	rc++;
	return rc;
}

// Implementation of MGFSurface projection.

//�Ȑ���ʂɖʒ��܂��̓x�N�g�����e���ċȐ����X�g�����߂�B
///���e�Ȑ��͖ʏ�̃p�����[�^�Ȑ���3�����Ȑ��Ƃ��Ă��ꂼ�ꏇ�ԂɁA
///vec_crv_uv, vec_crv�Ɋi�[�����B
///uv�Ȑ��̃g�������X��rc_zero()���A3�����Ȑ���line_zero()�����ꂼ��g�p���Ă���B
///get perpendicular or vector projection curve list.
///uv projection curves are put into vec_crv_uv(rc_zero() is used),
///3d projection curves are put into vec_crv(line_zero() is used) respectively.
//�߂�l�F
///		���e�Ȑ��̐�:		���e�Ȑ������܂���
///		0:			���e�Ȑ������܂�Ȃ�����
///		-1:			���������G���[
///		-2:			���������G���[�i�������Ȃ������j
///�ǋL�F����vec���^�����Ȃ�(null�j�Ƃ��A�ʒ����e����B
///Obtain the projected curve of a curve onto the surface.
///The direction of the projection is along the vector vec if the vec is not NULL,
///and normal to the surface if the vec is NULL.
///Output of 'project' is two kind of curves:
///one is general world coordinate curves('vec_crv'), and the other is (u,v) curves of
///the parameter space of the surfaces(vec_crv_uv).
///vec_crv_uv.size() is equal to vec_crv.size(). Let the size be n, then
/// (vec_crv_uv[i], vec_crv[i]) is one pair for 0<=i<n.
///Function's return value is:
/// >=0: number of curves obtained, <0 : Some error detected.
int MGCurve::project(
	const MGFSurface& surf,	//given surface.
	std::vector<UniqueCurve>& vec_crv_uv,	//uv projection curve will be appended.
	std::vector<UniqueCurve>& vec_crv,	//3d projection curve will be appended.
	const MGVector& vec	//projection vector.
						//if vec = NULL then calculate perpendicular project.
)const{
	return surf.projectbyApproximateAsLBRep(*this,vec_crv_uv,vec_crv,vec);
		//This is default process of a curve.
}

///�Ȑ���ʂɖʒ��܂��̓x�N�g�����e���ċȐ����X�g�����߂�B
///���e�Ȑ��͖ʏ�̃p�����[�^�Ȑ���3�����Ȑ��Ƃ��Ă��ꂼ�ꏇ�ԂɁA
///vec_crv_uv, vec_crv�Ɋi�[�����B
///uv�Ȑ��̃g�������X��rc_zero()���A3�����Ȑ���line_zero()�����ꂼ��g�p���Ă���B
///get perpendicular or vector projection curve list.
///uv projection curves are put into vec_crv_uv(rc_zero() is used),
///3d projection curves are put into vec_crv(line_zero() is used) respectively.
///�߂�l�F
///		���e�Ȑ��̐�:		���e�Ȑ������܂���
///		0:			���e�Ȑ������܂�Ȃ�����
///		-1:			���������G���[
///		-2:			���������G���[�i�������Ȃ������j
///�ǋL�F����vec���^�����Ȃ�(null�j�Ƃ��A�ʒ����e����B
///Obtain the projected curve of a curve onto the surface.
///The direction of the projection is along the vector vec if the vec is not NULL,
///and normal to the surface if the vec is NULL.
///Output of 'project' is two kind of curves:
///one is general world coordinate curves('vec_crv'), and the other is (u,v) curves of
///the parameter space of the surfaces(vec_crv_uv).
///vec_crv_uv.size() is equal to vec_crv.size(). Let the size be n, then
/// (vec_crv_uv[i], vec_crv[i]) is one pair for 0<=i<n.
///Function's return value is:
/// >=0: number of curves obtained, <0 : Some error detected.
int MGStraight::project(
	const MGFSurface& surf,	//given surface.
	std::vector<UniqueCurve>& vec_crv_uv,	//uv projection curve will be appended.
	std::vector<UniqueCurve>& vec_crv,	//3d projection curve will be appended.
	const MGVector& vec	//projection vector.
						//if vec = NULL then calculate perpendicular project.
)const{
	const MGPlane* plane=dynamic_cast<const MGPlane*>(&surf);
	if(plane){
		return plane->project(*this,vec_crv_uv,vec_crv,vec);
	}
	if(!vec.is_null())	//�x�N�g�����e
		return surf.projVector(*this, vec_crv_uv, vec_crv, vec);
	else
		return surf.projNormal(*this, vec_crv_uv, vec_crv);
}

///�Ȑ���ʂɖʒ��܂��̓x�N�g�����e���ċȐ����X�g�����߂�B
///���e�Ȑ��͖ʏ�̃p�����[�^�Ȑ���3�����Ȑ��Ƃ��Ă��ꂼ�ꏇ�ԂɁA
///vec_crv_uv, vec_crv�Ɋi�[�����B
///uv�Ȑ��̃g�������X��rc_zero()���A3�����Ȑ���line_zero()�����ꂼ��g�p���Ă���B
///get perpendicular or vector projection curve list.
///uv projection curves are put into vec_crv_uv(rc_zero() is used),
///3d projection curves are put into vec_crv(line_zero() is used) respectively.
///�߂�l�F
///		���e�Ȑ��̐�:		���e�Ȑ������܂���
///		0:			���e�Ȑ������܂�Ȃ�����
///		-1:			���������G���[
///		-2:			���������G���[�i�������Ȃ������j
///�ǋL�F����vec���^�����Ȃ�(null�j�Ƃ��A�ʒ����e����B
///Obtain the projected curve of a curve onto the surface.
///The direction of the projection is along the vector vec if the vec is not NULL,
///and normal to the surface if the vec is NULL.
///Output of 'project' is two kind of curves:
///one is general world coordinate curves('vec_crv'), and the other is (u,v) curves of
///the parameter space of the surfaces(vec_crv_uv).
///vec_crv_uv.size() is equal to vec_crv.size(). Let the size be n, then
/// (vec_crv_uv[i], vec_crv[i]) is one pair for 0<=i<n.
///Function's return value is:
/// >=0: number of curves obtained, <0 : Some error detected.
int MGLBRep::project(
	const MGFSurface& surf,	//given surface.
	std::vector<UniqueCurve>& vec_crv_uv,	//uv projection curve will be appended.
	std::vector<UniqueCurve>& vec_crv,	//3d projection curve will be appended.
	const MGVector& vec	//projection vector.
						//if vec = NULL then calculate perpendicular project.
)const{
	if(order()==2){
		const MGPlane* plane=dynamic_cast<const MGPlane*>(&surf);
		if(plane)
			return project(*plane,vec_crv_uv,vec_crv,vec);
	}
	return surf.projectbyRemovKnots(*this,vec_crv_uv,vec_crv,vec);
}
int MGLBRep::project(
	const MGPlane& plane,	//given surface.
	std::vector<UniqueCurve>& vec_crv_uv,	//uv projection curve will be appended.
	std::vector<UniqueCurve>& vec_crv,	//3d projection curve will be appended.
	const MGVector& vec	//projection vector.
						//if vec = NULL then calculate perpendicular project.
)const{
	if(order()==2){
		const MGKnotVector& t=knot_vector();
		int n=t.bdim();
		MGPosition Pi=eval(t(1));
		std::vector<UniqueCurve> vec_crv_uv2;	///<uv projection curve.
		std::vector<UniqueCurve> vec_crv2;	///<3d projection curve.
		for(int i=1; i<n; i++){
			MGPosition Pip1=eval(t(i+1));
			if(Pip1!=Pi){
				MGStraight sli(Pip1, Pi);
				plane.project(sli,vec_crv_uv2,vec_crv2,vec);
				if(i==1){
					vec_crv_uv.emplace_back(new MGLBRep(*(vec_crv_uv2.front()),0,1));
					vec_crv.emplace_back(new MGLBRep(*(vec_crv2.front()),0,1));
				}else{
					MGLBRep* lbiuv=static_cast<MGLBRep*>(vec_crv_uv.back().get());
					lbiuv->connect(0,2,MGLBRep(*(vec_crv_uv2.back()),0,1));
					MGLBRep* lbi=static_cast<MGLBRep*>(vec_crv.back().get());
					lbi->connect(0,2,MGLBRep(*(vec_crv2.back()),0,1));
				}
			}
			Pi=Pip1;
		}
		return (int)vec_crv.size();
	}

	return plane.projectbyRemovKnots(*this,vec_crv_uv,vec_crv,vec);
}

///�Ȑ���ʂɖʒ��܂��̓x�N�g�����e���ċȐ����X�g�����߂�B
///���e�Ȑ��͖ʏ�̃p�����[�^�Ȑ���3�����Ȑ��Ƃ��Ă��ꂼ�ꏇ�ԂɁA
///vec_crv_uv, vec_crv�Ɋi�[�����B
///uv�Ȑ��̃g�������X��rc_zero()���A3�����Ȑ���line_zero()�����ꂼ��g�p���Ă���B
///get perpendicular or vector projection curve list.
///uv projection curves are put into vec_crv_uv(rc_zero() is used),
///3d projection curves are put into vec_crv(line_zero() is used) respectively.
///�߂�l�F
///		���e�Ȑ��̐�:		���e�Ȑ������܂���
///		0:			���e�Ȑ������܂�Ȃ�����
///		-1:			���������G���[
///		-2:			���������G���[�i�������Ȃ������j
///�ǋL�F����vec���^�����Ȃ�(null�j�Ƃ��A�ʒ����e����B
///Obtain the projected curve of a curve onto the surface.
///The direction of the projection is along the vector vec if the vec is not NULL,
///and normal to the surface if the vec is NULL.
///Output of 'project' is two kind of curves:
///one is general world coordinate curves('vec_crv'), and the other is (u,v) curves of
///the parameter space of the surfaces(vec_crv_uv).
///vec_crv_uv.size() is equal to vec_crv.size(). Let the size be n, then
/// (vec_crv_uv[i], vec_crv[i]) is one pair for 0<=i<n.
///Function's return value is:
/// >=0: number of curves obtained, <0 : Some error detected.
int MGCompositeCurve::project(
	const MGFSurface& surf,	//given surface.
	std::vector<UniqueCurve>& vec_crv_uv,	//uv projection curve will be appended.
	std::vector<UniqueCurve>& vec_crv,	//3d projection curve will be appended.
	const MGVector& vec	//projection vector.
						//if vec = NULL then calculate perpendicular project.
)const{
	project_onto_surface(surf,vec_crv_uv,vec_crv,vec);
	return (int)vec_crv.size();
}

///�Ȑ���ʂɖʒ��܂��̓x�N�g�����e���ċȐ����X�g�����߂�B
///���e�Ȑ��͖ʏ�̃p�����[�^�Ȑ���3�����Ȑ��Ƃ��Ă��ꂼ�ꏇ�ԂɁA
///vec_crv_uv, vec_crv�Ɋi�[�����B
///uv�Ȑ��̃g�������X��rc_zero()���A3�����Ȑ���line_zero()�����ꂼ��g�p���Ă���B
///get perpendicular or vector projection curve list.
///uv projection curves are put into vec_crv_uv(rc_zero() is used),
///3d projection curves are put into vec_crv(line_zero() is used) respectively.
///�߂�l�F
///		���e�Ȑ��̐�:		���e�Ȑ������܂���
///		0:			���e�Ȑ������܂�Ȃ�����
///		-1:			���������G���[
///		-2:			���������G���[�i�������Ȃ������j
///�ǋL�F����vec���^�����Ȃ�(null�j�Ƃ��A�ʒ����e����B
///Obtain the projected curve of a curve onto the surface.
///The direction of the projection is along the vector vec if the vec is not NULL,
///and normal to the surface if the vec is NULL.
///Output of 'project' is two kind of curves:
///one is general world coordinate curves('vec_crv'), and the other is (u,v) curves of
///the parameter space of the surfaces(vec_crv_uv).
///vec_crv_uv.size() is equal to vec_crv.size(). Let the size be n, then
/// (vec_crv_uv[i], vec_crv[i]) is one pair for 0<=i<n.
///Function's return value is:
/// >=0: number of curves obtained, <0 : Some error detected.
int MGSurfCurve::project(
	const MGFSurface& surf,	//given surface.
	std::vector<UniqueCurve>& vec_crv_uv,	//uv projection curve will be appended.
	std::vector<UniqueCurve>& vec_crv,	//3d projection curve will be appended.
	const MGVector& vec	//projection vector.
						//if vec = NULL then calculate perpendicular project.
)const{
	if(base_surface()==&surf){
		MGLBRep* lbuv=new MGLBRep;
		m_curve.approximate_as_LBRep(*lbuv);
		vec_crv_uv.emplace_back(lbuv);
		MGLBRep* lbxyz=new MGLBRep;
		approximate_as_LBRep(*lbxyz);
		vec_crv.emplace_back(lbxyz);
	}else{
		surf.projectbyApproximateAsLBRep(*this,vec_crv_uv,vec_crv,vec);
	}
	return (int)vec_crv.size();
}

///�Ȑ���ʂɖʒ��܂��̓x�N�g�����e���ċȐ����X�g�����߂�B
///���e�Ȑ��͖ʏ�̃p�����[�^�Ȑ���3�����Ȑ��Ƃ��Ă��ꂼ�ꏇ�ԂɁA
///vec_crv_uv, vec_crv�Ɋi�[�����B
///uv�Ȑ��̃g�������X��rc_zero()���A3�����Ȑ���line_zero()�����ꂼ��g�p���Ă���B
///get perpendicular or vector projection curve list.
///uv projection curves are put into vec_crv_uv(rc_zero() is used),
///3d projection curves are put into vec_crv(line_zero() is used) respectively.
///�߂�l�F
///		���e�Ȑ��̐�:		���e�Ȑ������܂���
///		0:			���e�Ȑ������܂�Ȃ�����
///		-1:			���������G���[
///		-2:			���������G���[�i�������Ȃ������j
///�ǋL�F����vec���^�����Ȃ�(null�j�Ƃ��A�ʒ����e����B
///Obtain the projected curve of a curve onto the surface.
///The direction of the projection is along the vector vec if the vec is not NULL,
///and normal to the surface if the vec is NULL.
///Output of 'project' is two kind of curves:
///one is general world coordinate curves('vec_crv'), and the other is (u,v) curves of
///the parameter space of the surfaces(vec_crv_uv).
///vec_crv_uv.size() is equal to vec_crv.size(). Let the size be n, then
/// (vec_crv_uv[i], vec_crv[i]) is one pair for 0<=i<n.
///Function's return value is:
/// >=0: number of curves obtained, <0 : Some error detected.
int MGTrimmedCurve::project(
	const MGFSurface& surf,	//given surface.
	std::vector<UniqueCurve>& vec_crv_uv,	//uv projection curve will be appended.
	std::vector<UniqueCurve>& vec_crv,	//3d projection curve will be appended.
	const MGVector& vec	//projection vector.
						//if vec = NULL then calculate perpendicular project.
)const{
	if(is_same_range())
		return m_curve->project(surf,vec_crv_uv,vec_crv,vec);
	return surf.projectbyApproximateAsLBRep(*this,vec_crv_uv,vec_crv,vec);
}

///�^����ꂽ�Ȑ������g�ɖʒ��܂��̓x�N�g�����e���ċȐ����X�g�����߂�B����vec���^�����Ȃ��Ƃ��A�ʒ����e����B
///���e�Ȑ��͖ʏ�̃p�����[�^�Ȑ���3�����Ȑ��Ƃ��Ă��ꂼ�ꏇ�ԂɁAvec_crv_uv, vec_crv�Ɋi�[�����B
///uv�Ȑ��̃g�������X��rc_zero()���A3�����Ȑ���line_zero()�����ꂼ��g�p���Ă���B
///�߂�l�F
///		���e�Ȑ��̐�:		���e�Ȑ������܂���
///		0:			���e�Ȑ������܂�Ȃ�����
///		-1:			���������G���[
///		-2:			���������G���[�i�������Ȃ������j
///Obtain the projected curve of a curve onto the FSurface.
///The direction of the projection is along the vector vec if the vec is not
///NULL, and normal to the FSurface if the vec is NULL.
///Output of 'project' is two kind of curves:
///one is general world coordinate curves('vec_crv'), and the other is (u,v) curves of
///the parameter space of the FSurface(vec_crv_uv).
///vec_crv_uv.size() is equal to vec_crv.size(). Let the size be n, then
///(vec_crv_uv[i], vec_crv[i]) is one pair for 0<=i<n.
///Function's return value is:
/// >=0: number of curves obtained, <0 : Some error detected.
int MGFSurface::project(
	const MGCurve& crv,
	std::vector<UniqueCurve>& vec_crv_uv,
		//Projected curve(surface parameter (u,v) representation) will be appended.
	std::vector<UniqueCurve>& vec_crv,
		//Projected curve(world coordinate(x,y,z) representation) will be appended.
	const MGVector& vec
)const{
	//Invoke each curve project process.
	return crv.project(*this,vec_crv_uv,vec_crv,vec);
}

///�^����ꂽ�Ȑ������g�ɖʒ��܂��̓x�N�g�����e���ċȐ����X�g�����߂�B����vec���^�����Ȃ��Ƃ��A�ʒ����e����B
///���e�Ȑ���3�����Ȑ��Ƃ���vec_crv�Ɋi�[�����B
///uv�Ȑ��̃g�������X��line_zero()���g�p���Ă���B
///�߂�l�F
///		���e�Ȑ��̐�:		���e�Ȑ������܂���
///		0:			���e�Ȑ������܂�Ȃ�����
///		-1:			���������G���[
///		-2:			���������G���[�i�������Ȃ������j
///Obtain the projected curve of a curve onto the FSurface.
///The direction of the projection is along the vector vec if the vec is not
///NULL, and normal to the FSurface if the vec is NULL.
///Output of 'project' is general world coordinate curves('vec_crv')
///Function's return value is:
/// >=0: number of curves obtained, <0 : Some error detected.
int MGFSurface::project(
	const MGCurve& crv,			//given curve.
	std::vector<UniqueCurve>& vec_crv,
		//Projected curve(world coordinate(x,y,z) representation) will be appended.
	const MGVector& vec			//projection vector.
		//if vec = NULL then calculate perpendicular project.
)const{
	std::vector<UniqueCurve> vec_crv_uv;
	return crv.project(*this,vec_crv_uv,vec_crv,vec);
}

///�^����ꂽ�Ȑ������g�ɖʒ��܂��̓x�N�g�����e���ċȐ����X�g�����߂�B
///���e�Ȑ��͖ʏ�̃p�����[�^�Ȑ���3�����Ȑ��Ƃ��Ă��ꂼ�ꏇ�ԂɁA
///vec_crv_uv, vec_crv�Ɋi�[�����B
///uv�Ȑ��̃g�������X��rc_zero()���A3�����Ȑ���line_zero()�����ꂼ��g�p���Ă���B
///get perpendicular or vector projection curve list.
///uv projection curves are put into vec_crv_uv(rc_zero() is used),
///3d projection curves are put into vec_crv(line_zero() is used) respectively.
///�����F
///		const MGCurve&			crv,		(I/ )	given curve.
///		std::vector<UniqueCurve>&	vec_crv_uv,		( /O)	uv projection curve.
///		std::vector<UniqueCurve>&	vec_crv,		( /O)	3d projection curve.
///		const MGVector&			vec=mgNULL_VEC	(I/ )	projection vector.
///�߂�l�F
///		���e�Ȑ��̐�:		���e�Ȑ������܂���
///		0:			���e�Ȑ������܂�Ȃ�����
///		-1:			���������G���[
///		-2:			���������G���[�i�������Ȃ������j
///�ǋL�F����vec���^�����Ȃ�(null�j�Ƃ��A�ʒ����e����B
///Obtain the projected curve of a curve onto the surface.
///The direction of the projection is along the vector vec if the vec is not
///NULL, and normal to the surface if the vec is NULL.
///Output of 'project' is two kind of curves:
///one is general world coordinate curves('vec_crv'), and the other is (u,v) curves of
///the parameter space of the surfaces(vec_crv_uv).
///vec_crv_uv.size() is equal to vec_crv.size(). Let the size be n, then
/// (vec_crv_uv[i], vec_crv[i]) is one pair for 0<=i<n.
///Function's return value is:
/// >=0: number of curves obtained, <0 : Some error detected.
int MGPlane::project(
	const MGStraight& sl,	//given curve.
	std::vector<UniqueCurve>& vec_crv_uv,	//uv projection curve will be appended.
	std::vector<UniqueCurve>& vec_crv,	//3D projection curve will be appended.
	const MGVector& vec	//projection vector.
						//if vec = NULL then calculate perpendicular project.
)const{
	MGSTRAIGHT_TYPE kind=sl.straight_type();
	if(kind== MGSTRAIGHT_TYPE::MGSTRAIGHT_EMPTY)
		return 0;

	const MGPosition P=sl.start_point();
	if(vec.is_null()){
		MGPosition uvP; on(P,uvP);
		if(kind== MGSTRAIGHT_TYPE::MGSTRAIGHT_SEGMENT){
			const MGPosition Q=sl.end_point();
			MGPosition uvQ; on(Q,uvQ);
			vec_crv_uv.emplace_back(new MGStraight(uvQ,uvP));
			vec_crv.emplace_back(new MGStraight(eval(uvQ),eval(uvP)));
			return 1;
		}
		MGPosition dir(sl.direction());
		MGPosition uvdir; on(dir,uvdir);
		vec_crv_uv.emplace_back(new MGStraight(kind,uvdir,uvP));
		vec_crv.emplace_back(new MGStraight(kind,eval(uvdir),eval(uvP)));
	}else{
		MGStraight sl1(MGSTRAIGHT_TYPE::MGSTRAIGHT_UNLIMIT,vec,P);
		MGCSisect is1; relation(sl1,is1);
		const MGPosition& uvPprj=is1.param_surface();
		const MGPosition& Pprj=is1.point();
		if(kind== MGSTRAIGHT_TYPE::MGSTRAIGHT_SEGMENT){
			const MGPosition Q=sl.end_point();
			MGStraight sl2(MGSTRAIGHT_TYPE::MGSTRAIGHT_UNLIMIT,vec,Q);
			MGCSisect is2; relation(sl2,is2);
			const MGPosition& uvQprj=is2.param_surface();
			const MGPosition& Qprj=is2.point();
			vec_crv_uv.emplace_back(new MGStraight(uvQprj,uvPprj));
			vec_crv.emplace_back(new MGStraight(Qprj,Pprj));
			return 1;
		}
		MGPosition dir(sl.direction());
		MGStraight sl3(MGSTRAIGHT_TYPE::MGSTRAIGHT_UNLIMIT,vec,dir);
		MGCSisect is3; relation(sl3,is3);
		const MGPosition& uvDirprj=is3.param_surface();
		const MGPosition& Dirprj=is3.point();
		vec_crv_uv.emplace_back(new MGStraight(kind,uvDirprj,uvPprj));
		vec_crv.emplace_back(new MGStraight(kind,Dirprj,Dirprj));
	}
	return 1;
}

///�^����ꂽ�Ȑ������g�ɖʒ��܂��̓x�N�g�����e���ċȐ����X�g�����߂�B
///���e�Ȑ���3�����Ȑ��Ƃ���vec_crv�Ɋi�[�����B
///3�����Ȑ���tolerance��line_zero()���g�p���Ă���B
///get perpendicular or vector projection curve list.
///3d projection curves are put into vec_crv(line_zero() is used).
///�����F
///		const MGCurve&			crv,		(I/ )	given curve.
///		std::vector<UniqueCurve>&	vec_crv,		( /O)	3d projection curve.
///		const MGVector&			vec=mgNULL_VEC	(I/ )	projection vector.
///�߂�l�F
///		���e�Ȑ��̐�:		���e�Ȑ������܂���
///		0:			���e�Ȑ������܂�Ȃ�����
///		-1:			���������G���[
///		-2:			���������G���[�i�������Ȃ������j
///�ǋL�F����vec���^�����Ȃ�(null�j�Ƃ��A�ʒ����e����B
///Obtain the projected curve of a curve onto the surface.
///The direction of the projection is along the vector vec if the vec is not
///NULL, and normal to the surface if the vec is NULL.
///Output of 'project' is general world coordinate curves('vec_crv')
///Function's return value is:
/// >=0: number of curves obtained, <0 : Some error detected.
int MGPlane::project(
	const MGCurve& crv,	//given curve.
	std::vector<UniqueCurve>& vec_crv,//3d projection curve will be appended.
	const MGVector& vec	//projection vector.
						//if vec = NULL then calculate perpendicular project.
)const{
	std::vector<UniqueCurve> vec_crv_uv;
	return crv.project(*this,vec_crv_uv,vec_crv,vec);
}

//�g�������X�ɉ������ו������f�[�^�|�C���g�����߂�
void prjGetDataPoints(
	const MGFSurface& surf,
	const MGCurve& curve,
	const MGPosition& tuv0,//parameter range of curve to get the data point,
	const MGPosition& tuv1,//From tuv0 to tuv1.
	MGNDDArray& tau		//data points will be output.
){
	double t0=tuv0[0], t1=tuv1[0];
	MGInterval cspan(t0, t1);
	int n1=curve.divideNum(cspan);

	const MGSurface* srf=surf.get_surface_pointer();

	MGCurve* pcrv1=srf->parameter_curve(1,(tuv0[1]+tuv1[1])*.5);
	MGInterval cspan2(tuv0[1], tuv1[1]);
	int n2=pcrv1->divideNum(cspan2);
	delete pcrv1;

	MGCurve* pcrv2=srf->parameter_curve(0,(tuv0[2]+tuv1[2])*.5);
	MGInterval cspan3(tuv0[2], tuv1[2]);
	int n3=pcrv2->divideNum(cspan3);
	delete pcrv2;

	int n=n1;
	if(n<n2)
		n=n2;
	if(n<n3)
		n=n3;

	tau.resize(n+1);
	double delta=(t1-t0)/n;
	double t=t0;
	for(int i=0; i<n; i++, t+=delta){
		tau(i)=t;
	}
	tau(n)=t1;
}

//Get world position data at (u, v) and store the coordinates
//and the parameter in bp.
void prj_store_bp(
	const MGFSurface& f,
	int i,
	double u, double v,
	MGBPointSeq& bp
){
	MGPosition Psrf=f.eval(u, v);	//�n�_�����߂�
	bp.store_at(i,Psrf,0,0,3);
	bp(i,3)=u; bp(i,4)=v;
}

//�ʒ��ɓ��e�����_��ԋp����
//�߂�l�́A��_�܂��͖ʒ��_�����܂����Ƃ���1�A���܂�Ȃ������Ƃ���0��ԋp����
int project_normal(
	const MGFSurface& surf,
	const MGPosition& pos,
	const MGPosition& uv_guess,	//���ʃp�����[�^
	MGPosition& uv
){
	mgTolSetWCZero wczeroSet(MGTolerance::wc_zero()*0.5);//Set&save the error.
	int rc = surf.perp_point(pos, uv, &uv_guess);
	return rc;
}

//���e�����s���ĂT�����Ȑ����擾����
//curve�ɂ�����f�[�^�|�C���g�Ɠ��e�x�N�g�����瓊�e�_����쐬��
//�_�񂩂�T����(xyz,uv)�̋Ȑ�����쐬���ĕԋp����
void prj2OneCurve(
	const MGFSurface& f,
	const MGCurve& curve,	//target curve to prject.
	std::deque<MGPosition>& ranges,//start(ranges[0]) and end(ranges[1]) point parameter
							//of the curve and the face.
							//On return ranges will be so updated that processed rages[0] to [1]
							//are pop_front().
	MGLBRep*& crvProjected	//newed object will be returend when obtained.
							//When not obtained, null will be returned.
){
	crvProjected=0;
	MGPosition tuv0=ranges.front(); ranges.pop_front();
	MGPosition tuv1=ranges.front(); ranges.pop_front();

	MGNDDArray tau;//data point to get the projected points.
				//tau[0] is tuv[0][0], and tau[n-1]=tuv[1][0].
	prjGetDataPoints(f,curve,tuv0,tuv1,tau);

	int ntau = tau.length();
	int ntaum1=ntau-1;

	MGBPointSeq bp_add(ntau, 5);	//xyz,uv�����킹���a�\��[x,y,z,u,v]
	prj_store_bp(f,0,tuv0[1],tuv0[2],bp_add);//�n�_�����߂�

	const MGSurface* srf=f.get_surface_pointer();
	MGPosition uv, uv_guess(tuv0[1],tuv0[2]);
	int i=1;
	for(; i<ntaum1 ; i++){
		MGPosition Pcrv = curve.eval_position(tau(i));
		project_normal(f,Pcrv, uv_guess, uv);
		
		int periNum;
		srf->on_a_perimeter(uv(0), uv(1), periNum);
			//�ӂɂ̂��Ă����炻�̃p�����[�^���g�p����(uv��update����Ă���j
		uv_guess = uv;		//���ʓ_���X�V����
		prj_store_bp(f,i,uv[0],uv[1],bp_add);//�n�_�����߂�
	}
	prj_store_bp(f,i++,tuv1[1],tuv1[2],bp_add);//End point�����߂�
	bp_add.set_length(i);

	//�Ȑ����쐬���ăx�N�g���ɑ}������
	MGBox bxAddXYZ(3,bp_add.box());	//xyz�̃{�b�N�X�����߂�
	//�V���[�g�A�[�N�̃`�F�b�N������
	if(bxAddXYZ.len() <= MGTolerance::wc_zero())
		return;

	MGBPointSeq tmpBpXYZ(3,bp_add);
	MGNDDArray tauXYZ(tmpBpXYZ);
	crvProjected = new MGLBRep;
	crvProjected->buildByInterpolationDataPoints(tauXYZ, bp_add, 4, 4.);
	crvProjected->remove_knot(0,3);	//�ŏ��̂R����(XYZ)��line_zero�Ńm�b�g�폜
}

//Normalize the parameter of the input range[.][0], i.e. if the value is alomost equalt to
//a knot value, round the value into the knot value. This is usefull especially for
//start or end parameter value.
int normalize_check_isect(
	const MGCurve& curve,
	int retval,
	MGPosition range[2]
){
	//�p�����[�^�l�𐳋K������
	MGPosition& tuv0=range[0];
	MGPosition& tuv1=range[1];
	double ts=tuv0(0)=curve.param_normalize(tuv0[0]);
	double te=tuv1(0)=curve.param_normalize(tuv1[0]);
	if((te-ts)<=curve.param_error()){
		if(retval==1)
			return 0;	//���e�͈͂��덷�͈͂̏ꍇ��_�ƂȂ�̂�0��Ԃ�
		else
			return -2;
	}
	return retval;
}

class MGPrjBisect{
public:

MGPrjBisect(
	double ts,	//parameter range from ts to te.
	double te
);

//Virtual Destructor
virtual ~MGPrjBisect(){;};

//compare with the previous function value(the initial value is set
//by set_initial_t) and replace t with the previous one if necessary.
//The function's return value is the new parameter value.
virtual void compare_replace(
	double t	//parameter value to compare at.
)=0;

int solve(
	double tolerance//The tolerance to halt the bisection iteration.
);

protected:
	double m_ts, m_te;
		//the curve's parameter range from m_ts to m_te, or from m_te to m_ts.
		//note that m_ts<m_te does not always hold.
};

class MGProjBoundaryParam: public MGPrjBisect{
public:
	MGProjBoundaryParam(
		const MGFSurface& surf,
		const MGCurve& curve,
		const MGPosition& tuv,
		double ts, double te
	):MGPrjBisect(ts,te),m_surf(surf),m_curve(curve),m_uv(tuv[1],tuv[2]){
		if(tuv[0]==ts){
			m_te=ts;
			m_ts=te;
		}
	//m_te holds the curve parameter value whose point has
	//a perpendicular point ont the surface.
	}

	//compare with the previous function value(the initial value is set
	//by set_initial_t) and replace t with the previous one if necessary.
	//The function's return value is the new parameter value.
	virtual void compare_replace(
		double t	//parameter value to compare at.
	){
		MGPosition uv;
		MGPosition P=m_curve.eval(t);
		if(project_normal(m_surf,P,m_uv, uv)){
			m_uv=uv;
			m_te=t;
		}else{
			m_ts=t;
		}
	}

	MGPosition get_result(){return MGPosition(m_te,m_uv[0],m_uv[1]);};
	const MGFSurface& m_surf;
	const MGCurve& m_curve;
	MGPosition m_uv;
};

//curve��̓��e�\�ȃp�����[�^�Ɠ��e�s�\�ȃp�����[�^��^���āA
//�Ԃɂ��铊�e�\�ȋ��E�p�����[�^�l�����߂ĕԋp����
//�`�F�b�N�|�C���g�̈ړ��� < (�p�����[�^�͈�*rc_zero())�ɂȂ�ΏI���B
//Function's return value is MGPosition tuv, where
//tuv[0]=start point parameter of curve,
//(tuv[1], tuv[2])=the corresponding surface(this surface's) parameter (u,v)
//of the start point of the range.
MGPosition proj2GetParamIter(
	const MGFSurface& f,
	const MGCurve& curve,//target curve to project.
	const MGPosition& tuv,//tuv[0] is the curve's parameter value that has a perp point onto the
				//surface. and (tuv[1], tuv[2]) is the surface's parameter of the perp point.
				//tuv[0] is eithere ts or te.
	double ts,	//(ts, te) is curve's parameter range that indicates
	double te	//between ts and te there must be the boundary point.
				//That is one of the following situations occurs:
				//(1) there are no perp points at ts and there is a perp point at te,
				//(2) there is a perp point at ts and there are no perp points at te,
){
	//������
	const double crvError = curve.param_error()*2.;//, srfError = get_surface_pointer()->param_error();
	MGProjBoundaryParam solver(f,curve,tuv,ts,te);
	solver.solve(crvError);
	return solver.get_result();
}

//A virtual super class to solve non-linear equations by bicection methos.
MGPrjBisect::MGPrjBisect(
	double ts,	//parameter range from ts to te.
	double te
):m_ts(ts), m_te(te){}

//Compute the fn(t)'s parameter value that is the maxima 
//Function's return value will be the solution obtained.
int MGPrjBisect::solve(
	double tolerance//The tolerance to halt the bisection iteration.
){
	assert(tolerance>0.);
	int nrepition=0;
	double span=m_te-m_ts;
	tolerance+=MGTolerance::mach_zero()*2.;
	while(fabs(span)>=tolerance){//m_ts > m_te may happend.
		span*=.5; nrepition++;
	    double t=m_ts+span;
		compare_replace(t);
	}
	return nrepition;
}

//Obtain the projected curve of a curve onto the surface.
//The direction of the projection is along the vector vec if the vec is not
//NULL, and normal to the surface if the vec is NULL.
//Output of 'project' is two kind of curves:
//one is general world coordinate curves('vec_crv'), and the other is (u,v) curves of
//the parameter space of the surfaces(vec_crv_uv).
//vec_crv_uv.size() is equal to vec_crv.size(). Let the size be n, then
// (vec_crv_uv[i], vec_crv[i]) is one pair for 0<=i<n.
///Function's return value is:
/// >=0: number of curves obtained, <0 : Some error detected.
///Function's return value is:
/// >=0: number of curves obtained, <0 : Some error detected.
//�߂�l�F
//		���e�Ȑ��̐�:		���e�Ȑ������܂���
//		0:			���e�Ȑ������܂�Ȃ�����
//		-1:			���������G���[
//		-2:			���������G���[�i�������Ȃ������j
int MGFSurface::projectbyApproximateAsLBRep(
	const MGCurve& crv,
	std::vector<UniqueCurve>& vec_crv_uv,
		//Projected curve(surface parameter (u,v) representation) will be appended.
	std::vector<UniqueCurve>& vec_crv,
		//Projected curve(world coordinate(x,y,z) representation) will be appended.
	const MGVector& vec
)const{
	//������
	std::unique_ptr<MGLBRep> lb(new MGLBRep);
	crv.approximate_as_LBRep(*lb);
	if(!vec.is_null())	//�x�N�g�����e
		return projVector(*lb, vec_crv_uv, vec_crv, vec);
	else
		return projNormal(*lb, vec_crv_uv, vec_crv);
}

//Obtain the projected curve of a curve onto the surface.
//The direction of the projection is along the vector vec if the vec is not
//NULL, and normal to the surface if the vec is NULL.
//Output of 'project' is two kind of curves:
//one is general world coordinate curves('vec_crv'), and the other is (u,v) curves of
//the parameter space of the surfaces(vec_crv_uv).
//vec_crv_uv.size() is equal to vec_crv.size(). Let the size be n, then
// (vec_crv_uv[i], vec_crv[i]) is one pair for 0<=i<n.
//�߂�l�F
//		���e�Ȑ��̐�:		���e�Ȑ������܂���
//		0:			���e�Ȑ������܂�Ȃ�����
//		-1:			���������G���[
//		-2:			���������G���[�i�������Ȃ������j
///Function's return value is:
/// >=0: number of curves obtained, <0 : Some error detected.
int MGFSurface::projectbyRemovKnots(
	const MGCurve& crv,
	std::vector<UniqueCurve>& vec_crv_uv,
		//Projected curve(surface parameter (u,v) representation) will be appended.
	std::vector<UniqueCurve>& vec_crv,
		//Projected curve(world coordinate(x,y,z) representation) will be appended.
	const MGVector& vec
)const{
	//������
	std::unique_ptr<MGCurve> crvRemoveKnot(crv.clone());
	crvRemoveKnot->remove_knot();
	if(!vec.is_null())	//�x�N�g�����e
		return projVector(*crvRemoveKnot, vec_crv_uv, vec_crv, vec);
	else
		return projNormal(*crvRemoveKnot, vec_crv_uv, vec_crv);
}

//����������2�{��B�\���Ȑ���ڑ�����(������ނ̂Ƃ�)
MGLBRep* join2LBRep(const MGLBRep& crv1, const MGLBRep& crv2);

//�����������T��(xyzuv)��B�\���Ȑ����X�g��(x,y,z)��(u,v)��2��̋Ȑ��ɕ������A
//join_crvl_uv, _xyz�Ɋi�[����B�����A�A������Ȑ����A���ł���΁A�ڑ�����B
//join_crvl_uv,join_crv_xyz�ɐڑ������Ȑ����X�g������B
//Function's return value is the number of output curves.
void prj2Pieces(
	const MGLBRep& lb5,
	std::vector<UniqueCurve>& join_crvl_uv,
	std::vector<UniqueCurve>& join_crvl_xyz
){
	join_crvl_uv.emplace_back(new MGLBRep(2,lb5,0,3));
	join_crvl_xyz.emplace_back(new MGLBRep(3,lb5,0,0));
}
	
//�����������T��(xyzuv)��B�\���Ȑ����X�g��(x,y,z)��(u,v)��2��̋Ȑ��ɕ������A
//join_crvl_uv, _xyz�Ɋi�[����B�����A�A������Ȑ����A���ł���΁A�ڑ�����B
//join_crvl_uv,join_crv_xyz�ɐڑ������Ȑ����X�g������B
//Function's return value is the number of output curves.
void prjJoin(
	std::vector<UniqueLBRep>& crvl,
	std::vector<UniqueCurve>& join_crvl_uv,
	std::vector<UniqueCurve>& join_crvl_xyz
){
	size_t num = crvl.size();
	if(!num)
		return;

	if(num==1){	//�Ȑ����P�{�̂Ƃ��̏���
		prj2Pieces(*(crvl[0]),join_crvl_uv,join_crvl_xyz);
		return;
	}

	MGLBRep *cur_pcrv;
	MGLBRep *next_pcrv;
	cur_pcrv = crvl[0].release();
	for(size_t i=1; i<num; i++){
		next_pcrv = crvl[i].release();
		MGLBRep *pre_pcrv = join2LBRep(*cur_pcrv,*next_pcrv);
		if(pre_pcrv){//If successfully joined.
			delete cur_pcrv;
			delete next_pcrv;
			cur_pcrv = pre_pcrv;
		}else{
			prj2Pieces(*cur_pcrv,join_crvl_uv,join_crvl_xyz);
			delete cur_pcrv;
			cur_pcrv = next_pcrv;
		}
	}
	prj2Pieces(*cur_pcrv,join_crvl_uv,join_crvl_xyz);
	delete cur_pcrv;
}

///Define curve division number when a curve crv be projected onto this MGFSurface.
///The result is used in prj2GetParamRange().
int MGSurface::get_proj_divnum(const MGCurve& crv)const{
	return crv.intersect_dnum()+2;
}


//�x�N�g�����e�́A�J�[�u��܂�ŕ������čs���A��Őڑ�����
int MGFSurface::projVector(
	const MGCurve& crv,
	std::vector<UniqueCurve>& vec_crv_uv,
		//Projected curve(surface parameter (u,v) representation) will be appended.
	std::vector<UniqueCurve>& vec_crv,
		//Projected curve(world coordinate(x,y,z) representation) will be appended.
	const MGVector& vec
)const{
	std::vector<UniqueCurve> mult_crvl;
	std::vector<UniqueLBRep> temp_vec_crv;

	//�}���`�m�b�g�ɂȂ�Ƃ��͐؂��ē��e���čēx�ڑ�����
	int num_mult = crv.divide_multi(mult_crvl);
	for(int i=0; i<num_mult; i++)
		projVectorProc(*(mult_crvl[i]), temp_vec_crv, vec);
	prjJoin(temp_vec_crv, vec_crv_uv, vec_crv);
	return (int)vec_crv.size();
}

#define EXTEND_COEF 2.0		//�x�N�g�����e�̃X�C�[�v�ʂ̒��������߂�̂Ɏg�p
//Obtain the projected curve of a curve onto the surface.
//projVectorProc does not divide the curve at the c0 continuity, which is treated by
//project(). See projec().
//The direction of the projection is vec, which is not null.
//Output of 'projVectorProc' is curves of space dimension 5, which are (x,y,z,u,v).
//Here (x,y,z) is the world coordinates, and (u,v) is parameter of this surface.
int MGFSurface::projVectorProc(
	const MGCurve& crv,
	std::vector<UniqueLBRep>& crv_xyzuv_vector,//Projected curve of(x,y,z,u,v) will be appended.
	const MGVector& vec
)const{
	double extLeng;
	if(get_surface_pointer()->type() == MGSURFACE_TYPE::MGSURFACE_PLANE){
		MGPosition Ps=crv.start_point(), Pe=crv.end_point(), Pc=crv.center();
		MGPosition PsNear=closest(Ps); extLeng=PsNear.distance(Ps);
		MGPosition PeNear=closest(Pe); extLeng+=PeNear.distance(Pe);
		MGPosition PcNear=closest(Pc); extLeng+=PcNear.distance(Pc)+1.;
		//+1. to avoid the case that
		//all of (Ps, Pe, Pc) are on the plane and the sum of their distance is zero.
		extLeng *= EXTEND_COEF;	//���S�̈�
	}else{
		MGBox box(get_box() | crv.box());
		extLeng = (box.len()+1.) * EXTEND_COEF;	//���S�̈�
	}

	std::unique_ptr<MGSurface> apCutSrf(crv.sweep(vec, extLeng, -extLeng));
	if(!apCutSrf.get())
		return 0;

	double errsave=MGTolerance::set_wc_zero(MGTolerance::line_zero());
	MGSSisects isectl = isectFS(*apCutSrf);
	MGTolerance::set_wc_zero(errsave);

	const double srfError=param_error();
	int i;
	MGSSisects::iterator ssiter;
	for(ssiter = isectl.begin(), i = 0; ssiter != isectl.end(); ssiter++, i++){
		auto& isec=isectCast<MGSSisect>(ssiter);
		//�J�[�u�̌��������킹�Ă���std::vector<MGCurve>��push����
		std::unique_ptr<MGCurve> apCrv2d(isec.release_param1());
		std::unique_ptr<MGCurve> apCrv3d(isec.release_line());
		MGCurve &crvOth2d = isec.param2();
		assert(apCrv2d.get() && apCrv3d.get());
		//u�p�����[�^�������������n�_�Ƃ���(���J�[�u�̌����Ɠ����ɂ���)
		double spt = (crvOth2d.start_point())(0), ept = (crvOth2d.end_point())(0);
		if(spt > ept){
			apCrv2d->negate();
			apCrv3d->negate();
		}

		double	param_len = apCrv2d->box().len(),	//�V���[�g�A�[�N�̃`�F�b�N������
				world_len = apCrv3d->box().len();
		if(param_len < srfError || world_len < MGTolerance::wc_zero())continue;
		MGLBRep *pCrv2d = dynamic_cast<MGLBRep*>(apCrv2d.get());
		MGLBRep *pCrv3d = dynamic_cast<MGLBRep*>(apCrv3d.get());
		assert(pCrv2d && pCrv3d);
		int bdim = pCrv2d->bdim();
		MGBPointSeq bpXYZUV(5,pCrv3d->line_bcoef());
		const MGBPointSeq &bpUV = pCrv2d->line_bcoef();
		for(int i2=0; i2<bdim; i2++)
			bpXYZUV.store_at(i2,bpUV(i2),3,0);

		MGLBRep* lb=new MGLBRep;
		lb->buildLBRepFromMemberData(
			std::move(pCrv2d->knot_vector()),std::move(bpXYZUV));
		lb->copy_appearance(crv);
		crv_xyzuv_vector.emplace_back(lb);
	}
	return i;
}

///�J�[�u��܂�ŕ������čs���A��Őڑ�����
int MGFSurface::projNormal(
	const MGCurve& crv,
	std::vector<UniqueCurve>& vec_crv_uv,
	std::vector<UniqueCurve>& vec_crv
)const{
	std::vector<UniqueCurve> divided_curves;
	//�}���`�m�b�g�ɂȂ�Ƃ��͐؂��ē��e����
	int n=crv.divide_multi(divided_curves);
	std::vector<UniqueLBRep> crv_xyzuv_vector;//Projected curve of(x,y,z,u,v).
	for(int i=0; i<n; i++){
		projNormalProc(*(divided_curves[i]),crv_xyzuv_vector);
	}

	//���e�Ȑ���̏璷�m�b�g���폜���Đڑ�����
	prjJoin(crv_xyzuv_vector, vec_crv_uv, vec_crv);
	return (int)vec_crv.size();
}

//�X�^�[�g�p�����[�^��^�����e�\�ȃp�����[�^�͈͂�1�����擾����B
//�߂�l��	1:�͈͂����܂���(this�̍Ō�܂�)
//			0:�͈͂����܂�Ȃ�����(this�̍Ō�܂�)
//			-1:�͈͂����܂���(this�̓r���܂ŁA�܂��ʒ��͈͂����邩������Ȃ�)
//			-2;�͈͂����܂�Ȃ��������Athis�̓r���ŁA�܂��ʒ��͈͂����邩������Ȃ�
int prj2GetParamRange(
	const MGFSurface& f,
	const MGCurve& curve,
	int start_counter,	//input start counter of the curve parameter incrementation.
	MGPosition range[2],
		//range[0][0]=start point parameter of curve,
		//range[0][1-2]=the corresponding surface(this surface's) parameter (u,v)
		//of the start point of the range.
		//Regarding to range[1] : the same.
	int& next_counter	//Updated curve parameter incremental counter will be output.
){
	//������
	//int ndiv = curve.intersect_dnum()+2;
	int ndiv = f.get_proj_divnum(curve);

	//curve�𕪊����A���e�\�����ׂ�
	const double delta=curve.param_span()/ndiv;	//�`�F�b�N�|�C���g�̃p�����[�^�X�p��
	const double tend=curve.param_e();

	int t_is_on_surface;
	double t=curve.param_s()+delta*double(start_counter), t_pre;
	MGPosition uv;
	int i=start_counter;//counter for t's incremental number.
	if(t_is_on_surface=f.perp_one(curve.eval(t),uv)){
		range[0]=MGPosition(t,uv[0],uv[1]);
	}else{
		//Find the 1st t that is on surface.
		while(!t_is_on_surface && i<ndiv){
			t_pre=t;
			t+=delta; i++;
			t_is_on_surface=f.perp_one(curve.eval(t),uv);
		}
		if(t_is_on_surface){
			range[0]=proj2GetParamIter(f,curve,MGPosition(t,uv[0],uv[1]),t_pre,t);
		}else{
			if(f.perp_one(curve.eval(tend),uv)){
				range[0]=proj2GetParamIter(f,curve,MGPosition(tend,uv[0],uv[1]),tend-delta,tend);
				range[1]=MGPosition(tend,uv[0],uv[1]);
				return normalize_check_isect(curve,1,range);
			}
			return 0;
		}
	}

	//Here t_is_on_surface=true always holds.
	MGPosition uv_save;
	while(t_is_on_surface && t<tend && i<ndiv){
		uv_save=uv;
		t_pre=t;
		t+=delta; i++;
		t_is_on_surface=f.perp_point(curve.eval(t),uv,&uv_save);
	}

	int retval;
	if(!t_is_on_surface){
		range[1]=proj2GetParamIter(f,curve,MGPosition(t_pre,uv_save[0],uv_save[1]),t_pre,t);
		retval=-1;
	}else{
		uv_save=uv;
		if(f.perp_point(curve.eval(tend),uv,&uv_save)){
			range[1]=MGPosition(tend,uv[0],uv[1]);
		}else{
			range[1]=proj2GetParamIter(f,curve,MGPosition(t,uv[0],uv[1]),t,tend);
		}
		retval=1;
	}

	//�I������
	next_counter=i;		//����̌����J�n�p�����[�^��Ԃ�
	return normalize_check_isect(curve,retval,range);
}

//Obtain the projected curve of a curve onto the surface.
//projNormalProc does not divide the curve at the c0 continuity, which is treated by
//project(). See projec().
//The direction of the projection is normal to the surface.
//Output of 'projNormalProc' is curves of space dimension 5, which are (x,y,z,u,v).
//Here (x,y,z) is the world coordinates, and (u,v) is parameter of this surface.
void MGFSurface::projNormalProc(
	const MGCurve& crv,	//The target curve to project.
	std::vector<UniqueLBRep>& crv_xyzuv_vector//Projected curve of(x,y,z,u,v) will be appended.						
)const{
	//�p�����[�^�͈͂����߂�
	std::deque<MGPosition> ranges;//ranges[2*i] to ranges[2*i+1] is the projection range.
			//Let ranges[.]=tuv, then tuv[0] is the parameter of the curve, and (tuv[1], tuv[2])
			//is the parameter of this surface (u,v).

	double sparam = crv.param_s();
	int counter=0, rc;
	do{
		MGPosition tuv[2];
		rc = prj2GetParamRange(*this,crv,counter,tuv,counter);
		if(rc==1 || rc==-1){
			ranges.push_back(tuv[0]);
			ranges.push_back(tuv[1]);
		}
	}while(rc < 0);	//crv1�̓r���ł���ΌJ��Ԃ�

	//���e���s��
	while(ranges.size()>=2){
		//ranges will be pop_front by prj2OneCurve.
		MGLBRep* projectedCurve;
		prj2OneCurve(*this,crv,ranges,projectedCurve);

		//���e�Ȑ���̏璷�m�b�g���폜���Đڑ�����
		if(projectedCurve)
			crv_xyzuv_vector.emplace_back(projectedCurve);
	}
}
