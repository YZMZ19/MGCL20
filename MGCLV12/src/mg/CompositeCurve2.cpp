/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#include "StdAfx.h"
#include "mg/CCisects.h"
#include "mg/CompositeCurve.h"
#include "mg/LBRep.h"
#include "mg/Plane.h"
#include "mg/SurfCurve.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//
//Implements MGCompositeCurve Class.
//MGCompositeCurve is a composite of other leaf curves.
//Assumedly they are connected as C0 continuity. However, MGCompositeCurve
//does not check their continuity, but only put first or last as the user says
// (in connect_to_end or in connect_to_start), except two curves connecting
//are both MGLBRep or MGRLBRep. When the two connecting curves are both
//MGLBRep or MGRLBRep and they have more than C0 continuity,
//they are changed to one MGLBRep orMGRLBRep representation.
//Parameter ranges of the member curves are always continuous, from param_s() of the
//1st curve to param_e() of the last form MGCompositeCurve's paramter range.

//�֐����F
//		common
//�ړI�F
//		�^����ꂽ�Ȑ��Ǝ��g�̌�_�������͋��ʕ��������邩�ǂ������ׂ�B
//�����F
//		const MGCurve&			crv2,		(I/ )	�^������Ȑ�
//		std::vector<double>&	vec_param	( /O)	���ʕ����̃p�����[�^�͈�
//		MGCCisects&			isect		( /O)	��_
//				 4n�̔z��ŁAt(4*i+0),t(4*i+1)�����g�̃p�����[�^�͈�(t(4*i+0) < t(4*i+1))�A
//							 t(4*i+2),t(4*i+3)���^�Ȑ��̃p�����[�^�͈�(f(t(4*i+0))=f(t(4*i+2))
//�߂�l�F
//		3:��_�����ʕ��������܂���
//		2:��_�݂̂����܂���
//		1:���ʕ����݂̂����܂���
//		0:��_�����ʕ������Ȃ�����
//		-1:���ʃG�b�W�̎����v�Z�G���[
//		-2:���ʃG�b�W���S�ȏ㋁�܂���(�̂��Ă��Ȃ��ƌ��Ȃ�)
//�ǋL�F
//	�Ȑ������ʂ��ǂ����̌덷�ɂ�line_zero()�A���p�����[�^�͈͂̎����v�Z��
//	�덷�ɂ́A�p�����[�^�͈�*rc_zero()���g�p����
int MGCompositeCurve::common(
	const MGCurve& crv2,
	std::vector<double>& vecComSpan,
	MGCCisects& isect
)const{
	int obtainedIsect=0, obtainedCommon=0;
	const_iterator i=begin(), ie=end();
	for(; i!=ie; i++){
		std::vector<double> vecComSpan2;
		MGCCisects isect2;
		int obtained2=(**i).common(crv2,vecComSpan2,isect2);
		if(obtained2<0)
			return obtained2;
		size_t ncom=vecComSpan2.size();
		if(ncom){
			obtainedCommon=1;
			for(size_t j=0; j<ncom; j++)
				vecComSpan.push_back(vecComSpan2[j]);
		}
		size_t nisect=isect2.size();
		if(nisect){
			obtainedIsect=2;
			isect.append(std::move(isect2));
		}
	}
	return obtainedCommon+obtainedIsect;
}

//�֐����F
//		common
//�ړI�F
//		�^����ꂽ�Ȑ��Ǝ��g�̋��ʕ��������邩�ǂ������ׂ�B
//�����F
//		const MGCurve&			crv2,		(I/ )	�^������Ȑ�
//		std::vector<double>&	vec_param	( /O)	���ʕ����̃p�����[�^�͈�
//				 4n�̔z��ŁAt(4*i+0),t(4*i+1)�����g�̃p�����[�^�͈�(t(4*i+0) < t(4*i+1))�A
//							 t(4*i+2),t(4*i+3)���^�Ȑ��̃p�����[�^�͈�(f(t(4*i+0))=f(t(4*i+2))
//�߂�l�F
//		���ʕ����̐�:	���ʕ��������܂���
//		0:				���ʕ������Ȃ�����
//		-1:				���ʃG�b�W�̎����v�Z�G���[
//		-2:				���ʃG�b�W���S�ȏ㋁�܂���(�̂��Ă��Ȃ��ƌ��Ȃ�)
//�ǋL�F
//	�Ȑ������ʂ��ǂ����̌덷�ɂ�line_zero()���A�p�����[�^�͈͂̎����v�Z�̌덷�ɂ́A
//  �p�����[�^�͈�*rc_zero()���g�p����
int MGCompositeCurve::common(
	const MGCurve& crv2,
	std::vector<double>& vecComSpan
)const{
	const_iterator i=begin(), ie=end();
	for(; i!=ie; i++){
		std::vector<double> vecComSpan2;
		int obtained2=(**i).common(crv2,vecComSpan2);
		if(obtained2<0)
			return obtained2;
		size_t ncom=vecComSpan2.size();
		for(size_t j=0; j<ncom; j++)
			vecComSpan.push_back(vecComSpan2[j]);

	}
	return (int)(vecComSpan.size()/4);
}

//���I�t�Z�b�g�֐�
//�I�t�Z�b�g�����́A�@���������猩�ē��͋Ȑ��̐i�s���������𐳂Ƃ���B
//�@���x�N�g�����k���̏ꍇ�A�n�_�ɂ����ċȗ����S�����𐳂Ƃ���B�������A�ȗ����S�֋ȗ����a�ȏ�̃I�t�Z�b�g
//�͍s��Ȃ��B�g�������X��line_zero()���g�p���Ă���B�߂�l�́A�I�t�Z�b�g�Ȑ����X�g���ԋp�����B
//costant offset curve. if the norm_vector is given, the positive offset direction decide
//to left hand side from ahead, or the direction to center of curvature at start parameter.
//the offset value is less than radius of curvature. line_zero() is used.
std::vector<UniqueCurve> MGCompositeCurve::offset(
	double ofs_value,			//�I�t�Z�b�g��
	const MGVector& norm_vector	//�@���x�N�g��
)const{
	std::vector<UniqueCurve> retval(1);
	retval[0].reset(static_cast<MGCurve*>(0));
	return retval;
}

//�σI�t�Z�b�g�֐�
//�I�t�Z�b�g�ʂ͋�Ԏ���1�̐�B�\���ŗ^������B
//�I�t�Z�b�g�����́A�@���������猩�ē��͋Ȑ��̐i�s���������𐳂Ƃ���B
//�@���x�N�g�����k���̏ꍇ�A�n�_�ɂ����ċȗ����S�����𐳂Ƃ���B�������A�ȗ����S�֋ȗ����a�ȏ�̃I�t�Z�b�g
//�͍s��Ȃ��B�g�������X��line_zero()���g�p���Ă���B�߂�l�́A�I�t�Z�b�g�Ȑ����X�g���ԋp�����B
//valuable offset curve. if the norm_vector is given, the positive offset direction decide
//to left hand side from ahead, or the direction to center of curvature at start parameter.
//the offset value is less than radius of curvature. line_zero() is used.
std::vector<UniqueCurve> MGCompositeCurve::offset(
	const MGLBRep& ofs_value_lb,					//��Ԏ����P�̐�B�\���Ŏ������I�t�Z�b�g��
	const MGVector& norm_vector	//�@���x�N�g��
)const{
	std::vector<UniqueCurve> retval(1);
	retval[0].reset(static_cast<MGCurve*>(0));
	return retval;
}

//C2�A���Ȑ��̈��I�t�Z�b�g�֐�
//�I�t�Z�b�g�����́A�@���������猩�ē��͋Ȑ��̐i�s���������𐳂Ƃ���B
//�@���x�N�g�����k���̏ꍇ�A�n�_�ɂ����ċȗ����S�����𐳂Ƃ���B�������A�ȗ����S�֋ȗ����a�ȏ�̃I�t�Z�b�g
//�͍s��Ȃ��B�g�������X��line_zero()���g�p���Ă���B�߂�l�́A�I�t�Z�b�g�Ȑ����ԋp�����B
//costant offset curve of C2 continuous curve. if the norm_vector is given, the positive offset direction
//decide to left hand side from ahead, or the direction to center of curvature at start parameter.
//the offset value is less than radius of curvature. line_zero() is used.
MGLBRep MGCompositeCurve::offset_c2(
	double ofs_value,			//�I�t�Z�b�g��
	const MGVector& norm_vector	//�@���x�N�g��
)const{
	return MGLBRep();
}

//C2�A���Ȑ��̉σI�t�Z�b�g�֐�
//�I�t�Z�b�g�ʂ͋�Ԏ���1�̐�B�\���ŗ^������B
//�I�t�Z�b�g�����́A�@���������猩�ē��͋Ȑ��̐i�s���������𐳂Ƃ���B
//�@���x�N�g�����k���̏ꍇ�A�n�_�ɂ����ċȗ����S�����𐳂Ƃ���B�������A�ȗ����S�֋ȗ����a�ȏ�̃I�t�Z�b�g
//�͍s��Ȃ��B�g�������X��line_zero()���g�p���Ă���B�߂�l�́A�I�t�Z�b�g�Ȑ����ԋp�����B
//valuable offset curveof C2 continuous curve. if the norm_vector is given, the positive offset direction
//decide to left hand side from ahead, or the direction to center of curvature at start parameter.
//the offset value is less than radius of curvature. line_zero() is used.
MGLBRep MGCompositeCurve::offset_c2(
	const MGLBRep& ofs_value_lb,					//��Ԏ����P�̐�B�\���Ŏ������I�t�Z�b�g��
	const MGVector& norm_vector	//�@���x�N�g��
)const{
	return MGLBRep();
}

//Return sweep surface from crv
//Returned is a newed MGSurface, must be deleted.
MGSurface* MGCompositeCurve::sweep(
	const MGUnit_vector& uvec,			//Sweep Direction.
	double start_dist,					//distance to start edge.
	double end_dist			//distance to end edge.
)const{
	return new MGPlane();
}

//Obtain the projected curve of this onto the surface srf.
//The direction of the projection is along the vector vec if the vec is not
//NULL, and normal to the surface if the vec is NULL.
//Output of 'project' is two kind of curves:
//one is general world coordinate curves('vec_crv'), and the other is (u,v) curves of
//the parameter space of the surfaces(vec_crv_uv).
//vec_crv_uv.size() is equal to vec_crv.size(). Let the size be n, then
// (vec_crv_uv[i], vec_crv[i]) is one pair for 0<=i<n.
size_t MGCompositeCurve::project_onto_surface(
	const MGFSurface& srf,
	std::vector<UniqueCurve>& vec_crv_uv,
		//Projected curve(surface parameter (u,v) representation) will be appended.
	std::vector<UniqueCurve>& vec_crv,
		//Projected curve(world coordinate(x,y,z) representation) will be appended.
	const MGVector& vec
)const{
	const_iterator i=begin(), ie=end();
	for(; i!=ie; i++)
		srf.project(**i,vec_crv_uv,vec_crv,vec);
	return vec_crv.size();
}

///Approximate this curve by a polyline and output to lb2.
///The tolerance of the approximation is error.
void MGCompositeCurve::polygonize(
	double error,	///<tolerance allowed for the approximation
	MGLBRep& lb2	///<Obtained polyline will be output as an MGLBRep of order2.
)const{
	const_iterator i=begin(), ie=end();
	if(i!=ie)
		(**i++).polygonize(error,lb2);
	for(; i!=ie; i++){
		MGLBRep lbi;
		(**i).polygonize(error,lbi);
		lb2.connect(0,2,lbi);
	}
}

///Approximate this curve as a MGLBRep curve
///within the tolerance MGTolerance::line_zero().
///When parameter_normalization=0, reparameterization will not done, and
///the evaluation at the same parameter has the same values before and after
///of approximate_as_LBRep.
void MGCompositeCurve::approximate_as_LBRep(
	MGLBRep& lb,		///<Approximated obrep will be set.
	int ordr,		///<new order
	int parameter_normalization,
		//Indicates how the parameter normalization be done:
		//=0: no parameter normalization.
		//=1: normalize to range=(0., 1.);
		//=2: normalize to make the average length of the 1st derivative 
		//    is as equal to 1. as possible.
	bool neglectMulti///<Indicates if multiple knots be kept.
		///< true: multiplicity is removed.
		///< false: multiplicity is kept.
)const{
	const_iterator i=begin(), ie=end();
	if(i!=ie)
		(**i++).approximate_as_LBRep(lb,ordr,parameter_normalization,neglectMulti);
	for(; i!=ie; i++){
		MGLBRep lbi;
		(**i).approximate_as_LBRep(lbi,ordr,parameter_normalization,neglectMulti);
		int which;
		double ratio;
		int cntnty=lb.continuity(lbi,which,ratio);
		lb.connect(cntnty,2,lbi);
	}
}
