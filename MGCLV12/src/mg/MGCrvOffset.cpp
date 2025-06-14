/********************************************************************/
/* Copyright (c) 2023 System fugen G.K. and Yuzi Mizuno             */
/* All rights reserved.                                             */
/********************************************************************/
#include "StdAfx.h"
#include "mg/Tolerance.h"
#include "mg/Box.h"
#include "mg/BPointSeq.h"
#include "mg/Curve.h"
#include "mg/Ellipse.h"
#include "mg/LBRep.h"
#include "mg/RLBRep.h"
#include "mg/Straight.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

// Implementation of curve offset.


///�Ȑ����I�t�Z�b�g����̂ɏ\�����������m�b�g�x�N�g����ԋp����
///�I�t�Z�b�g�ʋȐ����l���ɓ���A�������̑������ɂ��킹�Ă���B
MGKnotVector offset_make_knotvector(const MGCurve* lb1, const MGLBRep* lb2){
	//���ƂȂ�m�b�g�x�N�g���𐶐�����
	const MGKnotVector& t1 = lb1->knot_vector();
	MGKnotVector tNew(t1, 4);	//�I�[�_�[4�̃m�b�g�x�N�g���ɍ��ւ���
	double error = lb1->param_error();
	for (int i = t1.order() - 1; i < t1.bdim(); i++) {
		double spara = t1(i), epara = t1(i + 1);	//�X�p���̎n�I�I�_
		if (epara - spara < error)
			continue;	//�}���`�m�b�g�̂Ƃ��̏���(RLBRep�̂�)

		//1�X�p���̕����������肷��
		MGInterval interval(spara, epara);
		int ndiv = lb1->divideNum(interval);//�I�t�Z�b�g�Ȑ��̕�����
		int tmp_ndiv = lb2 ?  lb2->divideNum(interval) : -1;//�I�t�Z�b�g�ʋȐ��̕�����
		if (tmp_ndiv > ndiv)
			ndiv = tmp_ndiv;//�������̑�������p����

		int id=tNew.locate(spara)+1;
		int ndivm1 = ndiv - 1;
		tNew.reserveInsertArea(ndivm1, id);

		double delta = (epara - spara) / ndiv;
		double tpara = spara + delta;
		for (int j = 0; j < ndivm1; j++, tpara += delta)
			tNew(id + j) = tpara;
	}
	return tNew;
}

///�Ȑ� crv��܂�(C0 continuous)�ŕ�������
size_t divide_multi_ofs(
	const MGCurve& crv,
	std::vector<UniqueCurve>& crv_list ///divided curves are output.
) {
	const MGKnotVector& t = crv.knot_vector();
	int k = t.order();
	int km1 = k - 1;
	int	start_index = km1, index = 0, multi = 1, vbdim = t.bdim();
	do {
		if (km1 <= 1)
			index = start_index + 1;
		else
			multi = t.locate_multi(start_index, km1, index);
		crv_list.emplace_back(crv.part(t(start_index), t(index)));
		start_index = index + multi - 1;
	} while (index < vbdim);	//���d�x��������Ȃ�������I���
	return crv_list.size();
}

//�I�t�Z�b�g�ʌŒ��C1�A���Ȑ��̃I�t�Z�b�g
//Curve offset funciont when offset value is constant.
//original curve is supposed to be C1 continuous everywhere.
std::unique_ptr<MGLBRep> C1CurveConstantOffset(
	const MGCurve& original,
	double  ofs_value,//�Œ�I�t�Z�b�g��
	bool pricipalNormal//true if offset direction is toward principal normal,
						//false if to binormal.
){
	MGKnotVector t = offset_make_knotvector(&original, nullptr);//�\�����������m�b�g�x�N�g�������߂�
	MGNDDArray tau; tau.buildByKnotVector(t);

	//����_�𐶐�����
	MGVector T, N, B, Nold, Bold;
	double curvature, torsion;
	original.Frenet_frame(tau(0), T, Nold, Bold, curvature, torsion);
		//save the 1st N, or B data.
	MGVector& dirOld = pricipalNormal ? Nold : Bold;
	MGVector& dir = pricipalNormal ? N : B;

	int n = tau.length();
	MGBPointSeq bp1(n, original.sdim());
	for (int i = 0; i < n; i++) {
		double taui = tau(i);
		original.Frenet_frame(taui, T, N, B, curvature, torsion);
		if (dir % dirOld < 0.)
			dir.negate();
		MGPosition pos = original.eval(taui) + dir * ofs_value;
		bp1.store_at(i, pos);
		dirOld = dir;
	}
	//�m�b�g�x�N�g�������߂ăI�t�Z�b�g�Ȑ��𐶐�����(���x�\���̋Ȑ��𐶐�����)
	std::unique_ptr<MGLBRep> newLB(new MGLBRep);
	if(newLB->buildByInterpolationDataPoints(tau, bp1))
		return nullptr;//if error detected.
	newLB->remove_knot();
	return newLB;
}

//�I�t�Z�b�g�ʉς�C1�A���Ȑ��̃I�t�Z�b�g
//Curve offset funciont when offset value is variable.
//original curve is supposed to be C1 continuous everywhere.
UniqueLBRep C1CurveVariableOffset(
	const MGCurve& crv,
	const MGLBRep& ofs_value_lb,	//Variable offset value input(LBRep of space dim =1)
	bool pricipalNormal//true if offset direction is principal normal,
					//false if binormal.
){
	MGKnotVector knotVector = offset_make_knotvector(&crv, &ofs_value_lb);	//�\�����������m�b�g�x�N�g�������߂�
	MGNDDArray tau;
	tau.buildByKnotVector(knotVector);

	MGVector T, N, B, Nold, Bold;
	double curvature, torsion;
	crv.Frenet_frame(tau(0), T, Nold, Bold, curvature, torsion);
	MGVector& dirOld = pricipalNormal ? Nold : Bold;
	MGVector& dir = pricipalNormal ? N : B;

	//Generate data of interpolation at tau.
	int n = tau.length();
	MGBPointSeq bp1(n, crv.sdim());
	for (int i = 0; i < n; i++) {
		double taui = tau(i);
		double delta = ofs_value_lb.eval_position(taui).ref(0);//offset value at tau i.
		crv.Frenet_frame(taui, T, N, B, curvature, torsion);
		if (dir % dirOld < 0.)
			dir.negate();

		MGPosition pos = crv.eval(taui) + dir * delta;
		bp1.store_at(i, pos);
		dirOld = dir;
	}

	//�m�b�g�x�N�g�������߂ăI�t�Z�b�g�Ȑ��𐶐�����
	UniqueLBRep offsetCrv(new MGLBRep);
	offsetCrv->setKnotVector(std::move(knotVector));
	if(offsetCrv->buildByInterpolationWithKTV(tau, bp1))	//���x�\���̋Ȑ��𐶐�����
		return nullptr;//When error detected in buildByInterpolationWithKTV().
	offsetCrv->remove_knot();
	return offsetCrv;
}

///Offset of costant deviation from this curve that is C1 cotinuous everywhere.
///The offset value ofs_value must be less than radius of curvature.
///The offset direction is Normal(obtained by Frenet_frame()) at each point.
UniqueLBRep MGCurve::offsetC1(
	double ofs_value, ///offset valuem may be negative.
	bool principalNormal/// true: Offset direction is to principal normal
						/// false: to binormal
) const {
	return C1CurveConstantOffset(*this, ofs_value, principalNormal);
}

//���I�t�Z�b�g�֐�
//�I�t�Z�b�g�����́A�@���������猩�ē��͋Ȑ��̐i�s���������𐳂Ƃ���B
//�@���x�N�g�����k���̏ꍇ�A�n�_�ɂ����ċȗ����S�����𐳂Ƃ���B
//�������A�ȗ����S�֋ȗ����a�ȏ�̃I�t�Z�b�g�͍s��Ȃ��B�g�������X��line_zero()���g�p���Ă���B
//�߂�l�́A�I�t�Z�b�g�Ȑ����X�g���ԋp�����B
std::vector<UniqueCurve> MGCurve::offset(
	double ofs_value,
	bool principalNormal/// true: Offset direction is to principal normal
						/// false: to binormal
)const{
	//�Ȑ���܂�(C0 continuity)�ŕ�������
	std::vector<UniqueCurve> crv_list, offseted;
	divide_multi_ofs(*this,crv_list);
	for (auto& crv : crv_list) {
		UniqueCurve crvOffset = C1CurveConstantOffset(*crv, ofs_value, principalNormal);
		if(crvOffset)
			offseted.push_back(std::move(crvOffset));
	}	
	return offseted;
}

//�σI�t�Z�b�g�֐�
//�I�t�Z�b�g�ʂ͋�Ԏ���1�̐�B�\���ŗ^������B
//�I�t�Z�b�g�����́A�@���������猩�ē��͋Ȑ��̐i�s���������𐳂Ƃ���B
//�@���x�N�g�����k���̏ꍇ�A�n�_�ɂ����ċȗ����S�����𐳂Ƃ���B
//�������A�ȗ����S�֋ȗ����a�ȏ�̃I�t�Z�b�g�͍s��Ȃ��B�g�������X��line_zero()���g�p���Ă���B
//�߂�l�́A�I�t�Z�b�g�Ȑ����X�g���ԋp�����B
std::vector<UniqueCurve> MGCurve::offset(
	const MGLBRep& ofs_value_lb,	//��Ԏ����P�̐�B�\���Ŏ������I�t�Z�b�g��
	bool principalNormal/// true: Offset direction is to principal normal
						/// false: to binormal
)const{
	//�Ȑ���܂�imultiple knot�j�ŕ�������
	std::vector<UniqueCurve> crv_list;
	divide_multi_ofs(*this, crv_list);

	std::vector<UniqueCurve> tmp_list;
	for (auto& crv : crv_list) {
		UniqueCurve crvOffset = C1CurveVariableOffset(*crv, ofs_value_lb, principalNormal);
		if (crvOffset)
			tmp_list.push_back(std::move(crvOffset));
	}

	return tmp_list;
}

/// Offset of constant deviation from this curve.
/// The offset value must be less than radius of curvature.
/// When this curve is not C1 continuous, this is divided into C1 curves,
/// and more than one offset curves are obtained.
/// line_zero() is used to approximate curves of the offset.
std::vector<UniqueCurve> MGRLBRep::offset(
	double ofs_value,
	bool principalNormal/// true: Offset direction is to principal normal
						/// false: to binormal
) const {
	std::vector<UniqueCurve> curves;
	curves.emplace_back(C1CurveConstantOffset(*this, ofs_value, principalNormal).release());
	return curves;
}

/// Offset of variable deviation from this curve.
/// When this curve is not C1 continuous, divided into C1 curves,
/// and more than one offset curves are obtained.
/// The direction of offset is toward the principal normal,
/// or to the direction to center of curvature.
/// line_zero() is used approximate the offset curve.
std::vector<UniqueCurve> MGRLBRep::offset(
	const MGLBRep& ofs_value_lb,			///<��Ԏ����P�̐�B�\���Ŏ������I�t�Z�b�g��
	bool principalNormal/// true: Offset direction is to principal normal
						/// false: to binormal
) const {
	std::vector<UniqueCurve> curves;
	curves.emplace_back(C1CurveVariableOffset(*this, ofs_value_lb, principalNormal).release());
	return curves;

}
