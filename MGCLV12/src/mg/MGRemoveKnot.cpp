/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#include "StdAfx.h"
#include "mg/Box.h"
#include "mg/BPointSeq.h"
#include "mg/Position.h"
#include "mg/LBRep.h"
#include "mg/RLBRep.h"
#include "mg/Straight.h"
#include "mg/SPointSeq.h"
#include "mg/SBRep.h"
#include "mg/RSBRep.h"
#include "mg/Tolerance.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

using namespace std;

// Implementation of knot remove.

//�m�b�g�폜�֐�(B�\���Ȑ��̂�)
//�g�������X��line_zero���g�p����B���̃m�b�g���ׂ������̂قǍ폜���₷��
//removal knot. line_zero tolerance is used.
void MGCurve::remove_knot(){;}

#define START_DELNUM 4000
//Remove knot if removed line has the difference less than line_zero();
//The difference is checked only for the space id of coef(.,j+k)
//for k=0, ..., nd-1.
void MGLBRep::remove_knot(int j, int nd){
	int k=order();
	int km1 = k-1;	//������
	int i, n=bdim();
	int mid=(km1+n)/2;
	if(mid<=km1) mid=k;
	else if(mid>START_DELNUM){mid=START_DELNUM;}
	double totalTol = 0.0;
	const double line0=MGTolerance::line_zero();

	int num_remained;
	for(i=n-1; i>=mid;){
		int ndel = remove_knot_one(line0,i,totalTol,num_remained,j,nd);
		i-=int(ndel+num_remained);
	}

	totalTol = 0.0;
	i=k;
//	int nt=int(bdim())-mid;
	int nt=int(bdim())-mid-1;if(nt<0) nt=0;
	while(int(bdim())-i>nt){
		remove_knot_one(line0,i,totalTol,num_remained,j,nd);
		i+=num_remained;
	}
}

//�m�b�g�폜�֐�(1�̃m�b�g)
//�߂�l�́A�폜�����m�b�g�̐�
//When snum!=0, tolerance of totalTol is checked only for coef(.,sid+j),
//of j=0, ..., snum-1. When snum=0, snum is set as sdim();
int MGLBRep::remove_knot_one(
	double line0,		//Tolerance allowed for the knot removal.
	int	id,			//�폜���悤�Ƃ���m�b�g�̔ԍ�
	double& totalTol,	//�덷���v
	int& num_knot,	//Remained knot number at knot(id) after removed.
	int sid,			//Space dimension start id of this LBRep's B-coef.
	int snum			//Num of space dimension for the totalTol tolerance check.
){
	int sd=sdim();
	if(!snum)
		snum=sd;
	MGPosition P(snum), Q(snum), R(snum);
	const MGBPointSeq& bcoef=line_bcoef();

	const int k = order();	//�I�[�_�[
	const int km1 = k - 1;	//����
	const int n=bdim();	//B�\������
	const double tau=knot(id);

	//Get the multiplicity of knot(id).
	//We need the end knot id of the multiplicity in 'id'.
	int m, idold=id;
	int multi=1;
	for(m=id+1; m<n && tau==knot(m); m++){
		multi++; id++;
	}
	for(m=idold-1; m>=k && tau==knot(m); m--)
		multi++;
	num_knot=multi;

	//const double line0=MGTolerance::line_zero();
	MGBPointSeq temp_bp(2 * km1 + 1, sd);

	int	last = id - multi;
	int first = id - km1;
	int ndel;				//�폜�����m�b�g�̐�
	if(k == 2){
	//�I�[�_�[2�̂Ƃ��̏���
		const int off = first - 1;
		bcoef.point(off,sid,snum,P);
		bcoef.point(last+1,sid,snum,Q);
		MGStraight st1(P, Q);
		bcoef.point(first,sid,snum,R);
		double tmpTol = st1.distance(R);
		totalTol += tmpTol;
		if(line0>0. && totalTol>line0)
			ndel=0;
			//�덷���v���g�������X�ȏ�̏ꍇ�m�b�g�폜���s��Ȃ�
		else ndel=1;
	}else{

	//�I�[�_�[3�ȏ�̂Ƃ��̏���
	for(ndel = 0; ndel<multi; ndel++){
		const int off = first - 1;
		temp_bp.store_at(0, coef(off));
		temp_bp.store_at(last + 1 - off, coef(last + 1));
		int i = first;	int j = last;
		int ii = 1;		int jj = last - off;
		while((j-i) > ndel){
			double ti=knot(i), tjmndel=knot(j - ndel);
			double alfi = (tau-ti) / (knot(i+k+ndel) - ti);
			double alfj = (tau - tjmndel) / (knot(j+k) - tjmndel);
			temp_bp.store_at(ii, (coef(i) - (1.0 - alfi) * temp_bp(ii-1)) / alfi);
			temp_bp.store_at(jj, (coef(j) - alfj * temp_bp(jj+1)) / (1.0 - alfj));
			i++;	j--;	ii++;	jj--;
		}

		if((j-i) < ndel){
			temp_bp.point(ii - 1,sid,snum,P);
			temp_bp.point(jj + 1,sid,snum,Q);
			totalTol += (P - Q).len();
		}else{
			double ti=knot(i);
			double alfi = (tau - ti) / (knot(i+k+ndel) - ti);
			temp_bp.point(ii+ndel+1,sid,snum,P);
			temp_bp.point(ii-1,sid,snum,Q);
			bcoef.point(i,sid,snum,R);
			totalTol += (R - (alfi*P + (1.0 - alfi)*Q)).len();
		}
		//�덷���v���g�������X�ȏ�̏ꍇ�m�b�g�폜���s��Ȃ�
		if(line0>0. && totalTol > line0) break;

		i = first;	j = last;
		while((j-i) > ndel){
			m_line_bcoef.store_at(i, temp_bp(i - off));
			m_line_bcoef.store_at(j, temp_bp(j - off));
			i++;	j--;
		}
		first--; last++;
	}

	}

	if(!ndel){totalTol=0.; return 0;}
	int l = id + 1;
	const int num_end_knot = n + km1;	//�ŏI�m�b�g�̔ԍ�
	for(; l <= num_end_knot; l++) knot(l-ndel) = knot(l);
	int j = (2*id - multi - km1)/2;//�ŏ��̃R���g���[���|�C���g id
	int i = j;
	for(l = 1; l < ndel; l++){
		if((l % 2) == 1)i++; else j--;
	}
	for(l = i+1; l<n; l++, j++) m_line_bcoef.store_at(j, coef(l));

	//�Ȑ����X�V����
	int newnbd=n - ndel;
	m_knot_vector.set_bdim(newnbd);
	m_line_bcoef.set_length(newnbd);
	num_knot=multi - ndel;
	if(num_knot) totalTol = 0.0;//�m�b�g�폜���s���Ȃ��Ƃ��덷���v���N���A����
	return ndel;
}

//�m�b�g�폜�֐�
void MGRLBRep::remove_knot(){
	double Pmax = 0.0, Wmin = 1.0;
	int sd=sdim();
	for(int i = 0;i < bdim(); i++){
		MGPosition pos = eval_position(knot(i));
		double P = pos.len();
		double W = coef(i, sd);		//�ŏI�������d�݂ł���
		if(P > Pmax) Pmax = P;
		if(W < Wmin) Wmin = W;
	}
	double save=MGTolerance::line_zero();
	double tol =  save* Wmin /(1 + Pmax);
	mgTolSetLineZero lineZeroSet(tol);
	m_line.remove_knot();
}

//�m�b�g�폜�֐�
//�g�������X��line_zero���g�p����B���̃m�b�g���ׂ������̂قǍ폜���₷��
void MGSBRep::remove_knot(){
	remove_knot_u();
	exchange_uv();
	remove_knot_u();
	exchange_uv();
}

//u�m�b�g���폜����
int MGSBRep::remove_knot_u(){
	int k = order_u(), n = bdim_u();//order and bdim along u.
	int km1=k-1;
	int mid=(km1+n)/2;
	if(mid<=km1)
		mid=k;
	int num_remained;

	int nDel = 0;		//�폜�����m�b�g�̐�
	double totalTol = 0.0;
	double line0=MGTolerance::line_zero();
	int i;
	for(i=n-1; i>=mid;){
		int nDelOne = remove_knot_u_one(line0,i, totalTol, num_remained);
		i-=(nDelOne+num_remained);
		nDel+=nDelOne;
	}
	i=k;
	int nt=bdim_u()-mid-1;
	if(nt<0) nt=0;

	totalTol=0.;
	while(int(bdim_u())-i>nt){
		int nDelOne=remove_knot_u_one(line0,i,totalTol,num_remained);
		i+=num_remained;
		nDel+=nDelOne;
	}
	return nDel;
}

//u�m�b�g���폜����
//�֐��̖߂�l�͍폜�����m�b�g�̐�
int MGSBRep::remove_knot_u_one(
	double line0,		//Tolerance allowed for the knot removal. 
						//When line0<=0., removal will be done uncoditionally.
	int	id,			//�폜���悤�Ƃ���m�b�g�̔ԍ�
	double& totalTol,	//�ŏ��͂O����͂���B���Ƃ�remove_knot_u_one���Ǘ�����B
	//�폜���J��Ԃ��Ă���Ƃ��͌덷�����Z�����B����ȏ�폜�ł��Ȃ����O�N���A�����B
	int& num_knot	//Remained knot number at knot(id) after removed.
){
	const int k = order_u();	//�I�[�_�[
	const int km1 = k - 1;		//����
	const int m = bdim_u(), n = bdim_v();
	const int sd=sdim();
	const double tau=knot_u(id);

	//Get the multiplicity of knot(id).
	//We need the end knot id of the multiplicity in 'id'.
	int a, multi=1, idold=id;
	for(a=id+1; a<n && tau==knot_u(a); a++){
		multi++; id++;
	}
	for(a=idold-1; a>=k && tau==knot_u(a); a--)
		multi++;
	num_knot=multi;

	int	last = id - multi;
	int first = id - km1;

	int l;
	int ndel=0;			//�폜�����m�b�g�̐�
	if(k==2){

//�I�[�_�[2�̂Ƃ��̏���
	const int off = first - 1;
	double maxTol = 0.0;
	int l;
	for(l=0; l<n; l++){
		MGStraight st1(coef(off, l), coef(last + 1, l));
		double tmpTol = st1.distance(coef(first, l));
		if(line0>0. && tmpTol>line0){ndel= 0; break;}//�m�b�g���폜�ł��Ȃ�����
		if(tmpTol>maxTol) maxTol = tmpTol;
	}
	if(l>=n){//When possible to delete knot.
		ndel=1;totalTol += maxTol;
	}

	}else{
//�I�[�_�[3�ȏ�̂Ƃ��̏���
	for(ndel = 0; ndel < multi; ndel++){

	MGSPointSeq temp(2*km1+1,n,sd);
	const int off = first - 1;
	int i, j, ii, jj;
	double maxTol = 0.0;
	for(l=0; l<n; l++){
		i = first;	j = last;
		ii = 1;		jj = last-off;

		temp.store_at(0, l, coef(off,l));
		temp.store_at(last+1-off,l,coef(last+1,l));
		double alfi, alfj, tmpTol = 0.0;
		while((j-i) > ndel){
			double ti=knot_u(i), tjmndel=knot_u(j-ndel);
			alfi = (tau-ti) / (knot_u(i+k+ndel) - ti);
			alfj = (tau-tjmndel) / (knot_u(j+k) - tjmndel);
			temp.store_at(ii,l,(coef(i,l) - (1.0 - alfi)*temp(ii-1,l)) / alfi);
			temp.store_at(jj, l, (coef(j,l) - alfj * temp(jj+1,l)) / (1.-alfj));
			i++;	j--;	ii++;	jj--;
		}
		if((j-i) < ndel){
			tmpTol = (temp(ii-1,l) - temp(jj+1,l)).len();
			if(line0>0. && tmpTol>line0){ maxTol=tmpTol; break;}
		}else{
			double ti=knot_u(i);
			alfi = (tau-ti) / (knot_u(i+k+ndel)-ti);
			tmpTol=(coef(i,l)-(alfi * temp(ii+ndel+1,l)+(1.-alfi)*temp(ii-1,l))).len();
			if(line0>0. && tmpTol>line0){ maxTol=tmpTol; break;}
		}
		if(tmpTol > maxTol) maxTol=tmpTol;
	}
	totalTol += maxTol;
	if(line0>0. && totalTol>line0) break;
	i = first;	j = last;
	while((j-i) > ndel){
		for(int vk=0; vk<n; vk++){
			m_surface_bcoef.store_at(i,vk,temp(i-off,vk));
			m_surface_bcoef.store_at(j,vk,temp(j-off,vk));
		}
		i++; j--;
	}
	first--;	last++;

	}
	}

	if(!ndel){totalTol=0.; return 0;}
	const int num_end_uknot = m + km1;		//�ŏIu�m�b�g�̔ԍ�
	for(l=id+1; l<=num_end_uknot; l++) knot_u(l-ndel)=knot_u(l);
	int j=(2*id-multi-km1)/2;	//�ŏ��̃R���g���[���|�C���g
	int i=j;
	for(l=1; l<ndel; l++){
		if((l % 2) == 1) i++; else j--;
	}
	for(l=i+1; l<m; l++){
		for(int vk = 0; vk < n; vk++) m_surface_bcoef.store_at(j,vk,coef(l,vk));
		j++;
	}
	//�Ȗʂ��X�V����
	int newnbd=m-ndel;
	m_uknot.set_bdim(newnbd);
	m_surface_bcoef.set_length(newnbd, n);
	num_knot=multi - ndel;
	if(num_knot) totalTol = 0.0;//�m�b�g�폜���s���Ȃ��Ƃ��덷���v���N���A����
	return ndel;
}

//u�m�b�g����폜����
//�֐��̖߂�l�͍폜�����m�b�g�̐�
int MGSBRep::remove_knot_one(
	double line0,		//Tolerance allowed for the knot removal. 
						//When line0<=0., removal will be done uncoditionally.
	int	id,	//�폜���悤�Ƃ���m�b�g�̔ԍ�
	double& tol,//�폜��̌덷���o�͂����
	bool u_knot	//�폜�Ώۂ��iu,v)�̂������knot vector������͂���
				//=true�̂Ƃ��Au-knot_vector���폜
){
	int num_knot;
	int delnum;
	if(u_knot) delnum=remove_knot_u_one(line0,id,tol,num_knot);
	else{
		exchange_uv();
		delnum=remove_knot_u_one(line0,id,tol,num_knot);
		exchange_uv();
	}
	return delnum;
}

//�m�b�g�폜�֐�
//�g�������X��line_zero���g�p����B���̃m�b�g���ׂ������̂قǍ폜���₷��
void MGRSBRep::remove_knot(){
	int i = 0, j = 0;
	double Pmax = 0.0, Wmin = 1.0;
	for(;i < bdim_u(); i++){
		for(j = 0; j < bdim_v(); j++){
			MGPosition pos = eval(knot_u(i), knot_v(j), 0);
			double P = pos.len();
			double W = coef(i, j, sdim());		//�ŏI�������d�݂ł���
			if(P>Pmax)
				Pmax = P;
			if(W<Wmin)
				Wmin = W;
		}
	}
	double save=MGTolerance::line_zero();
	double tol =  save* Wmin /(1 + Pmax);
	mgTolSetLineZero lineZeroSet(tol);
	m_surface.remove_knot();
}

#define INCNUM 20
///Rebuild this NURBS by reconstructing new knot configuration.
std::unique_ptr<MGRLBRep> MGRLBRep::rebuild_with_new_knot_configuration(
	double error,	//Error alowed to rebuild. If error<=0., MGTolerance::line_zero()
	                //will be employed.
	int parameter_normalization
		//Indicates how the parameter normalization be done:
		//=0: no parameter normalization.
		//=1: normalize to range=(0., 1.);
		//=2: normalize to make the average length of the 1st derivative 
		//    is as equal to 1. as possible.
)const{
	if(error<=0.)
		error=MGTolerance::line_zero();
	mgTolSetLineZero lineZeroSet(error);

	//����_�𐶐�����
	std::unique_ptr<MGRLBRep> rlb(new MGRLBRep);
	MGKnotVector& t=rlb->knot_vector();
	t=knot_vector();
	int len = t.length();
	t.change_number(len*INCNUM);
	rlb->buildByNewKnotVectorWithKTV(*this);
	rlb->remove_knot();

	if(parameter_normalization){
		rlb->change_range(0.,1.);
		if(parameter_normalization==2){
			MGVector v0=rlb->eval(0.,1);
			MGVector v1=rlb->eval(0.5,1);
			MGVector v2=rlb->eval(1.,1);
			double t2=(v0.len()+v1.len()+v2.len())/3.;
			rlb->change_range(0.,t2);
		}
	}
	return rlb;
}
