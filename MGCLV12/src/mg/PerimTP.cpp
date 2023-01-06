/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/

/** @file PerimTP.cpp
 *  Implementation for class mgPerimTP.
 */
#include "stdafx.h"
#include "mg/Tolerance.h"
#include "mg/FSurface.h"
#include "mg/TPOnFaceBoundary.h"
#include "mg/TPmaker.h"

namespace {
	// tp曲線のデータポイントサンプリング
	void sampling(
		const MGLBRep& tp,
		MGNDDArray& tau,//New sampling data points are appended.
		MGBPointSeq& ordinates//New sampling data are appended.
	){
		MGNDDArray tau2;
		tau2.buildByKnotVector(tp.knot_vector());
		int nold = tau2.length();
		tau2.change_number(nold * 8);

		int ntau = tau.length();
		int idStart= ntau ? tau2.locate(tau[ntau-1])+2 : 0;

		int ntau2 = tau2.length();
		int nNew = ntau + ntau2 - idStart;

		// Sampling(We reserve (nNew+1) area for the process after sampling).
		tau.reshape(nNew+1); tau.set_length(nNew);
		ordinates.reshape(nNew+1); ordinates.set_length(nNew);
		for (int j = ntau, i = idStart; j < nNew; j++, i++) {
			double tau2i= tau2(i);
			tau(j) = tau2i;
			ordinates.store_at(j, tp.eval(tau2i));
		}		
	}
}

mgPerimTP::mgPerimTP():m_perimeter(nullptr), m_tpmaker(nullptr), m_perimNum(-1){;}
void mgPerimTP::build(mgTPmaker* qsmain, int periNum, const MGLBRep * perim){
	m_tpmaker = qsmain;
	m_perimNum = periNum;
	m_perimeter = perim;
}

std::ostream& operator<<(std::ostream& os, const mgPerimTP& perimeter){
	os <<std::endl<<"mgPerimTP::";
	perimeter.toString(os);
	return os;
}
/// Returns the number of objects of class mgTPOnFaceBoundary this object holds.
size_t mgPerimTP::number_of_subtp() const { return m_tpVec.size(); }

//! Access to the first %mgTPOnFaceBoundary object in the perimeter.
mgPerimTP::MYELM& mgPerimTP::subtp_first(){ return m_tpVec.front(); }
const mgPerimTP::MYELM& mgPerimTP::subtp_first() const { return m_tpVec.front(); }

//! Access to the last %mgTPOnFaceBoundary object in the perimeter.
const mgPerimTP::MYELM& mgPerimTP::subtp_last() const { return m_tpVec.back(); }

const MGVector& mgPerimTP::unit_normal_start() const{
	if (m_perimNum<=1) {
		return m_tpmaker->m_cornerNormal[m_perimNum];
	}else{
		int next = (m_perimNum + 1) % 4;
		return m_tpmaker->m_cornerNormal[next];
	}
}

const MGVector& mgPerimTP::unit_normal_end() const {
	if (m_perimNum <= 1) {
		int next = m_perimNum + 1;
		return m_tpmaker->m_cornerNormal[next];
	}else{
		return m_tpmaker->m_cornerNormal[m_perimNum];
	}
}

bool mgPerimTP::need_negate(double cmnargs, double cmnarge) const{
	bool neg = m_perimNum<=1 ? false:true; // face's direction.
	if (cmnargs<cmnarge)
		neg = !neg;
	return neg;
}

// Computes the TP curve.
UniqueLBRep mgPerimTP::create_tp(){
	size_t n = m_tpVec.size();
	if(!n)
		return UniqueLBRep(nullptr);

	MGLBRep& subtpf = subtp_first()->tp();
	if(n==1 && subtpf.param_range()==m_perimeter->param_range()){
		return UniqueLBRep(new MGLBRep(std::move(subtpf)));
	}

	MGNDDArray tau;
	MGBPointSeq tpData(0, 3);

	// 始点方向にスキマがあるとき
	double ts = m_perimeter->param_s();
	if (!MGRZero(subtpf.param_s() - ts)) {
		tau.resize(1);tau(0)=ts;
		tpData.resize(1); tpData.store_at(0,unit_normal_start());
	}

	// 複数のTP曲線の細切れから通過点列をサンプリング
	for (auto& i: m_tpVec)
		sampling(i->tp(), tau, tpData);

	// 終点方向にスキマがあるとき
	double te = m_perimeter->param_e();
	const MGLBRep& subtpl = subtp_last()->tp();
	if (!MGRZero(subtpl.param_e() - te)) {
		tau.add_data(te);
		int nold = tpData.length();
		tpData.reshape(nold+1);
		tpData.store_at(nold,unit_normal_end());		
	}

	// tau完成
	assert(tau.length() == tpData.length());
	UniqueLBRep tp(new MGLBRep);
	tp->buildByInterpolationDataPoints(tau, tpData, 4, 10.);

	mgTolSetLineZero lineZeroSet(MGTolerance::rc_zero());
	tp->remove_knot();
	return tp;
}

//Test if this has no boudary conditions.
bool mgPerimTP::is_0B_perim() const{
	return m_tpVec.size()==0;
}

// Output the state of this object.
void mgPerimTP::toString(std::ostream& os) const{
	// Output the addresses of this and both adjacent perimeters.
	os<<"] THIS["<<this<<"]\n";
	os << "m_tpmaker="<< m_tpmaker<<", m_perimNum"<<m_perimNum<<std::endl;

	// Output positional data of both ends of the perimeter curve.
	os << "m_perimeter="<<*m_perimeter<<std::endl;

	// Output all of subtps if exist.
	os << "m_tpVec:: ";
	if(m_tpVec.empty()){
		os << " (No subtp)\n";
	}else{
		os << std::endl;
		for(auto& i:m_tpVec)
			os<<(*i)<< std::endl;
	}
}

// Creates the components of the tp curve as good as possible.
// (1) Set mgTPOnFaceBoundary in m_tpVec getting the common part of face boundary and m_perimeter;
// (2) Sort m_tpVec;
// (3) And return the size of m_tpVec.
size_t mgPerimTP::buidTpVec(const std::vector<const MGFSurface*>& faces){
	// Loop through all of outer boundary curves.
	for(auto& facei:faces){
		// Make the outer boundary curves a ring.
		std::vector<UniqueCurve> wcrvs = facei->outer_boundary();
		std::vector<UniqueCurve> pcrvs = facei->outer_boundary_param();
		for(size_t j=0, n= wcrvs.size(); j<n; j++){
			MGCurve& wcurvej = *wcrvs[j];
			MGCurve& pcurvej = *pcrvs[j];

			std::vector<double> spans(4);
			int nspan;
			{//To limit the update of MGTolerance::line_zero() in this scope.
				double lzero = MGTolerance::line_zero();
				mgTolSetLineZero lineZeroSet(lzero*8.);//Loosen the tolerance.
				nspan = perimeter().common(wcurvej, spans);
			}
			if(nspan <= 0)
				continue;

			double prm1s = spans[0], prm1e = spans[1], prm2s = spans[2], prm2e = spans[3];
			nspan *= 4;
			for(int k = 4; k < nspan; k+=4){
				if(prm1s > spans[k]){
					prm1s = spans[k];
					prm2s = spans[k+2];
				}
				if(prm1e < spans[k+1]){
					prm1e = spans[k+1];
					prm2e = spans[k+3];
				}
			}

			// Reject crossed common.
			MGUnit_vector v1 = perimeter().direction(prm1s);
			MGUnit_vector v2 = wcurvej.direction(prm2s);
			double inp = fabs(v1%v2);
			if(!MGAZero(inp - 1.))
				continue;

			MGInterval srange(prm1s, prm1e);
			bool neg = need_negate(prm2s, prm2e);
			m_tpVec.emplace_back(
				new mgTPOnFaceBoundary(srange, perimeter(),
					*facei,	prm2s, prm2e, wcurvej, pcurvej, neg));
		}
	}

	// Sort subtp curves in order to connect correctly them.
	sort(m_tpVec.begin(), m_tpVec.end(),
		 [](MYELM& obj1, MYELM& obj2){return obj1->param_sPeri() < obj2->param_sPeri();});

	return m_tpVec.size();
}
