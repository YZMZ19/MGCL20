/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#include "StdAfx.h"
#include "mg/SBRepVecTP.h"
#include "mg/LBRep.h"
#include "mg/Surface.h"
#include "mg/SBRep.h"
#include "mg/Interval.h"
#include "mg/Position.h"
#include "mg/Tolerance.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//
// Implements Tangent Plane Line B-Representation.
// A vector version of MGSBRepTP

////////Constructor////////

//Move Constructor
//all of the ownership of m_TP of tp2 will be
//transfered to this object.
MGSBRepVecTP::MGSBRepVecTP(MGSBRepVecTP&& vectp2){
	for(int i=0; i<4; i++){
		m_TP[i]=std::move(vectp2.m_TP[i]);
		m_prange[i]=std::move(vectp2.m_prange[i]);

		int i2=i*2;
		m_to_SE[i2]=vectp2.m_to_SE[i2];
		m_to_SE[i2+1]=vectp2.m_to_SE[i2+1];
	}
}

MGSBRepVecTP::MGSBRepVecTP(const MGSBRepVecTP& vectp2){
	for (int i = 0; i < 4; i++) {
		MYVEC& vectpi = m_TP[i];
		const MYVEC& vectp2i = vectp2.m_TP[i];
		vectpi.resize(vectp2i.size());
		int jj = 0;
		for (const auto& j : vectp2i) 
			vectpi[jj++] = UniqueLBRep(j->clone());

		m_prange[i] = vectp2.m_prange[i];
		int i2 = i * 2;
		m_to_SE[i2] = vectp2.m_to_SE[i2];
		m_to_SE[i2 + 1] = vectp2.m_to_SE[i2 + 1];
	}
}
MGSBRepVecTP::MGSBRepVecTP(const MGSBRepTP& tp2){
	for(int i=0; i<4; i++){
		if(!(tp2.specified(i))) continue;
		m_TP[i].resize(1);
		m_TP[i][0]=MYELM(new MGLBRep(tp2.TP(i)));
		m_prange[i]=MGInterval(tp2.TP(i).param_range());
		m_to_SE[i*2]=m_to_SE[i*2+1]=true;
	}
}

MGSBRepVecTP::MGSBRepVecTP(MGSBRepTP && tp2){
	for (int i = 0; i < 4; i++) {
		if (!(tp2.specified(i))) continue;
		m_TP[i].resize(1);
		m_TP[i][0] = MYELM(&tp2.TP(i));
		tp2.TP()[i] = nullptr;
		m_prange[i] = MGInterval(tp2.TP(i).param_range());
		m_to_SE[i * 2] = m_to_SE[i * 2 + 1] = true;
	}
}

//Move Assignment.
//all of the ownership of m_TP of tp2 will be
//transfered to this object.
MGSBRepVecTP& MGSBRepVecTP::operator=(MGSBRepVecTP&& vectp2){
	for(int i=0; i<4; i++){
		m_TP[i].clear();
		m_TP[i]=std::move(vectp2.m_TP[i]);
		m_prange[i]=std::move(vectp2.m_prange[i]);

		int i2=i*2;
		m_to_SE[i2]=vectp2.m_to_SE[i2];
		m_to_SE[i2+1]=vectp2.m_to_SE[i2+1];
	}
	return *this;
}

//////////////// Member Function ///////////////////////

//change the parameter range to [t0,t1] if t0<t1,
//                           to [t1,t0] if t1<t0.
//when t1<t0, the direction will be reversed.
void MGSBRepVecTP::change_range(
	bool along_u,	//the objective range is along u parameter or v.
	double t0, double t1
){
	double lennew=t1-t0;
	size_t j= along_u ? 0:1;
	for(size_t m = 0; m<2; m++, j += 2){
		MYVEC& tpvecm = m_TP[j];
		size_t n = tpvecm.size();
		if(!n)
			continue;

		double lenold = m_prange[j].length();
		double ostart = m_prange[j][0];
		double ratio = lennew/lenold;
		for(size_t i = 0; i<n; i++){
			MGLBRep& lb = *(tpvecm[i]);
			double t20 = t0+(lb.param_s()-ostart)*ratio;
			double t21 = t0+(lb.param_e()-ostart)*ratio;
			lb.change_range(t20, t21);
		}
		if(t0>t1)
			std::reverse(tpvecm.begin(), tpvecm.end());
		m_prange[j].change_range(t0, t1);
	}
}

void MGSBRepVecTP::change_range(
	bool along_u,	//the objective range is along u parameter or v.
	const MGInterval& prange
){
	change_range(along_u,prange[0].value(), prange[1].value());
}

//evaluate TP at the perimeter i's parameter t.
//Function's return value is:
//true if t was inside the parameter range of a tangent plane of m_TP[i].
//false if t was outside the parameter range of the m_TP[i] for all i.
bool MGSBRepVecTP::eval(
	int i,		//perimeter numeber.
	double t,		//parameter vaule of perimeter i.
	MGVector& normal//evaluated normal will be returned.
)const{
	assert(i<4);
	if(!specified(i) || !(m_prange[i].includes(t)))
		return false;

	size_t n=m_TP[i].size();
	for(size_t j=0; j<n; j++){
		MGLBRep& tp=*(m_TP[i][j]);
		if(tp.in_range(t)){
			normal=tp.eval(t);
			return true;
		}
	}
	return false;
}

// Compute the maximum (absolute) cos value of between vector deris[i](t) 
// and vector this->TP(i)(t) for i=0,1,2,3, where t is a common
// parameter of the data point obtained from deris[i]'s knot vector.
//Function's return value is the max out of cosmax[.].
double MGSBRepVecTP::get_perimeters_max_cos(
	const UniqueLBRep deris[4],
	double taumax[4],
	double cosmax[4]
)const{
	MGVector N(3), T(3);
	MGNDDArray tau;
	double max=0.;
	for(int i = 0; i < 4; i++){
		if(!specified(i)){
			taumax[i] = deris[i]->param_s();
			cosmax[i] = 0.;
			continue;
		}

		tau.buildByKnotVector(deris[i]->knot_vector());
		double taus=tau[0];
		double cmi=0.;
		double tmi = taus;
		if(eval(i,taus,T)){
			N = deris[i]->eval(taus);
			cmi = fabs(N.cangle(T));
			tmi = taus;
		}

		int ntau=tau.length();
		for(int j = 1; j < ntau; j++){
			double tauj=tau[j];
			double tmid = (tau[j-1]+tauj)*.5;
			double cm;
			if(eval(i,tmid,T)){
				N = deris[i]->eval(tmid);
				cm = fabs(N.cangle(T));
				if(cmi<cm){	cmi = cm; tmi = tmid;}
			}

			if(!eval(i,tauj,T)) continue;
			N = deris[i]->eval(tauj);
			cm = fabs(N.cangle(T));
			if(cmi<cm){cmi = cm; tmi = tauj;}
		}
		taumax[i] = tmi;
		cosmax[i] = cmi;
		if(max<cmi) max=cmi;
	}
	return max;
}

// Compute the maximum (absolute) sin value of between vector srf.normal(uv(t))
// and vector this->TP(i)(t) for i=0,1,2,3, where perim[i] is 
// the same as srf.perimeter_curve(i), and t is a common parameter
// of deris[i] and TP(i).
double MGSBRepVecTP::get_perimeters_max_sin(
	const MGSurface& srf,
	double         taumax[4],
	double         sinmax[4],
	bool*          evalf	//indicates perimeters to evalate if evalf!=null
			//When evalf[i] is true, perimeter i is evaluated for 0<=i<=3.
)const{
	MGVector N(3), T(3);
	MGPosition uv(2);
	MGNDDArray tau;
	double max=0.;
	for(int i = 0; i < 4; i++){
		if(!specified(i) || (evalf && !evalf[i])){
			taumax[i] = (i % 2 == 0) ? srf.param_s_u() : srf.param_s_v();
			sinmax[i] = 0.;
			continue;
		}

		int id = i%2;
		switch(i){
		case 0:
			uv(1) = srf.param_s_v();
			tau.buildByKnotVector(srf.knot_vector_u());
			break;
		case 1:
			uv(0) = srf.param_e_u();
			tau.buildByKnotVector(srf.knot_vector_v());
			break;
		case 2:
			uv(1) = srf.param_e_v();
			tau.buildByKnotVector(srf.knot_vector_u());
			break;
		case 3:
			uv(0) = srf.param_s_u();
			tau.buildByKnotVector(srf.knot_vector_v());
			break;
		};

		double taus=tau[0];
		uv(id) = taus;

		double cmi=0.;
		double tmi = taus;
		if(eval(i,taus,T)){
			N = srf.normal(uv);
			cmi = fabs(N.sangle(T));
			tmi = uv(id);
		}

		int ntau=tau.length();
		for(int j = 1; j < ntau; j++){
			double tauj=tau[j];
			double tmid = (tau[j-1]+tauj)*.5;
			uv(id) = tmid;
			if(eval(i,tmid,T)){
				N = srf.normal(uv);
				double cm = fabs(N.sangle(T));
				if(cmi<cm){cmi = cm; tmi = tmid;}
			}

			uv(id) = tauj;
			if(eval(i,tauj,T)){
				N = srf.normal(uv);
				double cm = fabs(N.sangle(T));
				if(cmi<cm){	cmi = cm; tmi = tauj;}
			}
		}
		taumax[i] = tmi;
		sinmax[i] = cmi;
		if(max<cmi) max=cmi;
	}
	return max;
}

//Set i-th perimeter's TP(std::vector version).
//vectp[i] must be newed objects, and all of the ownership are transferer to
//this instance.
void MGSBRepVecTP::set_TP(
	int i,					//perimeter number.
	std::vector<UniqueLBRep>&& vectp,
	const MGInterval& prange		//Whole perimeter's parameter range.
){
	assert(i<4);
	MYVEC& tpi=m_TP[i];
	tpi.clear();

	size_t n=vectp.size();
	if (!n) return;

	for(size_t j=0; j<n; j++)
		tpi.push_back(std::move(vectp[j]));
	m_prange[i]=prange;

	double plen=prange.length();
	m_to_SE[2*i]=MGREqual_base((double)prange[0],tpi[0]->param_s(),plen);
	m_to_SE[2*i+1]=MGREqual_base((double)prange[1],tpi[n-1]->param_e(),plen);
}

//Set i-th perimeter's TP as a null, as an unspecified one.
void MGSBRepVecTP::set_TP_null(int i){
	assert(i<4);
	m_TP[i].clear();
}

//Debug Function.
std::ostream& operator<<(
	std::ostream& ostrm, const MGSBRepVecTP& vectp
){
	ostrm<<"MGSBRepVecTP::"<<&vectp;
	for(int i=0; i<4; i++){
		ostrm<<std::endl;
		ostrm<<"Perimeter "<<i<<":";
		ostrm<<",m_prange="<<vectp.m_prange[i];
		ostrm<<",m_to_SE=["; if(vectp.m_to_SE[2*i])ostrm<<"1"; else ostrm<<"0";
		ostrm<<",";if(vectp.m_to_SE[2*i+1])ostrm<<"1"; else ostrm<<"0";
		ostrm<<"]"<<std::endl;
		//ostrm<<vectp.m_TP[i];		
	}
	return ostrm;
}
