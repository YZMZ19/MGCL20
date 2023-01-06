/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/

/** @file TPOnFaceBoundary.cpp
 *  Implementation for class mgTPOnFaceBoundary.
 */
#include "stdafx.h"
#include "mg/FSurface.h"
#include "mg/TPOnFaceBoundary.h"

std::ostream& operator<<(std::ostream& os, const mgTPOnFaceBoundary& subtp){
	subtp.toString(os); return os;
}

// The constructor.
mgTPOnFaceBoundary::mgTPOnFaceBoundary(
	const MGInterval& perimRange,//!< parameter range of the part of the perimeter curve 
							//! (==the subtp curve's paramter range).
	const MGLBRep& perimcrv,//!< perimeter curve of the target MGSBRep to build.

	const MGFSurface& face, //!< MGSurface or MGFace.
	double ts, double te,   //!< parameter range of the part of boundary curve of face.
		  //ts>te when bndcrv and perimcrv have the opposite direction.
	const MGCurve&    bndcrv,//!< face's boundary curve that has common part with perimcrv. 
	const MGCurve&    paramcrv,

	bool face_negated  //!< whether face normal is necessary to negate for the tp data.
):m_periPart(perimcrv, perimRange),
m_f(face), m_boundaryS(ts),m_boundaryE(te),
m_boundaryParam(paramcrv),
m_necessaryToNegate(face_negated){
	build();
}

// Output the state of this object.
void mgTPOnFaceBoundary::toString(std::ostream& os) const{
	os <<"mgTPOnFaceBoundary : "<<this<<" m_periPart::"<< m_periPart<<std::endl;
	os <<"m_boundaryS, m_boundaryE::("<<m_boundaryS<<","<<m_boundaryE<<")"<<std::endl;
	os <<"m_boundaryParam::("<<m_boundaryParam<<std::endl;
	os <<"m_necessaryToNegate::"<<std::boolalpha <<m_necessaryToNegate<<std::noboolalpha<<std::endl;
}

/// Returns the unit normal vector at the paramete uv of m_f.
/// Note that the direction may be negated according to m_necessaryToNegate.
MGVector mgTPOnFaceBoundary::unit_normal_at(const MGPosition& uv) const{
	MGVector dir = m_f.unit_normal(uv);
	if(m_necessaryToNegate)
		dir.negate();
	return dir;
}

void mgTPOnFaceBoundary::getStoreTPdata(
	double tPeri, //parameter of m_periPart to get normal of m_f.
	double tBoundary,//parameter of m_boundaryParam to get normal of m_f.
	int i,//position of tpd to store
	MGBPointSeq& tpd
){
	MGPosition P = m_periPart.eval(tPeri);
	MGPosition uv, uvguess = m_boundaryParam.eval(tBoundary);
	if(!m_f.perp_guess(P, uvguess, uv))
		uv = uvguess;
	MGVector N = unit_normal_at(uv);
	tpd.store_at(i, N);
}

//! Build tangent plane line lbrep(m_tp).
void mgTPOnFaceBoundary::build(){
	// 1. Construct abscissa.
	MGNDDArray tau; tau.buildByKnotVector(m_periPart.knot_vector());

	// 2. Construct ordinate.
	// 2.1 Determine start and end of the parameter curve on m_f.
	int n = tau.length();
	int nm1 = n - 1;
	double tau0 = tau[0];
	double ratio = (m_boundaryE-m_boundaryS)/(tau[nm1]-tau0);

	// 2.2 Evaluate unit length normal vector on m_f as the ordinate.
	MGBPointSeq tpd(n, m_periPart.sdim());
	MGPosition uv = m_boundaryParam.eval(m_boundaryS);
	MGVector N = unit_normal_at(uv);
	tpd.store_at(0, N);
	for(int i=1; i <nm1; i++){
		double taui = tau[i];
		double s = m_boundaryS + ratio*(taui - tau0);
		getStoreTPdata(taui,s, i, tpd);
	}
	uv = m_boundaryParam.eval(m_boundaryE);
	N = unit_normal_at(uv);
	tpd.store_at(nm1, N);

	// 3. build m_tp
	MGKnotVector& t = m_tp.knot_vector();
	t = MGKnotVector(tau, 4);
	m_tp.buildByInterpolationWithKTV(tau, tpd);
}
