#include "StdAfx.h"
#include "mg/Tolerance.h"
#include "mg/Straight.h"
#include "mg/SPointSeq.h"
#include "mg/Transf.h"
#include "topo/Loop.h"
#include "Tl2/TL2parameter.h"
#include "Tl2/TL2Triangles.h"
#include "Tl2/TL2LPlines.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/****************************************************************/
/*   Copyright (c) 2021 by System fugen G.K.                    */
/*                       All rights reserved.                   */
/****************************************************************/

//Get point(x,y,z,xn,yn,zn) from uv and surface..
//Here (x,y,z) is the position data, and (xn,yn,zn) is the unit normal at (x,y,z).
MGPosition xyz(const MGSurface& srf, const MGPosition& uv, bool need_normal) {
	int sdim = need_normal ? 6 : 3;
	MGPosition xyzPN(sdim);
	if (need_normal)
		xyzPN.store_at(3, srf.unit_normal(uv));
	xyzPN.store_at(0, srf.eval(uv));
	return  xyzPN;
}

/// <summary>
/// tessellation of a equal side 4-sided rectangle.
/// Let n=nPLine(idMax), nOpo=nPLine(idOpo), and m=nPLine(idPre), 
/// then -1<=nPline(idNxt)-m<=1, and 0<=(n-nOpo)<=1.
/// </summary>
/// <param name="idMax">
/// edge id whose point number is equal to or greater than the oposite one,
///  and the difference must be less than or equal to 1.
/// </param>
void mgTL2LPlines::makeStripsE4(
	int idMax	/// edge id whose point number is equal to or greater than the oposite one,
				/// and the difference must be less than or equal to 1.
)const {
	int nPoints = numberOfPoints();
	if (nPoints == 4) {
		makeFan4Points(idMax);
	}else if (nPoints == 5) {
		int eids[5], np[5];
		analyzePointNumber(eids, np);
		makeFan5Points(eids);
	} else {
		int idOpo = (idMax + 2) % 4, idNxt = (idMax + 1) % 4, idPre = (idMax + 3) % 4;
		int n = nPLine(idMax), nOpo = nPLine(idOpo);
		int ndif = n - nOpo;
		assert(0 <= ndif && ndif <= 1);//The difference must be less than 2.

		int m = nPLine(idNxt), md = nPLine(idPre);
		int mdif = m - md;
		assert(-1 <= mdif && mdif <= 1);//The difference must be at most 1.
		int idE4Sub2 = -1;   //-1 means makeStripsE4Sub2 can not be applied to this.
		int idMakeStrip = -1;//-1 means makeStrip can not be applied to this.
		if (ndif || mdif) {
			if (m == 2) {
				if (md == 2) //case that m==2 && m==2 && ndif==1.
					idMakeStrip = idMax;
				else //case that m==2 && md==3.
					idE4Sub2 = idNxt;
			} else if (md == 2) //case that m==3 && md==2.
				idE4Sub2 = idPre;
			else if (n == 2) //case that n=2, nOpo=2, m>=3, md>=3.
				idMakeStrip = idNxt;
			else if (nOpo == 2)//case that n==3 && nOpo==2.
				idE4Sub2 = idOpo;
		}

		if (idMakeStrip != -1) 
			makeStrip(idMakeStrip);
		else if (idE4Sub2 != -1) 
			makeStripsE4Sub2(idE4Sub2);
		else
			makeStripsE4Sub(idMax);//n, nOpo, m, md>=3(all are greater than 2).
	}
}

//Proprietry subprogram of makeStripsE4, handles the case when one of the edge's point number is 2.
void mgTL2LPlines::makeStripsE4Sub2(
	int id2	/// edge id whose point number is 2. The opposite edge's number must be3.
			//m or md must be >=3, and the difference must be at most 1.
)const {
	int idOpo = (id2 + 2) % 4, idNxt = (id2 + 1) % 4, idPre = (id2 + 3) % 4;
	int m = nPLine(id2);
	assert(m==2 && nPLine(idOpo) == 3);
	int n = nPLine(idNxt), nd = nPLine(idPre);
	assert(-1 <= n - nd && n - nd <= 1);//The difference must be at most 1.

	const mgTL2LPline& e = m_plines[id2];
	const mgTL2LPline& eOpo = m_plines[idOpo];
	const mgTL2LPline& eNxt = m_plines[idNxt];
	const mgTL2LPline& ePre = m_plines[idPre];

	bool nextIsGreaterThanPre = (n > nd);
	int nmin = nextIsGreaterThanPre ? nd : n;
	int n1 = (nmin + 1) / 2;//Number of points of the smaller point number 
	int n2 = nmin - n1 + 1;
	
	int n1Nxt = n1 + 1, n1Pre = n1;
	if (!nextIsGreaterThanPre)
		std::swap(n1Nxt, n1Pre);

	int iPre = nd - n1Pre;
	int n1Nxtm1 = n1Nxt - 1;
	std::shared_ptr <mgTL2Polyline> line=
		std::make_shared< mgTL2Polyline>(*this, mgTLEdgePoint(idNxt, n1Nxtm1), mgTLEdgePoint(idPre, iPre));
	mgTL2LPline lp(line);

	mgTL2LPlines lines1(m_stack,m_triangles);
	lines1[0] = e;                        lines1[1] = mgTL2LPline(eNxt, 0, n1Nxt);
	lines1[2] = mgTL2LPline(lp).negate(); lines1[3] = mgTL2LPline(ePre, iPre, n1Pre);
	lines1.makeStrip(nextIsGreaterThanPre ? 1 : 3);

	mgTL2LPlines lines2(m_stack, m_triangles);
	lines2[0] = lp;   lines2[1] = mgTL2LPline(eNxt, n1Nxtm1, n2);
	lines2[2] = eOpo; lines2[3] = mgTL2LPline(ePre, 0, n2);
	lines2.makeStripsE4Sub(1);
}

//Proprietry subprogram of makeStripsE4, handles edges of any number of points.
void mgTL2LPlines::makeStripsE4Sub(
	int idMax	/// edge id whose point number is equal to or greater than the oposite one,
				/// and the difference must be less than or equal to 1.
)const{
	int idOpo = (idMax + 2) % 4, idNxt = (idMax + 1) % 4, idPre = (idMax + 3) % 4;
	int n = nPLine(idMax), nOpo = nPLine(idOpo);
	assert(n - nOpo <= 1);//The difference must be less than 2.
	int m = nPLine(idNxt), md = nPLine(idPre);
	assert(-1 <= m - md && m - md <= 1);//The difference must be at most 1.

	int mMax=std::max(m,md);
	int nm1 = n - 1, mMm1 = mMax - 1;
	MGSPointSeq Suv(n, mMax, 2);

//Evaluate n(idMax) and nOpo(idOpo) edges
	const mgTL2LPline& eMax = m_plines[idMax];
	const mgTL2LPline& eOpo = m_plines[idOpo];
	const mgTL2LPline& eNxt = m_plines[idNxt];
	const mgTL2LPline& ePre = m_plines[idPre];

	//idMax, idOpo edge evaluation.
	double dnm1=double(nm1);
	for (int i = 0; i < n; i++){
		double ti= double(i) / dnm1;
		Suv.store_at(i, 0, eMax.eval(ti));
		Suv.store_at(i, mMm1, eOpo.eval(1.-ti));
	}

	//idNxt and idPre  edges Evaluation.
	double dmMm1 = double(mMm1);
	for (int j = 1; j < mMm1; j++) {
		double tj = double(j) / dmMm1;
		Suv.store_at(nm1, j, eNxt.eval(tj));
		Suv.store_at(0, j, ePre.eval(1.-tj));
	}

	MGPosition P0 = eMax.uv(0), P1 = eNxt.uv(0),
			   P2 = eOpo.uv(0), P3 = ePre.uv(0);

	//Evaluate inner points.
	//1. along ePre and eNxt edges.
	for (int i = 1; i < nm1; i++) {
		MGPosition Pi0 = Suv(i,0), PimMm1 = Suv(i,mMm1);
		MGTransf ePreToi(P0, P3, Pi0, PimMm1);
		MGTransf eNxtToi(P1, P2, Pi0, PimMm1);
		double ri = double(i)/dnm1, rnm1mi = double(nm1-i)/dnm1;
		for (int j = 1; j < mMm1; j++) {
			MGPosition P= Suv(0,j)*ePreToi*rnm1mi + Suv(nm1,j)*eNxtToi*ri;
			Suv.store_at(i, j, P);
		}
	}
	//2. along n(eMax) and nOpo(eOpo) edges.
	for (int j = 1; j < mMm1; j++) {
		MGPosition P0j = Suv(0,j), Pnm1j = Suv(nm1,j);
		MGTransf eMaxToj(P0, P1, P0j, Pnm1j);
		MGTransf eOpoToj(P3, P2, P0j, Pnm1j);
		double rj = double(j) / double(mMm1), rmm1mj = double(mMm1 - j) / double(mMm1);
		for (int i = 1; i < nm1; i++) {
			MGPosition P= Suv(i,0)*eMaxToj*rmm1mj + Suv(i, mMm1)*eOpoToj*rj;
			MGPosition Q = Suv(i, j);
			Suv.store_at(i, j, (P+Q)*.5);
		}
	}
	makeStrips(idMax,Suv);
}

/// <summary>
/// tessellation of 4-sided rectangle.
/// Let n=nPLine(idMin), nOpo=nPLine(idOpo), and m=nPLine(idPre), 
/// then n>=2, m>=2, nPline(idNxt)=m.
/// </summary>
void mgTL2LPlines::tessellateEOpo(
	int idMin	/// edge id whose neibor edges have the same number of points.
				/// and nd(the number of points of the opposite edge)>=n;
)const {
	int idOpo = (idMin + 2) % 4, idNxt = (idMin + 1) % 4, idPre = (idMin + 3) % 4;
	int n = nPLine(idMin), nOpo = nPLine(idOpo), m = nPLine(idNxt);
	assert(nOpo >= n && nPLine(idPre) == m);
	if (nOpo - 1 <= n) {//when nOpo=n or nOpo-1=n;
		makeStripsE4(idOpo);
		return;
	}

	const mgTL2LPline& eMin = m_plines[idMin];
	const mgTL2LPline& eNxt = m_plines[idNxt];
	const mgTL2LPline& ePre = m_plines[idPre];
	const mgTL2LPline& eOpo = m_plines[idOpo];

	int mm1 = m - 1;
	int d = (nOpo - n) / mm1; if (d == 0) d = 1;
	double dmm1 = double(mm1);
	mgTL2LPline eMinR(eMin, n - 1, -n);

	mgTL2LPline lpGreater = eOpo, lpSmaller;
	int n2 = nOpo;
	for (int i = mm1 - 1, j = 1; j <= mm1; j++, i--) {
		if (i) {
			n2 -= d; if (n2 < n) n2 = n;
			mgTL2LPPoint P03i(eNxt, i), P12j(ePre, j);
			getMixedLPline(eOpo, double(i) / dmm1, eMinR, P03i, P12j, lpSmaller, n2);
		} else {
			lpSmaller = eMinR;
		}

		mgTL2LPlines lines(m_stack, m_triangles);
		lines[0] = lpGreater;                       lines[1] = mgTL2LPline(ePre, j-1, 2);
		lines[2] = mgTL2LPline(lpSmaller).negate(); lines[3] = mgTL2LPline(eNxt, i, 2);
		lines.tessellate2m2n(0);

		lpGreater = lpSmaller;
	}
}

/// <summary>
/// tessellation of 4-sided rectangle.
/// Let me=nPLine(id), m=nPLine(idOpo), nd=nPLine(idPre), n=nPLine(idNxt),
/// ndif=abs(nd-n), and mdif=(me-m),
/// then me>m, (mdif>=2 or ndif>=2), and min(n,nd)>2.
/// </summary>
void mgTL2LPlines::tessellate4General(
	int id //edge id to slide along the neighbor edges.
)const {
	int idOpo = (id + 2) % 4, idNxt = (id + 1) % 4, idPre = (id + 3) % 4;
	const mgTL2LPline& e = m_plines[id];
	const mgTL2LPline& eNxt = m_plines[idNxt];
	const mgTL2LPline& ePre = m_plines[idPre];
	const mgTL2LPline& eOpo = m_plines[idOpo];

	int m = nPLine(idOpo), me = nPLine(id), n = nPLine(idNxt), nd = nPLine(idPre);
	int mdif = me - m, ndif = std::abs(nd - n);
	assert(mdif >= 1 && ndif >= 1 && (mdif >= 2 || ndif >= 2));

	int ndm1 = nd - 1;
	bool nextIsGreaterThanPre = (n > nd);
	int nmin = nextIsGreaterThanPre ? nd : n; assert(nmin > 2);
	int mm1 = m - 1, nminm1 = nmin - 1;
	double dnminm1 = double(nminm1);

	mgTL2LPline eOpoR(eOpo, mm1, -m);//eOpoR has the reversed direction of eOpo

	int dm = mdif / nminm1, dmRemainder = mdif % nminm1;
	int dmp1 = dm + 1;
	int dn = ndif / nminm1, dnRemainder = ndif % nminm1;
	int iPre = ndm1, iNxt = 0;//id of ePre and eNxt.

	mgTL2LPline lpGreater = e, lpSmaller;
	int m2 = me;//the number of lpGreater.
	int id2 = nextIsGreaterThanPre ? 3 : 1;// edge id whose number of points is 2(of lines below).
	for (int i = 1; i <= nminm1; i++) {
		int nPre = dn + 3, nNxt = 2;
		if (i > dnRemainder)
			nPre--;
		if (nextIsGreaterThanPre)
			std::swap(nPre, nNxt);
		int nPrem1 = nPre - 1, nNxtm1 = nNxt - 1;
		iPre -= nPrem1;
		iNxt += nNxtm1;

		int deltam = (i <= dmRemainder) ? dmp1 : dm;
		m2 -= deltam; assert(m2 >= m);
		if (i < nminm1) {
			mgTL2LPPoint Ppre(ePre, iPre), Pnxt(eNxt, iNxt);
			getMixedLPline(eOpoR, double(i) / dnminm1, e, Ppre, Pnxt, lpSmaller,m2);
		}else
			lpSmaller = eOpoR;

		std::unique_ptr<mgTL2LPlines>
			linesP(new mgTL2LPlines(m_stack, m_triangles));
		mgTL2LPlines& lines = *linesP;
		lines[0] = lpGreater; lines[1] = mgTL2LPline(eNxt, iNxt - nNxtm1, nNxt);
		lines[2] = mgTL2LPline(lpSmaller).negate(); lines[3] = mgTL2LPline(ePre, iPre, nPre);
		if(!lines.makeStrips2(id2))
			m_stack.push(std::move(linesP));

		lpGreater = lpSmaller;
	}
}

/// <summary>
/// make strips from 4-sided rectangle if makeStripsE4 is, tessellate2m2n, or
/// tessellate22mn,  applicable to this.
/// Function's retulrn value is true if tessellation was applied.
/// Let n=nPLine(id2), nOpo=nPLine(idOpo), m=nPLine(idPre), and md=nPLine(idNxt),
/// then n=2, nOpo>=2, and m, md>2.
/// </summary>
/// <param name="id2">
/// edge id whose point number is 2.
/// </param>
bool mgTL2LPlines::makeStrips2(
	int id2	// edge id whose point number is 2.
){
	int idOpo = (id2 + 2) % 4, idNxt = (id2 + 1) % 4, idPre = (id2 + 3) % 4;
	int n = nPLine(id2), nOpo = nPLine(idOpo);
	assert(n == 2);

	bool processed = true;
	int m = nPLine(idNxt), md = nPLine(idPre);
	int mmmd = m - md;
	int idMax = mmmd >= 0 ? idNxt : idPre;
	if (nOpo <= 3 && -1 <= mmmd && mmmd <= 1)
		makeStripsE4(idMax);
	else if (nOpo == 2)
		tessellate2m2n(idMax);
	else if (m == 2)
		tessellate22mn(id2);
	else if (md == 2)
		tessellate22mn(idPre);
	else 
		processed = false;

	return processed;
}

/// <summary>
/// Tessellate if makeStripsE4 or tessellateEOpo is applicable.
/// Function's return value is true if tessellated using makeStripsE4 or tessellateEOpo,
/// false if not.
/// </summary>
/// <param name="id2">
/// edge id whose point number is 2.
/// </param>
bool mgTL2LPlines::t4makeStripsE4OrtessellateEOpo(
	const int eids[5],
	const int np[MAX_LINE_NUM]//input analyzePointNumber's output.
)const{
	const int& eidMax = eids[4];//Maximum vertex number and the id of m_plines[i] will be stored.

	int eidPre = (eidMax + 3) % 4, eidNxt = (eidMax + 1) % 4, eidOpo = (eidMax + 2) % 4;
	int dif1 = np[eidMax] - np[eidOpo], dif2 = np[eidPre] - np[eidNxt];
	if (dif1 <= 1 && -1 <= dif2 && dif2 <= 1) {
		makeStripsE4(eidMax);
		return true;
	}

	int idMin = -1;
	if (dif1 == 0)
		idMin = dif2 < 0 ? eidPre : eidNxt;
	else if (dif2 == 0)
		idMin = eidOpo;
	if (idMin >= 0) {
		tessellateEOpo(idMin);
		return true;
	}
	return false;
}

/// <summary>
/// Tessellation of 4 side rectangle that has 2 or 3 2 point edge.
/// </summary>
void mgTL2LPlines::t4Side22(
	int numEdgeTwoPoint,//number of edges whose points are only 2.
	int eidMin,//edge whose point number is minimum.
	int eidMax //edge whose point number is maximum.
)const {
	assert(getNumberOfLines() == 4 && numEdgeTwoPoint >= 2);
	if (numEdgeTwoPoint == 3) {//Case that 3 of the edges have only 2 vertices.
		makeFan222n(eidMax);
	}else {	//numEdgeTwoPoint==2
		if (nPLine((eidMin+1)%4) == 2)//When next edge's point number is also 2.
			tessellate22mn(eidMin);
		else {
			int eidPre = (eidMin + 3) % 4;
			if (nPLine(eidPre) == 2)//When pre edge's point number is also 2.
				tessellate22mn(eidPre);
			else
				tessellate2m2n(eidMax);
		}
	}
}

/// <summary>
/// make strips from 4-sided rectangle if makeStripsE4 is applicable to this.
/// Function's retulrn value is true if tessellation was applied.
/// Let n=nPLine(id2), nOpo=nPLine(idOpo), m=nPLine(idPre), and md=nPLine(idNxt),
/// then n=2, nOpo>=2, and m, md>2.
/// </summary>
/// <param name="id2">
/// edge id whose point number is 2.
/// </param>
bool mgTL2LPlines::t4GeneralProcess()const {
	int eids[5], np[MAX_LINE_NUM];//np[i] is m_lines[i]'s number of points.
	analyzePointNumber(eids, np);
	if (t4makeStripsE4OrtessellateEOpo(eids,np))
		return true;

	int numEdgeTwoPoint = eids[0];//number of edges whose points are only 2.
	int eidMin = eids[2];//minimum number of vertices and the edge id will be stored.
	int eidMax = eids[4];//Maximum vertex number and the id of m_plines[i] will be stored.

	if (numEdgeTwoPoint >= 2) {
		t4Side22(numEdgeTwoPoint, eidMin, eidMax);
		return true;
	}

	return false;
}

///Make strip data from n*m points Suv and this boundary m_plines.
///The number of vertices of the idMax(n) is greater or equal to the one of the opposite,
///and the differecne must be at most 1. Furthermore, 
///Neighbor edges difference of idMax must be at most 1.
void mgTL2LPlines::makeStrips(
	int idMax,//edge id
	const MGSPointSeq& Suv //surface (u,v) date of n*mMax dimension, mMax=max(m,md).
)const {
	int idOpo = (idMax + 2) % 4, idNxt = (idMax + 1) % 4;
	const mgTL2LPline& eMax = m_plines[idMax];
	const mgTL2LPline& eNxt = m_plines[idNxt];
	const mgTL2LPline& ePre = m_plines[(idMax + 3) % 4];
	const mgTL2LPline& eOpo = m_plines[idOpo];
	int nOpo = nPLine(idOpo), n = nPLine(idMax);assert(n >= nOpo && n - nOpo <= 1);
	int nm1 = n - 1, nm2 = n - 2, npn = n + n;
	int npnm1 = n + nm1;

	int m = nPLine(idNxt), md = ePre.number_of_points(); assert(-1<=m-md && m-md<=1);
	int mm1=m-1, mm2 = m - 2, mm3=m-3;
	int mMn=std::min(m,md), mMxm1= std::max(m, md) -1;

	bool uvRequired = m_triangles.is_uv(), normal_is_required = m_triangles.need_normal();
	auto PonEdge =	[uvRequired,normal_is_required](const mgTL2LPline& edge, int i) {
		return uvRequired ? edge.uv(i) : edge.xyz(i, normal_is_required); 
	};
	const MGSurface& srf = eMax.TL2param().get_surface();
	auto Pinner = [uvRequired, normal_is_required, &Suv, &srf](int i, int j) {
		return uvRequired ? Suv(i,j) : xyz(srf, Suv(i,j), normal_is_required);
	};

	bool equalBothSides = (m == md && n == nOpo);
	int mLast = equalBothSides ? mMn - 2 : mMn - 3;
	int mdm1 = md - 1;
	for (int j = 0; j <= mLast; j++) {
		mgTL2Triangle& strip = *(new mgTL2Triangle(npn, mgTESTRIANG::mgTESTRIANG_STRIP));
		int jp1 = j + 1;
		int k = 0;
		strip[k++] = PonEdge(eNxt, j);
		strip[k++] = PonEdge(eNxt, jp1);
		for (int i = nm2; i >= 1; i--) {
			strip[k++] = j ? Pinner(i, j) : PonEdge(eMax,i);//for eMax edge.
			strip[k++] = (jp1==mMxm1) ? PonEdge(eOpo,nm1-i): Pinner(i, jp1);//for eOpo edge.
		}
		strip[k++] = PonEdge(ePre,mdm1-j);
		strip[k++] = PonEdge(ePre,mdm1 - jp1);
		m_triangles.push_back(&strip);
	}
	if (equalBothSides)
		return;

	if (m == md) {//and nOpo==n-1.
		mgTL2Triangle& strip2 = *(new mgTL2Triangle(npnm1, mgTESTRIANG::mgTESTRIANG_STRIP));
		int k = 0;
		strip2[k++] = PonEdge(eNxt,mm2);
		strip2[k++] = PonEdge(eNxt,mm1);
		for (int i = nm2, j = 1; i>=1; i--, j++) {
			strip2[k++] = Pinner(i, mm2);
			strip2[k++] = PonEdge(eOpo,j);
		}
		strip2[k] = PonEdge(ePre,1);
		m_triangles.push_back(&strip2);
		return;
	}

	//Here  m!=md and (n==nOpo or nOpo=n-1).
	int nstrip2= npnm1, nstrip3= npnm1;//default value of case 3.
	if (n == nOpo) {
		if (m > md) {//case 2
			nstrip2 = npn;
		} else {//case 4
			nstrip3 = npn;
		}
	} else {
		if (m > md) {//case 1
			nstrip2 = npn;
			nstrip3 = n + nm2;
		}
	}
	mgTL2Triangle& strip2 = *(new mgTL2Triangle(nstrip2, mgTESTRIANG::mgTESTRIANG_STRIP));
	mgTL2Triangle& strip3 = *(new mgTL2Triangle(nstrip3, mgTESTRIANG::mgTESTRIANG_STRIP));
	int k2 =0, k3 = 0;
	if (m > md) {//m=md+1.
	//case 1, 2 : m=md+1, n=nOpo+1 or n=nOpo.
		strip2[k2++] = PonEdge(eNxt, mm3);
		strip3[k3++] = strip2[k2++] = PonEdge(eNxt, mm2);
		strip3[k3++] = PonEdge(eOpo, 0);
		for (int i = nm2, j = 1; i >= 1; i--, j++) {
			strip2[k2++] = mm3 ? Pinner(i, mm3):PonEdge(eMax,i);
			strip3[k3++] = strip2[k2++] = Pinner(i, mm2);
			strip3[k3++] = PonEdge(eOpo, j);
		}
		strip2[k2++] = PonEdge(ePre, 1);
		strip2[k2] = PonEdge(ePre, 0);
		if (n == nOpo)//if case 2 
			strip3[k3] = PonEdge(ePre, 0);
	} else {//md=m+1.
	//case 3, 4 : md=m+1, n==nOpo or n=nOpo+1.
		if (n == nOpo) {//case 4
			k2 = nstrip2 - 1;
			strip3[k3++] = strip2[k2--] = PonEdge(eNxt, mm2);
			strip3[k3++] = PonEdge(eOpo, 0);
			int mdm2 = md - 2, mdm3 = md - 3;
			for (int i = nm2, j = 1; i >= 1; i--, j++) {
				strip2[k2--] = mdm3 ? Pinner(i, mdm3): PonEdge(eMax, i);
				strip3[k3++] = strip2[k2--] = Pinner(i, mdm2);
				strip3[k3++] = PonEdge(eOpo, j);
			}
			strip2[k2--] = PonEdge(ePre, 2);
			strip3[k3++] = strip2[k2--] = PonEdge(ePre, 1);
			strip3[k3++] = PonEdge(ePre, 0);
		}
		else {// case 3
			strip3[k3++] = strip2[k2++] = PonEdge(eNxt, mm2);
			strip3[k3++] = PonEdge(eOpo, 0);
			int mdm2 = md - 2, mdm3 = md - 3;
			for (int i = nm2, j = 1; i >= 1; i--, j++) {
				strip3[k3++] = strip2[k2++] = Pinner(i, mdm2);
				strip2[k2++] = mdm3 ? Pinner(i, mdm3): PonEdge(eMax, i);
				strip3[k3++] = PonEdge(eOpo, j);
			}
			strip3[k3++] = strip2[k2++] = PonEdge(ePre, 1);
			strip2[k2++] = PonEdge(ePre, 2);
		}
	}
	m_triangles.push_back(&strip2);
	m_triangles.push_back(&strip3);
}

//Mix line1 and line2 by line1d*r1+line2d*(1-r1). Here,  line1d is transformed line1
//as line1's start point becomes P0, and the end point, P1.
//line2d is the same.
//Generated mgTL2Polyline's point number is nNew.If nNew is 0, line1's number is employed.
void getMixedLPline(
	const mgTL2LPline& line1, double r1,
	const mgTL2LPline& line2,
	const MGPosition& P0, const MGPosition& P1,
	mgTL2LPline& lp, int nNew
) {
	int n1 = line1.number_of_points(), n2 = line2.number_of_points();
	MGTransf l1ToP0P1(line1.uv(0), line1.uv(n1 - 1), P0, P1);
	MGTransf l2ToP0P1(line2.uv(0), line2.uv(n2 - 1), P0, P1);

	if (nNew) n1 = nNew;//update n1.
	int n1m1 = n1 - 1;

	std::shared_ptr<mgTL2Polyline> poly=std::make_shared< mgTL2Polyline>(line1.TL2param(), n1);//new line.
	MGBPointSeq& bp = poly->line_bcoef();
	MGKnotVector& t = poly->knot_vector();
	double dn1m1 = double(n1m1);
	t(0) = t(1) = 0.;
	bp.store_at(0, P0);
	double r2 = 1. - r1;
	for (int i = 1; i < n1m1; i++) {
		double di=t(i + 1) = double(i);
		double diOverDn1m1 = di / dn1m1;
		MGPosition Q1 = line1.eval(diOverDn1m1) * l1ToP0P1 * r1;
		MGPosition Q2 = line2.eval(diOverDn1m1) * l2ToP0P1 * r2;
		bp.store_at(i, Q1 + Q2);
	}
	bp.store_at(n1m1, P1);
	t(n1) = t(n1 + 1) = double(n1m1);
	lp = mgTL2LPline(poly);
}

void getMixedLPline(
	const mgTL2LPline& line1, double r1,
	const mgTL2LPline& line2,
	const MGPosition& Ps,//Start point.
	const mgTL2LPPoint& Pe, //End point.
	mgTL2LPline& lp, int nNew
) {
	getMixedLPline(line1, r1, line2, Ps, Pe.uv(),lp,nNew);
	std::shared_ptr< mgTL2Polyline>& line = lp.sharedLine();
	line->setBoundaryID(Pe);
}

void getMixedLPline(
	const mgTL2LPline& line1, double r1,
	const mgTL2LPline& line2,
	const mgTL2LPPoint& Ps, //Start point.
	const mgTL2LPPoint& Pe, //End point.
	mgTL2LPline& lp, int nNew
) {
	getMixedLPline(line1, r1, line2, Ps.uv(), Pe.uv(),lp, nNew);
	std::shared_ptr< mgTL2Polyline>& line = lp.sharedLine();
	line->setBoundaryID(Ps, false);
	line->setBoundaryID(Pe);
}
