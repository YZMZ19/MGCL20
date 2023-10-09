/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#include "StdAfx.h"
#include "mg/Tolerance.h"
#include "mg/Box.h"
#include "mg/Position.h"
#include "mg/Transf.h"
#include "mg/Unit_vector.h"
#include "mg/CParam_list.h"
#include "mg/CCisects.h"
#include "mg/CSisects.h"
#include "mg/CSisect.h"
#include "mg/Position_list.h"
#include "mg/Straight.h"
#include "mg/LBRep.h"
#include "mg/RLBRep.h"
#include "mg/Ellipse.h"
#include "mg/SurfCurve.h"
#include "mg/BSumCurve.h"
#include "mg/Plane.h"
#include "mg/Sphere.h"
#include "mg/Cylinder.h"
#include "mg/SBRep.h"
#include "mg/RSBRep.h"
#include "mg/BSumSurf.h"
#include "topo/Face.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//MGStraighy Class Implementation
//


//
// メンバ関数
//

///Compute the closest point parameter value of this curve from a point.
double MGStraight::closest(const MGPosition& point)const{
	double t=perp_param(point);//t is the parameter value of the closest point on this.
	double error=MGTolerance::rc_zero()*m_direction.len();
	mgTolSetWCZero wczeroSet(error);//Set&save the error.
	t=range(t);
	return t;
}

///Compute the closest point parameter value pair of this curve and curve2.
///MGPosition P of the function return contains this and curve2's parameter
///as:     P(0)=this curve's parameter, P(1)=curve2's parameter value.
MGPosition MGStraight::closest(const MGCurve& curve2)const{
	const MGStraight* sl2=dynamic_cast<const MGStraight*>(&curve2);
	if(sl2)
		return closestSL(*sl2);

	MGUnit_vector sldir=direction();
	const MGPosition p_on_crv2=curve2.center();
	MGPosition origin=eval_position_unlimitted(perp_param(p_on_crv2));

	MGMatrix M; M.set_axis(sldir,2);
	std::unique_ptr<MGCurve> crv2dP(curve2.clone());
	MGCurve& crv2d=*crv2dP;
	crv2d-=origin;
	crv2d*=M;

	MGPosition param(2);
	double* tout=param.data();
	double& otherParam=tout[1]=crv2d.closest2D(MGDefault::origin_2D());
	double& thisParam=tout[0]=perp_param(curve2.eval(tout[1]));
	
	MGVector V=curve2.eval(otherParam)-eval(thisParam);
	double d=V%V;

	double paramS=param_s(), paramE=param_e();
	MGPosition SPoint=eval(paramS);
	double otherParam2=curve2.closest(SPoint);
	MGVector VS=curve2.eval(otherParam2)-SPoint;
	double dS=VS%VS;
	if(dS<d){
		d=dS;
		otherParam=otherParam2;
		thisParam=paramS;
	}
	MGPosition EPoint=eval(paramE);
	double otherParam3=curve2.closest(EPoint);
	MGVector VE=curve2.eval(otherParam3)-EPoint;
	double dE=VE%VE;
	if(dE<d){
		otherParam=otherParam3;
		thisParam=paramE;
	}

	return param;
}

// 自身と与えられた点との距離を返す。
double MGStraight::distance(
	const MGPosition& p
)const{
	if(straight_type()== MGSTRAIGHT_TYPE::MGSTRAIGHT_EMPTY)
		return 0.0;

	MGVector v1(p,root_point());  // 基点から与点へのベクトル。
	double v1len=v1.len();
	double dlen=m_direction.len();
	double hs = (v1%m_direction)/dlen;
						// 与点から直線への垂線の足と基点との距離。
	double t=hs/dlen;	//parameter of the straight.

	double d;
	if(t<m_sparam) d=(p-start_point()).len();
	else if(t>m_endparam) d=(p-end_point()).len();
	else{
		d= v1len*v1len - hs*hs; if(d<0.) d=0.;
		d= sqrt(d);	// 点と無限直線の最短距離。
	}
	// 最短距離を返す。
	return d;
}

///Compute the closest point parameter value pair of this MGStraight and straight2.
///MGPosition P of the function return contains this and straight2's parameter
///as:     P(0)=this MGStraight's parameter, P(1)=straight2's parameter value.
MGPosition MGStraight::closestSL(const MGStraight& straight2)const{
	// 与えられた直線との関係を調べる。
	MGCCisect is;
	MGPSRELATION rl=relation(straight2,is);
	MGPosition st(is.param1(), is.param2());
		//Function's return value: [0] this straight's parameter,
        //[1]: straight2's parameter.

	double& s=st(0);
	double& t=st(1);
	switch(rl){
	case MGPSRELATION::MGPSREL_ISECT:
	case MGPSRELATION::MGPSREL_PARALLEL:
	case MGPSRELATION::MGPSREL_COIN:
		break;

	default://MGPSREL_VIRTUAL_ISECT, or MGPSREL_TORSION,
			//and s_in_range=false or t_in_range=false.
		bool s_in_range=in_range(s);
		bool t_in_range=straight2.in_range(t);
		if(s_in_range && t_in_range)
			break;

		//case that !s_in_range or !t_in_range holds.
		double distance1=-1., distance2=-1;
		double s1,s2, t1,t2;
		if(!s_in_range){
			s1=range(s);
			MGPosition P=eval(s1); 
			t1=straight2.closest(P);
			MGVector V=P-straight2.eval(t1);
			distance1=V%V;
		}
		if(!t_in_range){
			t2=straight2.range(t);
			MGPosition Q=straight2.eval(t2); 
			s2=closest(Q);
			MGVector V=Q-eval(s2);
			distance2=V%V;
		}

		if(s_in_range){
			//case of !t_in_range
			s=s2;
			t=t2;
		}else if(t_in_range){
			//case of !s_in_range
			s=s1;
			t=t1;
		}else{
			//!s_in_range and !t_in_range
			if(distance1<distance2){
				s=s1;
				t=t1;
			}else{
				s=s2;
				t=t2;
			}
		}
	}
	return st;
}

// 自身と与えられた直線との距離を返す。
double MGStraight::distance(const MGStraight& straight2) const{
	if(straight_type()== MGSTRAIGHT_TYPE::MGSTRAIGHT_EMPTY ||
		straight2.straight_type()== MGSTRAIGHT_TYPE::MGSTRAIGHT_EMPTY)
		return 0.;

	MGPosition st=closestSL(straight2);
	MGPosition P1=eval_position_unlimitted(st[0]);
	MGPosition P2=straight2.eval_position_unlimitted(st[1]);
	return P1.distance(P2);
}

////////////isect with a curve.

// Straight と Curve の交点を求める。
MGCCisects MGStraight::isect(const MGCurve& curve)const{
	const MGStraight* sl2 = dynamic_cast<const MGStraight*>(&curve);
	if(sl2)
		return isect(*sl2);

	MGCCisects list=curve.isect(*this);
	list.exchange12();
	return list;
}

// Straight と Straight の交点を求める。
MGCCisects MGStraight::isect(const MGStraight& st) const{
	std::unique_ptr<MGCCisect> p(new MGCCisect);
	MGCCisects list(this, &st);
	MGPSRELATION rel=relation(st,*p);
	if(rel== MGPSRELATION::MGPSREL_ISECT || rel== MGPSRELATION::MGPSREL_COIN )
		list.push_back(std::move(p));
	return list;
}

//Compute intersections with MGRLBRep curve2.
MGCCisects MGStraight::isect(const MGRLBRep& curve2)const{
	MGCCisects list=curve2.isect(*this);
	list.exchange12();
	return list;
}

//Intersection with a ellipse.
MGCCisects MGStraight::isect(const MGEllipse& curve2)const{
	MGCCisects list=curve2.isect(*this);
	list.exchange12();
	return list;
}

//Intersection with a MGBSumCurve.
MGCCisects MGStraight::isect(const MGBSumCurve& curve2)const{
	MGCCisects list=curve2.isect(*this);
	list.exchange12();
	return list;
}

//isect with SurfCurve whose m_curve is not a MGTrimmedCurve of MGCompositeCurve.
MGCCisects MGStraight::isect_with_noCompoSC(const MGSurfCurve& curve2)const{
	MGCCisects list=curve2.isect_noCompo(*this);
	list.exchange12();
	return list;
}

//Compute intersections with MGLBRep curve2 that does not have C0 continuity in it.
MGCCisects MGStraight::isect_withC1LB(const MGLBRep& curve2)const{
	MGCCisects list=curve2.C1isect(*this);
	list.exchange12();
	return list;
}

MGCSisects MGStraight::isect(const MGSurface & srf) const{
	return srf.isectSl(*this);
}
MGCSisects MGStraight::isect(const MGPlane & srf) const{
	return srf.isectSl(*this);
}
MGCSisects MGStraight::isect(const MGFace& f)const{
	MGCSisects list;
	const MGSurface* srf=f.surface();
	if(!srf)
		return list;

	const MGBox& sbx=f.box();
	if(!sbx.crossing(*this))
		return list;

	list=srf->isectSl(*this,f.box_param());
	MGCSisects::iterator	i=list.begin(), iend=list.end(), i1;
	while(i!=iend){
		i1=i; i1++;
		auto& csi=isectCast<MGCSisect>(i);
		if(!f.in_range(csi.param_surface()))
			list.removeAt(i);
		i=i1;
	}
	return list;
}

//Compute intersection point of 1D sub curve of original curve.
//Parameter values of intersection point will be returned.
MGCParam_list MGStraight::intersect_1D(						
	double f,			// Coordinate value
	int coordinate	// Coordinate kind of the data f(from 0).
)const{
	MGCParam_list list(this);

	double u=direction().ref(coordinate);
	if(MGMZero(u)){
		if(MGAEqual(f,m_root_point(coordinate)))
			list.append(m_sparam.value());
	}else{
		double p=f-root_point().ref(coordinate);
		double t=p/u;
		if(in_range(t)) list.append(t);
	}
	return list;
}

//isect2D returns parameter values of this(t) and l2(s)
// of the intersection point of both 2D straight lines.
// This and l2 are treated as infinite lines.
//Function's return value is:
// true if two lines are parallel(or one of the directin is zero)
// false if intersection was obtained.
bool MGStraight::isect2D(const MGStraight& l2,double& t,double& s)const{
	double x1, a1, y1, b1, x2, a2, y2, b2;
	x1=root_point().ref(0); a1=m_direction.ref(0);
	y1=root_point().ref(1); b1=m_direction.ref(1);
	x2=l2.root_point().ref(0); a2=l2.m_direction.ref(0);
	y2=l2.root_point().ref(1); b2=l2.m_direction.ref(1);
	double det = a1*b2-b1*a2;
	if(MGMZero(det))
		return true;

	double y1my2=y1-y2; double x2mx1=x2-x1;
	t=(a2*y1my2+b2*x2mx1)/det;
	s=(a1*y1my2+b1*x2mx1)/det;
	return false;
}

// 点が直線上にあるかを試験する。直線上にあれば，その点のパラメーター値を，
//　直線上になくても最近傍点のパラメーター値を返す。
bool MGStraight::on(
	const MGPosition& p,	 // 指定点       
	double& d				 // パラメータ値 
)const{
	bool on; double t;
	if(perp_point(p, t))
		on=MGAZero(MGVector(p,eval(t)).len());
	else
		on=false;
	d=range(t);
	return on;
}

// 直線が平面上にあるか調べる。（平面上ならばtrue）
bool MGStraight::on(
	const MGPlane& pl	// Plane
)const{
	MGPosition uv(2);
	return direction().orthogonal(pl.normal()) && pl.on(root_point(), uv);
}

// 直線上の与えられたポイントにおけるパラメータ値を返す。
// If input point is not on the curve, return the nearest point on the
// curve.
double MGStraight::param(const MGPosition& p)const{
	// 垂線の足を求める。
	double d2;
	perp_point(p,d2);

	// パラメータ値を返す。
	return range(d2);
}

// 与点から直線への垂線の足のパラメータ値を返却する。
//Return the foot of the  straight line that is perpendicular to this line.
//Function's return value is parameter value of this straight line,
//may ***NOT*** be in_range.
double MGStraight::perp_param(
	const MGPosition& point		// 与点
)const{
	MGVector v2(point,root_point());// 始点から与点へのベクトル			
	double lensqr=m_direction.len(); lensqr*=lensqr;
	return v2%m_direction/lensqr;	// 与点から直線への垂線の足と始点との距離
}

// 与えられたポイントから曲線への垂線の足、そのパラメータ値を返却する。
// Function's return value is if point is obtained(1) or not(0)
int MGStraight::perp_point(
	const MGPosition& point,	// 指定点
	double& d1,					// 垂線の足のパラメータ値
	const double* d2			// guess parameter value of d1.
)const{
	d1=perp_param(point);	
	return in_range(d1);
}
	
// 与ポイントから直線へ下ろした垂線の足の，直線のパラメータ値を
// すべて求める。
MGCParam_list MGStraight::perps(
	const MGPosition& point	// 与ポイント
)const{
	MGCParam_list tlist(this);
	double l =perp_param(point);	
	if(in_range(l)) tlist.append(l);
	return tlist;
}

//Compute all the perpendicular points of this curve and the second one.
//That is, if f(s) and g(t) are the points of the two curves f and g,
//then obtains points where the following conditions are satisfied:
//  fs*(f-g)=0.    gt*(g-f)=0.
//Here fs and gt are 1st derivatives at s and t of f and g.
//MGPosition P in the MGPosition_list contains this and crv2's parameter
//as:     P(0)=this curve's parameter, P(1)=crv2's parameter value.
MGPosition_list MGStraight::perps(
	const MGCurve& crv2		//The second curve
)const{
	MGPosition_list list=crv2.perps(*this);
	return MGPosition_list(list,1,0);	
}
MGPosition_list MGStraight::perps(
	const MGStraight& l2		//The second curve
)const{
	if(MGMZero(m_direction.sangle(l2.m_direction)))// 平行
		return relation_parallel(l2);
	MGPosition P;
	perp_guess(1.,0.,l2,1.,0.,0.,0.,P);
	MGPosition_list list;
	if(in_range(P.ref(0)) && l2.in_range(P.ref(1)))
		list.append(P);
	return list;
}

MGPosition_list MGStraight::perps(
	const MGRLBRep& crv2	//The second curve
)const{
	MGPosition_list list=crv2.perps(*this);
	return MGPosition_list(list,1,0);	
}

MGPosition_list MGStraight::perps(
	const MGEllipse& crv2	//The second curve
)const{
	MGPosition_list list=crv2.perps(*this);
	return MGPosition_list(list,1,0);	
}

MGPosition_list MGStraight::perps(
	const MGSurfCurve& crv2	//The second curve
)const{
	MGPosition_list list=crv2.perps(*this);
	return MGPosition_list(list,1,0);	
}

MGPosition_list MGStraight::perps(
	const MGBSumCurve& crv2	//The second curve
)const{
	MGPosition_list list=crv2.perps(*this);
	return MGPosition_list(list,1,0);	
}

//Compute two straight lines relationship of parallel or coincidence.
//Parallelness of the two is assumed.
//Obtain the relationship when two lines coinside.
//ip of a intersection or nearest point will be returned.
//When this and sl2 do not coincide, MGPSREL_PARALLEL will be returned as
//the function's return value.
MGPSRELATION MGStraight::relation_coincide(
	const MGStraight& sl2,
	MGCCisect& ip
)const{
	MGPSRELATION rl;
	if(straight_type()== MGSTRAIGHT_TYPE::MGSTRAIGHT_SEGMENT)
		return relation_coincide1(sl2,ip);
	else if(sl2.straight_type()== MGSTRAIGHT_TYPE::MGSTRAIGHT_SEGMENT){
		rl=sl2.relation_coincide1(*this,ip);
		ip.exchange12();
		return rl;
	}

	double& t1=ip.param1();
	double& t2=ip.param2();
	MGPosition P1, P2;
	if(straight_type()== MGSTRAIGHT_TYPE::MGSTRAIGHT_UNLIMIT){
		P2=sl2.center();
		t2=sl2.perp_param(P2);
		t1=perp_param(P2);
		P1=eval_position_unlimitted(t1);
	}else if(sl2.straight_type()== MGSTRAIGHT_TYPE::MGSTRAIGHT_UNLIMIT){
		P1=center();
		t1=perp_param(P1);
		t2=sl2.perp_param(P1);
		P2=sl2.eval_position_unlimitted(t2);
	}else{
	//Following are the cases that this and sl2 are both MGSTRAIGHT_HALF_LIMIT.
		if(m_sparam.finite()){
			t1=param_s();
			P1=start_point();
		}else{
			t1=param_e();
			P1=end_point();
		}

		double s1, s2;
		s2=sl2.perp_param(P1);
		if(sl2.in_range(s2)){
			t2=s2;
			P2=sl2.eval_position_unlimitted(s2);
		}else{
			if(sl2.m_sparam.finite()){
				t2=sl2.param_s();
				P2=sl2.start_point();
			}else{
				t2=sl2.param_e();
				P2=sl2.end_point();
			}
			s1=perp_param(P2);
			if(in_range(s1)){
				t1=s1;
				P1=eval_position_unlimitted(t1);
			}
		}
	}

	ip.point()=(P1+P2)*.5;
	MGVector dif=P1-P2;
	if(dif%dif<=MGTolerance::wc_zero_sqr())
		return MGPSRELATION::MGPSREL_COIN;
	else
		return MGPSRELATION::MGPSREL_PARALLEL;
}

//relation_coincide when this is MGSTRAIGHT_SEGMENTt.
MGPSRELATION MGStraight::relation_coincide1(
	const MGStraight& sl2,
	MGCCisect& ip
)const{
	double& t1=ip.param1();
	double& t2=ip.param2();

	MGPosition Ps=start_point();
	MGPosition Pe=end_point();
	double t1s=sl2.perp_param(Ps);
	double t1e=sl2.perp_param(Pe);
	MGInterval t1I(t1s); t1I.expand(t1e);//Parameter range of this in sl2 parameter.
	bool t1s_is_low=(t1s<t1e);//Since t1s may be greater than t1e.

	MGInterval t2I=sl2.param_range();
	if(t1I.high()<=t2I.low()){
		t2=sl2.param_s();
		t1=t1s_is_low ? param_e():param_s();
	}else if(t2I.high()<=t1I.low()){
		t2=sl2.param_e();
		t1=t1s_is_low ? param_s():param_e();
	}else{
		MGInterval comI(t1I); comI&=t2I;
		t2=comI.mid_point();
		t1=perp_param(sl2.eval_position_unlimitted(t2));
	}
	MGPosition P1=eval_position_unlimitted(t1);
	MGPosition P2=sl2.eval_position_unlimitted(t2);
	ip.point()=(P1+P2)*.5;
	MGVector dif=P1-P2;
	if(dif%dif<=MGTolerance::wc_zero_sqr())
		return MGPSRELATION::MGPSREL_COIN;
	else
		return MGPSRELATION::MGPSREL_PARALLEL;
}

//Return two straight line's relationship.
//Whe two are parallel, MGPSREL_PARALLEL or MGPSREL_COIN will be returned.
//MGPSREL_PARALLEL takes place when no overlapping part exists, and 
//MGPSREL_COIN when some parts overlaps.
//MGPSREL_ISECT when they have an intersection.
//MGPSREL_VIRTUAL_ISECT when there is an intersection at extended part of the lines.
//MGPSREL_TORSION when the two do not lie on a same plane.
MGPSRELATION MGStraight::relation(
	const MGStraight& sl2,
	MGCCisect& ip
		//When MGPSREL_PARALLEL or MGPSREL_COIN, a pair of parameter values of the nearest
		//point of the two will be returned.
		//When MGPSREL_ISECT or MGPSREL_VIRTUAL_ISECT, the intersection point parameter values
		//will be returned.
		//When MGPSREL_TORSION, MGPSREL_VIRTUAL_ISECT intersection point parameter value after
		//transformed to lie on the same plane will be returned.
)const{	
	double& t1=ip.param1();
	double& t2=ip.param2();
	if(straight_type()== MGSTRAIGHT_TYPE::MGSTRAIGHT_EMPTY 
		|| sl2.straight_type()== MGSTRAIGHT_TYPE::MGSTRAIGHT_EMPTY){
		t1=param_s();
		t2=sl2.param_s();
	}else{
		// ２直線の関係のチェック、交点なしかどうか。
		if(sl2.sdim()==2 && sdim()==2){
			if(isect2D(sl2, t1, t2))//If parallel.
				return relation_coincide(sl2,ip);
		}else if(MGMZero(m_direction.sangle(sl2.m_direction))){// 平行
			return relation_coincide(sl2,ip);
		}else{
			// 二直線から自身の直線を含む平面を作成する。
			MGUnit_vector N=m_direction*sl2.m_direction;
			//N is normal to both m_direction and s.m_direction.
			MGMatrix mat(3); mat.set_axis(N,2);
			MGStraight tsl1=(*this)*mat; //sl1 is 2D straightline on x-y plane.
			MGStraight tsl2=sl2*mat;		//sl2 is 2D straightline on x-y plane.
			if(tsl1.isect2D(tsl2,t1,t2))//If parallel.
				return relation_coincide(sl2,ip);
		}
	}

//Here intersection point t1 and t2 are obtained after this and sl2 are so transformed
//that two lines lie on the same plane.

	MGPosition P1=eval_position_unlimitted(t1), P2=sl2.eval_position_unlimitted(t2);
	ip.point()=(P1+P2)*.5;
	MGCCRELATION& ccrel=ip.rel(); 
	MGPSRELATION rl;

	// 交点が線分上か延長線上にあるかを調べる。
	MGVector dif=P1-P2;
	if(dif%dif<=MGTolerance::wc_zero_sqr()){
		if(in_range(t1) && sl2.in_range(t2)){
			ccrel= MGCCRELATION::MGCCREL_ISECT;
			rl= MGPSRELATION::MGPSREL_ISECT;
		}else{
			ccrel= MGCCRELATION::MGCCREL_UNKNOWN;
			rl= MGPSRELATION::MGPSREL_VIRTUAL_ISECT;
		}
	}else{
		ccrel= MGCCRELATION::MGCCREL_UNKNOWN;
		rl= MGPSRELATION::MGPSREL_TORSION;
	}
	return rl;
}

// 自身と与えられた平面の関係を調べる。
MGPSRELATION MGStraight::relation(
	const MGPlane& pl,
	MGCSisect& ip
)const{
	MGPSRELATION rl; MGCSRELATION csrel;
	// 直線の方向ベクトルと平面の法線ベクトルの内積を求める。
	double cross = pl.normal()%m_direction;
	if(MGMZero(cross)){	// 直線と平面が平行、又は直線が平面上
		if (MGAZero(pl.distance(root_point()))){	// 直線が平面上
			rl = MGPSRELATION::MGPSREL_COIN; csrel= MGCSRELATION::MGCSREL_COIN;
		}else{			// 直線と平面が平行
			rl = MGPSRELATION::MGPSREL_PARALLEL; csrel= MGCSRELATION::MGCSREL_UNKNOWN;
		}
		ip=MGCSisect(root_point(), 0., MGPosition(0.,0.),csrel);
	}else{				// 交差するとき
		// 直線上のパラメータ値を求める。
		double t=(pl.distance()-pl.normal()%root_point())/cross;
		// 交点を求める。
		MGPosition point = MGPosition(root_point() + t*m_direction);
		MGPosition uv(2);
		pl.on(point,uv);	//Compute plane's parameter value of point.
		ip = MGCSisect(point,t,uv, MGCSRELATION::MGCSREL_IN);

		// 交点が直線上にあるかしらべる。
		if(in_range(t))	 // 直線上にあるとき
			rl = MGPSRELATION::MGPSREL_ISECT;
		else			 // 直線上にないとき
			rl = MGPSRELATION::MGPSREL_VIRTUAL_ISECT;
	}
	return rl;
}

//Compute parallel range of two straight lines.
//Two straight line this and l2 must be parallel.
//Function's return value MGPosition_list list is:
// list.entries() is 0(when no paralle part) or 2(when there is a parallel
// part). When list.entries() is 2, let their entries be P1 and P2.
// Then from P1(0) to P2(0) is the range of this straight line.
// From P1(1) to P2(1) is the range of line2 straight line.
MGPosition_list MGStraight::relation_parallel(const MGStraight& l2)const{
	assert(m_direction.parallel(l2.m_direction));

	MGPosition_list list;
	double s1=perp_param(l2.start_point()), s2=perp_param(l2.end_point());
	double t1=l2.perp_param(start_point()), t2=l2.perp_param(end_point());
	//s1, s2: perp point on this from l2 start and end point
	//t1, t2: perp point on l2 from this start and end point
	MGPosition Ps(param_s(),t1), Pe(param_e(),t2);
	MGPosition Qs(s1,l2.param_s()), Qe(s2,l2.param_e());

	if(straight_type()== MGSTRAIGHT_TYPE::MGSTRAIGHT_EMPTY ||
		l2.straight_type()== MGSTRAIGHT_TYPE::MGSTRAIGHT_EMPTY) return list;

	if(straight_type()== MGSTRAIGHT_TYPE::MGSTRAIGHT_SEGMENT){		//When this is segment.
		if(l2.in_range(t1)){
			list.append(Ps);
			if(l2.in_range(t2)) list.append(Pe);
			else{
				if(in_range(s1)) list.append(Qs);
				else list.append(Qe);
			}
		}else if(l2.in_range(t2)){
			if(in_range(s1)) list.append(Qs);
			else list.append(Qe);
			list.append(Pe);
		}else if(t1*t2<0.){
			list.append(Qs); list.append(Qe);
		}
	}
	else if(straight_type()== MGSTRAIGHT_TYPE::MGSTRAIGHT_UNLIMIT){	//When this is unlimit.
		list.append(Qs); list.append(Qe);
	}else if(l2.straight_type()== MGSTRAIGHT_TYPE::MGSTRAIGHT_UNLIMIT){//When l2 is unlimt.
		list.append(Ps); list.append(Pe);
	}else{										//When this is half limit.
		if(l2.straight_type()== MGSTRAIGHT_TYPE::MGSTRAIGHT_SEGMENT){	//l2 is segment.
			if(in_range(s1)){
				list.append(Qs);
				if(in_range(s2)) list.append(Qe);
				else list.append(Ps);
			}else if(in_range(s2)){
				list.append(Ps); list.append(Qe);
			}
		}else{									//When both are half_limit.
			if(in_range(s1)){
				if(l2.in_range(t1)){
					list.append(Ps); list.append(Qs);
				}else{
					list.append(Qs); list.append(Qe);
				}
			}else if(l2.in_range(t1)){
				list.append(Ps); list.append(Pe);
			}
		}
	}
	return list;
}

//一定オフセット関数
//オフセット方向は、法線方向から見て入力曲線の進行方向左側を正とする。
//法線ベクトルがヌルの場合、始点において曲率中心方向を正とする。
//戻り値は、オフセット曲線が返却される。
std::vector<UniqueCurve> MGStraight::offset(
	double ofs_value,			//オフセット量
	bool principalNormal/// true: Offset direction is to principal normal
						/// false: to binormal
)const {
	MGVector T, N, B;
	double crvtr, torsn;
	Frenet_frame(param_s(), T, N, B, crvtr, torsn);
	MGVector& dir = principalNormal ? N : B;

	std::unique_ptr<MGCurve> offsetSl(new MGStraight(*this));
	*offsetSl += dir * ofs_value;
	std::vector<UniqueCurve> voffset;
	voffset.emplace_back(offsetSl.release());
	return voffset;
}
