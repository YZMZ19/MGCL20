/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#include "StdAfx.h"
#include "mg/Box.h"
#include "mg/Vector.h"
#include "mg/Position.h"
#include "mg/Straight.h"
#include "mg/Ellipse.h"
#include "mg/LBRep.h"
#include "mg/RLBRep.h"
#include "mg/BSumCurve.h"
#include "mg/CompositeCurve.h"
#include "mg/CCisects.h"
#include "mg/CSisects.h"
#include "mg/CParam_list.h"
#include "mg/Position_list.h"
#include "mg/SurfCurve.h"
#include "mg/nlbit.h"
#include "mg/Plane.h"
#include "mg/Sphere.h"
#include "mg/Cylinder.h"
#include "mg/SBRep.h"
#include "mg/RSBRep.h"
#include "mg/BSumSurf.h"
#include "mg/Tolerance.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//Default constructor
MGSurfCurve::MGSurfCurve(const MGSurface& srf, const MGCurve& crv)
:MGCurve(crv),m_surface(&srf),m_curve(crv, crv.param_range()){
		invalidateBox();
}

//Default constructor
MGSurfCurve::MGSurfCurve(
	const MGFSurface& srf, const MGCurve& crv
):MGCurve(crv),m_surface(srf.get_surface_pointer()),
m_curve(crv, crv.param_range()){
	invalidateBox();
}

///Approximate this curve as a MGLBRep curve
///within the tolerance MGTolerance::line_zero().
///When parameter_normalization=0, reparameterization will not done, and
///the evaluation at the same parameter has the same values before and after
///of approximate_as_LBRep.
void MGSurfCurve::approximate_as_LBRep(
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
	const MGCompositeCurve* compo=base_composite();
	if(compo){
		MGLBRep lb2;
		bool CareRange=!(m_curve.is_same_range());
		MGLBRep* lbp=CareRange ? &lb2:&lb;
		int pnormal2=CareRange ? 0:parameter_normalization;
		MGInterval prng=m_curve.param_range();
		int ncurves=compo->number_of_curves();
		for(int i=0; i<ncurves; i++){
			MGLBRep lbi;
			MGLBRep* lbip=i ? &lbi:lbp;
			MGSurfCurve ci(*m_surface,compo->curve(i));
			ci.approximate_as_LBRep(*lbip,ordr,pnormal2,neglectMulti);
			if(CareRange && (i==0 || i==ncurves-1))
				lbip->limit(prng);
			if(i){
				int which;
				double ratio;
				int cn=lbp->continuity(lbi,which,ratio);
				if(cn<0)
					cn=0;
				lbp->connect(cn,2,lbi);
			}
		}
		if(CareRange)
			lbp->approximate_as_LBRep(lb,ordr,parameter_normalization,neglectMulti);
	}else{
		int knew=ordr;
		if(!knew)
			knew=4;
		MGCurve::approximate_as_LBRep(lb,knew,parameter_normalization,neglectMulti);
	}
}

//Test if m_curve is MGCompositeCurve. If composite, return
//the pointer. If not, return null.
const MGCompositeCurve* MGSurfCurve::base_composite()const{
	const MGCompositeCurve* compo
		=dynamic_cast<const MGCompositeCurve*>(m_curve.base_curve());
	return compo;
}

//Changing this object's space dimension.
void MGSurfCurve::change_dimension(
	int dim,		// new space dimension
	int start1, 		// Destination order of new object.
	int start2		// Source order of this object.
){ 
	assert(false);//This cannot be used.
}

//Return minimum box that includes whole of the curve.
//曲線部分を囲むボックスを返す。
void MGSurfCurve::compute_box(MGBox& bx)const{
	bx=box_limitted(param_range());
}

//Return minimum box that includes the curve of parameter interval.
// 入力のパラメータ範囲の曲線部分を囲むボックスを返す。
MGBox MGSurfCurve::box_limitted(
	const MGInterval& l // Parameter Range of the curve.
)const{
	//uv-curveを囲むボックス(uv-interval)を内包する
	//m_surfaceのボックス。
	return m_surface->box_limitted(m_curve.box_limitted(l));
}

//Construct new curve object by copying to newed area.
//User must delete this copied object by "delete".
MGSurfCurve* MGSurfCurve::clone()const{
	return new MGSurfCurve(*this);
}

// Evaluate n'th derivative data. n=0 means positional data evaluation.
// とりあえず3次微分まで展開して計算しています。
// 続きは後で考えましょう。
MGVector MGSurfCurve::eval(
	double t,	// Parameter value.
	int nderiv,	// Order of Derivative.
	int left	//Left continuous(left=true) or right continuous(left=false).
)const{
	if(nderiv>=4)
		return MGVector(sdim(), 0.);

	MGPosition uv= m_curve.eval(t, 0, left);//母曲面のuvパラメータを求める。
	if(nderiv<=0)
		return m_surface->eval(uv, 0, 0);

	MGVector dfdt = m_curve.eval(t, 1, left);
	double du = dfdt[0];
	double dv = dfdt[1];
	if(nderiv==1)
		return m_surface->eval(uv, 1, 0)*du	+ m_surface->eval(uv, 0, 1)*dv;

	MGVector d2fdt2 = m_curve.eval(t, 2, left);
	double du2 = d2fdt2[0];
	double dv2 = d2fdt2[1];
	if(nderiv==2)
		return m_surface->eval(uv, 2, 0) * du*du
		+ m_surface->eval(uv, 1, 1) * du*dv *2.
		+ m_surface->eval(uv, 1, 0) * du2
		+ m_surface->eval(uv, 0, 1) * dv2
		+ m_surface->eval(uv, 0, 2) * dv*dv;

	MGVector d3fdt3 = m_curve.eval(t, 3, left);
	return m_surface->eval(uv, 3, 0) * du*du*du
		+ m_surface->eval(uv, 0, 3) * dv*dv*dv
		+ m_surface->eval(uv, 2, 1) * du*du*dv *3.
		+ m_surface->eval(uv, 1, 2) * du*dv*dv *3.
		+ m_surface->eval(uv, 2, 0) * du2*du * 3.
		+ m_surface->eval(uv, 0, 2) * dv2*dv * 3.
		+ m_surface->eval(uv, 1, 1) *(du2*dv+ dv2*du)*3.
		+ m_surface->eval(uv, 1, 0) * d3fdt3[0]
		+ m_surface->eval(uv, 0, 1) * d3fdt3[1];
}

//Extrapolate this curve by an (approximate) chord length.
//The extrapolation is C2 continuous.
void MGSurfCurve::extend(
	double length,	//approximate chord length to extend. 
	bool start		//Flag of which point to extend, start or end point of the line.
					//If start is true extend on the start point.
){
	double ts=param_s(), te=param_e();
	double t=start ? ts:te;
	MGVector V1=eval(t,1);
	MGVector V2=m_curve.eval(t,1);
	double tan_ratio=V1.len()/V2.len();
	length/=tan_ratio;
	if(start)
		length*=-1.;

	m_curve.extend(length,start);
	double prm[4]={m_surface->param_s_u(),m_surface->param_s_v(),
		m_surface->param_e_u(),m_surface->param_e_v()};
	for(int i=0; i<4; i++){
		MGCParam_list plist=m_curve.isect_1D(prm[i],i%2);
		if(plist.size()){
			if(start)
				m_curve.limit(MGInterval(plist.back(),te));
			else
				m_curve.limit(MGInterval(ts,plist.front()));
		}
	}
}

//Test if input parameter value is inside parameter range of the line.
bool MGSurfCurve::in_range(double t) const{
	return m_curve.in_range(t);
}

//Provide divide number of curve span for function intersect.
int MGSurfCurve::intersect_dnum() const{
	int n1=m_curve.intersect_dnum(), n2,n3;
	const MGRSBRep* surf_is_MGRSBRep=dynamic_cast<const MGRSBRep*>(m_surface);

	const MGBox& uvbox=m_curve.box();
	double u0=uvbox[0].low_point(), u1=uvbox[0].high_point();
	const MGKnotVector& tu=m_surface->knot_vector_u();
	int ordr=tu.order();
	int i1=tu.locate(u0), i2=tu.locate(u1,1);
	if(i2<i1)
		i2=i1;
	if(surf_is_MGRSBRep)
		n2=(i2-i1+1)*ordr;
	else
		n2=(i2-i1+1)*(ordr-1);

	double v0=uvbox[1].low_point(), v1=uvbox[1].high_point();
	const MGKnotVector& tv=m_surface->knot_vector_v();
	ordr=tv.order();
	i1=tv.locate(v0); i2=tv.locate(v1,1);
	if(i2<i1)
		i2=i1;
	if(surf_is_MGRSBRep)
		n3=(i2-i1+1)*ordr;
	else
		n3=(i2-i1+1)*(ordr-1);

	//uvカーブのintersect_dnum() + Surfaceのintersect_dnum()のうち
	//最大のもので代用する。
	if(n1<n2){
		if(n2>n3)
			return n2;
		else
			return n3;
	}else if(n3>n1)
		return n3;
	else
		return n1;
}

//Intersection with a curve.
MGCCisects MGSurfCurve::isect(const MGCurve& curve2) const{
	MGCCisects list(this, &curve2);
	isect_of_each(curve2,list);
	return list;
}
MGCCisects MGSurfCurve::isect(const MGStraight& curve2) const{
	MGCCisects list(this, &curve2);
	isect_of_each(curve2, list);
	return list;
}
MGCCisects MGSurfCurve::isect(const MGLBRep& curve2) const{
	MGCCisects list(this, &curve2);
	isect_of_each(curve2, list);
	return list;
}
MGCCisects MGSurfCurve::isect(const MGSurfCurve& curve2) const{
	MGCCisects list(this, &curve2);
	isect_of_each(curve2, list);
	return list;
}

//isect of this SurfCurve whose m_curve is not a MGTrimmedCurve of MGCompositeCurve.
MGCCisects MGSurfCurve::isect_noCompo(const MGCurve& curve2)const{
	MGCCisects list=curve2.isect_with_noCompoSC(*this);
	list.exchange12();
	return list;
}

//isect of this SurfCurve whose m_curve is not a MGTrimmedCurve of MGCompositeCurve.
MGCCisects MGSurfCurve::isect_noCompo(const MGSurfCurve& curve2)const{
	MGCCisects list(this, &curve2);
	if(curve2.base_surface()==base_surface()){
		MGCCisects pisects=m_curve.isect(curve2.m_curve);
		if(pisects.size()==0)
			return list;
		MGCCisects::iterator i=pisects.begin(), ie=pisects.end();
		for(; i!=ie; i++){
			MGCCisect& isi=isectCast<MGCCisect>(i);
			double t1=isi.param1(), t2=isi.param2();
			MGPosition P=eval(t1);
			if(P!=curve2.eval(t2)){
				MGPosition t12(t1,t2);
				perp_guess(1.,-1.,curve2,1.,-1.,t1,t2,t12);
				t1=t12[0]; t2=t12[1];
				P=eval(t1);
			}
			list.emplace_back(new MGCCisect(P,t1,t2));
		}
		return list;
	}

	const MGCompositeCurve* compo2=curve2.base_composite();
	if(compo2){
		double ts=curve2.m_curve.param_s(), te=curve2.m_curve.param_e();
		MGCompositeCurve::const_iterator
			icurve=compo2->begin(), endcurve=compo2->end();
		const MGSurface& srf2=*(curve2.m_surface);
		for(; icurve!=endcurve; icurve++){
			const MGCurve& crvi=**icurve;
			if(crvi.param_e()<=ts || te<=crvi.param_s())
				continue;
			MGSurfCurve scrv2(srf2,crvi);
			MGCCisects listi=intersect(scrv2);
			MGCCisects::iterator j=listi.begin(), je=listi.end();
			for(; j!=je; j++){
				auto& isecj=isectCast<MGCCisect>(j);
				if(in_range(isecj.param1()))
					list.emplace_back(j->release());
			}
		}
		return list;
	}

	return intersect(curve2);
}

//isect of this SurfCurve whose m_curve is not a MGTrimmedCurve of MGCompositeCurve.
MGCCisects MGSurfCurve::isect_noCompo(const MGStraight& sline) const{
	MGCCisects list(this, &sline);

	double tmid=(param_s()+param_e())*0.5;
	MGVector N=m_surface->normal(m_curve.eval(tmid));
	MGPlane plane(sline.direction(), N, sline.root_point());
	//The plane that includes sline.
	MGCSisects cslist=isect(plane);//直線を含む平面と面上線との交点を求める。

	MGPosition p; double t1,t2;
	MGCSisects::iterator i;
	for(i=cslist.begin(); i!=cslist.end(); i++){
		auto& iseci=isectCast<MGCSisect>(i);
		t1=iseci.param_curve(); //交点のSurfCurve上のパラメータ値を求める。
		p=iseci.point();//交点の座標値を求める(eval(t1)よりはこちらのほうが良い？)
		if(sline.on(p,t2)) //交点が直線sline上に乗っていれば面上線と直線の交点になる。
			list.append(p,t1,t2);
	}
	return list;
}

//isect of each elements of this m_curve,
//if m_curve is MGTrimmedCurve of a MGCompositeCurve.
void MGSurfCurve::isect_of_each(
	const MGCurve& curve2,	//The isect objective curve.
	MGCCisects& list	//Obtained isect will be appended.
)const{
	const MGCompositeCurve* compo=base_composite();
	if(!compo){
		list=isect_noCompo(curve2);
		return;
	}

	MGCompositeCurve::const_iterator
		icurve=compo->begin(), endcurve=compo->end();
	double ts=m_curve.param_s(), te=m_curve.param_e();
	for(; icurve!=endcurve; endcurve++){
		const MGCurve& crvi=**icurve;
		if(crvi.param_e()<=ts || te<=crvi.param_s())
			continue;
		MGSurfCurve scrv(*m_surface,crvi);
		MGCCisects listi=scrv.isect_noCompo(curve2);
		MGCCisects::iterator j=listi.begin(), je=listi.end();
		for(; j!=je; j++){
			auto& isecj=isectCast<MGCCisect>(j);
			if(in_range(isecj.param1()))
				list.emplace_back(j->release());
		}
	}
}

//Intersection with a Surface
MGCSisects MGSurfCurve::isect(const MGSurface& surf)const{
	MGCSisects list(this,&surf);
	isect_of_each(surf,list);
	return list;
}
MGCSisects MGSurfCurve::isect(const MGPlane& surf)const{
	MGCSisects list(this, &surf);
	isect_of_each(surf, list);
	return list;
}

//isect of this SurfCurve whose m_curve is not a MGTrimmedCurve of MGCompositeCurve.
MGCSisects MGSurfCurve::isect_noCompo(const MGSurface& surf)const{
	return surf.isect_with_noCompoSC(*this);
}

class MGSurfCurvePlaneDist{
	const MGSurfCurve* m_curve;
	const MGPlane* m_plane;
public:
	MGSurfCurvePlaneDist(const MGSurfCurve* scrv, const MGPlane* plane)
		:m_curve(scrv), m_plane(plane){;};
	double operator()(double t)const{
		return m_plane->distance(m_curve->eval(t));};
};

// 面上線と平面の交点を求める。
MGCSisects MGSurfCurve::isect_noCompo(const MGPlane& pl) const{
	int ndiv=intersect_dnum();
	double t0=param_s(), t1=param_e();
	double delta=(t1-t0)/double(ndiv);

	double error=MGTolerance::wc_zero();
	MGCSisects list(this,&pl);

	MGPosition Ppre, Paft;
	double tpre=t1,taft=t0;
	double dpre=pl.distance(Ppre=eval(tpre));
	double daft=pl.distance(Paft=eval(taft));
	if(fabs(daft)<error*.5) list.append(Paft,t0,pl.uv(Paft));
	if(fabs(dpre)<error*.5) list.append(Ppre,t1,pl.uv(Ppre));

	//Prepare for mgNlbit.
	MGSurfCurvePlaneDist pdist(this,&pl);

	//Iterate by checking singned distance dpre and daft from the plane.
	//When dpre and daft have different signs, an intersection point must 
	//lie between tpre and taft.
	int i=0;
	while(i<=ndiv){
		dpre=daft; tpre=taft; Ppre=Paft;
		while(fabs(dpre)<=error){
			list.append(Ppre,tpre,pl.uv(Ppre));
			if(i>=ndiv) break;
			i++; tpre=t0+delta*double(i);
			Ppre=eval(tpre); dpre=pl.distance(Ppre);
		}
		if(i>=ndiv) break;
		i++; taft=t0+delta*double(i);
		Paft=eval(taft); daft=pl.distance(Paft);
		if(fabs(daft)<=error) continue;
		else if(dpre*daft<0.){
		//Now there exists a solution between tpre and taft.
			int ier;
			double x=mgNlbit(pdist, tpre,taft, error, 20, ier);
			MGPosition Pt=eval(x);
			list.append(Pt,x,pl.uv(Pt));
		}
	}

	return list;
}

//isect of each elements of this m_curve,
//if m_curve is MGTrimmedCurve of a MGCompositeCurve.
void MGSurfCurve::isect_of_each(
	const MGSurface& surf,	//The isect objective surface.
	MGCSisects& list	//Obtained isect will be appended.
)const{
	const MGCompositeCurve* compo=base_composite();
	if(!compo){
		list=surf.isect_with_noCompoSC(*this);
		return;
	}

	MGCompositeCurve::const_iterator
		icurve=compo->begin(), endcurve=compo->end();
	double ts=m_curve.param_s(), te=m_curve.param_e();
	for(; icurve!=endcurve; endcurve++){
		const MGCurve& crvi=**icurve;
		if(crvi.param_e()<=ts || te<=crvi.param_s())
			continue;
		MGSurfCurve scrv(*m_surface,crvi);
		MGCSisects listi=scrv.isect_noCompo(surf);
		MGCSisects::iterator j=listi.begin(), je=listi.end();
		for(; j!=je; j++){
			auto& isecj=isectCast<MGCSisect>(j);
			if(in_range(isecj.param_curve()))
				list.emplace_back(j->release());
		}
	}
}

//Access to i-th element of knot.
double MGSurfCurve::knot(int i) const{ return m_curve.knot(i);}

// Return ending parameter value.
double MGSurfCurve::param_e() const{
	// uvカーブのパラメータが面上線のパラメータになる。
	return m_curve.param_e();
}

//Normalize parameter value t to the nearest knot if their distance is
//within tolerance.
double MGSurfCurve::param_normalize(double t) const{
	// uvカーブのパラメータが面上線のパラメータになる。
	return m_curve.param_normalize(t);
}

//Return parameter range of the curve(パラメータ範囲を返す)
MGInterval MGSurfCurve::param_range() const{
	// uvカーブのパラメータが面上線のパラメータになる。
	return m_curve.param_range();
}

// Return starting parameter value.
double MGSurfCurve::param_s() const{
	// uvカーブのパラメータが面上線のパラメータになる。
	return m_curve.param_s();
}
	
//Compute part of this curve from parameter t1 to t2.
//Returned is the pointer to newed object, and so should be deleted
//by calling program, or memory leaked.
MGSurfCurve* MGSurfCurve::part(double t1, double t2, int multiple) const{
	assert(t2-t1>param_error());
	MGSurfCurve* ret = new MGSurfCurve(*this);
	if(multiple){
		MGCurve* crv=m_curve.part(t1,t2,multiple);
		ret->m_curve=MGTrimmedCurve(*crv,crv->param_range());
		delete crv;
	}else{
		ret->m_curve.limit(MGInterval(t1,t2));
	}
	ret->invalidateBox();
	return ret;
}

//Round t into curve's parameter range.
// 入力パラメータをパラメータ範囲でまるめて返却する。
double MGSurfCurve::range(double t) const{
	return m_curve.range(t);
}

//Return space dimension
int MGSurfCurve::sdim() const{
	return m_surface->sdim();
}

//Return curve type(曲線のタイプを返す)
MGCURVE_TYPE MGSurfCurve::type() const{
	return MGCURVE_TYPE::MGCURVE_SURFACE;
}

//Operator overload(演算子多重定義)

//Assignment.
//When the leaf object of this and crv2 are not equal, this assignment
//does nothing.
MGSurfCurve& MGSurfCurve::operator=(const MGGel& gel2){
	const MGSurfCurve* gel2_is_this=dynamic_cast<const MGSurfCurve*>(&gel2);
	if(gel2_is_this)
		operator=(*gel2_is_this);
	return *this;
}
MGSurfCurve& MGSurfCurve::operator=(MGGel&& gel2){
	MGSurfCurve* gel2_is_this=dynamic_cast<MGSurfCurve*>(&gel2);
	if(gel2_is_this)
		operator=(std::move(*gel2_is_this));
	return *this;
}

//Logical operator overload(論理演算子多重定義)
//Test if two curves are equal.
// 与曲線と自身が等しいかの比較判定を行う。
// 等しくなるのは自分自身と、自分自身から作ったコピーのみ。
bool MGSurfCurve::operator==(const MGSurfCurve& crv)const{
	return (m_surface==crv.m_surface) && (m_curve==crv.m_curve);
}
bool MGSurfCurve::operator<(const MGSurfCurve& gel2)const{
	if(m_surface==gel2.m_surface)
		return m_curve<gel2.m_curve;
    return (*m_surface)<*(gel2.m_surface);
}
bool MGSurfCurve::operator==(const MGGel& gel2)const{
	const MGSurfCurve* gel2_is_this=dynamic_cast<const MGSurfCurve*>(&gel2);
	if(gel2_is_this)
		return operator==(*gel2_is_this);
	return false;
}
bool MGSurfCurve::operator<(const MGGel& gel2)const{
	const MGSurfCurve* gel2_is_this=dynamic_cast<const MGSurfCurve*>(&gel2);
	if(gel2_is_this)
		return operator<(*gel2_is_this);
	return identify_type() < gel2.identify_type();
}

//Obtain parameter value if this curve is negated by "negate()".
double MGSurfCurve::negate_param(double t)const{
	assert(false);
	return 0.0;
}

MGPosition MGSurfCurve::negate_param(const MGPosition& t)const{
	assert(false);
	return MGPosition();
}

//Test if given point is on the curve or not. If yes, return parameter
//value of the curve. Even if not, return nearest point's parameter.
// 指定点が自身上にあるかを調べる。曲線上にあれば，そのパラメーター値を，
// なくても最近傍点のパラメータ値を返す。
// Function's return value is >0 if the point is on the curve,
// and 0 if the point is not on the curve.
/*bool MGSurfCurve::on(
	const MGPosition& P,//Point(指定点)
	double&	t			//Parameter of the curve(パラメータ)
	) const{
	double t0=param_s(), t1=param_e();
	MGVector V(eval(t0)-P);
	double lensqrMin0=V%V, tmin=t0;
	double lensqrMin=lensqrMin0;
	V=eval(t1)-P;
	double lensqr=V%V;
	if(lensqr<lensqrMin){
		lensqrMin=lensqr; tmin=t1;
	}

	int ndiv=intersect_dnum();
	double tg=t0;
	double delta=(t1-t0)/ndiv;

	for(int i=1; i<ndiv; i++){
		tg+=delta;
		V=eval(tg)-P;
		lensqr=V%V;
		if(lensqr<lensqrMin){lensqrMin=lensqr; tmin=tg;}
	}
	if(perp_guess(1.,0.,P,tmin,t)){
		V=eval(t)-P;
		lensqr=V%V;
		if(lensqr<=lensqrMin) lensqrMin=lensqr;
		else t=tmin;
	}else{
		//Try again by increasing ndiv.
		ndiv*=3;
		delta=(t1-t0)/ndiv;
		tg=t0;
		lensqrMin=lensqrMin0;
		for(i=1; i<ndiv; i++){
			tg+=delta;
			V=eval(tg)-P;
			lensqr=V%V;
			if(lensqr<lensqrMin){lensqrMin=lensqr; tmin=tg;}
		}
		if(perp_guess(1.,0.,P,tmin,t)){
			V=eval(t)-P;
			lensqr=V%V;
			if(lensqr<=lensqrMin) lensqrMin=lensqr;
			else t=tmin;
		}else t=tmin;
	}
	return lensqrMin<=MGTolerance::wc_zero_sqr();
}*/

//Compute all the perpendicular points of this curve and the second one.
//That is, if f(s) and g(t) are the points of the two curves f and g,
//then obtains points where the following conditions are satisfied:
//  fs*(f-g)=0.    gt*(g-f)=0.
//Here fs and gt are 1st derivatives at s and t of f and g.
//MGPosition P in the MGPosition_list contains this and crv's parameter
//as:     P(0)=this curve's parameter, P(1)=crv's parameter value.
MGPosition_list MGSurfCurve::perps(
	const MGCurve& crv2		//The second curve
)const{
	MGPosition_list list;
	perps_of_each(crv2,list);
	return list;	
}
MGPosition_list MGSurfCurve::perps(
	const MGStraight& crv2		//The second curve
)const{
	MGPosition_list list;
	perps_of_each(crv2,list);
	return list;	
}
MGPosition_list MGSurfCurve::perps(
	const MGRLBRep& crv2		//The second curve
)const{
	MGPosition_list list;
	perps_of_each(crv2,list);
	return list;	
}
MGPosition_list MGSurfCurve::perps(
	const MGEllipse& crv2		//The second curve
)const{
	MGPosition_list list;
	perps_of_each(crv2,list);
	return list;	
}
MGPosition_list MGSurfCurve::perps(
	const MGLBRep& crv2		//The second curve
)const{
	MGPosition_list list;
	perps_of_each(crv2,list);
	return list;	
}
MGPosition_list MGSurfCurve::perps(
	const MGSurfCurve& crv2		//The second curve
)const{
	MGPosition_list list;
	perps_of_each(crv2,list);
	return list;	
}
MGPosition_list MGSurfCurve::perps(
	const MGBSumCurve& crv2		//The second curve
)const{
	MGPosition_list list;
	perps_of_each(crv2,list);
	return list;	
}

//isect of this SurfCurve whose m_curve is not a MGTrimmedCurve of MGCompositeCurve.
MGPosition_list MGSurfCurve::perps_noCompo(const MGCurve& curve2)const{
	MGPosition_list list=curve2.perps_with_noCompoSC(*this);
	return MGPosition_list(list,1,0);	
}

//perpendicular points of each elements of this m_curve,
//if m_curve is MGTrimmedCurve of a MGCompositeCurve.
void MGSurfCurve::perps_of_each(
	const MGCurve& curve2,	//The perps objective curve.
	MGPosition_list& list	//Obtained perpendicular points will be appended.
)const{
	const MGCompositeCurve* compo=base_composite();
	if(!compo){
		list=perps_noCompo(curve2);
		return;
	}

	MGCompositeCurve::const_iterator
		icurve=compo->begin(), endcurve=compo->end();
	double ts=m_curve.param_s(), te=m_curve.param_e();
	for(; icurve!=endcurve; endcurve++){
		const MGCurve& crvi=**icurve;
		if(crvi.param_e()<=ts || te<=crvi.param_s())
			continue;
		MGSurfCurve scrv(*m_surface,crvi);
		MGPosition_list listi=scrv.perps_noCompo(curve2);
		MGPosition_list::iterator j=listi.begin(), je=listi.end();
		for(; j!=je; j++){
			if(in_range(j->ref(0)))
				list.append(*j);
		}
	}
}

///Divide this curve at the designated knot multiplicity point.
///Function's return value is the number of the curves after divided.
int MGSurfCurve::divide_multi(
	std::vector<UniqueCurve>& crv_list,	//divided curves are appended.
	int multiplicity	///<designates the multiplicity of the knot to divide at,
						///<When multiplicity<=0, order()-1 is assumed,
						///<When multiplicity>=order(), order() is assumed.
)const{
	std::vector<UniqueCurve> clist2;
	int ncurve=m_curve.divide_multi(clist2,multiplicity);
	for(int i=0; i<ncurve; i++){
		crv_list.emplace_back(new MGLBRep(MGSurfCurve(*m_surface,*(clist2[i]))));
	}
	return ncurve;
}

//Exchange ordering of the coordinates.
//Exchange coordinates (i) and (j).
void MGSurfCurve::coordinate_exchange(int i, int j){
	assert(false);
}

//copy as a newed curve. The new curve will be MGLBRep or MGRLBRep.
//When original curve was a MGRLBRep, the new curve will be a MGRLBRep.
//Otherwise,  the new curve will be a MGLBRep.
//Returned object must be deleted.
MGCurve* MGSurfCurve::copy_as_nurbs()const{return new MGLBRep(*this);}

//Construct new curve object by changing
//the original object's space dimension.
//User must delete this copied object by "delete".
MGCurve* MGSurfCurve::copy_change_dimension(
	int sdim,			// new space dimension
	int start1, 		// Destination order of new line.
	int start2		// Source order of this line.
)const{
	MGLBRep tmp(*this);
	return tmp.copy_change_dimension(sdim,start1,start2);
}

//Update this by limiting the parameter range of the curve.
// 自身に指定したパラメータ範囲のｌｉｍｉｔをつける。
void MGSurfCurve::limit(const MGInterval& I){
	m_curve.limit(I);
	invalidateBox();
}

//Negate the curve direction(曲線の方向を反転する)
void MGSurfCurve::negate(){
	assert(false);
}

//Unlimit parameter range of the curve(limitをはずす)
MGCurve& MGSurfCurve::unlimit(){
	m_curve.unlimit();
	invalidateBox();
	return *this;
}

//Unlimit parameter range of the curve to the end point direction
//(終点方向にlimitをはずす)
MGCurve& MGSurfCurve::unlimit_end(){
	m_curve.unlimit_end();
	invalidateBox();
	return *this;
}

//Unlimit parameter range of the curve to the start point direction
//(始点方向にlimitをはずす)
MGCurve& MGSurfCurve::unlimit_start(){
	m_curve.unlimit_start();
	invalidateBox();
	return *this;
}

//Update the curve by translation.
// 与ベクトルだけ曲線を平行移動して自身とする。
MGSurfCurve& MGSurfCurve::operator+= (const MGVector& v){
//	assert(false);
	return *this;
}

//Update the curve by translation.
// 与ベクトルだけ曲線をマイナス方向に平行移動して自身とする。
MGSurfCurve& MGSurfCurve::operator-= (const MGVector& v){
//	assert(false);
	return *this;
}

//Update the curve by multiplying scale.
// 与えられたスケールを曲線にかける。
MGSurfCurve& MGSurfCurve::operator*= (double scale){
//	assert(false);
	return *this;
}

//Update the curve by transformation of matrix.
// 与えられた変換で直線の変換を行い自身の直線とする。
MGSurfCurve& MGSurfCurve::operator*= (const MGMatrix& mat){
//	assert(false);
	return *this;
}

//Update the curve by transformation of transf.
// 与えられた変換で曲線のトランスフォームを行い自身とする。
MGSurfCurve& MGSurfCurve::operator*= (const MGTransf& tr){
//	assert(false);
	return *this;
}

MGSurfCurve MGSurfCurve::operator+ (const MGVector& v) const{
	assert(false); return *this;
}
MGSurfCurve operator+ (const MGVector& v, const MGSurfCurve& cv2){
	assert(false); return cv2;
}
MGSurfCurve MGSurfCurve::operator- (const MGVector& v) const{
	assert(false); return *this;
}
MGSurfCurve MGSurfCurve::operator* (double scale) const{
	assert(false); return *this;
}
MGSurfCurve operator* (double scale, const MGSurfCurve& cv2){
	assert(false); return cv2;
}
MGSurfCurve MGSurfCurve::operator* (const MGMatrix& mat) const{
	assert(false); return *this;
}
MGSurfCurve MGSurfCurve::operator* (const MGTransf& tr) const{
	assert(false); return *this;
}

void add_knots(
	int istart,int n,const MGKnotVector& st, double st0, double st1,
	//istart=id of st, n=num of knots found between c.eval(ts) and c.eval(te)
	//st0=u or v value of c.eval(ts), st1=u or v value of c.eval(te).
	double ts, double te, MGNDDArray& tau
	//ts - te is a knot span of tau where knots be inserted.
){
	double span=te-ts;
	int iend;
	if(n<0){
		tau*=-1.;
		double tsave=ts;
		ts=te*-1.;
		te=tsave*-1.;
		iend=istart;
		istart+=n;
		double save=st0;
		st0=st1;
		st1=save;
	}else
		iend=istart+n;

	double stspan=st1-st0;
	for(int i=istart+1; i<iend; i++){
		double stnex=(i-1>istart) ? st[i-1]:st0;
		stnex+=st[i];
		stnex+=st[i+1];
		stnex+=(i+2<iend) ? st[i+2]:st1;
		stnex/=4.;
		double delta=span*(stnex-st0)/stspan;
		double tnex=ts+delta;
		tau.add_data(tnex);
	}
	if(n<0)
		tau*=-1.;
}

#define INCREASE_NUM 10
///Get data points for approximate_as_LBRep2.
void MGSurfCurve::data_points_for_approximate_as_LBRep2(
	int is, int ie,//approximation parameter range, from knot_vector()[is] to [ie].
	MGKnotVector& t,//New knot configuration will be output.
				//t's order is input. other information of t will be updated.
	MGNDDArray& tau,//Data point for t will be output.
	bool neglectMulti///<Indicates if multiple knots be kept.
		///< true: multiplicity is removed.
		///< false: multiplicity is kept.
)const{
	const MGKnotVector& told=knot_vector();
	int kold=told.order();
	int nnew=ie-is+kold-1;
	MGKnotVector ttmp(is-kold+1,nnew,told);
	tau.buildByKnotVector(ttmp);
	MGNDDArray tau2(tau);

	const MGSurface& srf=*(base_surface());
	const MGCurve& c=*(base_curve());

	const MGKnotVector& tu=srf.knot_vector_u();
	double uerr=tu.param_error();
	const MGKnotVector& tv=srf.knot_vector_v();
	double verr=tv.param_error();
	double err=ttmp.param_error();
	double te=tau2(0);
	MGPosition uvim1=c.eval(te);
	for(int i=1; i<nnew; i++){
		double	ts=te;			//スパンの始点
		te=tau2(i);		//スパンの終点
		if(te-ts <= err)
			continue;	//マルチノットのときの処理

		MGPosition uvi=c.eval(te);
		double u0=uvim1[0], u1=uvi[0];
		double v0=uvim1[1], v1=uvi[1];
		int iuim1=tu.locate(u0);
		int ivim1=tv.locate(v0);
		int iui=tu.locate(u1,1);
		int ivi=tv.locate(v1,1);
		int nu=iui-iuim1, nv=ivi-ivim1;
		uvim1=uvi;
		int nua=abs(nu), nva=abs(nv);
		if(nua<=1 && nva<=1)
			continue;

		if(nua>=nva){
			if(fabs(u1-u0)<=uerr)
				continue;
			add_knots(iuim1,nu,tu,u0,u1,ts,te,tau);
		}else{
			if(fabs(v1-v0)<=verr)
				continue;
			add_knots(ivim1,nv,tv,v0,v1,ts,te,tau);
		}
	}
	int nIncreased=tau.length();
	tau.change_number(nIncreased*INCREASE_NUM);
	tau.remove_too_near();
	int knew=t.order();
	if(knew>tau.length())
		knew=tau.length();
	t=MGKnotVector(tau,knew);
}
