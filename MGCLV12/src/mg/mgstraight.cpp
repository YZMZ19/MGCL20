/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#include "StdAfx.h"
#include "mg/Box.h"
#include "mg/Position.h"
#include "mg/Transf.h"
#include "mg/Unit_vector.h"
#include "mg/CParam_list.h"
#include "mg/Straight.h"
#include "mg/LBRep.h"
#include "mg/RLBRep.h"
#include "mg/Plane.h"
#include "mg/Tolerance.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//MGStraighy Class Implementation
//

//
// Constructor
	
//Void constructor
MGStraight::MGStraight()
: MGCurve()	,m_sparam(0.0),m_endparam(-1.),m_knotV(nullptr){;}

//Copy constructor.
MGStraight::MGStraight(const MGStraight& sl)
:MGCurve(sl), m_direction(sl.m_direction), m_endparam(sl.m_endparam),
m_knotV(nullptr), m_root_point(sl.m_root_point), m_sparam(sl.m_sparam){;}

//Move constructor.
MGStraight::MGStraight(MGStraight&& sl)
:MGCurve(std::move(sl)), m_direction(std::move(sl.m_direction))
, m_sparam(sl.m_sparam), m_endparam(sl.m_endparam)
, m_knotV(sl.m_knotV), m_root_point(std::move(sl.m_root_point)){
	sl.m_knotV=nullptr;
}

//Straight specifying all the member data. All of the data are employed as 
//member data of this straight.
MGStraight::MGStraight (
	const MGEReal& endparam,	//end parameter value
	const MGEReal& sparam,		//start parameter value
	const MGVector& direction,	//Direction, which will be the direction of this.
	const MGPosition& Origin	//Origin
):m_direction(direction), m_endparam(endparam),
m_knotV(nullptr), m_root_point(Origin), m_sparam(sparam){;}

//Straight specifying all the member data. All of the data are employed as 
//member data of this straight.
MGStraight::MGStraight (
	double endparam,			//end parameter value
	double sparam,				//start parameter value
	const MGVector& direction,	//Direction, which will be the direction of this.
	const MGPosition& Origin	//Origin
):m_direction(direction), m_endparam(endparam),
m_knotV(nullptr), m_root_point(Origin), m_sparam(sparam){;}

// 始点、方向ベクトル、直線のタイプを指定して直線を生成する。
MGStraight::MGStraight (
	MGSTRAIGHT_TYPE t,
	const MGVector& v,
	const MGPosition& p
):MGCurve(),m_root_point(p),m_direction(v.normalize())
,m_sparam(0.0),m_endparam(MGINFINITE_TYPE::MGINFINITE_PLUS),m_knotV(nullptr){
	if (t== MGSTRAIGHT_TYPE::MGSTRAIGHT_SEGMENT) m_endparam=v.len();
	else if(t== MGSTRAIGHT_TYPE::MGSTRAIGHT_EMPTY)
		m_endparam=-1.;
	else if(t== MGSTRAIGHT_TYPE::MGSTRAIGHT_UNLIMIT)
		m_sparam=MGEReal(MGINFINITE_TYPE::MGINFINITE_MINUS);

	int dimv=v.sdim(), dimp=p.sdim();
	if(dimv<dimp) m_direction=MGVector(dimp,v);
	else if(dimv>dimp) m_root_point=MGVector(dimv,p);
}

// ２点から直線を生成する。
//Parameter value of the start point is set to be 0.
MGStraight::MGStraight(
	const MGPosition& end,		//End point.
	const MGPosition& start		//Start point.
):MGCurve(),m_root_point(start),m_sparam(0.0),m_knotV(nullptr){
	// ２点から単位ベクトルを作成する。
	MGVector v(end,start);
	m_direction = v.normalize();

	// 単位方向ベクトルをもつので
	// 作成したベクトルの長さが終点のパラメータ値
	m_endparam=v.len();
	if(m_root_point.sdim()<m_direction.sdim())
		m_root_point=MGVector(m_direction.sdim(),start);
}

//始終点の座標、パラメータ値から直線を生成する。
//MGSTRAIGHT_SEGMENT straight from two points.
//Start point is start and end point is end.
//In this version, can specify start and end parameter values.
MGStraight::MGStraight (
	const MGPosition& Pe,const MGPosition& Ps, 
	const double Te,	const double Ts
):MGCurve(),m_sparam(Ts),m_endparam(Te),m_knotV(nullptr){
	assert(Ts<Te);
	// ２点から方向ベクトルを作成する。
	MGVector v(Pe,Ps);
	m_direction = v/(Te-Ts);
	m_root_point = Ps - Ts * m_direction;
}

// 始点、単位方向ベクトル、終点のパラメータ値を指定して直線を生成する。
MGStraight::MGStraight (
	const MGUnit_vector & v,
	double d,
	const MGPosition& p
):MGCurve(),m_direction(v),m_sparam(0.0),m_endparam(d),
	m_root_point(p),m_knotV(nullptr){
	int dimv=v.sdim(), dimp=p.sdim();
	if(dimv<dimp) m_direction=MGVector(dimp,v);
	else if(dimv>dimp) m_root_point=MGVector(dimv,p);
	if (m_endparam < 0.0) {
		m_direction = -m_direction;
		m_endparam = -m_endparam;
	}
}
	
//Construct the infinite straight line that is a perpendicular bisect
//of the two point P1 and P2 and that is normal to the vector N.
//The line's direction is N*(P2-P1).
//N is the normal of the plane P1, P2, and the constructed line lie on.
MGStraight::MGStraight (
	const MGPosition& P1,	//point 1.
	const MGPosition& P2,	//point 2.
	const MGVector& N
):MGCurve(),m_root_point((P1+P2)*.5),m_direction(N*(P2-P1))
,m_sparam(MGINFINITE_TYPE::MGINFINITE_MINUS),
m_endparam(MGINFINITE_TYPE::MGINFINITE_PLUS),m_knotV(nullptr){
}

//Construct Straight Line copying original line. Able to change
//space dimension and ordering of axis.
MGStraight::MGStraight (
	int dim,				  //New space dimension.
	const MGStraight& line2,  //Original line.
	int start1, 			  //store order of new line.
	int start2)	 		  //source order of original line.
:MGCurve(line2), m_root_point(dim,line2.m_root_point,start1, start2)
	, m_direction(dim, line2.m_direction, start1,start2)
	,m_sparam(line2.m_sparam), m_endparam(line2.m_endparam),m_knotV(0){
	invalidateBox();
	save_length_zero();	
}

///Construct the unlimitted straight that pass through the point uv,
///and the direction is the middle vector of (-v0, v1).
///All of v0, v1, uv, and this straight are objects of space dimension 2.
MGStraight::MGStraight(
	const MGVector& v0,	//a vector whose end is uv.
	const MGVector& v1,	//a vector whose start is uv.
	const MGPosition& uv//origin of sl.
):MGCurve(),m_knotV(nullptr){
	MGVector v2=v0*-1.;
	double angHalf=v1.angle2pai(v2, MGDefault::z_unit_vector());
	angHalf*=.5;
	MGMatrix rotate; rotate.set_rotate_2D(angHalf);
	MGVector vmid=v1*rotate;
	set_straight(MGSTRAIGHT_TYPE::MGSTRAIGHT_UNLIMIT,vmid,uv);
}

///Approximate this curve as a MGLBRep curve
///within the tolerance MGTolerance::line_zero().
///When parameter_normalization=0, reparameterization will not done, and
///the evaluation at the same parameter has the same values before and after
///of approximate_as_LBRep.
void MGStraight::approximate_as_LBRep(
	MGLBRep& lb,		///<Approximated lbrep will be set.
	int ordr,		///<new order. When this is MGLBRep, if ordr=0,
					///ordr=order() will be assumed, else ordr=4 is assumed.
					///For MGStraight, ordr=2 is always employed.
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
	MGPosition Ps=start_point(), Pe=end_point();
	double ts=0., te;
	if(parameter_normalization==0){
		ts=param_s(), te=param_e();
	}else if(parameter_normalization==1){
		te=1.;
	}else{
		te=(Pe-Ps).len();
	}

	ordr=2;//Straight's order is always 2.
	int sd=sdim(), bdim=2;//Straight's B-rep dimension is always 2.

	//Set the knot vector data.
	MGKnotVector& t=lb.knot_vector();
	t.size_change(ordr,bdim);
	t(0)=t(1)=ts;
	t(2)=t(3)=te;

	//Set the B-coefficients data.
	MGBPointSeq& bp=lb.line_bcoef();
	bp.resize(bdim,sd);
	bp.store_at(0,Ps);
	bp.store_at(1,Pe);
}

// 指定線分を囲むボックスを返却する。
MGBox MGStraight::box_limitted(const MGInterval& l) const{
	// 直線範囲のパラメータを作成して入力Intervalとの積をとる。
	MGStraight sl(*this); sl.limit(l);
	return sl.box();
}

//Changing this object's space dimension.
void MGStraight::change_dimension(
	int dim,		// new space dimension
	int start1, 		// Destination order of new object.
	int start2) 		// Source order of this object.
{
	m_root_point=MGPosition(dim,m_root_point,start1, start2);
	m_direction=MGVector(dim, m_direction, start1,start2);
	save_length_zero();	
	invalidateBox();
}

//Change parameter range, be able to change the direction by providing
//s1 greater than s2.
void MGStraight::change_range(
	double s1,		//Parameter value for the start of the original. 
	double s2)		//Parameter value for the end of the original. 
{
	if(m_knotV){delete m_knotV; m_knotV=0;}
	if(s1>s2){
		negate();
		double save=s1; s1=s2; s2=save;
	}
	if(m_sparam.minus_infinite()&&m_endparam.plus_infinite())
		return;

	double t1=m_sparam, t2=m_endparam;
	if(m_sparam.minus_infinite()){
		m_root_point+=(t2-s2)*m_direction;
		m_endparam=s2;
		return;
	}
	if(m_endparam.plus_infinite()){
		m_root_point+=(t1-s1)*m_direction;
		m_sparam=s1;
		return;
	}
	double s2ms1=s2-s1, t2mt1=t2-t1;
	if(MGMZero(s2ms1)){
		if(MGMZero(t2mt1)){
			m_root_point+=(t2-s2)*m_direction;
			m_sparam=s1; m_endparam=s2;
		}
		return;
	}
	double ratio=t2mt1/s2ms1;
	m_root_point += m_direction*(t1-s1*ratio);
	m_direction *= ratio;
	m_sparam=s1; m_endparam=s2;
}

//Compute box of the whole line.
void MGStraight::compute_box(MGBox& bx) const{
	bx.set_null();
	MGInterval p_limit=param_range();
	if(p_limit.empty())
		return;

	// Intervalで示されるboxを作成する。
	if(p_limit.finite()){
		bx=MGBox(eval_position(p_limit.low_point()),
			eval_position(p_limit.high_point()));
		return;
	}

	int dim=sdim();
	bx=MGBox(dim);
	if (p_limit.finite_below()&&p_limit.infinite_above()) {
		for(int i=0; i<dim; i++){
			if(MGRZero2(m_direction.ref(i),m_direction.len()))
				bx(i)=MGInterval(start_point().ref(i),start_point().ref(i));
			else if(m_direction.ref(i)>0.)
				bx(i)=
				MGInterval(MGINTERVAL_TYPE::MGINTERVAL_FINITE_BELOW,start_point().ref(i));
			else
				bx(i)=
				MGInterval(MGINTERVAL_TYPE::MGINTERVAL_FINITE_ABOVE,start_point().ref(i));
		}
	} else if (p_limit.finite_above()&&p_limit.infinite_below()) {
		for(int i=0; i<dim; i++){
			if(MGRZero2(m_direction.ref(i),m_direction.len()))
				bx(i)=MGInterval(root_point().ref(i),root_point().ref(i));
			else if(m_direction.ref(i)>0.)
				bx(i)=
				MGInterval(MGINTERVAL_TYPE::MGINTERVAL_FINITE_ABOVE,root_point().ref(i));
			else
				bx(i)=
				MGInterval(MGINTERVAL_TYPE::MGINTERVAL_FINITE_BELOW,root_point().ref(i));
		}
	} else {
		for(int i=0; i<dim; i++){
			if(MGRZero2(m_direction.ref(i),m_direction.len())){
				MGInterval Ii=MGInterval(root_point().ref(i),root_point().ref(i));
				bx(i)=Ii;
			} else{
				MGInterval Ii=MGInterval(MGINTERVAL_TYPE::MGINTERVAL_INFINITE);
				bx(i)=Ii;
			}
		}
	} 
}

//Exchange ordering of the coordinates.
//Exchange coordinates (i) and (j).
void MGStraight::coordinate_exchange(int i, int j){
	assert(i<sdim() && j<sdim());
	m_direction.swap(i,j); m_root_point.swap(i,j);
	invalidateBox();
}

//Construct new curve object by copying to newed area.
//User must delete this copied object by "delete".
MGStraight* MGStraight::clone() const{return new MGStraight(*this);}

///Convert this curve to Bezier curve.
///If this is MGLBRep or MGStraight, the shape is exactly the same
///as the original. Otherwise, this is apporoximated by MGLBRep.
void MGStraight::convert_to_Bezier(MGLBRep& bezier)const{
	bezier.copy_appearance(*this);
	MGKnotVector& t=bezier.knot_vector();
	t.size_change(4,4);
	t(0)=t(1)=t(2)=t(3)=0.;
	t(4)=t(5)=t(6)=t(7)=1.;

	MGPosition P0=start_point();
	MGPosition P3=end_point();
	MGBPointSeq& bp=bezier.line_bcoef();
	MGVector P0toP3(P3,P0);
	P0toP3/=3.;
	MGPosition P1=P0+P0toP3;
	MGPosition P2=P1+P0toP3;
	bp.resize(4,P0.sdim());
	bp.store_at(0,P0);
	bp.store_at(1,P1);
	bp.store_at(2,P2);
	bp.store_at(3,P3);
}

//copy as a newed curve. The new curve will be MGLBRep or MGRLBRep.
//When original curve was a MGRLBRep, the new curve will be a MGRLBRep.
//Otherwise,  the new curve will be a MGLBRep.
//Returned object must be deleted.
MGCurve* MGStraight::copy_as_nurbs() const{
	return new MGLBRep(*this);
}

//Construct new curve object by changing
//the original object's space dimension.
//User must delete this copied object by "delete".
MGStraight* MGStraight::copy_change_dimension(
	int sdim,		// new space dimension
	int start1, 		// Destination order of new line.
	int start2 		// Source order of this line.
)const{
	return new MGStraight(sdim,*this,start1,start2);
}

//Compute curvilinear integral of the 1st two coordinates.
//This integral can be used to compute area sorounded by the curve.
//線積分を求める。
double MGStraight::curvilinear_integral(double t1, double t2) const{
	t1=range(t1); t2=range(t2);
	return (t2-t1)*(root_point().ref(0)*m_direction.ref(1)
				  - root_point().ref(1)*m_direction.ref(0));
}

// 終点の座標値を返却する。
MGPosition MGStraight::end_point() const{
	return eval_position(m_endparam.value());
}

// Evaluate n'th derivative data. n=0 means positional data evaluation.
MGVector MGStraight::eval(
	double t,		// Parameter value.
	int nderiv,	// Order of Derivative.
	int left		//Left continuous(left=true)
					//or right continuous(left=false).
)const{
	if(nderiv==0)      return eval_position(t);
	else if(nderiv==1) return m_direction;
	else               return MGVector(sdim(),MGVector (0.0,0.0));
}

// パラメータ値を与えて位置、一次微分値、二次微分値を求める。
void MGStraight::eval_all (double d,
	MGPosition & p,				   // 位置      
	MGVector & v1,				   // 一次微分値
	MGVector & v2				   // 二次微分値
)const{
	p = eval_position (d);			// 位置
	v1 = m_direction;				// 一次微分値
	v2 = MGVector(sdim(),0.0);		// 二次微分値
}

// 直線上の与えられたパラメータ値における一次微分値を求める。
MGVector MGStraight::eval_deriv(double d) const{
	return m_direction;
}

/// 与えられたパラメータ値に相当する直線上の点を返却する。
///Evaluate positional data at parameter t.
///Evaluation is done after t is limitted within the parameter range.
MGPosition MGStraight::eval_position(double d)const{
	return m_root_point+range(d)*m_direction;
}


/// 与えられたパラメータ値に相当する直線上の点を返却する。
///Evaluate positional data at parameter t.
///Evaluation is done without limitting the parameter value t within the parameter range.
MGPosition MGStraight::eval_position_unlimitted(double d)const{
	return m_root_point+d*m_direction;
}

//Extrapolate this curve by an (approximate) chord length.
//The extrapolation is C2 continuous.
void MGStraight::extend(
	double length,	//approximate chord length to extend. 
	bool start		//Flag of which point to extend, start or end point of the line.
					//If start is true extend on the start point.
){
	assert(length>0.);

	const MGVector& dir=direction();
	double tlen=length/dir.len();
	if(start){
		if(infinite_below())
			return;
		m_sparam-=tlen;
	}else{
		if(infinite_above())
			return;
		m_endparam+=tlen;
	}
	if(m_knotV){
		delete m_knotV;
		m_knotV=0;
	}
	invalidateBox();
}

//Test if input parameter value is inside parameter range of the line.
bool MGStraight::in_range(double t)const{
	double error=MGTolerance::rc_zero()*m_direction.len();
				//error is maximum distance of 
				// two points that are assumed to be the same.
	mgTolSetWCZero wczeroSet(error);//Set&save the error.
	return (m_sparam<=t && t<=m_endparam);
}

//Test if this cure is co-planar with the 2nd curve curve2.
//MGPlane expression will be out to plane if this is co-planar.
//Function's return value is true if co-planar.
bool MGStraight::is_coplanar(const MGCurve& curve2, MGPlane& plane)const{
	const MGStraight* sl2=dynamic_cast<const MGStraight*>(&curve2);
	if(sl2){
		//When both are straight.
		if(direction().parallel(sl2->direction())){
			plane=MGPlane(*sl2,root_point());
			return true;
		}else{
			MGUnit_vector normal=direction()*sl2->direction();
			plane=MGPlane(normal,root_point());
			return plane.on(sl2->root_point());
		}
	}

	MGStraight line2;
	MGPosition point;
	int plkind=-1;
	const MGLBRep* lb=dynamic_cast<const MGLBRep*>(&curve2);
	if(lb)
		plkind=lb->planar(plane,line2,point);
	else{
		const MGRLBRep* rlb=dynamic_cast<const MGRLBRep*>(&curve2);
		if(rlb)
			plkind=lb->planar(plane,line2,point);
	}
	if(plkind==0) return false;
	if(plkind==1){
		plane=MGPlane(*this,point);
		return true;
	}else if(plkind==2){
		if(!direction().parallel(line2.direction())) return false;
		plane=MGPlane(line2,root_point());
		return true;
	}else if(plkind==3){
		return plane.on(*this);
	}

	//When curve2 is neither MGLBRep nor MGRLBRep.
	if(!curve2.is_planar(plane)) return false;
	return plane.on(*this);
}

//Test if the input parameter t is the start point parameter or not.
bool MGStraight::is_startpoint_parameter(double t)const{
	if(m_sparam.minus_infinite())
		return false;

	if(m_endparam.plus_infinite()){
		return m_sparam==t;
	}

	return MGREqual_base(param_s(),t,param_span());
}

//Test if the input parameter t is the start point parameter or not.
bool MGStraight::is_endpoint_parameter(double t)const{
	if(m_endparam.plus_infinite())
		return false;

	if(m_sparam.minus_infinite()){
		return m_endparam==t;
	}

	return MGREqual_base(param_s(),t,param_span());
}

//Test if this cure is linear or not, that is, is straight or not.
//MGStraight expression will be out to straight if this is linear or not.
//Function's return value is true if linear.
bool MGStraight::is_linear(MGStraight& straight)const{
	straight=(*this);
	return true;
}

//Test if this cure is planar or not.
//MGPlane expression will be out to plane if this is planar.
//Function's return value is true if planar.
bool MGStraight::is_planar(MGPlane& plane)const{
	const MGVector& sldir=direction();
	double xx=sldir%mgX_UVEC; xx*=xx;
	double yy=sldir%mgY_UVEC; yy*=yy;
	double zz=sldir%mgZ_UVEC; zz*=zz;
	MGVector dir2;
	if(xx<=yy){
		if(xx<=zz)
			dir2=mgX_UVEC;
		else
			dir2=mgZ_UVEC;
	}else{
		if(yy<=zz)
			dir2=mgY_UVEC;
		else
			dir2=mgZ_UVEC;
	}
	plane=MGPlane(sldir,dir2,root_point());
	return true;
}


//Access to i-th element of knot.
//i=0, 1 and returns start or end parameter value of the straight.
double MGStraight::knot(int i) const{
	assert(i<=3);
	if(i<=1)
		return param_s();
	else
		return param_e();
}

//Returns the knot vector of the curve.
const MGKnotVector& MGStraight::knot_vector() const{
	if(!m_knotV){
		m_knotV=new MGKnotVector(2,2);
		(*m_knotV)(0)=(*m_knotV)(1)=param_s();
		(*m_knotV)(2)=(*m_knotV)(3)=param_e();
	}
	return *m_knotV;
}
MGKnotVector& MGStraight::knot_vector(){
	if(!m_knotV){
		m_knotV=new MGKnotVector(2,2);
		(*m_knotV)(0)=(*m_knotV)(1)=param_s();
		(*m_knotV)(2)=(*m_knotV)(3)=param_e();
	}
	return *m_knotV;
}
// パラメータ値が昇順であたえられた場合正の値で、降順の場合
// 負の値でパラメータ間の直線の代数的距離を返却する。
double MGStraight::length (double d1, double d2) const {
	return (range(d2) - range(d1))*m_direction.len();
}

// 自身の直線が有界の場合、その直線の距離を返却する。
// 非有界のとき－２を返却する。
double MGStraight::length () const {
	double d = -2.0;
	if (straight_type() == MGSTRAIGHT_TYPE::MGSTRAIGHT_SEGMENT)
		d=(m_endparam.value()-m_sparam.value())*m_direction.len();
	else if (straight_type() == MGSTRAIGHT_TYPE::MGSTRAIGHT_EMPTY) d = 0.0;
	return d;
}

// 直線上の与えられたパラメータで示される点から指定距離はなれた点
// のパラメータ値をかえす。
double MGStraight::length_param (double t,double len) const {
	return range(t+len/m_direction.len());
}

// 自身の直線に指定されたｌｉｍｉｔを付与する。
void MGStraight::limit(const MGInterval& l) {
	if(m_knotV){delete m_knotV; m_knotV=0;}
	MGInterval p_limit = param_range(); // 直線範囲のパラメータを作成
	p_limit &= l;						// 指定Intervalとの積をとる
	m_sparam=p_limit.low();
	m_endparam=p_limit.high();
	invalidateBox();
}

//Compute sub straight line limitted by an box.
//This box's coordinates consist of world coordinates.
void MGStraight::limit(const MGBox& box){
	MGInterval prange=param_range();
	MGCParam_list list;
	double bound, t;
	int dim=sdim();
	for(int i=0; i<dim; i++){
		bound=box.ref(i).low_point();
		list=isect_1D(bound, i);
		if(list.entries()){
			t=list.first();
			if(prange.includes(t)){
				if(m_direction.ref(i)>=0.)
					prange.set_low_point(t);
				else
					prange.set_high_point(t);
			}
		}
		bound=box.ref(i).high_point();
		list=isect_1D(bound, i);
		if(list.entries()){
			t=list.first();
			if(prange.includes(t)){
				if(m_direction.ref(i)>=0.)
					prange.set_high_point(t);
				else
					prange.set_low_point(t);
			}
		}
	}
	limit(prange);
}

//Compute nearest point on the line to the origin.
MGPosition MGStraight::nearest_to_origin() const{	
	double lensqr=m_direction.len(); lensqr*=lensqr;
	return MGPosition
		(root_point()-m_direction*((root_point()%m_direction)/lensqr));
}

// 直線の方向を反転する(方向ベクトルを逆向きにする）。
// 始終点があるときは入れ換える。
void MGStraight::negate(){
	if(m_knotV){delete m_knotV; m_knotV=0;}
	MGEReal save=m_sparam;
	m_sparam=-m_endparam; m_endparam=-save;
	m_direction = -m_direction;
}

//Obtain parameter value if this curve is negated by "negate()".
double MGStraight::negate_param(double t)const{
	return -t;
}

//Obtain so transformed 1D curve expression of this curve that
//f(t)={sum(xi(t)*g[i]) for i=0(x), 1(y), 2(z)}-g[3], where f(t) is the output
//of oneD and xi(t) is i-th coordinate expression of this curve.
//This is used to compute intersections with a plane g[4].
std::unique_ptr<MGCurve> MGStraight::oneD(
	const double g[4]			//Plane expression(a,b,c,d) where ax+by+cz=d.
)const{
	MGStraight* sl=new MGStraight(*this);
	MGVector G(3,g);
	MGPosition one(1);
	one(0)=G%m_direction; sl->m_direction=one;
	one(0)=G%m_root_point-g[3]; sl->m_root_point=one;
	return std::unique_ptr<MGCurve>(sl);
}

//Obtain parameter space error.
double MGStraight::param_error()const{
	MGEReal len=m_endparam-m_sparam;
	if(len.finite())
		return len.value()*MGTolerance::rc_zero();
	else{
		double vlen=m_direction.len();
		return vlen*MGTolerance::rc_zero();
	}
}

//Normalize parameter value t to the nearest knot if their distance is
//within tolerance. For straight, the knots are start and end points.
double MGStraight::param_normalize(double t) const{
	double plen;
	if(straight_type()== MGSTRAIGHT_TYPE::MGSTRAIGHT_SEGMENT) plen=param_span();
	else plen=1.;
	double tnew;
	if(straight_type()== MGSTRAIGHT_TYPE::MGSTRAIGHT_SEGMENT
		|| straight_type()== MGSTRAIGHT_TYPE::MGSTRAIGHT_HALF_LIMIT){
		tnew=param_s();
		if(MGRZero2(t-tnew,plen)) return tnew;
	}
	if(straight_type()== MGSTRAIGHT_TYPE::MGSTRAIGHT_SEGMENT){
		tnew=param_e();
		if(MGRZero2(t-tnew,plen)) return tnew;
	}
	return t;
}

//Return parameter range of the curve(パラメータ範囲を返す)
MGInterval MGStraight::param_range()const{
	return MGInterval(m_sparam, m_endparam);
}

//Compute part of this curve from parameter t1 to t2.
//Returned is the pointer to newed object, and so should be deleted
//by calling program, or memory leaked.
MGStraight* MGStraight::part(double t1, double t2, int multiple)const{
	assert(t2>=t1);
	MGStraight* sl=new MGStraight(*this);
	sl->limit(MGInterval(t1,t2));
	return sl;
}

// 入力パラメータをパラメータ範囲で丸めて返却する。
double MGStraight::range(double d) const {
	MGInterval limit = param_range();
	if(straight_type()== MGSTRAIGHT_TYPE::MGSTRAIGHT_EMPTY)
		d=0.0;
	else if(limit>=d)
		d=m_sparam;
	else if(limit<=d)
		d=m_endparam;
	return d;
}

//Function to avoid m_direction.len()=zero.
//The line's m_direction is set as a unit vector and m_endparam
//is set to zero.
void MGStraight::save_length_zero(){
	double dlen=m_direction.len();
	if(MGMZero(dlen)){
		m_direction=m_direction.normalize();
		m_endparam=m_sparam+dlen;
		if(m_knotV){delete m_knotV; m_knotV=0;}
	}
}

//Return space dimension
int MGStraight::sdim() const{	
	int dim1=m_root_point.sdim();
	if(dim1==0) return 0;
	int dim2=m_direction.sdim();
	if(dim1<dim2) dim1=dim2;
	return dim1;
}

// 直線のタイプ，方向ベクトル，始点を指定して直線を生成する。
//Straight from straight line type, direction vector, and an origin.
//Construct a straight and replce this with it.
//This fuction does not convert input vec to a unit vector.
//If you like the conversion, use MGStraight() constructor.
MGStraight& MGStraight::set_straight(
	MGSTRAIGHT_TYPE type,				//Type
	const MGVector& vec,				//Direction
	const MGPosition& Q					//Origin
){
	if(m_knotV){delete m_knotV; m_knotV=0;}
	m_root_point=Q;
	m_direction=vec;
	m_sparam=0.;
	m_endparam=MGEReal(MGINFINITE_TYPE::MGINFINITE_PLUS);

	if (type== MGSTRAIGHT_TYPE::MGSTRAIGHT_SEGMENT)
		m_endparam=1.;
	else if(type== MGSTRAIGHT_TYPE::MGSTRAIGHT_EMPTY)
		m_endparam=-1.;
	else if(type== MGSTRAIGHT_TYPE::MGSTRAIGHT_UNLIMIT)
		m_sparam=MGEReal(MGINFINITE_TYPE::MGINFINITE_MINUS);
	int dimv=vec.sdim(), dimp=Q.sdim();
	if(dimv<dimp) m_direction=MGVector(dimp,vec);
	else if(dimv>dimp) m_root_point=MGPosition(dimv,Q);
	save_length_zero();
	invalidateBox();
	return *this;
}

// 自身の直線の始点を返却する。
//Return start(root) point of the straight.
MGPosition MGStraight::start_point() const{
	return eval_position(m_sparam.value());
}

// 直線のタイプを返却する。
//Return the straight line type.
MGSTRAIGHT_TYPE MGStraight::straight_type() const{
	if(m_sparam>m_endparam) return MGSTRAIGHT_TYPE::MGSTRAIGHT_EMPTY;
	else if(m_sparam.finite()){
		if(m_endparam.finite()) return MGSTRAIGHT_TYPE::MGSTRAIGHT_SEGMENT;
		else return MGSTRAIGHT_TYPE::MGSTRAIGHT_HALF_LIMIT;
	}else{
		if(m_endparam.finite()) return MGSTRAIGHT_TYPE::MGSTRAIGHT_HALF_LIMIT;
		else return MGSTRAIGHT_TYPE::MGSTRAIGHT_UNLIMIT;
	}
}

// 自身の直線からｌｉｍｉｔを取り除く。
MGCurve& MGStraight::unlimit(){
	m_sparam=MGEReal(MGINFINITE_TYPE::MGINFINITE_MINUS);
	m_endparam=MGEReal(MGINFINITE_TYPE::MGINFINITE_PLUS);
	invalidateBox();
	return *this;
}

//Unlimit parameter range of the curve to the end point direction
//(終点方向にlimitをはずす)
MGCurve& MGStraight::unlimit_end(){
	m_endparam=MGEReal(MGINFINITE_TYPE::MGINFINITE_PLUS);
	invalidateBox();
	return *this;
}

//Unlimit parameter range of the curve to the start point direction
//(始点方向にlimitをはずす)
MGCurve& MGStraight::unlimit_start(){
	m_sparam=MGEReal(MGINFINITE_TYPE::MGINFINITE_MINUS);
	invalidateBox();
	return *this;
}

//Update the root point of this straight.
void MGStraight::update_root(const MGPosition& rootP){
	m_root_point=rootP;
	invalidateBox();
}

//
// 演算子の多重定義
//

//Copy Assignment.
MGStraight& MGStraight::operator=(const MGStraight& sl2){
	if(this==&sl2)
		return *this;

	MGCurve::operator=(sl2);
	m_direction=sl2.m_direction;
	m_endparam=sl2.m_endparam;
	if(m_knotV){
		delete m_knotV; m_knotV=nullptr;
	}
	m_root_point=sl2.m_root_point;
	m_sparam=sl2.m_sparam;
	return *this;
}

//Move assignment.
MGStraight& MGStraight::operator=(MGStraight&& sl2){
	MGCurve::operator=(std::move(sl2));
	m_direction=std::move(sl2.m_direction);
	m_endparam=sl2.m_endparam;
	if(m_knotV){
		delete m_knotV; m_knotV=sl2.m_knotV;
		sl2.m_knotV=nullptr;
	}
	m_root_point=std::move(sl2.m_root_point);
	m_sparam=sl2.m_sparam;
	return *this;
}

//When the leaf object of this and crv2 are not equal, this assignment
//does nothing.
MGStraight& MGStraight::operator=(const MGGel& gel2){
	const MGStraight* gel2_is_this=dynamic_cast<const MGStraight*>(&gel2);
	if(gel2_is_this)
		operator=(*gel2_is_this);
	return *this;
}

//When the leaf object of this and crv2 are not equal, this assignment
//does nothing.
MGStraight& MGStraight::operator=(MGGel&& gel2){
	MGStraight* gel2_is_this=dynamic_cast<MGStraight*>(&gel2);
	if(gel2_is_this)
		operator=(std::move(*gel2_is_this));
	return *this;
}

// 直線の平行移動を行いオブジェクトを生成する。
MGStraight MGStraight::operator+(const MGVector& vec) const{
	MGStraight sl2 = *this;
	sl2 += vec;
	return sl2;
}
MGStraight operator+(const MGVector& v, const MGStraight& sl){
	return sl+v;
}

// 直線の平行移動を行い自身の直線とする。
MGStraight& MGStraight::operator+=(const MGVector & v){
	m_root_point += v;
	m_box+=v;
	return *this;
}

// 直線の逆方向に平行移動を行いオブジェクトを生成する。
MGStraight MGStraight::operator-(const MGVector & v)const{
	MGStraight sl2 = *this;
	sl2 -= v;
	return sl2;
}

//直線の逆方向に平行移動を行い自身の直線とする。
MGStraight& MGStraight::operator-= (const MGVector & v){
	m_root_point -= v;
	m_box-=v;
	return *this;
}

// 与えられたスケールで直線の変換を行いオブジェクトを生成する。
//Generate scaled straight.
MGStraight MGStraight::operator*(double scale) const{
	MGStraight sl(*this);
	sl *=scale;
	return sl;
}

// 与えられたスケールで直線の変換を行いオブジェクトを生成する。
//Generate scaled straight.
MGStraight operator*(double scale, const MGStraight& sl){
	return sl*scale;
}

// 与えられたスケールで直線の変換を行い自身の直線とする。
//Scaling the straight.
MGStraight& MGStraight::operator*=(double scale){
	m_root_point *=scale;
	m_direction *=scale;
	save_length_zero();
	invalidateBox();
	return *this;
}

// 与えられた変換で直線の変換を行いオブジェクトを生成する。
MGStraight MGStraight::operator*(const MGMatrix& t)const{
	MGStraight sl2 = *this;
	sl2 *= t;
	return sl2;
}

// 与えられた変換で直線の変換を行い自身の直線とする。
MGStraight& MGStraight::operator*=( const MGMatrix& t ){
	m_root_point *= t;
	m_direction *=t;
	save_length_zero();
	invalidateBox();
	return *this;
}

// 与えられた変換で直線のトランスフォームを行いオブジェクトを生成する。
MGStraight MGStraight::operator*(const MGTransf& t)const{
	MGStraight sl2 = *this;
	sl2 *= t;
	return sl2;
}

// 与えられた変換で直線のトランスフォームを行い自身の直線とする。
MGStraight& MGStraight::operator*=(const MGTransf& t){
	m_root_point *= t;
	m_direction *=t.affine();
	save_length_zero();
	invalidateBox();
	return *this;
}

//
// 論理演算子の多重定義
bool MGStraight::operator==(const MGStraight& sl2)const{
	if(sdim()==0 && sl2.sdim()==0)
		return 1;
	 return (m_root_point == sl2.m_root_point 
		&& m_direction == sl2.m_direction
		&& m_sparam==sl2.m_sparam
		&& m_endparam==sl2.m_endparam);
}

bool MGStraight::operator<(const MGStraight& gel2)const{
	return m_root_point.len()<gel2.m_root_point.len();
}
bool MGStraight::operator==(const MGGel& gel2)const{
	const MGStraight* gel2_is_this=dynamic_cast<const MGStraight*>(&gel2);
	if(gel2_is_this)
		return operator==(*gel2_is_this);
	return false;
}
bool MGStraight::operator<(const MGGel& gel2)const{
	const MGStraight* gel2_is_this=dynamic_cast<const MGStraight*>(&gel2);
	if(gel2_is_this)
		return operator<(*gel2_is_this);
	return identify_type() < gel2.identify_type();
}
