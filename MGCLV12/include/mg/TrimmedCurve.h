/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#ifndef _MGTrimmedCurve_HH_
#define _MGTrimmedCurve_HH_

#include <vector>
#include "mg/Default.h"
#include "mg/Interval.h"
#include "mg/Box.h"
#include "mg/Position.h"
#include "mg/Position_list.h"
#include "mg/LBRep.h"
#include "mg/RLBRep.h"
#include "mg/Ellipse.h"
#include "mg/Straight.h"
#include "mg/BSumCurve.h"

//Define MGTrimmedCurve Class.

class MGInterval;
class MGBox;
class MGVector;
class MGUnit_vector;
class MGPosition;
class MGTransf;
class MGCParam_list;
class MGStraight;
class MGEllipse;
class MGCompositeCurve;
class MGKnotVector;
class MGLBRep;
class MGRLBRep;
class MGSurface;
class MGCCisects;
class MGCSisects;
class MGIfstream;
class MGOfstream;
/** @file */
/** @addtogroup GEO
 *  @{
 */

///MGTrimmedCurve is a part of an original curve of a limitted parameter range.

///MGTrimmedCurve is a part of original curve that has limitted parameter range.
///MGTrimmedCurve is a temporal curve, and does not have update functions.
class MG_DLL_DECLR MGTrimmedCurve:public MGCurve{

public:

MG_DLL_DECLR friend MGTrimmedCurve operator+ (const MGVector& v, const MGTrimmedCurve& cv2);
MG_DLL_DECLR friend MGTrimmedCurve operator* (double scale, const MGTrimmedCurve& cv2);


////////Special member functions/////////
MGTrimmedCurve();
~MGTrimmedCurve()=default;
MGTrimmedCurve(const MGTrimmedCurve&);///<Copy constructor.
MGTrimmedCurve& operator= (const MGTrimmedCurve&);///<Copy assignment.
MGTrimmedCurve(MGTrimmedCurve&&)=default;		///<Move constructor.
MGTrimmedCurve& operator= (MGTrimmedCurve&&)=default;///<Move assignment.
MGTrimmedCurve(const MGCurve& crv, const MGInterval range);

///subcurve of the input curve crv from t1 to t2.
///t1 may be >t2. In this case, the curve will be from t2 to t1.
MGTrimmedCurve(const MGCurve& crv, double t1, double t2);

///Assignment.
///When the leaf object of this and crv2 are not equal, this assignment
///does nothing.
MGTrimmedCurve& operator=(const MGGel& gel2);
MGTrimmedCurve& operator=(MGGel&& gel2);

///Comparison of two curves.
bool is_same_curve(const MGCurve& curve2)const;
bool operator==(const MGCompositeCurve& gel2)const;

bool operator==(const MGTrimmedCurve& gel2)const;
std::partial_ordering operator<=>(const MGTrimmedCurve& gel2)const;

//gel2 must be the same class as this.
bool equal_test(const MGGel& gel2)const override;

//gel2 must be the same class as this.
std::partial_ordering ordering_test(const MGGel& gel2)const override;

///Approximate this curve as a MGLBRep curve
///within the tolerance MGTolerance::line_zero().
///When parameter_normalization=0, reparameterization will not done, and
///the evaluation at the same parameter has the same values before and after
///of approximate_as_LBRep.
void approximate_as_LBRep(
	MGLBRep& lb,///<Approximated obrep will be set.
	int ordr=0,	///<new order. When this is MGLBRep, if ordr=0,
				///<ordr=order() will be assumed, else ordr=4 is assumed.
	int parameter_normalization=0,
		///<Indicates how the parameter normalization be done:
		///<  =0: no parameter normalization.
		///<  =1: normalize to range=(0., 1.);
		///<  =2: normalize to make the average length of the 1st derivative 
		///<      is as equal to 1. as possible.
	bool neglectMulti=false///<Indicates if multiple knots be kept.
		///< true: multiplicity is removed.
		///< false: multiplicity is kept.
)const override;

///Returns B-Rep Dimension.
int bdim() const;

///Return minimum box that includes the curve of parameter interval.
/// 入力のパラメータ範囲の曲線部分を囲むボックスを返す。
MGBox box_limitted(
	const MGInterval& ///< Parameter Range of the curve.
) const;

///Changing this object's space dimension.
void change_dimension(
	int dim,		///< new space dimension
	int start1, 	///< Destination order of new object.
	int start2 	///< Source order of this object.
){
	assert(false);//This cannot be used.
}
///Change parameter range, be able to change the direction by providing
///t1 greater than t2.
void change_range(
	double t1,	///<Parameter value for the start of original. 
	double t2	///<Parameter value for the end of original.
){
	assert(false);//The use is prohibitted.
}

///Construct new curve object by copying to newed area.
///User must delete this copied object by "delete".
MGCurve* clone() const;

///与えられた曲線と自身の共通部分があるかどうか調べる。

///引数：
///		const MGCurve&			curve2,		(I/ )	与えられる曲線
///		std::vector<double>&	vecComSpan	( /O)	共通部分のパラメータ範囲
///		 4nの配列で、vecComSpan(4*i+0),vecComSpan(4*i+1)が自身のパラメータ範囲
///					(vecComSpan(4*i+0) < vecComSpan(4*i+1))、
///				 vecComSpan(4*i+2),vecComSpan(4*i+3)がcurve2のパラメータ範囲
///		MGCCisects&			isect		( /O)	交点
///戻り値：
///		3:交点も共通部分も求まった
///		2:交点のみが求まった
///		1:共通部分のみが求まった
///		0:交点も共通部分もなかった
///		-1:共通エッジの収束計算エラー
///		-2:共通エッジが４個以上求まった(のっていないと見なす)
///追記：
///	曲線が共通かどうかの誤差にはline_zero()、をパラメータ範囲の収束計算の
///	誤差には、パラメータ範囲*rc_zero()を使用した
int common(
	const MGCurve& curve2,
	std::vector<double>& vecComSpan,
	MGCCisects& isect
) const;

///与えられた曲線と自身の共通部分があるかどうか調べる。

///引数：
///		const MGCurve&			curve2,		(I/ )	与えられる曲線
///		std::vector<double>&	vecComSpan	( /O)	共通部分のパラメータ範囲
///		 4nの配列で、vecComSpan(4*i+0),vecComSpan(4*i+1)が自身のパラメータ範囲
///					(vecComSpan(4*i+0) < vecComSpan(4*i+1))、
///				 vecComSpan(4*i+2),vecComSpan(4*i+3)がcurve2のパラメータ範囲
///戻り値：
///		共通部分の数:	共通部分が求まった
///		0:				共通部分がなかった
///		-1:				共通エッジの収束計算エラー
///		-2:				共通エッジが４個以上求まった(のっていないと見なす)
///追記：
///	曲線が共通かどうかの誤差にはline_zero()を、パラメータ範囲の収束計算の誤差には、
///  パラメータ範囲*rc_zero()を使用した
int common(
	const MGCurve& curve2,
	std::vector<double>& vecComSpan
) const;

///Compute the box from the scratch.
void compute_box(MGBox& bx) const;

///copy as a newed curve.

///The new curve will be MGLBRep or MGRLBRep.
///When original curve was a MGRLBRep, the new curve will be a MGRLBRep.
///Otherwise,  the new curve will be a MGLBRep.
///Returned object must be deleted.
MGCurve* copy_as_nurbs() const override;

///Construct new curve object by changing
///the original object's space dimension.
///User must delete this copied object by "delete".
MGCurve* copy_change_dimension(
	int sdim,			///< new space dimension
	int start1=0, 		///< Destination order of new line.
	int start2=0 		///< Source order of this line.
)const;

///Construct new curve object by copying to newed area,
///and limitting the parameter range to prange.
///Returned is newed object and must be deleted.
///Returned curve is not TrimmedCurve, but a part copy of the original
///curve.
MGCurve* copy_limitted(const MGInterval& prange) const;
MGCurve* copy_limitted() const;

/// Return parameter curve's pointer
const MGCurve* base_curve() const {return m_curve;};

///Compute curvilinear integral of the 1st two coordinates.

///This integral can be used to compute area sorounded by the curve.
///(線積分）を求める。
///curvilinear_integral from t1 to t2 can be obtained by
///Integral of (x*dy-y*dx) about t, where curve is expressed by
///f(t)=(x(t),y(t)), dx=dx/dt, and dy=dy/dt.
double curvilinear_integral(double t1, double t2) const;

///Divide this curve at the designated knot multiplicity point.
///Function's return value is the number of the curves after divided.
int divide_multi(
	std::vector<UniqueCurve>& crv_list,///<divided curves are appended.
	int multiplicity=-1	///<designates the multiplicity of the knot to divide at,
						///<When multiplicity<=0, order()-1 is assumed,
						///<When multiplicity>=order(), order() is assumed.
)const override;

///Draw using vbo.
void drawSE(
	mgVBO& vbo,///<The target graphic object.
	double t0,			///<Start parameter value of the curve.
	double t1			///<End parameter value of the curve,
						///<Draw will be performed from t0 to t1.
)const;

/// Evaluate n'th derivative data. n=0 means positional data evaluation.
MGVector eval(
	double t,			///< Parameter value.
	int nderiv=0,	///< Order of Derivative.
	int left=0			///<Left continuous(left=true)
						///<or right continuous(left=false).
)const{return m_curve->eval(t,nderiv,left);};

///Extrapolate this curve by an (approximate) chord length.
///The extrapolation is C2 continuous.
void extend(
	double length,	///<approximate chord length to extend. 
	bool start=false///<Flag of which point to extend, start or end point of the line,
					///<If start is true extend on the start point.
);

/// Return This object's typeID
long identify_type() const;

///Test if input parameter value is inside parameter range of the line.
bool in_range(double t) const{return m_range.includes(t);}

///Provide divide number of curve span for function intersect.
int intersect_dnum() const;

///Test if this cure is planar or not.
///MGPlane expression will be out to plane if this is planar.
///Function's return value is true if planar.
bool is_planar(MGPlane& plane)const;

///<Test if this parameter range is the same as the original m_curve's.
bool is_same_range()const{return m_sameRange;};

///Intersection of Curve and Surface
MGCCisects isect(const MGCurve& curve2) const;
MGCCisects isect(const MGStraight& curve2) const;
MGCCisects isect(const MGLBRep& curve2) const;
MGCCisects isect(const MGSurfCurve& curve2) const;

///Intersection with a Surface
MGCSisects isect(const MGSurface& surf) const override;
MGCSisects isect(const MGPlane& surf) const override;

///Access to i-th element of knot
double knot(int i) const;
	
///Returns the knot vector.
const MGKnotVector& knot_vector() const;

///Cmpute curve length of the interval.

///If t1 is greater than t2, return negative value.
/// 与えられたパラメータ値間の曲線の長さを返す。
/// パラメータが昇順で与えられたときは正値、降順のときは負値を返す。
double length(double t1, double t2) const;

///Inverse function of length. Compute the point that is away from
///the point t by length len.
/// lengthの逆関数。指定パラメータtで示される点から指定距離len
/// 曲線上に沿って離れた点を示すパラメータ値を返す。
double length_param( double t, double len) const;

///Update this by limiting the parameter range of the curve.
/// 自身に指定したパラメータ範囲のｌｉｍｉｔをつける。
void limit(const MGInterval& );

///Obtain parameter value if this curve is negated by "negate()".
double negate_param(double t)const{return m_curve->negate_param(t);};

///Returns the order.
int order() const{return m_curve->order();}

/// Output function.
std::ostream& toString(std::ostream&) const;

///IGES output function.
int out_to_IGES(
	MGIgesOfstream& igesfile,
	int SubordinateEntitySwitch=0
)const;

/// Return ending parameter value.
double param_e() const;

///Normalize parameter value t to the nearest knot if their distance is
///within tolerance.
double param_normalize(double t) const;

///Return parameter range of the curve(パラメータ範囲を返す)
MGInterval param_range() const{return m_range;}

/// Return starting parameter value.
double param_s() const;
	
///Compute part of this curve from parameter t1 to t2.

///Returned is the pointer to newed object, and so should be deleted
///by calling program, or memory leaked.
MGCurve* part(
	double t1,///<target parameter range from t1,
	double t2,///< to t2.
	int multiple=0	///<Indicates if start and end knot multiplicities
					///<are necessary. =0:unnecessary, !=0:necessary.
) const;
	
///Compute all foot points of the perpendicular line from point to the curve.

/// 与ポイントから曲線へ下ろした垂線の足の，曲線のパラメータ値を
/// すべて求める。
MGCParam_list perps(
	const MGPosition& P		///<Point(指定点)
)const;

///Compute all the perpendicular points of this curve and the second one.

///If f(s) and g(t) are the points of the two curves f and g,
///then obtains points where the following conditions are satisfied:
///  fs*(f-g)=0.    gt*(g-f)=0.
///Here fs and gt are 1st derivatives at s and t of f and g.
///MGPosition P in the MGPosition_list contains this and crv's parameter
///as:     P(0)=this curve's parameter, P(1)=crv2's parameter value.
MGPosition_list perps(const MGCurve& crv2)const{return perps_in_range(crv2);};
MGPosition_list perps(const MGStraight& crv2)const{return perps_in_range(crv2);};
MGPosition_list perps(const MGRLBRep& crv2)const{return perps_in_range(crv2);};
MGPosition_list perps(const MGEllipse& crv2)const{return perps_in_range(crv2);};
MGPosition_list perps(const MGLBRep& crv2)const{return perps_in_range(crv2);};
MGPosition_list perps(const MGSurfCurve& crv2)const;
MGPosition_list perps(const MGBSumCurve& crv2)const{return perps_in_range(crv2);};


///Obtain the projected curve of a curve onto the surface.

///The direction of the projection is along the vector vec if the vec is not NULL,
///and normal to the surface if the vec is NULL.
///Output of 'project' is two kind of curves:
///one is general world coordinate curves('vec_crv'), and the other is (u,v) curves of
///the parameter space of the surfaces(vec_crv_uv).
///vec_crv_uv.size() is equal to vec_crv.size(). Let the size be n, then
/// (vec_crv_uv[i], vec_crv[i]) is one pair for 0<=i<n.
///Function's return value is:
/// >=0: number of curves obtained, <0 : Some error detected.
int project(
	const MGFSurface& surf,	//given surface.
	std::vector<UniqueCurve>& vec_crv_uv,	//uv projection curve will be appended.
	std::vector<UniqueCurve>& vec_crv,	//3d projection curve will be appended.
	const MGVector& vec	//projection vector.
						//if vec = NULL then calculate perpendicular project.
)const;

///Round t into curve's parameter range.
double range(double t) const;

///Return space dimension
int sdim() const{return m_curve->sdim();};

///Return sweep surface from crv.

///Returned is a newed MGSurface, must be deleted.
///The sweep surface is defined as:
///This curve(say c(t)) is the rail and the straight line segments from
///C(t)+start_dist*uvec to C(t)+end_dist*uvec are the generatrix.
MGSurface* sweep(
	const MGUnit_vector& uvec,	///<Sweep Direction.
	double start_dist,			///<distance to start edge.
	double end_dist				///<distance to end edge.
)const;				

///Return curve type.
MGCURVE_TYPE type() const{return MGCURVE_TYPE::MGCURVE_TRIMMED;}

///Unlimit parameter range of the curve.
MGCurve& unlimit();

///Unlimit parameter range of the curve to the end point direction
MGCurve& unlimit_end();

///Unlimit parameter range of the curve to the start point direction
MGCurve& unlimit_start();

///Get the name of the class.
std::string whoami()const{return "TrimmedCurve";};

protected:

///Compute intersection point of 1D sub curve of original curve.
///Parameter values of intersection point will be returned.
MGCParam_list intersect_1D(						
	double f,			///< Coordinate value
	int coordinate=0	///< Coordinate kind of the data f(from 0).
)const;	

///Obtain so transformed 1D curve expression of this curve that
///f(t)={sum(xi(t)*g[i]) for i=0(x), 1(y), 2(z)}-g[3], where f(t) is the output
///of oneD and xi(t) is i-th coordinate expression of this curve.
///This is used to compute intersections with a plane g[4].
std::unique_ptr<MGCurve> oneD(
	const double g[4]			///<Plane expression(a,b,c,d) where ax+by+cz=d.
) const;

///メンバデータを読み込む関数
/// 戻り値boolは正常に読み出しが出来ればtrue、失敗すればfalseになる
void ReadMembers(MGIfstream& buf);

///メンバデータを書き込む関数
/// 戻り値boolは正常に書き込みが出来ればtrue、失敗すればfalseになる
void WriteMembers(MGOfstream& buf) const;
	
private:

	const MGCurve* m_curve;	///<Curve pointer.
	mutable bool m_sameRange;///<Indicates if m_range is the same as the original m_curve's.
	MGInterval     m_range; ///<Parameter range of the above m_curve.
	mutable std::unique_ptr<MGKnotVector> m_knotV;///<When knot_vector() is invoked, the knot vector is set.

//////////////////////////////////////////////////////////////////////////////////////
///Following functions are defined to prohibit their use in TrimmedCurve.///
//////////////////////////////////////////////////////////////////////////////////////

///Exchange ordering of the coordinates.
///Exchange coordinates (i) and (j).
void coordinate_exchange(int i, int j);

///Intersection of Curve.
MGCCisects isect_in_range(const MGCurve& curve2)const;

///Compute intersections with MGLBRep curve2 that does not have C0 continuity in it.
MGCCisects isect_withC1LB(const MGLBRep& curve2)const;

///Intersection with a surface.
MGCSisects isect_in_range(const MGSurface& surf)const;

///Exclusive function for common.
///When this curve is a TrimmedCurve of MGComposite, narrow the parameter range
///by this m_range.
void narrow_into_range(
	const MGCurve& curve2,
	std::vector<double>& vecComSpan
)const;

///Negate the curve direction(曲線の方向を反転する)
void negate();

///Compute all the perpendicular points of this curve and the second one.
///That is, if f(s) and g(t) are the points of the two curves f and g,
///then obtains points where the following conditions are satisfied:
///  fs*(f-g)=0.    gt*(g-f)=0.
///Here fs and gt are 1st derivatives at s and t of f and g.
///MGPosition P in the MGPosition_list contains this and crv's parameter
///as:     P(0)=this curve's parameter, P(1)=crv's parameter value.
MGPosition_list perps_in_range(
	const MGCurve& crv2		///<The second curve
)const;

///ノット削除関数(B表現曲線のみ)
///トレランスはline_zeroを使用する。元のノットが細かいものほど削除しやすい
///戻り値は、削除したノットの数
///removal knot. line_zero tolerance is used.
void remove_knot();

///Object transformation.
MGTrimmedCurve& operator+=(const MGVector& v);
MGTrimmedCurve& operator-=(const MGVector& v);
MGTrimmedCurve& operator*=(double scale);
MGTrimmedCurve& operator*=(const MGMatrix& mat);
MGTrimmedCurve& operator*=(const MGTransf& tr);

///Transformation object construction(NOT ALLOWED TO USE).
MGTrimmedCurve operator+ (const MGVector& v) const;
MGTrimmedCurve operator- (const MGVector& v) const;
MGTrimmedCurve operator* (double scale) const;
MGTrimmedCurve operator* (const MGMatrix& mat) const;
MGTrimmedCurve operator* (const MGTransf& tr) const;

friend class MGSurfCurve;
friend class MGLBRep;

};

/** @} */ // end of GEO group
#endif
