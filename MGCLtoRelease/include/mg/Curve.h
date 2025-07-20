/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/

#ifndef _MGCurve_HH_20190406
#define _MGCurve_HH_20190406

#include <vector>
#include <memory>
#include "mg/MGCL.h"
#include "mg/Default.h"
#include "mg/Position.h"
#include "mg/Geometry.h"
#include "mg/KnotVector.h"
#include "mg/CCisects.h"
#include "mg/CSisects.h"
#include "mgGL/VBO.h"

//
//Define MGCurve Class.

class MGInterval;
class MGBox;
class MGVector;
class MGUnit_vector;
class MGPosition_list;
class MGTransf;
class MGCParam_list;
class MGPoint;
class MGCurve;
class MGStraight;
class MGEllipse;
class MGLBRep;
class MGRLBRep;
class MGSurfCurve;
class MGBSumCurve;
class MGTrimmedCurve;
class MGCompositeCurve;
class MGSurface;
class MGPlane;
class MGSphere;
class MGCylinder;
class MGSBRep;
class MGRSBRep;
class MGBSumSurf;
class MGCCisects;
class MGCSisects;
class MGIfstream;
class MGOfstream;
class MGFace;
class MGShell;
class MGCFisects;
class MGPPRep;
class MGCommonON;
class mgVBO;

/** @file */

/** @addtogroup GEO
 *  @{
 */

///@cond
#define GRAPH_LENGTH_DONOM 5.//For curvatureLengthDisplay().
///@endcond

///MGCurve is an abstract class which represents a whole curve.
class MG_DLL_DECLR MGCurve:public MGGeometry{

public:

//////////////Constructor///////////////

///Void constructor(初期化なしでオブジェクトを作成する。).
MGCurve();

///Copy constructor.
MGCurve(const MGCurve& curve);

/// Virtual Destructor.
virtual ~MGCurve();

//////////// Operator overload(演算子多重定義) ////////////

///Assignment.

///When the leaf object of this and geo2 are not equal, this assignment
///does nothing.
virtual MGCurve& operator=(const MGCurve& gel2){
	MGGeometry::operator=(gel2);return *this;
};

////////////Logical operator overload/////////

///Object transformation.
virtual MGCurve& operator+=(const MGVector& v)=0;
virtual MGCurve& operator-=(const MGVector& v)=0;
virtual MGCurve& operator*=(double scale)=0;
virtual MGCurve& operator*=(const MGMatrix& mat)=0;
virtual MGCurve& operator*=(const MGTransf& tr)=0;

//////////// Member Function ////////////

///Approximate this curve as a MGLBRep.

///Approximate this curve as a MGLBRep curve
///within the tolerance MGTolerance::line_zero().
///When parameter_normalization=0, reparameterization will not done, and
///the evaluation at the same parameter has the same values before and after
///of approximate_as_LBRep.
virtual void approximate_as_LBRep(
	MGLBRep& lb,///<Approximated lbrep will be set.
	int ordr=0,	///<new order. When this is MGLBRep, if ordr=0,
				///<ordr=order() will be assumed, else ordr=4 is assumed.
	int parameter_normalization=0,
		///<Indicates how the parameter normalization be done:
		///<   =0: no parameter normalization.
		///<   =1: normalize to range=(0., 1.);
		///<   =2: normalize to make the average length of the 1st derivative 
		///<       is as equal to 1. as possible.
	bool neglectMulti=false///<Indicates if multiple knots be kept.
		///< true: multiplicity is removed.
		///< false: multiplicity is kept.
)const;

///Generate arrow data of the tangent at the parameter value t of the curve.

///data[0] is the origin, data[1] is top of the arrow,
///data[2], [3] are two bottoms of arrowhead.
void arrow(double t,MGPosition data[4])const;

///Returns B-Rep Dimension.
virtual int bdim() const=0;

///Return minimum box that includes the curve of parameter interval.

/// 入力のパラメータ範囲の曲線部分を囲むボックスを返す。
virtual MGBox box_limitted(
	const MGInterval& ///< Parameter Range of the curve.
) const = 0;

///Obtain ceter coordinate of the geometry.
virtual MGPosition center() const;

///Obtain ceter parameter value of the geometry.
virtual MGPosition center_param() const;

///Changing this object's space dimension.
virtual void change_dimension(
	int sdim,			///< new space dimension
	int start1=0, 		///< Destination order of new object.
	int start2=0 		///< Source order of this object.
)=0;

///Change parameter range.

///Be able to change the direction by providing t1 greater than t2.
virtual void change_range(
	double t1,	///<Parameter value for the start of original. 
	double t2	///<Parameter value for the end of original. 
)=0;

///Construct new geometry object by copying to newed area.

///User must delete this copied object by "delete".
virtual MGCurve* clone()const=0;

///Compute the closest point parameter value of this curve from a point.
virtual double closest(const MGPosition& point) const;

///Compute the nearest point from input point on this curve's (x,y) 2D part.
virtual double closest2D(const MGPosition& point) const;

///Compute the closest point parameter value pair of this curve and curve2.

///MGPosition P of the function return contains this and curve2's parameter
///as:     P(0)=this curve's parameter, P(1)=curve2's parameter value.
virtual MGPosition closest(const MGCurve& curve2) const;

///Test if this curve is cn continuous.

///曲線がCn連続かどうか調べる
///LBRep以外はかならずtrueが返却される
bool cn_continuity(int n)const;

///Test if this has a common line part with the 2nd curve.

///目的：与えられた曲線と自身の共通部分があるかどうか調べる。
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
///	誤差には、パラメータ範囲*rc_zero()を使用する。
virtual int common(
	const MGCurve& curve2,
	std::vector<double>& vecComSpan,
	MGCCisects& isect
) const;

///Test if this has a common line part with the 2nd curve.

///関数名：common
///目的：与えられた曲線と自身の共通部分があるかどうか調べる。
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
///  パラメータ範囲*rc_zero()を使用する。
virtual int common(
	const MGCurve& curve2,
	std::vector<double>& vecComSpan
) const;

///Exchange ordering of the coordinates.

///Exchange coordinates (i) and (j).
virtual void coordinate_exchange(int i, int j)=0;

///copy as a newed curve.

///The new curve will be MGLBRep or MGRLBRep.
///When original curve was a MGRLBRep, the new curve will be a MGRLBRep.
///Otherwise,  the new curve will be a MGLBRep.
///Returned object must be deleted.
virtual MGCurve* copy_as_nurbs() const=0;

///Convert this curve to Bezier curve.

///If this is MGLBRep or MGStraight, the shape is exactly the same
///as the original. Otherwise, this is apporoximated by MGLBRep.
virtual void convert_to_Bezier(MGLBRep& bezier)const;

///Construct new curve object by changing the original object's space dimension.

///Returned is a newed object, 
///user must delete this copied object by "delete".
virtual MGCurve* copy_change_dimension(
	int sdim,			///< new space dimension
	int start1=0, 		///< Destination order of new line.
	int start2=0 		///< Source order of this line.
)const=0;

///Construct new curve object limitting the parameter range to prange.

///Construct new curve object by copying to newed area,
///and limitting the parameter range to prange.
///Returned is a newed object and must be deleted.
virtual MGCurve* copy_limitted(const MGInterval& prange) const;

///Return curvature at the given point.

///When the curve is 2D, curvature has sign. when 3D, curvature is
///always plus.
/// 与えられた点における曲線の曲率を返却する。
virtual double curvature( double ) const;

///Compute curvilinear integral of the 1st two coordinates.

///線積分を求める。
///This integral can be used to compute area sorounded by the curve.
///Second form is from param_s() to param_e();
///curvilinear_integral from t1 to t2 can be obtained by
///Integral of (x*dy-y*dx) about t, where curve is expressed by
///f(t)=(x(t),y(t)), dx=dx/dt, and dy=dy/dt.
virtual double curvilinear_integral(double t1, double t2) const;
virtual double curvilinear_integral() const{
	return curvilinear_integral(param_s(), param_e());
};

///Compute direction unit vector of the geometry.
MGUnit_vector direction(const MGPosition& param) const;

///Return tangent vector at the given point.

/// 与えられた点における曲線の接ベクトルを返す。
virtual MGUnit_vector direction(double) const;

///////display member function.

///Display direction arrows of this curve.
virtual void display_arrows(mgSysGL& sgl)const;

///Display break points of this curve.
virtual void display_break_points(mgSysGL& sgl)const;

///Display curvature of this curve.
virtual void display_curvatures(
	mgSysGL& sgl,///<Sgl to make curvature pictures in.
	int		density,///<densitiy of the graph.
	bool	use_radius,///<true:radius display, false:curvature display.
	double	scale=1.	///<scaling of the graph.
)const;

///Divide this curve at the designated knot multiplicity point.
///Function's return value is the number of the curves after divided.
virtual int divide_multi(
	std::vector<UniqueCurve>& crv_list,	//divided curves are appended.
	int multiplicity=-1	///<designates the multiplicity of the knot to divide at,
						///<When multiplicity<=0, order()-1 is assumed,
						///<When multiplicity>=order(), order() is assumed.
)const;

///get the a divide number for offset, intersection, or others.
virtual int divide_number() const{return divideNum(param_range());};

///Draw this curve into vbo, approximating with polyline.
virtual void drawSE(
	mgVBO& vbo,///<Target graphic object.
	double t0,			///<Start parameter value of the curve.
	double t1			///<End parameter value of the curve,
						///<Draw will be performed from t0 to t1.
)const;

///Draw this curve into vbo, approximating with polyline.
virtual void drawWire(
	mgVBO& vbo,///<Target graphic object.
	int line_density=1	///<line density to draw a surface in wire mode.
)const{drawSE(vbo,param_s(),param_e());};

///Return end point(終点を返却する)
virtual MGPosition end_point() const;

///@brief Evaluate n'th derivative data.
///n=0 means positional data evaluation.
virtual MGVector eval(
	double,				///< Parameter value.
	int nderiv=0,	///< Order of Derivative.
	int left=0			///<Left continuous(left=true)
						///<or right continuous(left=false).
) const = 0;

///Compute position, 1st and 2nd derivatives.
/// パラメータ値を与えて位置、一次微分値、二次微分値をもとめる。
virtual void eval_all(
	double,			///<Input parameter value(パラメータ値)
	MGPosition&,	///<Position(位置)
	MGVector&,		///<1st derivative(1次微分値)
	MGVector&		///<2nd derivative(2次微分値)
) const;

///Compute 1st derivative.
virtual MGVector eval_deriv(double) const;

///@brief Evaluate deviations of two curves(this and curve2) at npoint discrete points.
///(1)Search the common curve spans which have the distance within tolerance.
///(2)Compute the nearest points from npoint discrete points of this to curve2.
///Let sti=sts[i], then
///sti[0] is this curve's parameter value s, and sti[1] is the parameter value t
///of curve2 which is the nearest point from the point s.
///If this and curve2 have the minimum distance more than tolerance,
///sts.size()==1 and sts[0] is the minimum distance points of this and curve2.
void eval_discrete_deviation(
	const MGCurve& curve2,///<2nd target curve.
	std::vector<MGPosition>& sts,///<Parameter values of this and curve2 will be output.
			///<sts[i] is i-th parameter value.
	int npoint=20,		///<indicates how many discrete points be obtained.
	double tolerance=0.1///<tolerance to get two edge to compute deviation.
)const;

///Evaluate line data at data point tau.
virtual void eval_line(
	const MGNDDArray& tau,	///<Data points.
	MGBPointSeq& value		///<Values evaluated. value(i,.)=eval(tau[i]);
)const;

///Compute positional data.
virtual MGPosition eval_position(double) const;

/// Evaluate n'th derivative data. n=0 means positional data evaluation.
MGVector evaluate(
	const MGPosition& t,	///< Parameter value.
				///<t's space dimension is geometry's manifold dimension.
	const int* nderiv=0	///<Order of derivative of i-th parameter
				///<in nderiv[i].
				///<When nderiv=null, nderiv[i]=0 is assumed for all i.
) const;

///Extrapolate this curve by an (approximate) chord length.

///The extrapolation is C2 continuous.
virtual void extend(
	double length,	///<approximate chord length to extend. 
	bool start=false///<Flag of which point to extend, start or end point of the line.
					///<If start is true extend on the start point.
)=0;

///Compute Frenet_frame, curvature and torsion in 3D space.
virtual void Frenet_frame2(
	double t,			///<Input parameter value(パラメータ値)
	MGVector& V2,///<2nd derivative at t.
	MGVector& T,///<Tangent
	MGVector& N,///<Principal Normal
	MGVector& B	///<Binormal
)const;

///Compute Frenet_frame, curvature and torsion in 3D space.
virtual void Frenet_frame(
	double t,			///<Input parameter value(パラメータ値)
	MGVector& T,	///<Tangent(Unit vector)
	MGVector& N,	///<Principal Normal(Unit vector)
	MGVector& B,	///<Binormal(Unit vector)
	double& curvature,	///<Curvature is always >=0.
	double& torsion		///<Tortion
)const;

///Get average tangent length.

///The average means the average of start, mid, and end point's tangent.
double get_average_tangent_length()const;

///Extracts control points.

///Fucntion's return value is 
///true if control points was obtained, false if not.
virtual bool get_control_points(
	MGBPointSeq& cpoints	///<Control points will be output.
)const{return false;};

///Find C0 cotinuity parameter values o fhtis curve,
///and push back the parameter values to param.

///Obtained parameter values are pushed in ascending order.
void getParamsC0Continuity(std::vector<double>& param)const;

///Computes an approximate maximum curvature graph length.
///Function's return value is the length.
double curvatureLengthDisplay(
	bool use_radius = false//true if curvature radius is to output.
)const;

///Test if this curve has the same direction with curve2.

///Test at the point s(of this) and t(of curve2).
///Function's return value is true if they have the same direction.
///"same direction" means their tangent vectors have the angle less than 90 degree.
bool has_same_direction_at(
	double s,
	const MGCurve& curve2,
	double t
)const;

/// Return This object's typeID
virtual long identify_type() const=0;

///Test if input parameter value is inside parameter range of the line.
virtual bool in_range(double t) const;

///Test if input parameter value is inside parameter range of the line.
bool in_range(const MGPosition& t) const;

///Curve to curve intersection.

/// Curve と Curve の交点を求める。
MGCCisects intersect_brute_force(const MGCurve&) const;

///Curve to curve intersection.

///***Caution***intersect can be used only for finite curve, i.e.
///parameter range of the computation is only from param_s() to param_e().
///For example, intersect cannot be applied to infinite straight line.
virtual MGCCisects intersect(const MGCurve&) const;

///Compute the intersections of two objects.

///Provide divide number of curve span for function intersect.
virtual int intersect_dnum()const=0;

///intersections with a plane.
MGCSisects intersect_with_plane(const MGPlane& surf)const;


///Compute the intersections of two objects.
///Intersections are obtained from two objects, which are known using
///the MGisects::object1() and object2().
///****NOTE****
///When two objects' manifold dimension are the same, object1 is this object
///at the invocation of MGObject::intersection(), and object2 is the argument
///object.
///However, their manifold dimension are not the same, object1 is always
///the lower dimension's object and object2 is the higer dimension's object.
MGisects isect(const MGObject& obj2)const override;
virtual MGCCisects isect(const MGCurve& crv2)const=0;

///Intersection of Curve and other geometry.
virtual MGCCisects isect(const MGStraight& curve2)const;
virtual MGCCisects isect(const MGLBRep& curve2)const;
virtual MGCCisects isect(const MGSurfCurve& curve2)const;
MGCCisects isect(const MGTrimmedCurve& curve2)const;
MGCCisects isect(const MGCompositeCurve& curve2)const;

private:
MGCSisects isectSurf(const MGSurface&)const;

public:
virtual MGCSisects isect(const MGSurface&f)const;
virtual MGCSisects isect(const MGPlane& f)const;

virtual MGCSisects isect(const MGFace&)const;

///Compute intersection point of 1D sub curve of original curve.

///Parameter values of intersection point will be returned.
MGCParam_list isect_1D(						
	double f,			///< Coordinate value
	int coordinate=0	///< Coordinate kind of the data f(from 0).
) const;

///Test if this is a Bezier Curve.

///Functions's return value is true if Bezier, false if not.
///If input ordr>=2, order is also tested if this Bezier's order is the same as input order.
///If input ordr<=1, any ordr>=2 is allowed for Bezier curve.
///Bezier curve is defined as follows. Here t=knot_vector(), k is this LBRep's order,
///n=bdim(), and m=(n-k)/(k-1).
///(1) n=k+(k-1)*m.
///(2) t(0)=t(1)=,...,=t(k-1)=0
///(3) t(i)=t(i+1)=,...,=t(i+k-2)=j+1
///         for i=k, k+(k-1),...,k+j*(k-1) and j=0,...,m-1.
///(4) t(n)=t(n+1)=,...,=t(n+k-1)=m+1
virtual const MGLBRep* is_Bezier(int ordr=0)const;

///Test if this is a closed curve.
bool is_closed()const{return start_point()==end_point();};

///Terst if this is a closed curve, given the tolerance.
bool is_closedWithError(double err)const;

///Test if this cure is co-planar with the 2nd curve curve2.

///MGPlane expression will be out to plane if this is co-planar.
///Function's return value is true if co-planar.
virtual bool is_coplanar(const MGCurve& curve2, MGPlane& plane)const;

///Test if the input parameter t is the start point parameter or not.
virtual bool is_startpoint_parameter(double t)const;

///Test if the input parameter t is the start point parameter or not.
virtual bool is_endpoint_parameter(double t)const;

///Test if the vector from P to this->eval(t) is perpendicular.

///Perpendicular to the tangent of this curve at t.
bool is_perpendicular(const MGPosition& P, double t)const;

///Test if this cure is linear or not, that is, is straight or not.

///MGStraight expression will be out to straight even if this is linear or not.
///Function's return value is true if linear.
virtual bool is_linear(MGStraight& straight)const;

///Test if this cure is planar or not.

///MGPlane expression will be out to plane if this is planar.
///Function's return value is true if planar.
virtual bool is_planar(MGPlane& plane)const;

///Access to i-th element of knot
virtual double knot(int i) const=0;
	
///Returns the knot vector of the curve.
virtual const MGKnotVector& knot_vector() const=0;

///Returns the knot vector of the curve.
MGKnotVector& knot_vector();

///Cmpute curve length of the interval.

///If t1 is greater than t2, return negative value.
/// 与えられたパラメータ値間の曲線の長さを返す。
/// パラメータが昇順で与えられたときは正値、降順のときは負値を返す。
virtual double length(double t1, double t2) const;

///Compute the whole curve length.

///If the curve is infinite, return -1.
/// 自身の曲線が有界の場合、その曲線の距離を返却する。非有界の場
/// 合はー１を返却をする。
virtual double length() const {return length(param_s(), param_e());}

///Inverse function of length.

///Compute the point that is away from the point t by length len.
/// lengthの逆関数。指定パラメータtで示される点から指定距離len
/// 曲線上に沿って離れた点を示すパラメータ値を返す。
virtual double length_param( double t, double len) const;

///Update this by limiting the parameter range of the curve.
virtual void limit(const MGInterval& rng) = 0;
void limit(double t0, double t1);

///Return manifold dimension, 0:point, 1:curve, 2:surface.
int manifold_dimension() const{ return 1;};

///Negate the curve direction(曲線の方向を反転する).
virtual void negate() = 0;

///Obtain the parameter value to t when this curve is negated by "negate()".
virtual double negate_param(double t)const= 0;

///Transform the coordinates of boundary of this geometry so that
///new coordinate of boundary is the same coordinate as the new one of
///this geometry after negate() of this geometry is done.
///That is, boundary coordinates are parameter world of this geometry.
void negate_transform(MGGeometry& boundary)const;

/// Offset curve of constant deviation.  This must be C1 cotinuous everywhere.
/// Although the offset value ofs_value must be less than radius of curvature,
/// this is unchecked.
/// The offset direction is Normal(obtained by Frenet_frame()) at each point.
UniqueLBRep offsetC1(
	double ofs_value, ///offset value, may be negative.
	bool principalNormal = true /// true: Offset direction is to principal normal
								/// false: to binormal
) const;

/// Offset of constant deviation from this curve.
/// Although the offset value ofs_value must be less than radius of curvature,
/// this is unchecked.
/// The direction of offset is toward the principal normal,
/// or to the direction to center of curvature(if offset value is positive).
/// When this curve is not C1 continuous, this is divided into C1 curves,
/// and more than one offset curves are obtained.
/// line_zero() is used to approximate curves of the offset.
virtual std::vector<UniqueCurve> offset(
	double ofs_value, ///offset value, may be negative.
	bool principalNormal=true /// true: Offset direction is to principal normal
							  /// false: to binormal
) const;

/// Offset of variable deviation from this curve.
/// Although the offset value ofs_value must be less than radius of curvature,
/// this is unchecked.
/// The direction of offset is toward the principal normal,
/// or to the direction to center of curvature(if offset value is positive).
/// When this curve is not C1 continuous, divided into C1 curves,
/// and more than one offset curves are obtained.
/// line_zero() is used approximate the offset curve.
virtual std::vector<UniqueCurve> offset(
	const MGLBRep& ofs_value_lb,	///offset value, may be negative.
			///ofs_value_lb's parameter space must be the same as this.
			///ofs_value_lb's space dimension is supposed to be one.
	bool principalNormal = true /// true: Offset direction is to principal normal
								/// false: to binormal
) const;

///@cond
///get the a divide number for offset, intersection, or others.
virtual int divideNum(
	const MGInterval& interval	///<分割数を求めるパラメータ範囲
)const;
///@endcond

///Test if given point is on the curve or not.
///If given point is on the curve, return parameter
///value of the curve. Even if not, return nearest point's parameter t.
/// 指定点が自身上にあるかを調べる。曲線上にあれば，そのパラメーター値を，
/// なくても最近傍点のパラメータ値を返す。
/// Function's return value is >0 if the point is on the curve,
/// and 0 if the point is not on the curve.
virtual bool on(
	const MGPosition& point,///<point(指定点)
	double& t	///<Parameter of the curve(パラメータ) will be returned.
) const;

///Test if given point is on this geometry or not.
///If the point is on this geometry, return parameter
///value of the geometry. Even if not, return nearest point's parameter.
/// 指定点が自身上にあるかを調べる。曲線上にあれば，そのパラメーター値を，
/// なくても最近傍点のパラメータ値を返す。
/// Function's return value is >0 if the point is on the geometry,
/// and 0 if the point is not on the geometry.
bool on(
	const MGPosition& P,///<Point(指定点)
	MGPosition&	t		///<Parameter of the geometry(パラメータ)
) const;

///Returns the order.
virtual	int order() const=0;

///Compute parameter value of given point.

/// 自身の上の指定点を表すパラメータ値を返す。
/// If input point is not on the curve, return the nearest point on the
/// curve.
virtual double param(
	const MGPosition &	///<Point(指定点)
) const;

///Return ending parameter value.
virtual double param_e() const=0;

///Obtain parameter space error.
virtual double param_error() const;

///Normalize parameter value t to the nearest knot if their distance is within tolerance.
virtual double param_normalize(double t) const=0;

///Return parameter range of the curve(パラメータ範囲を返す).
virtual MGInterval param_range() const;

///Round the parameter t into this parameter range.
double param_round_into_range(double t)const;

///Return parameter range of the geometry(パラメータ範囲を返す).
MGBox parameter_range() const;

/// Return starting parameter value.
virtual double param_s() const=0;

/// Return starting or ending parameter value that is nearer to the param t.
double param_se(double t) const;

///Compute parameter span length 
virtual double param_span() const{return param_e()-param_s();};

///Compute part of this curve from parameter t1 to t2.

///Returned is the pointer to newed object, and so should be deleted
///by calling program, or memory leaked.
virtual MGCurve* part(
	double t1,///<Start prameter value.
	double t2,///<End prameter value.
	int multiple=0	///<Indicates if start and end knot multiplicities
					///<are necessary. =0:unnecessary, !=0:necessary.
) const=0;

///Return perpendicular point from a point P.

///Given guess starting paramter values.
///Function's return value is:
///   perp_guess=true if perpendicular points obtained,
///   perp_guess=false if perpendicular points not obtained,
virtual int perp_guess(
	double t0,	///<parameter range of this, start.
	double t1,  ///<End. (t0>=t1) indicates no range specified.
	const MGPosition& P,	///<Point(指定点)
	double tg,				///<Guess parameter values of this curve.
	double& t				///Output parameter
) const;

///Return perpendicular points of two curves.

///Given guess starting paramter values.
///Function's return value is:
///   perp_guess=true if perpendicular points obtained,
///   perp_guess=false if perpendicular points not obtained,
virtual int perp_guess(
	double s0,	///<parameter range of this, start.
	double s1,	///< End. <When s0>=s1, no limit for this parameter range.
	const MGCurve& curve2,	///<2nd curve.
	double t0,	///<parameter range of curve2, start.
	double t1,	///< End. When t0>=t1, no limit for curve2 parameter range.
	double sg,	///<Guess parameter values of the two curves, for this.
	double tg, 	///<curve2's parameter.
	MGPosition& st	///<perpendicular points' parameter values
					///<will be output.
	///<st(0): this curve's parameter, st(1):curve2's parameter.
) const;

///Compute a foot point of the perpendicular line from point p to the curve.

///If more than one points are found, return nearest one.
/// 指定点からの自身への垂線の足とパラメータ値を返す。
/// Function's return value is if point is obtained(1) or not(0)
virtual int perp_point(
	const MGPosition& p,	///<Point(指定点)
	double& t,				///<Parameter of the curve(パラメータ値)
	const double* g=0		///<guess parameter value of line
) const;

///Compute all the perpendicular points of this curve and the second one.

///If f(s) and g(t) are the points of the two curves f and g,
///then obtains points where the following conditions are satisfied:
///  fs*(f-g)=0.    gt*(g-f)=0.
///Here fs and gt are 1st derivatives at s and t of f and g.
///**** NOTE 1 ****
///perpendiculars is general function of perps, used in perps.
///General users should use function perps, not perpendiculars, since
///perps is optimized for each curve type.
///**** NOTE 2 ****
///perpendiculars can not be used for infinite parameter range curve.
///param_s() and param_e() of both curves must return their finite
///parameter range.
///MGPosition P in the MGPosition_list contains this and crv's parameter
///as:     P(0)=this curve's parameter, P(1)=crv's parameter value.
MGPosition_list perpendiculars(
	const MGCurve& crv		///<The second curve
) const;

///Compute all foot points of the perpendicular line from point to the curve.

/// 与ポイントから曲線へ下ろした垂線の足の，曲線のパラメータ値を
/// すべて求める。
virtual MGCParam_list perps(
	const MGPosition& P		///<Point(指定点)
) const;

///Compute all the perpendicular points of this curve and the second one.

///That is, if f(s) and g(t) are the points of the two curves f and g,
///then obtains points where the following conditions are satisfied:
///  fs*(f-g)=0.    gt*(g-f)=0.
///Here fs and gt are 1st derivatives at s and t of f and g.
///MGPosition P in the MGPosition_list contains this and crv's parameter
///as:     P(0)=this curve's parameter, P(1)=curve2's parameter value.
virtual MGPosition_list perps(const MGCurve& crv2)const=0;
virtual MGPosition_list perps(const MGStraight& crv2)const=0;
virtual MGPosition_list perps(const MGRLBRep& crv2)const;
virtual MGPosition_list perps(const MGEllipse& crv2)const;
virtual MGPosition_list perps(const MGLBRep& crv2)const;
virtual MGPosition_list perps(const MGSurfCurve& crv2)const;
virtual MGPosition_list perps(const MGBSumCurve& crv2)const;
MGPosition_list perps(const MGCompositeCurve& crv2)const;
MGPosition_list perps(const MGTrimmedCurve& crv2)const;

///Compute the parameter value of the closest point from the straight to this object.

///sl is the eye projection line whose direction is from yon to hither, and if
///sl had multiple intersection points, The closest point to the eye will be selected.
virtual MGPosition pick_closest(const MGStraight& sl)const;

///Approximate this curve by a polyline and output to lb2.

///The tolerance of the approximation is error.
virtual void polygonize(
	double error,	///<tolerance allowed for the approximation
	MGLBRep& lb2	///<Obtained polyline will be output as an MGLBRep of order2.
)const;


///Obtain the projected curve of a curve onto the surface.

//曲線を面に面直またはベクトル投影して曲線リストを求める。
///投影曲線は面上のパラメータ曲線と3次元曲線としてそれぞれ順番に、
///vec_crv_uv, vec_crvに格納される。
///uv曲線のトレランスはrc_zero()を、3次元曲線はline_zero()をそれぞれ使用している。
///get perpendicular or vector projection curve list.
///uv projection curves are put into vec_crv_uv(rc_zero() is used),
///3d projection curves are put into vec_crv(line_zero() is used) respectively.
//戻り値：
///		投影曲線の数:		投影曲線が求まった
///		0:			投影曲線が求まらなかった
///		-1:			内部処理エラー
///		-2:			収束処理エラー（収束しなかった）
///追記：引数vecが与えられない(null）とき、面直投影する。
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
virtual int project(
	const MGFSurface& surf,	//given surface.
	std::vector<UniqueCurve>& vec_crv_uv,	//uv projection curve will be appended.
	std::vector<UniqueCurve>& vec_crv,	//3d projection curve will be appended.
	const MGVector& vec	//projection vector.
						//if vec = NULL then calculate perpendicular project.
)const;

///Round t into curve's parameter range.

/// 入力パラメータをパラメータ範囲でまるめて返却する。
virtual double range(double t) const;

///Round t into geometry's parameter range.

/// 入力パラメータをパラメータ範囲でまるめて返却する。
///t's space dimension is geometry's manifold dimension.
MGPosition range(const MGPosition& t) const;

///Rebuild this curve.
std::unique_ptr<MGCurve> rebuild(
	int how_rebuild=1,
		///< intdicates how rebuild be done.
		///<  =0: no approximation(only parameter change)
		///<  =1: if this is rational spline(MGRLBRep), reconstructed with new knot configuration
		///<      as rational spline(MGRLBRep).
		///<      Otherwise approximated by non-rational spline(MGLBRep) with new knot configuration.
		///<  =2: approximated by non-rational spline(MGLBRep) with new knot configuration
		///<      if this is rational spline. If this is not rational spline, same as =1.
	int parameter_normalization=2,
		///< Indicates how the parameter normalization be done:
		///< =0: no parameter normalization.
		///< =1: normalize to range=(0., 1.);
		///< =2: normalize to make the average length of the 1st derivative 
		///<     is as equal to 1. as possible.
		///< =3: specify parameter range in param_range.
	double tol=-1.,	///<tolerance allowed for the approximation
		///< When tol<=0., MGTolerance::line_zero() will be employed.
	int ordr=0,	///<order of the new MGLBRep, >=4 is recommended.
		///< When order=0 is input, the original order is unchanged if this curve is
		///< MGLBRep or MGRLBRep. Otherwise order is set to 4.
	const double* param_range=0
		///<Input new paramter range, from param_range[0] to param_range[1].
		///<When param_range=0, parameter_normalization is set to 2 and the starting parameter is set to 0.
)const;

///Remove redundant knot, and reduce the b-rep dimension.

///ノット削除関数(B表現曲線のみ)
///トレランスはline_zeroを使用する。元のノットが細かいものほど削除しやすい
///Remove redundant knot, and reduce the b-rep dimension.
///The tolerance used is MGTolerance::line_zero().
virtual void remove_knot();

///Update curve by rotating around straight line.

/// 指定点を通る指定ベクトルを軸として回転させたものを自身とする。
virtual MGCurve& rotate_self(
	const MGVector& v,			///<Vector of the line to rotate around.
	double,						///<Angle of rotation.
	const MGPosition & = mgORIGIN	///<A point on the line to rotate around.
);

///Obtain polar-rotated curve of this.

///This curve's (x,y) are updated. No other coordinates are unchanged.
///The returned curve is always MGLBRep.
///Rotation is performed from angle1 to angle2, around angleBase.
///That is, when angle1=angle2, no change is done.
///When angle2 is angleBase, all the data will lie on the straight of from origin to
///(cos(angleBase), sin(angleBase)).
std::unique_ptr<MGLBRep> scalePolar(
	double angleBase,///<base angle.
	double angle1,	///< 1st angle.
	double angle2	///< 2nd angle.
)const;

///Return space dimension.
virtual int sdim() const =0 ;

///Return start point(始点を返却する).
virtual MGPosition start_point() const;

///Return sweep surface from crv.

///Returned is a newed MGSurface, must be deleted.
///The sweep surface is defined as:
///This curve(say c(t)) is the rail and the straight line segments from
///C(t)+start_dist*uvec to C(t)+end_dist*uvec are the generatrix.
virtual MGSurface* sweep(
	const MGUnit_vector& uvec,	///<Sweep Direction.
	double start_dist,			///<distance to start edge.
	double end_dist			///<distance to end edge.
) const =0;

///Return tangent point from a point P, given guess starting paramter tg.
///Searching is done only from t0 to t1.
virtual int tangent_guess(
	double t0,///< parameter range of this, starting.
	double t1,///< Ending. If t0>t1, whole range of this.
	const MGPosition& P,	///<Point(指定点)
	double tg,				///<Guess parameter values of the two curves
	double& t				///<Output parameter
)const;

///Trim the end part of this curve at the parameter t.

///The new curve range is [start_of_original, t]
///t must be inside this parameter rage, else does nothing.
void trim_end(double t);

///Trim the start part of this curve at the parameter t.

///The new curve range is [t,end_of_original]
///t must be inside this parameter rage, else does nothing.
void trim_start(double t);

///Trim the start part and end part of this curve at the parameter ts and te.

///The new curve range is [ts,te]
///Both ts and te must be inside this parameter rage.
void trim_start_and_end(double ts, double te);

///Return curve type(曲線のタイプを返す).
virtual MGCURVE_TYPE type() const =0;

///Unlimit parameter range of the curve(limitをはずす).
virtual MGCurve& unlimit() =0;

///Unlimit parameter range of the curve to the end point direction.
virtual MGCurve& unlimit_end() =0;

///Unlimit parameter range of the curve to the start point direction.
virtual MGCurve& unlimit_start() =0;

///Obtain polar coordinates system MGLBRep of this curve.

///This curve's (x,y) coordinates are changed polar coordinates system(r,theta)
///where r is the distance from origin and theta is the angel with x coordinate.
///The space dimension of this curve must be >=2;
///If this space dimension is lager than 2, the remaining coordinates are set unchanged
///to MGLBRep.
virtual std::unique_ptr<MGLBRep> PolarCoordinatesLBRep()const;

/// Output virtual function.
virtual std::ostream& toString(std::ostream&) const;

///Get the name of the class.
virtual std::string whoami()const{return "Curve";};

protected:

///メンバデータを読み出す関数.
virtual void ReadMembers(MGIfstream& buf);

///メンバデータを書き込む関数.
virtual void WriteMembers(MGOfstream& buf) const;

///Compute intersection point of 1D sub curve of original curve.

///Parameter values of intersection point will be returned.
virtual MGCParam_list intersect_1D(						
	double f,			///< Coordinate value
	int coordinate=0	///< Coordinate kind of the data f(from 0).
) const;	

///Approximate this curve as a MGLBRep curve from knot_vector[is] to [ie].

///This is an internal program of MGLBRep constructor.
void approximate_as_LBRep2(
	MGLBRep& lb,		///<Approximated LBRep will be set.
	int order,		///<new order
	int is,///<approximation parameter range, from knot_vector()[is].
	int ie,///<approximation parameter range, to knot_vector()[ie].
	bool neglectMulti=false///<Indicates if multiple knots be kept.
		///< true: multiplicity is removed.
		///< false: multiplicity is kept.
)const;

///Get data points for approximate_as_LBRep2.
virtual void data_points_for_approximate_as_LBRep2(
	int is,///<approximation parameter range, from knot_vector()[is].
	int ie,///< to [ie].
	MGKnotVector& t,///<New knot configuration will be output.
				///<t's order is input. other information of t will be updated.
	MGNDDArray& tau,///<Data point for t will be output.
	bool neglectMulti=false///<Indicates if multiple knots be kept.
		///< true: multiplicity is removed.
		///< false: multiplicity is kept.
)const;

///Obtain an extrapolated PP-Rep curve by the parameter value.
void extrapolated_pp(
	double tau,		///<The parameter value at the end of extended point,
					///<When tau<param_s(), extension will be done at the starting point,
					///<When tau>param_e(), extension will be done at the end point.
	double dk,     ///<Coefficient of how curvature should vary at the connecting point.
	MGPPRep& pp   ///<PP-rep will be output.
)const;

///Compute intersections with MGLBRep curve2 that does not have C0 continuity in it.
virtual MGCCisects isect_withC1LB(const MGLBRep& curve2)const;

///isect with SurfCurve whose m_curve is not a MGTrimmedCurve of MGCompositeCurve.
virtual MGCCisects isect_with_noCompoSC(const MGSurfCurve& curve2)const;

///Obtain transformed 1D curve expression of this curve.

///Obtain so transformed 1D curve expression of this curve that
///f(t)={sum(xi(t)*g[i]) for i=0(x), 1(y), 2(z)}-g[3], where f(t) is the output
///of oneD and xi(t) is i-th coordinate expression of this curve.
///This is used to compute intersections with a plane g[4].
virtual std::unique_ptr<MGCurve> oneD(
	const double g[4]			///<Plane expression(a,b,c,d) where ax+by+cz=d.
) const=0;

///Perpendicular points with C1 conitnuity LBRep.

///MGPosition P in the MGPosition_list contains this and crv's parameter
///as:     P(0)=this curve's parameter, P(1)=crv's parameter value.
virtual MGPosition_list perps_withC1LB(
   const MGLBRep& lbC1
)const;

///Perpendicular points of this to curve2.

///curve2 is a MGSurfcurve and the composite parameter curve
///must not a MGTrimmedCurve of MGCompositeCurve.
virtual MGPosition_list perps_with_noCompoSC(const MGSurfCurve& curve2)const;

///Perpendicular points with straight.

///MGPosition P in the MGPosition_list contains this and crv's parameter
///as:     P(0)=this curve's parameter, P(1)=crv's parameter value.
MGPosition_list perpsSl(
	const MGStraight& sl	///<The second curve
)const;

private:
	
///getMaxCurvatureLengthApprox() computes an approximate maximum curvature graph length.
///Function's return value is the length.
double getMaxCurvatureLengthApprox(
	bool use_radius=false//true if curvature radius is to output.
)const;

///thisのスタートパラメータを与え共通パラメータ範囲を1つだけ取得する
///sparamはスタートパラメータであり、検索終了時のパラメータを返す
///戻り値は	1:範囲が求まった(thisの最後まで)
///			0:範囲が求まらなかった(thisの最後まで)
///			-1:範囲が求まった(thisの途中まで、まだ共通範囲があるかもしれない)
int common_one_span(
	const MGCurve& curve2,
	MGCommonON SEon[2],///<data if curve2's start or end point is on this curve.
	MGCommonON &sparam, ///<Input start parameter of this curve to search starting point of the next
					///<common part.
					///<On return from common_one_span, the next startign parameter will be set.
	double span,	///<parameter span to advance the next on point check of this curve.
	double comSpan[4]
)const;

///this上のcurve2に乗っているパラメータonparamと乗っていないパラメータoffparam
///を与え、ちょうどはずれる点を求める。
///チェックポイントの移動量 < (パラメータ範囲*rc_zero())になれば終了。
///戻り値は　収束して求まったパラメータ値である。
void common_boundary(
	const MGCurve& curve2,
	const MGCommonON& onparam,
	const MGCommonON& offparam,
	double& param1, double& param2
)const;

friend class MGFace;
friend class MGFSurface;
friend class MGSurface;
friend class MGLBRep;
friend class MGTrimmedCurve;
friend class MGSurfCurve;
friend class MGCompositeCurve;
friend class MGBSumCurve;

};

///@cond
///The class for function object for mgGausp to compute the length() of the curve.
class MGCurveLengthDrive{
	const MGCurve* m_curve;
public:
	MGCurveLengthDrive(const MGCurve* curve):m_curve(curve){;};

	///Compute the length of the 1st derivative at the curve parameter t.
	double operator()(double t) const;
};

///The class for function object for mgDefint to compute the length_param() of the curve.
class MGCurveLenParamDrive {
	const MGCurve* m_curve;
	double m_len, m_ts;
public:
	MGCurveLenParamDrive(const MGCurve* curve, double len, double ts)
		:m_curve(curve), m_len(len), m_ts(ts){;};
	double operator()(double t)const;
};
///@endcond

namespace MGCL{

///Compute curvature in 3D space, ie, the value is not negative.
MG_DLL_DECLR double Curvature(
	const MGVector& v1,		///<First derivative.
	const MGVector& v2	///<Second derivative.
);

///Compute torsion.
MG_DLL_DECLR double Torsion(
	const MGVector& v1,		///<First derivative.
	const MGVector& v2,		///<Second derivative.
	const MGVector& v3	///<Third derivative.
);
	
///Generate arrow data from (root, vecx, vecy).
MG_DLL_DECLR void one_arrow(
	const MGPosition& root,	///<root of the arrow
	const MGVector& vecx,	///<the vector from the roo to the head of the arrrow
	const MGUnit_vector& vecy,///<vecy that is normal to the vector from root to head
	MGPosition& head,		///<head of the arrow will be returned.
	MGPosition& headtail1,	///<two tail of arrowhead line segments will be returned.
	MGPosition& headtail2	///<2nd tail.
);

/// @brief  Creates a curve that has weight.
/// @param  curve 曲線オブジェクト
///Returned object is a newed object. User must delete it.
MG_DLL_DECLR MGRLBRep* convert_to_rational(const MGCurve& curve);

};

/** @} */ // end of GEO group

#endif
