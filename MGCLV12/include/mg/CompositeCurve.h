/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#ifndef _MGCompositeCurve_HH_
#define _MGCompositeCurve_HH_

#include <deque>
#include "mg/Interval.h"
#include "mg/Position_list.h"
#include "mg/Straight.h"
#include "mg/RLBRep.h"
#include "mg/Ellipse.h"
#include "mg/LBRep.h"
#include "mg/SurfCurve.h"
#include "mg/BSumCurve.h"

class MGPlane;
class MGIgesDirectoryEntry;
class MGSurfCurve;
class MGStraight;
class MGCompositeCurve;
class mgVBO;

/** @file */
/** @addtogroup GEO
 *  @{
 */

///MGCompositeCurve is a composite of other leaf curves.

///Assumedly they are connected as C0 continuity. However, MGCompositeCurve
///does not check their continuity, but only put first or last as the user says
/// (in connect_to_end or in connect_to_start).
///Parameter ranges of the member curves are always continuous. param_s() of the
///1st curve to param_e() of the last form MGCompositeCurve's paramter range.
///number_of_curves() indicates the number of leaf curves(in other words,
///element curves to construct the MGCompositeCurve).
class MG_DLL_DECLR MGCompositeCurve:public MGCurve{

public:

///Iterator definition.
typedef std::deque<MGCurve*> container_type;
typedef container_type::iterator iterator;
typedef container_type::const_iterator const_iterator;
typedef container_type::reverse_iterator reverse_iterator;
typedef container_type::const_reverse_iterator const_reverse_iterator;

///Translation by a vector.
MG_DLL_DECLR friend MGCompositeCurve operator+ (const MGVector& v, const MGCompositeCurve& lb);

///Scaling transfromation.
MG_DLL_DECLR friend MGCompositeCurve operator* (double scale, const MGCompositeCurve&);

////////Special member functions/////////
MGCompositeCurve(){	;};
~MGCompositeCurve();
MGCompositeCurve(const MGCompositeCurve&);///Copy constructor.
MGCompositeCurve& operator= (const MGCompositeCurve&);///Copy assignment.
MGCompositeCurve(MGCompositeCurve&&);		///Move constructor.
MGCompositeCurve& operator= (MGCompositeCurve&&);///Move assignment.

///Constructor of one curve.
///crv is a newed object pointer and MGCompositeCurve takes the ownership.
explicit MGCompositeCurve(MGCurve* crv);


//////////// Operator overload(演算子多重定義) ////////////

///Assignment.
///When the leaf object of this and crv2 are not equal, this assignment
///does nothing.
MGCompositeCurve& operator=(const MGGel& gel2);
MGCompositeCurve& operator=(MGGel&& gel2);

///Transformation object construction
MGCompositeCurve operator+ (const MGVector& v) const;
MGCompositeCurve operator- (const MGVector& v) const;
MGCompositeCurve operator* (double scale) const;
MGCompositeCurve operator* (const MGMatrix& mat) const;
MGCompositeCurve operator* (const MGTransf& tr) const;

///Object transformation.
MGCompositeCurve& operator+=(const MGVector& v);
MGCompositeCurve& operator-=(const MGVector& v);
MGCompositeCurve& operator*=(double scale);
MGCompositeCurve& operator*=(const MGMatrix& mat);
MGCompositeCurve& operator*=(const MGTransf& tr);

///Comparison of two curves.
bool is_same_curve(const MGCurve& curve2)const;
bool operator==(const MGCompositeCurve& gel2)const;
bool operator==(const MGGel& gel2)const;
bool operator==(const MGTrimmedCurve& gel2)const;
bool operator<(const MGCompositeCurve& gel2)const;
bool operator<(const MGGel& gel2)const;

//////////// Member Function ////////////

///Approximate this curve as a MGLBRep curve
///within the tolerance MGTolerance::line_zero().
///When parameter_normalization=0, reparameterization will not done, and
///the evaluation at the same parameter has the same values before and after
///of approximate_as_LBRep.
void approximate_as_LBRep(
	MGLBRep& lb,	///<Approximated obrep will be set.
	int ordr=0,		///<new order. When this is MGLBRep, if ordr=0,
					///ordr=order() will be assumed, else ordr=4 is assumed.
	int parameter_normalization=0,	///<Indicates how the parameter normalization be done:
		///<  =0: no parameter normalization.
		///<  =1: normalize to range=(0., 1.);
		///<  =2: normalize to make the average length of the 1st derivative 
		///<    is as equal to 1. as possible.

	bool neglectMulti=false///<Indicates if multiple knots be kept.
		///< true: multiplicity is removed.
		///< false: multiplicity is kept.
)const override;

///Returns B-Rep Dimension.
///bdim of MGCompositeCurve is the sum of each member curves.
int bdim() const;

const_iterator begin()const{return m_composite.begin();};
iterator begin(){return m_composite.begin();};

///Return minimum box that includes the curve of parameter interval.
/// 入力のパラメータ範囲の曲線部分を囲むボックスを返す。
MGBox box_limitted(
	const MGInterval& ///< Parameter Range of the curve.
) const;

///Changing this object's space dimension.
void change_dimension(
	int sdim,			///< new space dimension
	int start1=0, 		///< Destination order of new object.
	int start2=0 		///< Source order of this object.
);

///Change parameter range, be able to change the direction by providing
///t1 greater than t2.
void change_range(
	double t0,	///<Parameter value for the start of original. 
	double t1	///<Parameter value for the end of original. 
);

///Compute the closest point parameter value of this curve from a point.
double closest(const MGPosition& point) const;

///Return minimum box that includes whole of the curve.
///曲線部分を囲むボックスを返す。
///Or'ed box of all the member curves.
void compute_box(MGBox& bx) const;

///関数名：
///		common
///目的：
///		与えられた曲線と自身の交点もしくは共通部分があるかどうか調べる。
///引数：
///		const MGCurve&			crv2,		(I/ )	与えられる曲線
///		std::vector<double>&	vec_param	( /O)	共通部分のパラメータ範囲
///		MGCCisects&			isect		( /O)	交点
///				 4nの配列で、t(4*i+0),t(4*i+1)が自身のパラメータ範囲(t(4*i+0) < t(4*i+1))、
///							 t(4*i+2),t(4*i+3)が与曲線のパラメータ範囲(f(t(4*i+0))=f(t(4*i+2))
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
	const MGCurve& crv2,
	std::vector<double>& vecComSpan,
	MGCCisects& isect
)const;

///関数名：
///		common
///目的：
///		与えられた曲線と自身の共通部分があるかどうか調べる。
///引数：
///		const MGCurve&			crv2,		(I/ )	与えられる曲線
///		std::vector<double>&	vec_param	( /O)	共通部分のパラメータ範囲
///				 4nの配列で、t(4*i+0),t(4*i+1)が自身のパラメータ範囲(t(4*i+0) < t(4*i+1))、
///							 t(4*i+2),t(4*i+3)が与曲線のパラメータ範囲(f(t(4*i+0))=f(t(4*i+2))
///戻り値：
///		共通部分の数:	共通部分が求まった
///		0:				共通部分がなかった
///		-1:				共通エッジの収束計算エラー
///		-2:				共通エッジが４個以上求まった(のっていないと見なす)
///追記：
///	曲線が共通かどうかの誤差にはline_zero()を、パラメータ範囲の収束計算の誤差には、
///  パラメータ範囲*rc_zero()を使用した
int common(
	const MGCurve& crv2,
	std::vector<double>& vecComSpan
)const;

///Connect the input curve to the end(to_end) or start(to_start) of this curve.
///(1) End(start) point of this curve is assumedly the same as the start(end) point
///of add_curve. However, connect_to_xxx does not check the continuity except when
///connecting two curves are both MGLBRep.
///(2) add_curve must be a newed object pointer and connect_to_end takes the ownership.
///(3) add_curve can be a MGCompositeCurve. In this case each curves in add_curve
///become members of this MGCompositeCurve.
///(4) When add_curve is a SurfCurve, or a TrimmedCurve, it is changed to non SurfCurve
///or non TrimmedCurve. Thus the original surface or curve can be deleted or modified.
///(5) When MGCompositeCurve was not an empty curve, the original part of
///the MGCompositeCurve's parameter range will not be changed by this connect_to_end.
///Instead add_curve's parameter range will be so modified that the magnitude of
///1st derivatives of the two curves are the same at the connecting point and
///the modified add_curve's parameter range is continuous to the original.
///(6) connect_to_end(start) will change add_curve's direction if necessary.
///Function's return value is the new parameter range of add_curve after added.
///connect() connects to either start or end of this depending the distance
///of the start or end points of this curve to the start or end points of add_curve.
///connect_to_end() connects either start or end points of add_curve to the end to this curve.
///connect_to_start() connects either start or end points of add_curve to the start to this curve.
MGInterval connect(MGCurve* add_curve);
MGInterval connect_to_end(MGCurve* add_curve);
MGInterval connect_to_start(MGCurve* add_curve);

///Exchange ordering of the coordinates.
///Exchange coordinates (i) and (j).
void coordinate_exchange(int i, int j);

///Construct new curve object by copying to newed area.
///User must delete this copied object by "delete".
MGCompositeCurve* clone() const;

///copy as a newed curve. The new curve will be MGLBRep or MGRLBRep.
///When original curve was a MGRLBRep, the new curve will be a MGRLBRep.
///Otherwise,  the new curve will be a MGLBRep.
///Returned object must be deleted.
MGCurve* copy_as_nurbs() const override;

///Construct new curve object by changing
///the original object's space dimension.
///User must delete this copied object by "delete".
MGCompositeCurve* copy_change_dimension(
	int sdim,			///< new space dimension
	int start1=0, 		///< Destination order of new line.
	int start2=0		///< Source order of this line.
)const;

///Get the i-th curve in this CompositeCurve.
const MGCurve& curve(int i)const{return *(m_composite[i]);};
MGCurve& curve(int i){return *(m_composite[i]);};

///Compute curvilinear integral of the 1st two coordinates from parameter t0
///to t1.
///This integral can be used to compute area sorrounded by the curve.
///The sum of all the member curve's curvilinear_integral.
double curvilinear_integral(double t1, double t2) const;

///Compute curvilinear integral of the 1st two coordinates.
///(All the parameter range of the curve.)
///This integral can be used to compute area sorrounded by the curve.
///The sum of all the member curve's curvilinear_integral.
double curvilinear_integral() const;

void display_break_points(mgSysGL& sgl)const;
void display_control_polygon(mgSysGL& sgl)const;
void display_curvatures(
	mgSysGL& sgl,	///<sgl to make pictures in.
	int		density,///<densitiy of the graph.
	bool	use_radius,///<true:radius display, false:curvature display.
	double	scale=1.	///<scaling of the graph.
)const;

///Divide this curve at the designated knot multiplicity point.
///Function's return value is the number of the curves after divided.
int divide_multi(
	std::vector<UniqueCurve>& crv_list,	//divided curves are appended.
	int multiplicity=-1	///<designates the multiplicity of the knot to divide at.
						///<When multiplicity<=0, order()-1 is assumed.
						///<When multiplicity>=order(), order() is assumed.
)const override;

void drawSE(
	mgVBO& vbo,///<Target graphic object.
	double t0,			///<Start parameter value of the curve.
	double t1			///<End parameter value of the curve.
						///<Draw will be performed from t0 to t1.
)const;

const_iterator end()const{return m_composite.end();};
iterator end(){return m_composite.end();};

/// Evaluate n'th derivative data. n=0 means positional data evaluation.
MGVector eval(
	double t,			///< Parameter value.
	int nderiv=0,	///< Order of Derivative.
	int left=0			///<Left continuous(left=true)
						///<or right continuous(left=false).
)const;

///Extrapolate this curve by an (approximate) chord length.
///The extrapolation is C2 continuous.
void extend(
	double length,	///<approximate chord length to extend. 
	bool start=false///<Flag of which point to extend, start or end point of the line.
					///<If start is true extend on the start point.
);

///Find which curve parameter range input t belongs to.
///Function's return value m is the id of m_composite which satisfy
///m_composite[m]<=t<m_composite[m+1];
///m_composite.size() should be >=1, or find will abort.
int find(double t, bool right_continuous=true) const;

/// Return This object's typeID
long identify_type() const;

///Test if this cure is planar or not.
///MGPlane expression will be output to plane if this is planar.
///Function's return value is true if planar.
bool is_planar(MGPlane& plane)const;

///Intersection of CompositeCurve and another curve crv2.
MGCCisects isect(const MGCurve& crv2)const{return isect_1by1(crv2);};
MGCCisects isect(const MGStraight& curve2)const override{return isect_1by1(curve2);};
MGCCisects isect(const MGLBRep& curve2)const override{return isect_1by1(curve2);};
MGCCisects isect(const MGSurfCurve& curve2)const override{return isect_1by1(curve2);};

///Intersection with a Surface
MGCSisects isect(const MGSurface& surf) const override;
MGCSisects isect(const MGPlane& surf) const override;

///Access to i-th element of knot.
double knot(int i) const;

///Returns the knot vector of the curve.
///This should not be used.
const MGKnotVector& knot_vector() const;

///Cmpute curve length of the interval.
///If t1 is greater than t2, return negative value.
/// パラメータが昇順で与えられたときは正値、降順のときは負値を返す。
///The sum of all the member curve's length.
double length(double t1, double t2) const;
double length()const{return MGCurve::length();}

///Update this by limiting the parameter range of the curve.
/// 自身に指定したパラメータ範囲のｌｉｍｉｔをつける。
void limit(const MGInterval& rng);

///Negate the curve direction(曲線の方向を反転する)
void negate();

///Obtain parameter value if this curve is negated by "negate()".
double negate_param(double t)const;

///Get the number of included curves in this CompositeCurve.
int number_of_curves()const{return int(m_composite.size());};

//****Offset of CompositeCurve is CompositeCurve of the offsets of each
//member curves********

///一定オフセット関数
///オフセット方向は、法線方向から見て入力曲線の進行方向左側を正とする。
///法線ベクトルがヌルの場合、始点において曲率中心方向を正とする。ただし、曲率中心へ曲率半径以上のオフセット
///は行わない。トレランスはline_zero()を使用している。戻り値は、オフセット曲線リストが返却される。
///costant offset curve. if the norm_vector is given, the positive offset direction decide
///to left hand side from ahead, or the direction to center of curvature at start parameter.
///the offset value is less than radius of curvature. line_zero() is used.
std::vector<UniqueCurve> offset(
	double ofs_value,							///<オフセット量
	const MGVector& norm_vector = mgNULL_VEC	///<法線ベクトル
) const;

///可変オフセット関数
///オフセット量は空間次元1の線B表現で与えられる。
///オフセット方向は、法線方向から見て入力曲線の進行方向左側を正とする。
///法線ベクトルがヌルの場合、始点において曲率中心方向を正とする。ただし、曲率中心へ曲率半径以上のオフセット
///は行わない。トレランスはline_zero()を使用している。戻り値は、オフセット曲線リストが返却される。
///valuable offset curve. if the norm_vector is given, the positive offset direction decide
///to left hand side from ahead, or the direction to center of curvature at start parameter.
///the offset value is less than radius of curvature. line_zero() is used.
std::vector<UniqueCurve> offset(
	const MGLBRep& ofs_value_lb,					///<空間次元１の線B表現で示したオフセット量
	const MGVector& norm_vector = mgNULL_VEC	///<法線ベクトル
)const;

///C2連続曲線の一定オフセット関数
///オフセット方向は、法線方向から見て入力曲線の進行方向左側を正とする。
///法線ベクトルがヌルの場合、始点において曲率中心方向を正とする。ただし、曲率中心へ曲率半径以上のオフセット
///は行わない。トレランスはline_zero()を使用している。戻り値は、オフセット曲線が返却される。
///costant offset curve of C2 continuous curve. if the norm_vector is given, the positive offset direction
///decide to left hand side from ahead, or the direction to center of curvature at start parameter.
///the offset value is less than radius of curvature. line_zero() is used.
MGLBRep offset_c2(
	double ofs_value,								///<オフセット量
	const MGVector& norm_vector = mgNULL_VEC	///<法線ベクトル
) const;

///C2連続曲線の可変オフセット関数
///オフセット量は空間次元1の線B表現で与えられる。
///オフセット方向は、法線方向から見て入力曲線の進行方向左側を正とする。
///法線ベクトルがヌルの場合、始点において曲率中心方向を正とする。ただし、曲率中心へ曲率半径以上のオフセット
///は行わない。トレランスはline_zero()を使用している。戻り値は、オフセット曲線が返却される。
///valuable offset curveof C2 continuous curve. if the norm_vector is given, the positive offset direction
///decide to left hand side from ahead, or the direction to center of curvature at start parameter.
///the offset value is less than radius of curvature. line_zero() is used.
MGLBRep offset_c2(
	const MGLBRep& ofs_value_lb,					///<空間次元１の線B表現で示したオフセット量
	const MGVector& norm_vector = mgNULL_VEC	///<法線ベクトル
) const;

///Test if given point is on the curve or not. If yes, return parameter
///value of the curve. Even if not, return nearest point's parameter.
/// 指定点が自身上にあるかを調べる。曲線上にあれば，そのパラメーター値を，
/// なくても最近傍点のパラメータ値を返す。
/// Function's return value is >0 if the point is on the curve,
/// and 0 if the point is not on the curve.
bool on(
	const MGPosition&,	///<Point(指定点)
	double&				///<Parameter of the curve(パラメータ)
)const;

///Returns the order.
///Returns the maximum order among the curves.
int order() const;

/// Return ending parameter value.
double param_e() const;

///Normalize parameter value t to the nearest knot if their distance is
///within tolerance.
double param_normalize(double t) const;

/// Return starting parameter value.
double param_s() const;
	
///Compute part of this curve from parameter t0 to t1.
///Returned is the pointer to newed object, and so should be deleted
///by calling program, or memory leaked.
MGCurve* part(
	double t0,///< Start paramter.
	double t1,///<End parameter.
	int multiple=0	///<Indicates if start and end knot multiplicities
					///<are necessary. =0:unnecessary, !=0:necessary.
) const;

///Return perpendicular point from a point P,
///given guess starting paramter values.
///Function's return value is:
///   perp_guess=true if perpendicular points obtained,
///   perp_guess=false if perpendicular points not obtained,
int perp_guess(
	double t0,///< parameter range of this.
	double t1,///< (t0>=t1) indicates no range specified.
	const MGPosition& P,	///<Point(指定点)
	double tg,				///<Guess parameter values of this curve.
	double& t				///<Output parameter
) const;

///Return perpendicular points of two curves,
///given guess starting paramter values.
///Function's return value is:
///   perp_guess=true if perpendicular points obtained,
///   perp_guess=false if perpendicular points not obtained,
int perp_guess(
	double s0,	///<parameter range of this.
	double s1,	///<When s0>=s1, no limit for this parameter range.
	const MGCurve& curve2,		///<2nd curve.
	double t0,	///<parameter range of curve2.
	double t1,	///<When t0>=t1, no limit for curve2 parameter range.
	double sg,	///<Guess parameter values of this curve.
	double tg,	///<Guess parameter values of curve2's.
	MGPosition& st	///<perpendicular points' parameter values
					///<will be output.
	///<st(0): this curve's parameter, st(1):curve2's parameter.
) const;

///Compute all the perpendicular points of this curve and the second one.
///That is, if f(s) and g(t) are the points of the two curves f and g,
///then obtains points where the following conditions are satisfied:
///  fs*(f-g)=0.    gt*(g-f)=0.
///Here fs and gt are 1st derivatives at s and t of f and g.
///MGPosition P in the MGPosition_list contains this and crv's parameter
///as:     P(0)=this curve's parameter, P(1)=crv2's parameter value.
MGPosition_list perps(const MGCurve& crv2)const{return perps_1by1(crv2);};
MGPosition_list perps(const MGStraight& crv2)const{return perps_1by1(crv2);};
MGPosition_list perps(const MGRLBRep& crv2)const{return perps_1by1(crv2);};
MGPosition_list perps(const MGEllipse& crv2)const{return perps_1by1(crv2);};
MGPosition_list perps(const MGLBRep& crv2)const{return perps_1by1(crv2);};
MGPosition_list perps(const MGSurfCurve& crv2)const;
MGPosition_list perps(const MGBSumCurve& crv2)const{return perps_1by1(crv2);};

///Compute all foot points of the perpendicular line from point to
///the curve.
/// 与ポイントから曲線へ下ろした垂線の足の，曲線のパラメータ値を
/// すべて求める。
MGCParam_list perps(
	const MGPosition& P		///Point(指定点)
) const;

///Approximate this curve by a polyline and output to lb2.
///The tolerance of the approximation is error.
void polygonize(
	double error,	///<tolerance allowed for the approximation
	MGLBRep& lb2	///<Obtained polyline will be output as an MGLBRep of order2.
)const;

///曲線を面に面直またはベクトル投影して曲線リストを求める。
///投影曲線は面上のパラメータ曲線と3次元曲線としてそれぞれ順番に、
///vec_crv_uv, vec_crvに格納される。
///uv曲線のトレランスはrc_zero()を、3次元曲線はline_zero()をそれぞれ使用している。
///get perpendicular or vector projection curve list.
///uv projection curves are put into vec_crv_uv(rc_zero() is used),
///3d projection curves are put into vec_crv(line_zero() is used) respectively.
///戻り値：
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
int project(
	const MGFSurface& surf,	//given surface.
	std::vector<UniqueCurve>& vec_crv_uv,	//uv projection curve will be appended.
	std::vector<UniqueCurve>& vec_crv,	//3d projection curve will be appended.
	const MGVector& vec	//projection vector.
						//if vec = NULL then calculate perpendicular project.
)const;

///Obtain the projected curve of this onto the surface srf.
///The direction of the projection is along the vector vec if the vec is not
///NULL, and normal to the surface if the vec is NULL.
///Output of 'project' is two kind of curves:
///one is general world coordinate curves('vec_crv'), and the other is (u,v) curves of
///the parameter space of the surfaces(vec_crv_uv).
///vec_crv_uv.size() is equal to vec_crv.size(). Let the size be n, then
/// (vec_crv_uv[i], vec_crv[i]) is one pair for 0<=i<n.
size_t project_onto_surface(
	const MGFSurface& srf,///<The target surface.
	std::vector<UniqueCurve>& vec_crv_uv,
		///<Projected curve(surface parameter (u,v) representation) will be appended.
	std::vector<UniqueCurve>& vec_crv,
		///<Projected curve(world coordinate(x,y,z) representation) will be appended.
	const MGVector& vec///<The direction of the projection.
)const;

///Release the pointer of the last curve.
///Returned will be the released MGCurve pointer.
MGCurve* release_back();

///Release the pointer of the 1st curve.
///Returned will be the released MGCurve pointer.
MGCurve* release_front();

///ノット削除関数(B表現曲線のみ)
///トレランスはline_zeroを使用する。元のノットが細かいものほど削除しやすい
///Remove redundant knot, and reduce the b-rep dimension.
///The tolerance used is MGTolerance::line_zero().
void remove_knot();

///Return the space dimension. It is the maximum space dimension of the
///member curves.
int sdim() const;

///Return sweep surface from crv
///Returned is a newed MGSurface, must be deleted.
///The sweep surface is defined as:
///This curve(say c(t)) is the rail and the straight line segments from
///C(t)+start_dist*uvec to C(t)+end_dist*uvec are the generatrix.
MGSurface* sweep(
	const MGUnit_vector& uvec,	///<Sweep Direction.
	double start_dist,			///<distance to start edge.
	double end_dist				///<distance to end edge.
) const;	

///Return curve type.
MGCURVE_TYPE type() const{return MGCURVE_TYPE::MGCURVE_COMPOSITE;};

///Unlimit parameter range of the curve(limitをはずす)
MGCurve& unlimit();

///Unlimit parameter range of the curve to the end point direction
MGCurve& unlimit_end();

///Unlimit parameter range of the curve to the start point direction
MGCurve& unlimit_start();

///Output function.
std::ostream& toString(std::ostream&) const;

///Output to IGES stream file.
int out_to_IGES(
	MGIgesOfstream& igesfile,
	int SubordinateEntitySwitch=0
)const;

protected:

///Compute intersection point of 1D sub curve of original curve.
///Parameter values of intersection point will be returned.
MGCParam_list intersect_1D(						
	double f,			///< Coordinate value
	int coordinate=0	///< Coordinate kind of the data f(from 0).
) const;	

///Obtain so transformed 1D curve expression of this curve that
///f(t)={sum(xi(t)*g[i]) for i=0(x), 1(y), 2(z)}-g[3], where f(t) is the output
///of oneD and xi(t) is i-th coordinate expression of this curve.
///This is used to compute intersections with a plane g[4].
std::unique_ptr<MGCurve> oneD(
	const double g[4]			///<Plane expression(a,b,c,d) where ax+by+cz=d.
) const;

///メンバデータを読み出す関数
/// 戻り値boolは正常に読み出しが出来ればtrue、失敗すればfalseになる
void ReadMembers(MGIfstream& buf);

///メンバデータを書き込む関数
/// 戻り値boolは正常に書き込みが出来ればtrue、失敗すればfalseになる
void WriteMembers(MGOfstream& buf) const;

///Get the name of the class.
std::string whoami()const{return "CompositeCurve";};

private:

	container_type m_composite;	///<m_composite consists of newed MGCurve pointer.

///Provide divide number of curve span for function intersect.
///****This is declared here to prohibit the use****.
int intersect_dnum() const;

///Intersection of CompositeCurve and another curve crv2.
MGCCisects isect_1by1(const MGCurve& crv2) const;

///Intersection with a surface. Obtains the isects by computing isects
///of each elemet curve.
MGCSisects isect_1by1(const MGSurface& surf) const;

///Compute all the perpendicular points of this curve and the second one.
///That is, if f(s) and g(t) are the points of the two curves f and g,
///then obtains points where the following conditions are satisfied:
///  fs*(f-g)=0.    gt*(g-f)=0.
///Here fs and gt are 1st derivatives at s and t of f and g.
///MGPosition P in the MGPosition_list contains this and crv's parameter
///as:     P(0)=this curve's parameter, P(1)=crv2's parameter value.
MGPosition_list perps_1by1(
	const MGCurve& crv2		///<The second curve
)const;

friend class MGHHisect;

};

/** @} */ // end of GEO group
#endif
