/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#ifndef _MGSurfCurve_HH_
#define _MGSurfCurve_HH_

#include "mg/Default.h"
#include "mg/Curve.h"
#include "mg/Surface.h"
#include "mg/TrimmedCurve.h"

//
//Define MGSurfCurve Class.

class MGInterval;
class MGBox;
class MGVector;
class MGUnit_vector;
class MGPosition;
class MGPosition_list;
class MGTransf;
class MGCParam_list;
class MGStraight;
class MGEllipse;
class MGRLBRep;
class MGRSBRep;
class MGCCisects;
class MGCSisects;
class MGIfstream;
class MGOfstream;
/** @file */
/** @addtogroup GEO
 *  @{
 */

///MGSurfCurve is a curve on a surface.

///Defined by a surface and its parameter space line
///represented by (u,v). Let the surface S(u,v), and the parameter space
///curve of 2D be f(t). Then MGSurfCurve is defined as S(f(t)).
///MGSurfCurve is a TEMPORAL curve, and does not have update functions.
class MG_DLL_DECLR MGSurfCurve:public MGCurve{

public:

///Translation.
MG_DLL_DECLR friend MGSurfCurve operator+ (const MGVector& v, const MGSurfCurve& cv2);

///Scaling.
MG_DLL_DECLR friend MGSurfCurve operator* (double scale, const MGSurfCurve& cv2);


////////Special member functions/////////
MGSurfCurve():MGCurve(),m_surface(0){;};
~MGSurfCurve()=default;
MGSurfCurve(const MGSurfCurve&)=default;///Copy constructor.
MGSurfCurve& operator= (const MGSurfCurve&)=default;///Copy assignment.
MGSurfCurve(MGSurfCurve&&)=default;		///Move constructor.
MGSurfCurve& operator= (MGSurfCurve&&)=default;///Move assignment.

///Construct from a surface and a curve.
MGSurfCurve(
	const MGSurface& srf, ///< Surface
	const MGCurve& crv	///< Original uv-curve
); 
MGSurfCurve(
	const MGFSurface& srf, ///< FSurface
	const MGCurve& crv	///< Original uv-curve
); 

////////////Operator overload(演算子多重定義)////////////

public:

///Assignment.
///When the leaf object of this and crv2 are not equal, this assignment
///does nothing.
MGSurfCurve& operator=(const MGGel& gel2);
MGSurfCurve& operator=(MGGel&& gel2);

///Comparison of two curves.
bool operator==(const MGSurfCurve& gel2)const;
bool operator==(const MGGel& gel2)const;
bool operator<(const MGSurfCurve& gel2)const;
bool operator<(const MGGel& gel2)const;

///Approximate this curve as a MGLBRep curve.

///Approximation is done within the tolerance MGTolerance::line_zero().
///When parameter_normalization=0, reparameterization will not done, and
///the evaluation at the same parameter has the same values before and after
///of approximate_as_LBRep.
void approximate_as_LBRep(
	MGLBRep& lb,	///<Approximated obrep will be set.
	int ordr=0,		///<new order. When this is MGLBRep, if ordr=0,
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

///Test if m_curve is MGCompositeCurve.
///If composite, return the pointer. If not, return null.
const MGCompositeCurve* base_composite()const;

/// Return parameter curve's pointer
const MGCurve* base_curve() const {return &m_curve;};

/// Return surface's pointer
const MGSurface* base_surface() const {return m_surface;};

///Returns B-Rep Dimension.
int bdim() const	{return m_curve.bdim();}

///Return minimum box that includes the curve of parameter interval.
/// 入力のパラメータ範囲の曲線部分を囲むボックスを返す。
MGBox box_limitted(
	const MGInterval & ///< Parameter Range of the curve.
) const;

///Return minimum box that includes whole of the curve.
///Compute the box from the scratch.
void compute_box(MGBox& bx) const;

///Changing this object's space dimension.
void change_dimension(
	int sdim,			///< new space dimension
	int start1=0, 		///< Destination order of new object.
	int start2=0		///< Source order of this object.
); 

///Change parameter range, be able to change the direction by providing
///t1 greater than t2.
void change_range(
	double t1,		///<Parameter value for the start of original. 
	double t2	///<Parameter value for the end of original.
){	assert(false);}

///Construct new curve object by copying to newed area.
///User must delete this copied object by "delete".
MGSurfCurve* clone() const;

///copy as a newed curve. The new curve will be MGLBRep or MGRLBRep.
///When original curve was a MGRLBRep, the new curve will be a MGRLBRep.
///Otherwise,  the new curve will be a MGLBRep.
///Returned object must be deleted.
MGCurve* copy_as_nurbs() const override;

///Construct new curve object by changing
///the original object's space dimension.
///User must delete this copied object by "delete".
///This function returned MGLBRep*(order = 4)
MGCurve* copy_change_dimension(
	int sdim,			///< new space dimension
	int start1=0, 		///< Destination order of new line.
	int start2=0 		///< Source order of this line.
)const;

///Divide this curve at the designated knot multiplicity point.
///Function's return value is the number of the curves after divided.
int divide_multi(
	std::vector<UniqueCurve>& crv_list,	//divided curves are appended.
	int multiplicity=-1	///<designates the multiplicity of the knot to divide at,
						///<When multiplicity<=0, order()-1 is assumed,
						///<When multiplicity>=order(), order() is assumed.
)const override;

/// Evaluate n'th derivative data. n=0 means positional data evaluation.
MGVector eval(
	double,				///< Parameter value.
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

/// Return This object's typeID
long identify_type() const;

///Test if input parameter value is inside parameter range of the line.
bool in_range(double t) const;

///Provide divide number of curve span for function intersect.
int intersect_dnum() const;

///Intersection of Curve
MGCCisects isect(const MGCurve& curve2) const override;
MGCCisects isect(const MGStraight& curve2) const override;
MGCCisects isect(const MGLBRep& curve2) const override;
MGCCisects isect(const MGSurfCurve& curve2) const override;


///Intersection with a Surface
MGCSisects isect(const MGSurface& surf) const override;
MGCSisects isect(const MGPlane& surf) const override;

///Access to i-th element of knot.
double knot(int i) const;
	
///Returns the knot vector.
const MGKnotVector& knot_vector() const{return m_curve.knot_vector();}

///Update this by limiting the parameter range of the curve.
/// 自身に指定したパラメータ範囲のｌｉｍｉｔをつける。
void limit(const MGInterval& );

///Returns the order.
int order() const{return m_curve.order();}

/// Return ending parameter value.
double param_e() const;

///Normalize parameter value t to the nearest knot if their distance is
///within tolerance.
double param_normalize(double t) const;

///Return parameter range of the curve(パラメータ範囲を返す)
MGInterval param_range() const;

/// Return starting parameter value.
double param_s() const;

///Compute part of this curve from parameter t1 to t2.
///Returned is the pointer to newed object, and so should be deleted
///by calling program, or memory leaked.
MGSurfCurve* part(
	double t1,///< parameter range from t1,
	double t2,///< to t2.
	int multiple=0	///<Indicates if start and end knot multiplicities
					///<are necessary. =0:unnecessary, !=0:necessary.
) const;

///Compute all the perpendicular points of this curve and the second one.
///That is, if f(s) and g(t) are the points of the two curves f and g,
///then obtains points where the following conditions are satisfied:
///  fs*(f-g)=0.    gt*(g-f)=0.
///Here fs and gt are 1st derivatives at s and t of f and g.
///MGPosition P in the MGPosition_list contains this and crv's parameter
///as:     P(0)=this curve's parameter, P(1)=curve2's parameter value.
MGPosition_list perps(const MGCurve& crv2)const;
MGPosition_list perps(const MGStraight& crv2)const;
MGPosition_list perps(const MGRLBRep& crv2)const;
MGPosition_list perps(const MGEllipse& crv2)const;
MGPosition_list perps(const MGLBRep& crv2)const;
MGPosition_list perps(const MGSurfCurve& crv2)const;
MGPosition_list perps(const MGBSumCurve& crv2)const;


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

///Round t into curve's parameter range.
/// 入力パラメータをパラメータ範囲でまるめて返却する。
double range(double t) const;

///Return space dimension
int sdim() const ;

///Return sweep surface from crv
///Returned is a newed MGSurface, must be deleted.
///The sweep surface is defined as:
///This curve(say c(t)) is the rail and the straight line segments from
///C(t)+start_dist*uvec to C(t)+end_dist*uvec are the generatrix.
MGSurface* sweep(
	const MGUnit_vector& uvec,	///<Sweep Direction.
	double start_dist,	///<distance to start edge.
	double end_dist		///<distance to end edge.
)const;

///Unlimit parameter range of the curve(limitをはずす)
MGCurve& unlimit();

///Unlimit parameter range of the curve to the end point direction
MGCurve& unlimit_end();

///Unlimit parameter range of the curve to the start point direction
MGCurve& unlimit_start();

///Return curve type.
MGCURVE_TYPE type() const;

/// Output function.
std::ostream& toString(std::ostream&) const;

///IGES output function. PD126(approximated as a NURBS line).
///Function's return value is the directory entry id created.
int out_to_IGES(
	MGIgesOfstream& igesfile,
	int SubordinateEntitySwitch=0
)const;

protected:

///Obtain so transformed 1D curve expression of this curve that
///f(t)={sum(xi(t)*g[i]) for i=0(x), 1(y), 2(z)}-g[3], where f(t) is the output
///of oneD and xi(t) is i-th coordinate expression of this curve.
///This is used to compute intersections with a plane g[4].
///For SurfCurve this is not allowed to use.
std::unique_ptr<MGCurve> oneD(
	const double g[4]	///<Plane expression(a,b,c,d) where ax+by+cz=d.
)const{assert(false); std::unique_ptr<MGCurve> tmpPtr; return tmpPtr;}

///メンバデータを読み出す関数
/// 戻り値boolは正常に読み出しが出来ればtrue、失敗すればfalseになる
void ReadMembers(MGIfstream& buf);

///メンバデータを書き込む関数
/// 戻り値boolは正常に書き込みが出来ればtrue、失敗すればfalseになる
void WriteMembers(MGOfstream& buf) const;

///Get the name of the class.
std::string whoami()const{return "SurfCurve";};

protected:
	
///Get data points for approximate_as_LBRep2.
void data_points_for_approximate_as_LBRep2(
	int is,///<approximation parameter range, from knot_vector()[is].
	int ie,///< to [ie].
	MGKnotVector& t,///<New knot configuration will be output.
				///<t's order is input. other information of t will be updated.
	MGNDDArray& tau,///<Data point for t will be output.
	bool neglectMulti=false///<Indicates if multiple knots be kept.
		///< true: multiplicity is removed.
		///< false: multiplicity is kept.
)const;

private:

////////////Member Data//////////

	const MGSurface* m_surface;
	MGTrimmedCurve m_curve;

//////////Following operators are declared to prohibit the use,
//////////since these cannot be provided as member functions.

///Exchange ordering of the coordinates.
///Exchange coordinates (i) and (j).
void coordinate_exchange(int i, int j);

///isect of this SurfCurve whose m_curve is not a MGTrimmedCurve of MGCompositeCurve.
MGCCisects isect_noCompo(const MGCurve& curve2)const;

///isect of this SurfCurve whose m_curve is not a MGTrimmedCurve of MGCompositeCurve.
MGCCisects isect_noCompo(const MGStraight& curve2)const;

///isect of this SurfCurve whose m_curve is not a MGTrimmedCurve of MGCompositeCurve.
MGCCisects isect_noCompo(const MGSurfCurve& curve2)const;

///isect of this SurfCurve whose m_curve is not a MGTrimmedCurve of MGCompositeCurve.
MGCSisects isect_noCompo(const MGSurface& surf)const;

///isect of this SurfCurve whose m_curve is not a MGTrimmedCurve of MGCompositeCurve.
MGCSisects isect_noCompo(const MGPlane& surf)const;

///isect of each elements of this m_curve,
///if m_curve is MGTrimmedCurve of a MGCompositeCurve.
void isect_of_each(
	const MGCurve& curve2,	///<The isect objective curve.
	MGCCisects& list	///<Obtained isect will be appended.
)const;

///isect of each elements of this m_curve,
///if m_curve is MGTrimmedCurve of a MGCompositeCurve.
void isect_of_each(
	const MGSurface& surf,	///<The isect objective surface.
	MGCSisects& list	///<Obtained isect will be appended.
)const;

///Negate the curve direction(曲線の方向を反転する)
void negate();

///Obtain parameter value if this curve is negated by "negate()".
double negate_param(double t)const;
MGPosition negate_param(const MGPosition& t)const;

///perps of this SurfCurve whose m_curve is not a MGTrimmedCurve of MGCompositeCurve.
MGPosition_list perps_noCompo(const MGCurve& curve2)const;

///perpendicular points of each elements of this m_curve,
///if m_curve is MGTrimmedCurve of a MGCompositeCurve.
void perps_of_each(
	const MGCurve& curve2,	///<The perps objective curve.
	MGPosition_list& list	///<Obtained perpendicular points will be appended.
)const;

///Object transformation.
MGSurfCurve& operator+=(const MGVector& v);
MGSurfCurve& operator-=(const MGVector& v);
MGSurfCurve& operator*=(double scale);
MGSurfCurve& operator*=(const MGMatrix& mat);
MGSurfCurve& operator*=(const MGTransf& tr);

///Transformation object construction(NOT ALLOWED TO USE).
MGSurfCurve operator+ (const MGVector& v) const;
MGSurfCurve operator- (const MGVector& v) const;
MGSurfCurve operator* (double scale) const;
MGSurfCurve operator* (const MGMatrix& mat) const;
MGSurfCurve operator* (const MGTransf& tr) const;

friend class MGLBRep;
friend class MGStraight;
friend class MGPlane;
};

/** @} */ // end of GEO group
#endif
