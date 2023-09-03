/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#ifndef _MGStraight_HH_
#define _MGStraight_HH_

#include "mg/EReal.h"
#include "mg/Position.h"
#include "mg/Unit_vector.h"
#include "mg/Curve.h"

// MGStraight.h
// Definition of MGStraight class

class MGBox;
class MGInterval;
class MGVector;
class MGTransf;
class MGCCisect;
class MGCSisect;
class MGRLBRep;
class MGEllipse;
class MGLBRep;
class MGSurfCurve;
class MGBSumCurve;
class MGCompositeCurve;
class MGTrimmedCurve;
class MGKnotVector;
class MGPlane;
class MGSurface;
class MGCCisects;
class MGPosition_list;
class MGIfstream;
class MGOfstream;
/** @file */
/** @addtogroup GEO
 *  @{
 */

/// MGStraight is a curve of any space dimension, represent a straight line.

/// Parameterization of the MGStraight is as below using parameter t:
/// point(t) = m_root_point + t*m_direction.
///MGStraight can be a line segment that has start and end points.
///Let t be the parameter of the straight line(i.e. length from m_root_point),
///then straight lien f(t) can be expressed as:
///  f(t)=m_root_point+m_direction*t.
class MG_DLL_DECLR MGStraight: public MGCurve{
friend class MGCylinder;

////////////Member data. メンバデータ//////////
/// Parameterization of the MGStraight is as below using parameter t:
/// point(t) = m_root_point + t*m_direction.

MGPosition	    m_root_point;///<Root point. 基点
MGVector		m_direction;///<Direction vector. 直線の方向ベクトル
MGEReal			m_sparam;	///<Start point's parameter
MGEReal		    m_endparam;	///<End point's parameter value. 
mutable MGKnotVector* m_knotV;///<When knot_vector() is invoked, the knot vector is set.

public:

///Translation.
MG_DLL_DECLR friend MGStraight operator+ (const MGVector& v, const MGStraight& sl);

///Scaling.
MG_DLL_DECLR friend MGStraight operator* (double scale, const MGStraight&);


////////Special member functions/////////
MGStraight();	///void constructor.
~MGStraight(){if(m_knotV) delete m_knotV;};	///Destructor.
MGStraight(const MGStraight& sl);///Copy constructor.
MGStraight& operator=(const MGStraight& sl);///assignment.
MGStraight(MGStraight&& sl);		///Move constructor.
MGStraight& operator= (MGStraight&& rhs);///Move assignment.


///Straight specifying all the member data. All of the data are employed as 
///member data of this straight.
explicit MGStraight (
	const MGEReal& endparam,	///<end parameter value
	const MGEReal& sparam,		///<start parameter value
	const MGVector& direction,	///<Direction, which will be the direction of this.
	const MGPosition& Origin	///<Origin
);
explicit MGStraight(
	double endparam,			///<end parameter value
	double sparam,				///<start parameter value
	const MGVector& direction,	///<Direction, which will be the direction of this.
	const MGPosition& Origin	///<Origin
);

///Straight from straight line type, direction vector, and an origin.

///This constrcutor converts input vec to a unit vector.
///If you do not like the conversion, use set_straight().
MGStraight(
	MGSTRAIGHT_TYPE type,	///<Type
		///<Straight line type(直線の種類)
		///<enum MGSTRAIGHT_TYPE {
		///<	MGSTRAIGHT_EMPTY 		//Empty. 空
		///<	,MGSTRAIGHT_SEGMENT		//Line segment. 線分
		///<	,MGSTRAIGHT_HALF_LIMIT	//Half unlimit. 半直線
		///<	,MGSTRAIGHT_UNLIMIT		//Unlimit line for both direction. 無限直線
		///<};
	const MGVector& vec=mgX_UVEC,	///<Direction
	const MGPosition & = mgORIGIN	///<Origin
);

///MGSTRAIGHT_SEGMENT straight from two points.

///Start point is start and end point is end.
///Parameter value of the start point is set to be 0.
MGStraight(
	const MGPosition& end,		///<End point.
	const MGPosition& start		///<Start point.
);

///MGSTRAIGHT_SEGMENT straight from two points.
///Start point is start and end point is end.
///In this version, can specify start and end parameter values.
MGStraight(
	const MGPosition& endP,		///<End point.
	const MGPosition& startP,	///<Start point.
	const double endT,			///<End Parameter.
	const double startT = 0.0	///<Start Parameter.
);

///Straight line from direction vector, end point parameter, and an origin.
MGStraight(
	const MGUnit_vector& v,	///<Unit direction vector
	double d,				///<Parameter value of end point
	const MGPosition& p= mgORIGIN///<Origin
);
	
///Construct the infinite straight line that is a perpendicular bisect
///of the two point P1 and P2.

///Constructed straight is normal to the vector N.
///The line's direction is N*(P2-P1).
///N is the normal of the plane that P1, P2, and the constructed line lie on.
explicit MGStraight(
	const MGPosition& P1,	///<point 1.
	const MGPosition& P2,	///<point 2.
	const MGVector& N///<Normal.
);

///Construct the unlimitted straight that pass through the point uv,
///and the direction is the middle vector of (-v0, v1).

///All of v0, v1, uv, and this straight are objects of space dimension 2.
explicit MGStraight(
	const MGVector& v0,	//a vector whose end is uv.
	const MGVector& v1,	//a vector whose start is uv.
	const MGPosition& uv//origin of sl.
);


///Construct Straight Line copying original line. Able to change
///space dimension and ordering of axis.
MGStraight(
	int dim,				///<New space dimension.
	const MGStraight& linne2,///<Original line.
	int start1=0, 			///<Destination order of new line.
	int start2=0 			///<Source order of original line.
);

///Assignment.
///When the leaf object of this and crv2 are not equal, this assignment
///does nothing.
MGStraight& operator=(const MGGel& gel2);
MGStraight& operator=(MGGel&& gel2);

///Transformation
MGStraight operator+ (const MGVector& vec)const;
MGStraight operator- (const MGVector& vec)const;
MGStraight operator* (double scale)const;
MGStraight operator* (const MGMatrix&)const;
MGStraight operator* (const MGTransf&)const;

///Object transformation.
MGStraight& operator+=(const MGVector& v);
MGStraight& operator-=(const MGVector& v);
MGStraight& operator*=(double scale);
MGStraight& operator*=(const MGMatrix& mat);
MGStraight& operator*=(const MGTransf& tr);

///comparison
bool operator==(const MGStraight& sl2)const;
std::partial_ordering operator<=>(const MGStraight& gel2)const;

//gel2 must be the same class as this.
bool equal_test(const MGGel& gel2)const;

//gel2 must be the same class as this.
std::partial_ordering ordering_test(const MGGel& gel2)const;

///Approximate this curve as a MGLBRep curve.

///Approximated within the tolerance MGTolerance::line_zero().
///When parameter_normalization=0, reparameterization will not done, and
///the evaluation at the same parameter has the same values before and after
///of approximate_as_LBRep.
void approximate_as_LBRep(
	MGLBRep& lb,	///<Approximated lbrep will be set.
	int ordr=0,		///<new order. When this is MGLBRep, if ordr=0,
					///<ordr=order() will be assumed, else ordr=4 is assumed.
					///<For MGStraight, ordr=2 is always employed.
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
int bdim()const{return 2;};

///Minimum box that includes the line limitted by interval l.
MGBox box_limitted(const MGInterval& l)const;

///Changing this object's space dimension.
void change_dimension(
	int sdim,			///< new space dimension
	int start1=0, 		///< Destination order of new object.
	int start2=0 		///< Source order of this object.
);

///Change parameter range, be able to change the direction by providing
///t1 greater than t2.
void change_range(
	double t1,		///<Parameter value for the start of the original. 
	double t2		///<Parameter value for the end of the original. 
);

///Compute the box of the whole of the curve.
///Compute the box from the scratch.
void compute_box(MGBox& bx) const;

///Exchange ordering of the coordinates.
///Exchange coordinates (i) and (j).
void coordinate_exchange(int i, int j);

///Construct new curve object by copying to newed area.

///User must delete this copied object by "delete".
MGStraight* clone()const;

///Compute the closest point parameter value of this curve from a point.
double closest(const MGPosition& point) const;

///Compute the closest point parameter value pair of this curve and curve2.

///MGPosition P of the function return contains this and curve2's parameter
///as:     P(0)=this curve's parameter, P(1)=curve2's parameter value.
MGPosition closest(const MGCurve& curve2) const;

///Convert this curve to Bezier curve.

///If this is MGLBRep or MGStraight, the shape is exactly the same
///as the original. Otherwise, this is apporoximated by MGLBRep.
void convert_to_Bezier(MGLBRep& bezier)const;

///Copy as a newed curve.

///The new curve will be MGLBRep.
///Returned object must be deleted.
MGCurve* copy_as_nurbs()const override;

///Construct new curve object by changing the original object's space dimension.

///User must delete this copied object by "delete".
MGStraight* copy_change_dimension(
	int sdim,			///< new space dimension
	int start1=0, 		///< Destination order of new line.
	int start2=0 		///< Source order of this line.
)const;

///Return curvature, i.e. return 0.
double curvature( double ) const{return 0.;};

///Compute curvilinear integral of the 1st two coordinates.
///This integral can be used to compute area sorounded by the curve.
double curvilinear_integral(double t1, double t2)const;

///Return direction vector of the line at a parameter.
MGUnit_vector direction(double)const{return m_direction;};

///Return direction vector of the line.
const MGVector& direction()const{return m_direction;};

///Return direction vector length of the straight line.
double direction_len()const{return m_direction.len();};

///Return the ditance between a position.
double distance( const MGPosition& = mgORIGIN )const;

///Return distance between two striaght lines.
double distance( const MGStraight& )const;

///Draw as a wire.
void drawWire(
	mgVBO& vbo,///<The target vbo element.
	int line_density=1///<line density to draw surface in wire mode.
)const;

///Draw using vbo.
void drawSE(
	mgVBO& vbo,///<The target vbo element.
	double t0,			///<Start parameter value of the curve.
	double t1			///<End parameter value of the curve,
						///<Draw will be performed from t0 to t1.
)const;

///Return end parameter value.
const MGEReal& end()const{ return m_endparam;};

///Return end point parameter, valid only when type=MGSTRAIGHT_SEGMENT.
double end_param()const{return param_e();};

///Return end point coordinate, valid only when type=MGSTRAIGHT_SEGMENT.
MGPosition end_point()const;

/// Evaluate n'th derivative data. nderiv=0 means
/// positional data evaluation.
MGVector eval(
	double t,			///< Parameter value.
	int nderiv=0,	///< Order of Derivative.
	int left=0			///<Left continuous(left=true)
						///<or right continuous(left=false).
)const;

///Evaluate positional data, 1st derivative, and 2nd derivative
///at parameter t. 
	void eval_all (
	double t,		///<Parameter of the straight.
	MGPosition&,	///<Positional data will be returned.
	MGVector&,		///<1st derivative will be returned.
	MGVector&		///<2nd derivative will be returned.
) const;

///Evaluate 1st derivative at parameter t.
MGVector eval_deriv(double t)const;

///Evaluate positional data at parameter t.

///Evaluation is done after t is limitted within the parameter range.
MGPosition eval_position(double t)const;

///Evaluate positional data at parameter t.

///Evaluation is done without limitting the parameter value t within the parameter range.
MGPosition eval_position_unlimitted(double t)const;

///Extrapolate this curve by an (approximate) chord length.

///The extrapolation is C2 continuous.
void extend(
	double length,	///<approximate chord length to extend. 
	bool start=false///<Flag of which point to extend, start or end point of the line,
					///<If start is true extend on the start point.
);

///Test if thie straight is a finite line.
bool finite()const{return m_sparam.finite() && m_endparam.finite();};

/// Return This object's typeID
long identify_type()const;

///Test if infinite_above(infinite to larger parameter direction).
bool infinite_above()const{return m_endparam.plus_infinite();};

///Test if infinite_below(infinite to smaller parameter direction).
bool infinite_below()const{return m_sparam.minus_infinite();};

///Test if input parameter value is inside parameter range of the line.
bool in_range(double t)const;

///Provide divide number of curve span for function intersect.
int intersect_dnum() const{return 2;}

///Test if this cure is co-planar with the 2nd curve curve2.

///MGPlane expression will be out to plane if this is co-planar.
///Function's return value is true if co-planar.
bool is_coplanar(const MGCurve& curve2, MGPlane& plane)const;

///Test if the input parameter t is the start point parameter or not.
bool is_startpoint_parameter(double t)const;

///Test if the input parameter t is the start point parameter or not.
bool is_endpoint_parameter(double t)const;

///Test if this cure is linear or not, that is, is straight or not.

///MGStraight expression will be out to straight if this is linear or not.
///Function's return value is true if linear.
bool is_linear(MGStraight& straight)const;

///Test if this cure is planar or not.

///MGPlane expression will be out to plane if this is planar.
///Function's return value is true if planar.
bool is_planar(MGPlane& plane)const;

///Intersection of straight and a curve.
MGCCisects isect(const MGCurve&)const;
MGCCisects isect(const MGStraight& curve2)const;
MGCCisects isect(const MGRLBRep& curve2)const;
MGCCisects isect(const MGEllipse& curve2)const;
MGCCisects isect(const MGBSumCurve& curve2)const;

MGCSisects isect(const MGSurface& srf)const;
MGCSisects isect(const MGPlane& srf)const;
MGCSisects isect(const MGFace& f)const;

///Access to i-th element of the knot of this.

///i=0, 1 and returns start or end parameter value of the straight.
double knot(int i) const;

///Returns the knot vector of the curve.
const MGKnotVector& knot_vector() const;
MGKnotVector& knot_vector();

///Return curve length from parameter t1 to t2. If t1>t2, the result is
///minus value.
double length(double t1, double t2) const;

///Return the whole curve length.

///When type=MGSTRAIGHT_SEGMENT, return the whole length of straight
///line segment. Else, return -1.
double length() const;

///Return the parameter of the line away from point t by length len
///along the straight.
double length_param(double t, double len) const;

///Compute sub straight line limitted by an interval.
void limit(const MGInterval& );

///Compute sub straight line limitted by an box.
///This box's coordinates consist of world coordinates.
void limit(const MGBox& box);

///Compute nearest point on the line to the origin.
MGPosition nearest_to_origin() const;

///Negate the direction of the curve.
void negate();

///Obtain parameter value if this curve is negated by "negate()".
double negate_param(double t)const;

///一定オフセット関数

///オフセット方向は、法線方向から見て入力曲線の進行方向左側を正とする。
///法線ベクトルがヌルの場合、空間次元のエレメント番号が最も大きいエレメントが1の単位ベクトルを正とする。
///トレランスはline_zero()を使用する。戻り値は、オフセット曲線リストが返却される。
///costant offset curve. if the norm_vector is given, the positive offset direction decide
///to left hand side from ahead, or the MGUnit_vector() direction.
///line_zero() is used. return value is number of offset curve.
std::vector<UniqueCurve> offset(
	double ofs_value,						///<オフセット量
	const MGVector& norm_vector = mgNULL_VEC///<法線ベクトル
)const override;

///Test if input point is on the line or not.

///Even if the point is not on the line, return nearest point's parameter value.
///Function's return value is true if input point is on the line,
/// and  false if the point is not on the line.
bool on(
	const MGPosition&,	///<Input a point. 指定点
	double&				///<Parameter value of the curve will be returned,
						///<パラメータ値
)const;

///Test if the straight is on a plane or not.

///Function's return value is true if on the line,
/// and  false if not on the line.
bool on(
	const MGPlane&		///< Plane
)const;

///Returns the order.
int order() const{return 2;};

///Return parameter of the straight of input point.

/// If input point is not on the curve, return the nearest point's
///parameter on the curve.
double param (
	const MGPosition &
) const;

/// Return ending parameter value.
double param_e() const{return m_endparam.value();};

///Obtain parameter space error.
double param_error() const;

///Normalize parameter value t to the nearest knot if their distance is
///within tolerance.
///For straight, the knots are start and end points.
double param_normalize(double t) const;

///Return parameter range of the curve.
MGInterval param_range() const;

/// Return starting parameter value.
double param_s() const{return m_sparam.value();};
	
///Compute part of this curve from parameter t1 to t2.

///Returned is the pointer to newed object, and so should be deleted
///by calling program, or memory leaked.
MGStraight* part(
	double t1, ///<from parameter.
	double t2, ///<To parameter.
	int multiple=0	///<Indicates if start and end knot multiplicities
					///<are necessary. =0:unnecessary, !=0:necessary.
) const;

///Return the foot of the straight line that is perpendicular from a point P.

///Function's return value is parameter value of this straight line,
///may not be in_range.
double perp_param(
	const MGPosition&	P ///< 与点
) const;

///Return the foot of the straight line that is perpendicular from a point P.
///Function's return value is if point is obtained(1) or not(0).
int perp_point(
	const MGPosition& P,	///<point
	double& d1,				///<parameter value will be returned.
	const double* d2=NULL	///<guess parameter value, dummy, not used.
) const;
	
///Return all the foots of the  straight lines that is perpendicular
///to the line. Actually only one point for straight.
MGCParam_list perps(
	const MGPosition& point
)const;

///Compute all the perpendicular points of this curve and the second one.

///That is, if f(s) and g(t) are the points of the two curves f and g,
///then obtains points where the following conditions are satisfied:
///  fs*(f-g)=0.    gt*(g-f)=0.
///Here fs and gt are 1st derivatives at s and t of f and g.
///MGPosition P in the MGPosition_list contains this and crv's parameter
///as:     P(0)=this curve's parameter, P(1)=crv2's parameter value.
MGPosition_list perps(const MGCurve& crv2)const;
MGPosition_list perps(const MGStraight& crv2)const;
MGPosition_list perps(const MGRLBRep& crv2)const;
MGPosition_list perps(const MGEllipse& crv2)const;
MGPosition_list perps(const MGSurfCurve& crv2)const;
MGPosition_list perps(const MGBSumCurve& crv2)const;

///Obtain the projected curve of this onto the surface.

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
)const override;

///Round the input parameter value t into the parameter range of the line.
double range(double) const;

///Return two straight line's relationship.

///Whe two are parallel, MGPSREL_PARALLEL or MGPSREL_COIN will be returned.
///MGPSREL_PARALLEL takes place when no overlapping part exists, and 
///MGPSREL_COIN when some parts overlaps.
///MGPSREL_ISECT when they have an intersection.
///MGPSREL_VIRTUAL_ISECT when there is an intersection at extended part of the lines.
///MGPSREL_TORSION when the two do not lie on a same plane.
MGPSRELATION relation(
	const MGStraight& sl2,///<2nd straight line.
	MGCCisect& ip
		///<When MGPSREL_PARALLEL or MGPSREL_COIN, a pair of parameter values of the nearest
		///<point of the two will be returned.
		///<When MGPSREL_ISECT or MGPSREL_VIRTUAL_ISECT, the intersection point parameter values
		///<will be returned.
		///<When MGPSREL_TORSION, MGPSREL_VIRTUAL_ISECT intersection point parameter value after
		///<transformed to lie on the same plane will be returned.
) const;
	
///Return relationship with a plane.
MGPSRELATION relation(
	const MGPlane &,	///<Plane.
	MGCSisect &		///<Intersection point will be returned if exist.
) const;

const MGPosition& root_point() const{return m_root_point;};

///Return space dimension
int sdim() const;

///Straight from straight line type, direction vector, and an origin.
///Construct a straight and replce this with it.
///This fuction does not convert input vec to a unit vector.
///If you like the conversion, use MGStraight() constructor.
MGStraight& set_straight(
	MGSTRAIGHT_TYPE type,	///<Type
		///<Straight line type(直線の種類)
		///<enum MGSTRAIGHT_TYPE {
		///<	MGSTRAIGHT_EMPTY 		//Empty. 空
		///<	,MGSTRAIGHT_SEGMENT		//Line segment. 線分
		///<	,MGSTRAIGHT_HALF_LIMIT	//Half unlimit. 半直線
		///<	,MGSTRAIGHT_UNLIMIT		//Unlimit line for both direction. 無限直線
		///<};
	const MGVector& vec= mgX_UVEC,				///<Direction
	const MGPosition & = mgORIGIN			///<Origin
);

///Return straight line's direction.
const MGVector& sl_direction()const{return m_direction;}

///Return start point parameter value.
const MGEReal& start() const{ return m_sparam;};

///Return start(root) point of the straight.
MGPosition start_point() const;

///Return the straight line type.

///Straight line type
///enum MGSTRAIGHT_TYPE {
///	MGSTRAIGHT_EMPTY 		//Empty.
///	,MGSTRAIGHT_SEGMENT		//Line segment.
///	,MGSTRAIGHT_HALF_LIMIT	//Half unlimit.
///	,MGSTRAIGHT_UNLIMIT		//Unlimit line for both direction.
///};
MGSTRAIGHT_TYPE straight_type() const;

///Return sweep surface of this.

///Returned is a newed MGSurface, must be deleted.
///The sweep surface is defined as:
///This curve(say c(t)) is the rail and the straight line segments from
///C(t)+start_dist*uvec to C(t)+end_dist*uvec are the generatrix.
MGSurface* sweep(
	const MGUnit_vector& uvec,	///<Sweep Direction.
	double start_dist,			///<distance to start edge.
	double end_dist				///<distance to end edge.
)const;

///Return curve type, i.e. MGCURVE_STRAIGHT.
MGCURVE_TYPE type() const{return MGCURVE_TYPE::MGCURVE_STRAIGHT;};

///Unlimit the parameter range, i.e. change to infinite striahgt line
///for both direction.
MGCurve& unlimit();

///Unlimit parameter range of the curve to the end point direction.
MGCurve& unlimit_end();

///Unlimit parameter range of the curve to the start point direction
MGCurve& unlimit_start();

//Update the root point of this straight.
void update_root(const MGPosition& rootP);

/// String Output function.
std::ostream& toString(std::ostream &) const;

///Iges output.
int out_to_IGES(
	MGIgesOfstream& igesfile,
	int SubordinateEntitySwitch=0
)const;

///Get the name of the class.
std::string whoami()const{return "Straight";};

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
void ReadMembers(MGIfstream& buf);
	
///メンバデータを書き込む関数
void WriteMembers(MGOfstream& buf) const;

private:

///Compute the closest point parameter value pair of this MGStraight and straight2.
///MGPosition P of the function return contains this and straight2's parameter
///as:     P(0)=this MGStraight's parameter, P(1)=straight2's parameter value.
MGPosition closestSL(const MGStraight& straight2)const;

///isect2D returns parameter values of this(t) and l2(s)
/// of the intersection point of both 2D straight lines.
/// This and l2 are treated as infinite lines.
///Function's return value is:
/// true if two lines are parallel(or one of the directin is zero)
/// false if intersection was obtained.
bool isect2D(const MGStraight& l2, double& t,double& s) const;

///Compute intersections with MGLBRep curve2 that does not have C0 continuity in it.
MGCCisects isect_withC1LB(const MGLBRep& curve2)const;

///isect with SurfCurve whose m_curve is not a MGTrimmedCurve of MGCompositeCurve.
MGCCisects isect_with_noCompoSC(const MGSurfCurve& curve2)const;

///Compute parallel range of two straight lines.
///Two straight line this and lines must be parallel.
///Function's return value MGPosition_list list is:
/// list.entries() is 0(when no paralle part) or 2(when there is a parallel
/// part). When list.entries() is 2, let their entries be P1 and P2.
/// Then from P1(0) to P2(0) is the range of this straight line.
/// From P1(1) to P2(1) is the range of line2 straight line.
MGPosition_list relation_parallel(const MGStraight& line2) const;

///Function to avoid m_direction.len()=zero.
///The line's m_direction is set as a unit vector
///and m_endparam is set to zero.
void save_length_zero();

//Compute two straight lines relationship of parallel or coincidence.
//Parallelness of the two is assumed.
//Obtain the relationship when two lines coinside.
//ip of a intersection or nearest point will be returned.
//When this and s do not coincide, MGPSREL_PARALLEL will be returned as
//the function's return value.
MGPSRELATION relation_coincide(const MGStraight& s,	MGCCisect& ip)const;

//relation_coincide when this is MGSTRAIGHT_SEGMENTt.
MGPSRELATION relation_coincide1(
	const MGStraight& sl2,
	MGCCisect& ip
)const;


};

/** @} */ // end of GEO group
#endif
