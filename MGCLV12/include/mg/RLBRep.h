/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#ifndef _MGRLBRep_HH_
#define _MGRLBRep_HH_

#include "mg/LBRep.h"

// MGRLBRep.h
//

// Forward Declaration
class  MGPosition;
class  MGKnotArray;
class  MGCParam_list;
class  MGPosition_list;
class  MGIfstream;
class  MGOfstream;

/** @file */
/** @addtogroup GEO
 *  @{
 */

/// Defines Rational Line B-Representation.

/// This NURBS is a homogeneous form, i.e., B-Coefficients have
/// weight included values. 
/// When usual NURBS form is (xi, yi, zi, wi) ,
/// MGRLBRep form is (xi*wi, yi*wi, zi*wi, wi) for i=0,..., n-1.
class MG_DLL_DECLR MGRLBRep: public MGCurve {

public:

///Vector translation.
MG_DLL_DECLR friend MGRLBRep operator+ (const MGVector& v, const MGRLBRep& lb);

///Scaling.
MG_DLL_DECLR friend MGRLBRep operator* (double scale, const MGRLBRep&);

///Default constructor.
MGRLBRep(){ ; };

///Construct Line NURBS, providing all the member data.
///***** This is the fundamental constructor(when homogeneous=1).*****
MGRLBRep(
	const MGKnotVector& t,		///<Knot Vector.
	const MGBPointSeq& bcoef,	///<Line B-Coef, each of coefficients
		///<includes weight multiplied when homogeneous=true(1),
		///<and not includes when homogeneous =false.
		///<Maximum space dimension id of bcoef is for weight of the rational.
	int homogeneous=1///<Indicates if bcoef includes weight as above.
);

///Construct Line NURBS, providing all the member data.
MGRLBRep(
	const MGKnotVector& t,	///<Knot Vector.
	const MGBPointSeq& bcoef,///<Line B-Coef, each of coefficients does not include weights data.
	const std::vector<double>& weights///<Weights.
);

/// Construct ellipse NURBS.
explicit MGRLBRep(const MGEllipse& ellipse);///Original ellipse.

///Construct 2D ellipse RLBRep, whose center is origin.
///The ellipse is expressed as below using parameter t.
/// x(t)=a*cos(t),  y(t)=b*sin(t),   angle1<=t<=angle2
MGRLBRep(
	double a, double b,
	double angle1, double angle2
);

/// Construct a conic section NURBS.
///This conic is defined by ths start and end point, and each tangent,
///and mid-point of the conic.
MGRLBRep(
	const MGPosition& P0,///<Start point
	const MGVector& T0,	///<Start point's tangent
	const MGPosition& P,///<Mid point of the conic section
	const MGPosition& P2,///<End point
	const MGVector& T2	///<End point's tangent
);

///**** Conversion Constructor.****

/// Convert from Non ratoinal form to Rational form.
///When homogeneous==true(non zero), brep is homogeneous form MGLBRep.
///When homogeneous==false(zero), brep is ordinary MGLBRep and
///will be converted to MGRLBRep. That is, weight=1 elements will be
///added as last space dimension element.
///***** This is the fundamental constructor. *****
explicit MGRLBRep(
const MGLBRep& brep,///<Original LBRep. This can be ordinary LBRep, or 
	///<homogeneous form of MGRLBRep. When homogeneous form,
	///<the last space dimension elements are weights.
int homogeneous=0	///<true(non zero): homogeneous form,
					///<false(zero):ordinary LBRep.
);

/// Construct a Line NURBS by changing space dimension and ordering of
///coordinates.
MGRLBRep(
	int dim,			///< New space dimension.
	const MGRLBRep& lbrep,///< Original Line B-rep.
	int start1=0, 		///< Destination order of new line.
	int start2=0 		///< Source order of original line.
);

//////////// Operator overload ////////////

///Assignment.
///When the leaf object of this and crv2 are not equal, this assignment
///does nothing.
MGRLBRep& operator=(const MGGel& gel2);
MGRLBRep& operator=(MGGel&& gel2);

///Transformation object construction
MGRLBRep operator+ (const MGVector& ) const;
MGRLBRep operator- (const MGVector& ) const;
MGRLBRep operator* (double) const;
MGRLBRep operator* (const MGMatrix& ) const;
MGRLBRep operator* (const MGTransf& ) const;

///Object transformation.
MGRLBRep& operator+=(const MGVector& v);
MGRLBRep& operator-=(const MGVector& v);
MGRLBRep& operator*=(double scale);
MGRLBRep& operator*=(const MGMatrix& mat);
MGRLBRep& operator*=(const MGTransf& tr);

///comparison
bool operator==(const MGLBRep& gel2)const;
bool operator==(const MGRLBRep& gel2)const;
std::partial_ordering operator<=>(const MGRLBRep& gel2)const;

//gel2 must be the same class as this.
bool equal_test(const MGGel& gel2)const override;

//gel2 must be the same class as this.
std::partial_ordering ordering_test(const MGGel& gel2)const override;

//////////// Member Function ////////////

///Set(update) the knot vector, KV=MGKnotVector.
///This is move operation conforming.
template<class KV>
void setKnotVector(
	KV&& t	///<knot vector(MGKnotVector)
){
	m_line.setKnotVector(std::forward<KV>(t));
};

///Approximate this curve as a MGLBRep curve
///within the tolerance MGTolerance::line_zero().
///When parameter_normalization=0, reparameterization will not done, and
///the evaluation at the same parameter has the same values before and after
///of approximate_as_LBRep.
void approximate_as_LBRep(
	MGLBRep& lb,	///<Approximated obrep will be set.
	int ordr=0,		///<new order. When this is MGLBRep, if ordr=0,
					///ordr=order() will be assumed, else ordr=4 is assumed.
	int parameter_normalization=0,
		///<Indicates how the parameter normalization be done:
		///<  =0: no parameter normalization.
		///<  =1: normalize to range=(0., 1.);
		///<  =2: normalize to make the average length of the 1st derivative 
		///     is as equal to 1. as possible.
	bool neglectMulti=false///<Indicates if multiple knots be kept.
		///< true: multiplicity is removed.
		///< false: multiplicity is kept.
)const override;

///Approximate oldBrep by the MGRLBRep of new knot configuration.
///The new knot vector is input from the member data.
///Use setKnotVector() before buildByNewKnotVectorWithKTV() to set the knotvector.
///The new approximated MGRLBRep is replaced with this.
///The parameter range of t must be inside the one of this.
///The constructed MGRLBRep is an approximation
///of this with the parameter range from t.param_s() to t.param_e();
int buildByNewKnotVectorWithKTV(
	const MGRLBRep& oldBrep//The original B-Rep.
){
	invalidateBox();copy_appearance(oldBrep);
	return oldBrep.m_line.rebuildByNewKnotVector(m_line.m_knot_vector,m_line);
};

///Approximate this, given by a new knot configuration newT.
///The new approximated MGRLBRep is output into newRLB.
///The parameter range of newT must be inside the one of this.
///newRLB is an approximation of this with the parameter range
///from newT.param_s() to newT.param_e().
///Function's return value is error code:
/// =0 successfully rebuilt, !=0: input newT is invalid.
///When error!=0 is returned, this is copied int newLB.
int rebuildByNewKnotVector(
	const MGKnotVector& newT,//New knot configuration is input.
	MGRLBRep& newRLB//Rebuilt MGRLBRep is output, can be this.
)const{
	newRLB.invalidateBox();
	newRLB.copy_appearance(*this);
	return m_line.rebuildByNewKnotVector(newT, newRLB.m_line);
};

///Gets new B-Rep by subdividing the original one into a part.
///New one is exactly the same as the original except that it is partial.
///id1 and id2 are id's of this knot_vector, and indicate the parameter range
///as from t[id1] to t[id2]. Here t=this->knot_vector().
///shrinkToKnots() employs the partial knot vector of t and this B-coefficients.
///And so, knot multiplicity of start and end of the new knot vector is not guaranteed.
///It depends on the original one.
void shrinkToKnots(
	int id1,///< start id of this knot vector.
	int id2,///< End id of this knot vector.
	MGRLBRep& newRLB//Shrinked MGRLBRep is output, can be this.
)const;

///Gets new B-Rep by computing a part of the original. New one is exactly
///the same as the original except that it is partial.
///If multiple==true(!=0), knot(i)=t1 and knot(n+i)=t2 for i=0,..., k-1
///are guaranteed. Here, n=bdim() and k=order().
///Both t1 and t2 must be inside te range of this.
void shrinkToParameters(
	double t1,///<New parameter range. t1 must be less than t2.
	double t2,///< End of the parameter range. 
	MGRLBRep& newRLB,//Shrinked MGRLBRep is output, can be this.
	int multiple=0///<Indicates if start and end knot multiplicities
					///<are necessary. =0:unnecessary, !=0:necessary.
)const;

///Gets a new B-Rep by adding knots to an original B-Rep.
void addKnots(const MGKnotArray& knots);///<Knots to add.

///Returns NURBS Dimension.
int bdim() const{return m_line.bdim();}
	
/// 入力のパラメータ範囲の曲線部分を囲むボックスを返す。
///Return minimum box that includes the partial line.
MGBox box_limitted(const MGInterval& l) const;

///Return minimum box that includes the whole line.
MGBox box_unlimit() const;

///Change parameter range, be able to change the direction by providing
///t1 greater than t2.
void change_range(
	double t1,			///<Parameter value for the start of original. 
	double t2			///<Parameter value for the end of original. 
){
	m_line.change_range(t1,t2);
	invalidateBox();
};

///Changing this object's space dimension.
void change_dimension(
	int dim,			///< new space dimension
	int start1=0, 		///< Destination order of new object.
	int start2=0 		///< Source order of this object.
);

///Change order of the NURBS. When new order is greater than the original,
///new B-rep is guaranteed to be the same line as the original. However,
///if new order is less than the original one, new line is not the same
///in general.
void change_order(
	int order		///<New order number. 
){
	m_line.change_order(order);
};

///Access to (i,j)th element of coef
///( left-hand side version)
double& coef(int i, int j){invalidateBox(); return m_line.coef(i,j);};		

///Access to (i,j)th element of coef
///(right hand side version)
double coef(int i, int j) const{return m_line.coef(i,j);};

///Extract (i,j)elements for 0<=j<=sdim() as a vector
///of (sdim()+1) space dimension. The last elemnt is weight.
MGVector coef(int i) const{return m_line.coef(i);};

///Returns a pointer to the line b-coef data.
const double* coef_data(int i=0, int j=0) const
{return m_line.coef_data(i,j);};

///Compute the box of the whole of the curve.
///Compute the box from the scratch.
void compute_box(MGBox& bx) const;

///Connect brep2 to this brep to make one B-Representation.
///This parameter range will not be changed, instead brep2's range
///will be so changed that brep2 has the same 1st derivative magnitude
///as the original this brep's at the connecting point(start or end point of
///this).
///continuity and which can be obtained using the fucntion continuity().
void connect(
	int continuity,	///<continuity. must be>=0.
	int which,		///<which point of this to which of brep2.
				///< =0: start of this and start of brep2.
				///< =1: start of this and end of brep2.
				///< =2: end of this and start of brep2.
				///< =3: end of this and end of brep2.
	const MGRLBRep& brep2	///<B-Rep 2.
);

///Compute continuity with brep2.
/// Function's return value is:
/// -1: G(-1) continuity, i.e. two lines are discontinuous.
///  0: G0 continuity, i.e. two lines are connected,
///     but tangents are discontinuous
///  1: C1 continuity, i.e. 1st deriv's  are also continuous,
///     when weights are so arranged.
///  2: C2 continuity, i.e. 2nd deriv's  are also continuous,
///     when weights are so arranged.
int continuity(
	const MGRLBRep& brep2,///<The 2nd target RLBRep.
	int& which,	///<Indicates which point of this is connected
				///< to which of brep2, is meaingfull when continuity>=0,
				///< =0: start of this to start of brep2,
				///< =1: start of this to end of brep2,
				///< =2: end of this to start of brep2,
				///< =3: end of this to end of brep2.
	double& ratio ///< Ratio of 1st derivatives of the two line will
			///< be returned,
			///< ratio= d2/d1, where d1=1st deriv of this and d2=of brep2
) const;

///Exchange ordering of the coordinates.
///Exchange coordinates (i) and (j).
void coordinate_exchange(int i, int j);

///Construct new curve object by copying to newed area.
///User must delete this copied object by "delete".
MGRLBRep* clone() const;

///copy as a newed curve. The new curve will be MGLBRep or MGRLBRep.
///When original curve was a MGRLBRep, the new curve will be a MGRLBRep.
///Otherwise,  the new curve will be a MGLBRep.
///Returned object must be deleted.
MGCurve* copy_as_nurbs() const override{return clone();};

///Construct new curve object by changing
///the original object's space dimension.
///User must delete this copied object by "delete".
MGRLBRep* copy_change_dimension(
	int sdim,			///< new space dimension
	int start1=0, 		///< Destination order of new line.
	int start2=0 		///< Source order of this line.
)const;

///Construct new curve object by copying to newed area,
///and limitting the parameter range to prange.
///Returned is newed object and must be deleted.
MGCurve* copy_limitted(const MGInterval& prange) const;

///Divide this curve at the designated knot multiplicity point.
///Function's return value is the number of the curves after divided.
int divide_multi(
	std::vector<UniqueCurve>& crv_list,	//divided curves are appended.
	int multiplicity=-1	///<designates the multiplicity of the knot to divide at,
						///<When multiplicity<=0, order()-1 is assumed,
						///<When multiplicity>=order(), order() is assumed.
) const override;

///Display control polygons using mgVBO::MGDrawPointSeq()
void display_control_polygon(mgSysGL& sgl)const;

///Draw this line's 1st and 2nd coordinates in 2D space.

///Done using drawing function moveto( , ) and lineto( , ).
///wind[] is the window of the screen to draw the line in.
///Clipping will be performed about the wind[].
///(wind[0], wind[1]) is the center coordinates of the window.
///wind[2] is width and wind[3] is hight of the window. When wind[2]<=0,
///no clipping is performed. Even when wind[2]<=0, wind[3] is necessary 
///to input to specify the resolution of the line. In this case,
///wind[0] and wind[1] are not referended.
///ynum is the resolution of the line, is the number of
///straight line segments for the curve length of wind[3](height of window).
///***draw_2D does not perform box including judment, always performs clipping
///operation and draws the line. Users must do obvious box inclusion test
///if maximum drawing performance is necessary.
void draw_2D(
	void (*moveto)(int, int),///< Move to function.
	void (*lineto)(int, int),///< Line to function.
	const double wind[4],	///<window box to draw the line in.
	int ynum			///<Resolution of the line.
)const;

///Draw this line's 1st and 2nd coordinates in 2D space.
void draw_2D(
	void (*moveto)(float, float),///< Move to function.
	void (*lineto)(float, float),///< Line to function.
	const double wind[4],	///<window box to draw the line in.
	int ynum			///<Resolution of the line.
)const;

///Draw this line's 1st and 2nd coordinates in 2D space.
void draw_2D(
	void (*moveto)(double, double),///< Move to function.
	void (*lineto)(double, double),///< Line to function.
	const double wind[4],	///<window box to draw the line in.
	int ynum			///<Resolution of the line.
)const;

///Draw this line's coordinate'th coordinate in 2D space as 1D.

///Drawn by (t, LBRep(coordinate)) when t_is_x is true, 
///or as ( LBRep(coordinate),t)  when t_is_x is false,  
///Here t is the parameter of the LBRep.
///using drawing function moveto(int, int) and lineto(int,int).
///The other behaviours are the same as draw_2D.
void draw_1D(
	void (*moveto)(int, int),///<Move to function.
	void (*lineto)(int, int),///<Line to function.
	int coordinate,		///<id of coordinate, that is =0:x, =1:y, and so on.
	bool t_is_x,		///<=true:t is x coordinate, and false:t is y.
	const double wind[4],///<window box to draw the line in.
	int ynum		///<Resolution of the line.
)const;

///Draw this line's coordinate'th coordinate in 2D space as 1D.
void draw_1D(
	void (*moveto)(float, float),///<Move to function.
	void (*lineto)(float, float),///<Line to function.
	int coordinate,		///<id of coordinate, that is =0:x, =1:y, and so on.
	bool t_is_x,			///<=true:t is x coordinate, and false:t is y.
	const double wind[4],	///<window box to draw the line in.
	int ynum		///<Resolution of the line.
)const;

///Draw this line's coordinate'th coordinate in 2D space as 1D.
void draw_1D(
	void (*moveto)(double, double),///<Move to function.
	void (*lineto)(double, double),///<Line to function.
	int coordinate,		///<id of coordinate, that is =0:x, =1:y, and so on.
	bool t_is_x,			///<=true:t is x coordinate, and false:t is y.
	const double wind[4],	///<window box to draw the line in.
	int ynum	///<Resolution of the line.
)const;

///Draw this in 3D using vbo.
void drawSE(
	mgVBO& vbo,///<The target graphic object.
	double t0,			///<Start parameter value of the curve.
	double t1			///<End parameter value of the curve.
						///<Draw will be performed from t0 to t1.
)const;

/// Evaluate right continuous n'th derivative data.
/// nderiv=0 means positional data evaluation.
MGVector eval(
	double t,		///< Parameter value.
	int nderiv=0,///< Order of Derivative.
	int left=0		///<Left continuous(left=true)
					///<or right continuous(left=false).
) const;

///Compute position, 1st and 2nd derivatives.
void eval_all(
	double tau,			///<Input parameter value(パラメータ値)
	MGPosition& P,		///<Position(位置)
	MGVector& V1,		///<1st derivative(1次微分値)
	MGVector& V2		///<2nd derivative(2次微分値)
) const;

///Evaluate all of i'th derivative data for 0<=i<=nderiv.

///Output will be put on deriv[j+i*sdim()]
///for 0<=i<=nderiv and 0<=j<sdim(), i.e. 
///deriv[j+i*sdim()] is i-th derivative data for 0<=j<sdim(). 
void eval_all(
	double tau,		///< Parameter value to evaluate.
	int nderiv,	///< Order of Derivative.
	double* deriv,	///< Output area of size (nderiv+1)*sdim().
	int left=0		///<Left continuous(left=true)
					///<or right continuous(left=false).
) const;

///Extrapolate this curve by an (approximate) chord length.
///The extrapolation is C2 continuous.
void extend(
	double length,	///<approximate chord length to extend. 
	bool start=false///<Flag of which point to extend, start or end point of the line.
					///<If start is true extend on the start point.
);

///Modify the original NURBS by extrapolating the specified perimeter.
///The extrapolation is C2 continuous if the order >=4.
void extend(
	int start,			///<Flag of start or end poit of the line,
						///<If start is true extend on the start point.
	double length,		///<chord length to extend. 
	double dk=0.         ///<Coefficient of how curvature should vary at
///<    extrapolation start point. When dk=0, curvature keeps same, i.e.,
///<    dK/dS=0. When dk=1, curvature becomes zero at length extrapolated point,
///<    i.e. dK/dS=-K/length at extrapolation start point,
///<    (S=parameter of arc length, K=Curvature at start point)
///<    That is, when dk reaches to 1 from 0, curve changes to flat.
);

///Extrapolate the curve by the parameter value.
void extend_with_parameter(
	double tau,	///<The parameter value at the end of extended point,
				///<When tau<param_s(), extension will be done at the starting point,
				///<When tau>param_e(), extension will be done at the end point.
	double dk   ///<Coefficient of how curvature should vary at the connecting point.
				///<See extend();
);

///Extracts control points.

///Fucntion's return value is 
///true if control points was obtained, false if not.
bool get_control_points(
	MGBPointSeq& cpoints	///<Control points will be output.
)const;

///Return homogeneous Line B-Representation of the rational B-Spline.
const MGLBRep& homogeneous() const {return m_line;}
MGLBRep& homogeneous(){invalidateBox();return m_line;}

/// Return This object's typeID
long identify_type() const;

///Provide divide number of curve span for function intersect.
int intersect_dnum() const;

///Test if this cure is co-planar with the 2nd curve curve2.

///MGPlane expression will be out to plane if this is co-planar.
///Function's return value is true if co-planar.
bool is_coplanar(const MGCurve& curve2, MGPlane& plane)const;

///Test if this cure is planar or not.

///MGPlane expression will be out to plane if this is planar.
///Function's return value is true if planar.
bool is_planar(MGPlane& plane)const;

///Intersection point of spline and curve.
MGCCisects isect(const MGCurve&) const;
MGCCisects isect(const MGStraight& curve2)const;

///Intersection of Spline and Surface

///Intersection of MGRLBRep and Surface
MGCSisects isect(const MGSurface& surf) const;
MGCSisects isect(const MGPlane& surf) const;

///Compute intersection point of 2D sub NURBS of original B-rep.

///Parameter values of this at intersection points will be returned.
///Straight line sl and this(RLBRep) will be projected to 2D plane of
///coordinate kind (coordinate, coordinate+1), then intersection will
///be computed. For example when sl and this are 3 dimension (x,y,z),
///and coodinate =2, 2D data (z,x) are extracted from sl and this, then
///intersection will be performed.
///**Note** MGStraight sl is treated as infinite straight line,
///even if it is finite.
MGCParam_list isect_2D(
	const MGStraight& sl,///< Straight line.
	int coordinate=0	///< Coordinate kind of 2D sub space.
) const;	

///Compute intersection points of 3D sub NURBS of original B-rep.

///Parameter values of thisat intersection points will be returned.
///This(RLBRep) will be projected to 3D plane of coordinate kind 
///(coordinate, coordinate+1, coordinate+2), then intersection will
///be computed. This is valid only when sdim()>=4. For example when
///pl and this are 4 dimension (x,y,z,p), and coodinate =1,
///3D data (y,z,p) are extracted from pl and this, then
///intersection will be performed.
MGCParam_list isect_3D(
	const MGPlane& pl,	///< Plane.
	int coordinate=0	///< Coordinate kind of 3D sub space.
) const;	

///Access to i-th element of knot(left-hand side version).
double& knot(int i){return m_line.knot(i);}			

///Access to i-th element of knot(right hand side version).
double knot(int i) const{return m_line.knot(i);}

///Returns a pointer to the knot vector data.
const double* knot_data() const{return m_line.knot_data();}		

///Returns the knot vector. RHS version.
const MGKnotVector& knot_vector() const{return m_line.knot_vector();}

///Returns the knot vector. LHS version.
MGKnotVector& knot_vector(){return m_line.knot_vector();}	

///Get the sub interval line of the original line.
void limit(const MGInterval& itvl);

///Returns the B-coef's. RHS version.
const MGBPointSeq& line_bcoef() const{return m_line.line_bcoef();};

///Returns the B-coef's. LHS version.
MGBPointSeq& line_bcoef(){invalidateBox(); return m_line.line_bcoef();};

///Change direction of the line.
void negate(){m_line.negate();};

///Obtain parameter value if this curve is negated by "negate()".
double negate_param(double t)const{return m_line.negate_param(t);};

///Return non_homogeneous B-Coefficients with weights of
///the rational B-Spline.
///This MGBPointSeq includes weights.
MGBPointSeq non_homogeneous_bcoef() const;

///Test if this is actually non_rational, i.e. , all of the weights are
///same values. If non_rational return true, else false.
int non_rational() const;

///Returns the order.
int order() const{return m_line.order();}

/// Return ending parameter value.
double param_e() const{return m_line.param_e();};

///Normalize parameter value t to the nearest knot if their distance is
///within tolerance.
double param_normalize(double t) const{return m_line.param_normalize(t);}

/// Return starting parameter value.
double param_s() const{return m_line.param_s();};

///Compute part of this curve from parameter t1 to t2.
///Returned is the pointer to newed object, and so should be deleted
///by calling program, or memory leaked.
MGRLBRep* part(
	double t1,///< From parameter.
	double t2,///< To parameter.
	int multiple=0	///<Indicates if start and end knot multiplicities
					///<are necessary. =0:unnecessary, !=0:necessary.
) const;

///Compute all the perpendicular points of this curve and the second one.

///If f(s) and g(t) are the points of the two curves f and g,
///then obtains points where the following conditions are satisfied:
///  fs*(f-g)=0.    gt*(g-f)=0.
///Here fs and gt are 1st derivatives at s and t of f and g.
///MGPosition P in the MGPosition_list contains this and crv's parameter
///as:     P(0)=this curve's parameter, P(1)=crv2's parameter value.
MGPosition_list perps(const MGCurve& crv2)const;
MGPosition_list perps(const MGStraight& crv2)const;

///Check if the line B-rep is planar.

///Funtion's return value is;
/// 0: Not planar, nor a point, nor straight line.
/// 1: NURBS is a point.		2: NURBS is a straight line.
/// 3: NURBS is planar.
int planar(
	MGPlane& plane	///<When Brep is not straight line nor a point,
		///< plane is returned. Even when not planar(return value is 0), plane nearest is returned.
	, MGStraight& line		///<When Brep is a line, line is returned.
	, MGPosition& point		///<When Brep is a point, point is returned.
)const;

///Change the NURBS by decreasing B-Rep dimension by ndec.

///Result is an approximation of the origimal NURBS.
///Return value is error flag(=0: successfully reduced, !=0:failure)
int reduce(
	int ndec			///<Number of B-Rep dimension to decrease 
){
	invalidateBox();
	return m_line.reduce(ndec);
};

///Change an original NURBS to new one with subdivided knot configuration.
///Knots t must be subdivided knots.
void refine(
		const MGKnotVector& t	///<knot vector
){
	m_line.refine(t);
	invalidateBox();
};

///Rebuild this NURBS by reconstructing new knot configuration.
std::unique_ptr<MGRLBRep> rebuild_with_new_knot_configuration(
	double error,	//Error alowed to rebuild. If error<=0., MGTolerance::line_zero()
	                //will be employed.
	int parameter_normalization
		//Indicates how the parameter normalization be done:
		//=0: no parameter normalization.
		//=1: normalize to range=(0., 1.);
		//=2: normalize to make the average length of the 1st derivative 
		//    is as equal to 1. as possible.
)const;

///removal knot. line_zero tolerance is used.
void remove_knot();

///Returns the space dimension.
int sdim() const{return m_line.sdim()-1;};

///Return sweep surface from crv
///Returned is a newed MGSurface, must be deleted.
///The sweep surface is defined as:
///This curve(say c(t)) is the rail and the straight line segments from
///C(t)+start_dist*uvec to C(t)+end_dist*uvec are the generatrix.
MGSurface* sweep(
	const MGUnit_vector& uvec,			///<Sweep Direction.
	double start_dist,					///<distance to start edge.
	double end_dist				///<distance to end edge.
) const;

///Return the curve type.
MGCURVE_TYPE type() const{return MGCURVE_TYPE::MGCURVE_RSPLINE;}

/// unlimit this line.
MGCurve& unlimit(){return *this;};

///Unlimit parameter range of the curve to the end point direction
MGCurve& unlimit_end(){return *this;};

///Unlimit parameter range of the curve to the start point direction
MGCurve& unlimit_start(){return *this;};

///Debug Function
public:

///IGES output function. PD126.
int out_to_IGES(
	MGIgesOfstream& igesfile,
	int SubordinateEntitySwitch=0
)const;

///String output.
std::ostream& toString(std::ostream&) const;

/// Offset of constant deviation from this curve.
/// The offset value must be less than radius of curvature.
/// When this curve is not C1 continuous, this is divided into C1 curves,
/// and more than one offset curves are obtained.
/// line_zero() is used to approximate curves of the offset.
std::vector<UniqueCurve> offset(
	double ofs_value,
	bool principalNormal = true /// true: Offset direction is to principal normal
								/// false: to binormal
) const override;

/// Offset of variable deviation from this curve.
/// When this curve is not C1 continuous, divided into C1 curves,
/// and more than one offset curves are obtained.
/// The direction of offset is toward the principal normal,
/// or to the direction to center of curvature.
/// line_zero() is used approximate the offset curve.
std::vector<UniqueCurve> offset(
	const MGLBRep& ofs_value_lb,			///<空間次元１の線B表現で示したオフセット量
	bool principalNormal = true /// true: Offset direction is to principal normal
								/// false: to binormal
) const override;

protected:

///Compute intersection point of 1D sub NURBS of original B-rep.
///Parameter values of this at intersection points will be returned.
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

///Get the name of the class.
std::string whoami()const{return "RLBRep";};

private:

////////////Member Data//////////////

	MGLBRep  m_line;			///< Maximum space dimension id is for weights.

///Draw this by converting straight line segments.
void drawgl(
	mgVBO& vbo,///<The target graphic object.
	double tstart, double tend	///<start and end parameter value of this.
)const;

///Compute intersection points of n-Dimensional sub NURBS and n-dimension
///plane that passes through the origin.
///Parameter values of this at intersection points will be returned.
///MGVector N is the normal vector of the plane.
MGCParam_list isect_nD(						
	const MGVector& N,		///< Normal of n-dimension plane.
	int dimension,		///< Number of dimension.
	int coordinate		///< Coordinate kind of n-dimensional sub NURBS.
) const;

///Internal function for draw_1D, _2D.
///See draw_2D for the comments of move, line, wind, and ynum.
void draw_all2D(
	int kfunc,		///<Kind of function move and line,
		///<1:move(int,int), 2:move(float, float), otherwise:move(double,double).
	int (*moveto)(...), int (*lineto)(...),
	const double wind[4], ///<window box to draw the line in.
	int ynum	///<Resolution of the line.
)const;

///Internal function for draw_1D, _2D.
///See draw_2D for the comments of move, line, wind, and ynum.
void draw_all1D(
	int coordinate,	///<indicates if draw_1D(>=0) or draw_2D(<0)
					///<and coordinate kind if draw_1D.
	bool t_is_x,	///<=true:t is x coordinate, and false:t is y.
	int kfunc,		///<Kind of function move and line,
		///<1:move(int,int), 2:move(float, float), otherwise:move(double,double).
	int (*moveto)(...), int (*lineto)(...),
	const double wind[4], ///<window box to draw the line in.
	int ynum	///<Resolution of the line.
)const;

///Split conic RLBRep at the mid point of i-th span.
///This RLBRep mmust be conic section.
void split_conic(int i);

///2本のB表現曲線を接続する(同じ種類のとき)
///MGCurve* join(const MGCurve& crv1) const;

};

///@cond
///Function to compute control point P1 and weight w1 of rational form of
///an ellipse segment. Pi and Ti are points and tangents of start and end
///for i=0,2. P is mid point of the ellipse.
///Function's output is if obtained(!=0:true) or not(=0:false).
///When obtained, =1:as finite control point, =2:as infinite.
///When T0, T2, P0, and P2 are not in one plane, function return 0.
///
///(P0,1.) (P1,w1) (P2,1.) constitute the ellipse control polygon
///of order 3 in homogeneous form.
///
///See "The NURBS Book" of W.Tiller and L.Piegl publised by Springer.
MG_DLL_DECLR int MGRLBRep_ellipse_weight
	(const MGPosition& P0, const MGVector& T0,
	 const MGPosition& P,
	 const MGPosition& P2, const MGVector& T2,
	 MGPosition& P1, double& w1
);
///@endcond

/** @} */ // end of GEO group
#endif
