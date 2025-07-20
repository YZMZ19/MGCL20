/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#pragma once

#include "mg/MGCL.h"
#include "mg/Vector.h"
#include "mg/KnotVector.h"
#include "mg/BPointSeq.h"
#include "mg/Curve.h"

// MGLBRep.h
// Forward Declaration
class  MGInterval;
class  MGNDDArray;
class  MGPosition;
class  MGOscuCircle;
class  MGKnotArray;
class  MGCParam_list;
class  MGPosition_list;
class  MGLBRepEndC;
class  MGRLBRep;
class  MGLBRep;
class  MGPPRep;
class  MGSBRepTP;
class  MGIfstream;
class  MGOfstream;

// Defines Line B-Representation.
/** @file */
/** @addtogroup GEO
 *  @{
 */

///MGLBRep is a class for B-SPline representation.

///For a general NURBS(non uniform rational B-Spline) is MGRLBRep.
///LBRep abbrebiates Line B-Representation. MGLBRep consists of a knot vector(MGKnotVector)
///and a control polygon(MGBPointSeq) whose B-representaiton dimension are the same.
class MG_DLL_DECLR MGLBRep : public MGCurve{
///Member Data
private:
MGKnotVector m_knot_vector;	// Knot Vector.
MGBPointSeq  m_line_bcoef;	// Line B-Coef.

public:

///translation by a vector.
MG_DLL_DECLR friend MGLBRep operator+ (const MGVector& v, const MGLBRep& lb);

///scaling by a scalar.
MG_DLL_DECLR friend MGLBRep operator* (double scale, const MGLBRep&);

public:

////////Special member functions/////////
virtual ~MGLBRep()=default;
MGLBRep(const MGLBRep&)=default;
MGLBRep(MGLBRep&&)=default;
MGLBRep& operator=(const MGLBRep&)=default;
MGLBRep& operator=(MGLBRep&&)noexcept =default;

///Default(dummy) constructor.
///Dummy constructor that specifies area length.
///The contents of knot vector and control vertices are garbage.
MGLBRep(
	int bdim=0,	///< b-rep dimension of the lbrep.
	int order=0,	///<order of the lbrep.
	int sdim=0	///<space dimension of the lbrep
):MGCurve(), m_line_bcoef(bdim, sdim), m_knot_vector(order, bdim){;};

// ****  Approximation Constructor ****

///Construct Curve B-Rep.
///This is an approximation, and the tolerance is MGTolerance::line_zero().
explicit MGLBRep(
	const MGCurve& crv,	///<Original Curve.
	int order=0,///<Order. When order=0 is input, and crv was a MGLBRep,
		///<the original order will be used. Otherwise(order=0 and crv was not an MGLBRep)
		///<order will be set to 4.
	int parameter_normalization=0,
		///<Indicates how the parameter normalization be done:
		///<  =0: no parameter normalization.
		///<  =1: normalize to range=(0., 1.);
		///<  =2: normalize to make the average length of the 1st derivative 
		///<    is as equal to 1. as possible.
	bool neglectMulti=false///<Indicates if multiple knots be kept.
		///< true: multiplicity is removed.
		///< false: multiplicity is kept.
);

// **** 3.Conversion Constructor.****

///Convert PP-Rep to B-rep.
MGLBRep(const MGPPRep& pprep);

///Construct a Line B-Rep by changing space dimension and ordering of
///coordinates.
MGLBRep(
	int dim,				///< New space dimension.
	const MGLBRep& lbrep,	///< Original Line B-rep.
	int start1=0, 		///< Destination order of new line.
	int start2=0 		///< Source order of original line.
);

////////Operator overload////////

///Assignment.
///When the leaf object of this and crv2 are not equal, this assignment
///does nothing.
MGLBRep& operator=(const MGGel& gel2);
MGLBRep& operator=(MGGel&& gel2);

///Transformation object construction
MGLBRep operator+ (const MGVector& ) const;
MGLBRep operator- (const MGVector& ) const;
MGLBRep operator* (double) const;
MGLBRep operator* (const MGMatrix& ) const;
MGLBRep operator* (const MGTransf& ) const;

///Object transformation.
MGLBRep& operator+=(const MGVector& v);
MGLBRep& operator-=(const MGVector& v);
MGLBRep& operator*=(double scale);
MGLBRep& operator*=(const MGMatrix& mat);
MGLBRep& operator*=(const MGTransf& tr);

///comparison
bool operator==(const MGRLBRep& gel2)const;

bool operator==(const MGLBRep& gel2)const;
std::partial_ordering operator<=>(const MGLBRep& gel2)const;

//gel2 must be the same class as this.
bool equal_test(const MGGel& gel2)const override;

//gel2 must be the same class as this.
std::partial_ordering ordering_test(const MGGel& gel2)const override;

/////////Member Function////////

///Set(update) the knot vector, KV=MGKnotVector.
///This is move operation conforming.
template<class KV>
void setKnotVector(
	KV&& t	///<knot vector(MGKnotVector)
){
	m_knot_vector=std::forward<KV>(t);
};

///Construct this MGLBRep, given all the necessary member data,
///CV=MGBPointSeq, KV=MGKnotVector.
///This is move operation conforming. When one of vertex or t
///is movable, use this form after void constructor.
template<class KV, class CV>
void buildLBRepFromMemberData(
	KV&& t,	///<knot vector(MGKnotVector)
	CV&& vertex///<Control Vertex of Line BRep(MGBPointSeq).
){
	assert(t.bdim()==vertex.length());
	m_line_bcoef=std::forward<CV>(vertex);
	m_knot_vector=std::forward<KV>(t);
	invalidateBox();
};

///Construct Line B-rep by intepolation from Point data only.
///If circular is true, start and end points must be the same,
///(that is points[0]==points[n-1]) and the MGLBRep constructed is 
///smoothly connected at the point.
///If circular is true, order will be always 4, and input order is neglected.
///The knot vector is generated from data points of the input points'
///chord length starting 0.
///Function's return value is 0:successfully built, !=0:too near data is input.
///When error is returned, B-Rep of polyline of points is constructed that has uniform
///knot vector.
int buildByInterpolation(
	const MGBPointSeq& points,	///<Point seq data
	int order=4,		///<Order
	bool circular=false	///<Circular flag. If true, circular.
);

///Construct Line B-rep of a specified order, given data point abscissa and
///the ordinates.
///The knot vector is generated from data points tau.
///Function's return value is 0:successfully built,
///!=0:tau is illegal(too near data is input).
///knot vector.
int buildByInterpolationDataPoints(
	const MGNDDArray& tau,		///<Data point abscissa
	const MGBPointSeq& points,	///<Point seq data
	int order=4,			///<order
	double ratio=-1.	///<Maximum of data point ratio of pre and after spans.
			///< Let d(i)=tau[i]-tau[i-1], then if d(i)/d(i-1)>ratio or
			///< d(i-1)/d(i)>ratio, either tau[i] or tau[i-1] will be removed.
			///< This is done to prevent control polygon computation error.
			///< When ratio<0. no data point removal will be done.
);

///Construct Line B-rep of order by interpolation from Point data
///and end condition. (tau(i), value(i,.)) for 0<=i<=n(the length of value).
///For the start and end point, tau does not have multiplicity. However,
///if tau has multiplicity at inner point, this means 1st derivative data is
///provided for the associated value, i.e.:
///If tau has multiplicity 2 as tau(i)=tau(i+1), value(i,.) is 1st derivative
///at tau(i) and value(i+1,.) is positional data at tau(i)(=tau(i+1)).
///If tau has multiplicity 3 as tau(i)=tau(i+1)=tau(i+2),
///value(i,.) is 1st derivative at tau(i)- ,
///value(i+1,.) is positional data at tau(i)(=tau(i+1)), 
///value(i+2,.) is 1st derivative at tau(i)+.
///Maximum multiplicity allowed is 3.
int buildByInterpolationEC(
	const MGLBRepEndC& begin,	///<Begin end condition
	const MGLBRepEndC& end,		///<End end conditoion
	const MGNDDArray& tau,		///<Data point abscissa
	const MGBPointSeq& value,	///<Data point ordinate
	int order=4
);

///Construct Line B-rep by interpolation, given the knot vector in member m_knot_vector,
///and (tau, points) as arguments.

///Before buildByInterpolationWithKTV(), use setKnotVector to construct this knot vector.
///Function's return value is error code of blgint_(of cskernel) if not zero.
///When error is returned, B-Rep of polyline of points is constructed that has uniform
///knot vector.
int buildByInterpolationWithKTV(
	const MGNDDArray& tau,	//Data point abscissa
	const MGBPointSeq& points	//Point seq data
);

///Construct Line B-rep of any order by interpolation from Point data
///with end condition as arguments and the knot vector as the member data.

/// (tau(i), value(i,.)) for 0<= i <n(n is tau.length()),
///tau(i) and knot vector t must satisfy Shoenberg's variation diminishing
///constraint.
///Before buildByInterpolationECWithKTV(), use setKnotVector to construct this knot vector.
///For the start and end point, tau does not have multiplicity. However,
///if tau has multiplicity at inner point, this means 1st derivative data is
///provided for the associated value, i.e.:
///If tau has multiplicity 2 as tau(i)=tau(i+1), value(i,.) is 1st derivative
///at tau(i) and value(i+1,.) is positional data at tau(i)(=tau(i+1)).
///If tau has multiplicity 3 as tau(i)=tau(i+1)=tau(i+2),
///value(i,.) is 1st derivative at tau(i)- ,
///value(i+1,.) is positional data at tau(i)(=tau(i+1)), 
///value(i+2,.) is 1st derivative at tau(i)+.
///Maximum multiplicity allowed is 3.
///Function's return value is error code:
/// =0: successfully built.
/// !=0: error occured and the MGLBRep is built as uniform BSpline of order 2.
int buildByInterpolationECWithKTV(
	const MGLBRepEndC& begin,	///<Begin end condition
	const MGLBRepEndC& end,		///<End end conditoion
	const MGNDDArray& tau,		///<Data point abscissa
	const MGBPointSeq& value	///<Data point ordinate
);

///Construct this MGLBRep, given MGPPRep and movable knot vector KV.
template<class KV>
void buildLBRepFromPPRep(
	KV&& t,	///<knot vector
	const MGPPRep& pp	//PP-rep
){
	m_knot_vector=std::forward<KV>(t);
	buildLBfromPPRep(pp);
};

///Construct Line B-rep of any order number by least square approximation
///from Point data with approximation weights and knot vector of B-Rep.

///Knot vector is input from the member data, use setKnotVector()
///before buildByL2ApproximateWithKTV() to set the knotvector.
///weight[i] is for points(i,.). For detail information of the approximation
///method, see "A Practical Guide to Splines by Carl de Boor" Springer-Verlag.
void buildByL2ApproximateWithKTV(
	const MGNDDArray& tau,		///<Data point abscissa
	const MGBPointSeq& points,	///<Point seq data
	const double* weight		///<Weights for each points 
);

///Construct 3D B-Rep by mixing two 2D B-Rep.

///The two 2D B-Rep's directions and start and end points must be the same.
///Second 2D B-Rep can be girth representaion, ie,
///Let brep1 is f(t)=(f1(t),f2(t)) and brep2 is g(s)=(g1(s),g2(s)), where
///f1,f2 are two coordinates, g1 is parameter t of f(t), 
///g2(s) is the missing coordinate of f(t). Given parameter s,
/// ( f1(g1(s)), f2(g1(s)), g2(s)) is a 3D space point.
void buildByMixing2DLBRep(
	int coordinate1,	///<Missing oordinate kind of the brep1 
							///< 0:x, 1:y, 2:z.
	const MGLBRep& brep1,	///<Original 2D B-Rep1. Coordinates are
		///<(y,z), (z,x), (x,y) according to coordinate1.
	int coordinate2,	///<Missing coordinate kind of the brep2.
							///< 0:x, 1:y, 2:z, and 3:girth rep.
	const MGLBRep& brep2	///<Original 2D B-Rep2. Coordinates are
		///<(y,z), (z,x), (x,y) and (t, g2) according to coordinate2.
		///<t is parameter of brep1 and g2 is x, y, or z according to coordinate1.
);

///Construct Line B-rep of order 4 from point and point-kind followed by
///osculating circle data.

///point_kind[i] is point kind of the point points(i,.):
/// =0:G2 point, =1:G0 point, =2:G1 point.
///If two consecutive points are 1 or 2,
///the span is a straight line. point_kind 2 is a start of G2 curve.
///If two straight line span meet at points(i), osculating circle can be generated
///at this point by providing circle data at this point in circle.
void buildFromPointKindCircles(
	const MGLBRepEndC& begin,	///<Begin end condition
	const MGLBRepEndC& end,		///<End end conditoion
	const MGBPointSeq& points,	///<Point seq data
	const int* point_kind,		///<Point kind of above point.
	const MGOscuCircle& circle	///<Provides osculating circle data.
);

///Build line B-Rep by Schoenberg and Reinsch smoothing function.

///The end conditions are supposed to be free, given
///data points (tau,y), weights dy at data points, and a deviation.
///If dy[i] gets larger, deviation at tau(i) gets larger.
///n can be any number greater than or equal to 2.
///***End conditions are free end condition.***
void buildSRSmoothedLB_of_FreeEnd(
	const MGNDDArray& tau,	///<Data point abscissa
	const MGBPointSeq& y,	///<Data point ordinates.
	const double* dy,///<dy[i] is the weights  at tau[i] for i=0,..., tau.length()-1.
	double deviation,///<if dev_is_sum is true,
		///<deviation is the upper bound of Sum(((points(i)-pout(i))/dp[i])**2.
		///<if dev_is_sum is false, deviation is max_deviation of each point at tau[i],i.e.,
		///<dev_is_sum=true: deviation>=Sum(((points(i)-pout(i))/dp[i])**2),
		///<dev_is_sum=false:deviation>=Max((points(i)-pout(i))**2),
		///<for i=0,...,n-1. Here pout(i) is the this->eval(tau(i)).
	bool dev_is_sum=false///<See deviation.
);

///Approximate oldBrep by the MGLBRep of new knot configuration.

///The new knot vector is input from this member data.
///Use setKnotVector() before buildByNewKnotVectorWithKTV() to set the knotvector.
///The new approximated MGLBRep is replaced with this.
///The parameter range of t must be inside the one of this.
///The constructed MGLBRep is an approximation of oldBrep
///with the parameter range from t.param_s() to t.param_e();
int buildByNewKnotVectorWithKTV(
	const MGLBRep& oldBrep//The original B-Rep.
){
	return oldBrep.rebuildByNewKnotVector(m_knot_vector, *this);
};

///Build line B-Rep by Schoenberg and Reinsch smoothing function, given
///1st derivatives info at the ends.
///Data points (tau,y), weights dy at data points, and a mean deviation deviation
///are input.
///If dy[i] gets larger, deviation at tau(i) gets larger.
///n can be any number greater than or equal to 2.
void buildSRSmoothedLB_of_1stDeriv(
	const MGLBRepEndC& begin,///<Begin end condition
	const MGLBRepEndC& end,	///<End end conditoion.
		///<begin.cond() and end.cond() must be MGENDC_1D or MGENDC_12D.
	const MGNDDArray& tau,	///<Data point abscissa
	const MGBPointSeq& y,	///<Data point ordinates.
	const double* dy,///< dy[i] is the weight  at tau[i] for i=0,..., tau.length()-1.
	double deviation,///<if dev_is_sum is true,
		///<deviation is the upper bound of Sum(((points(i)-pout(i))/dp[i])**2.
		///<if dev_is_sum is false, deviation is max_deviation of each point at tau[i],i.e.,
		///<dev_is_sum=true: deviation>=Sum(((points(i)-pout(i))/dp[i])**2),
		///<dev_is_sum=false:deviation>=Max((points(i)-pout(i))**2),
		///<for i=0,...,n-1. Here pout(i) is the this->eval(tau(i)).
	bool dev_is_sum = false///<See deviation.
);

///Approximate this, given the new knot configuration newT.

///The new approximated MGLBRep is output into newLB.
///The parameter range of newT must be inside the one of this.
///newLB is an approximation of this with the parameter range
///from newT.param_s() to newT.param_e().
///Function's return value is error code:
/// =0 successfully rebuilt, !=0: input newT is invalid.
///When error!=0 is returned, this is copied int newLB.
int rebuildByNewKnotVector(
	const MGKnotVector& newT,//New knot configuration is input.
	MGLBRep& newLB//Rebuilt MGLBRep is output, can be this.
)const;

///Approximate this curve as a MGLBRep.

/// Approximation is done within the tolerance MGTolerance::line_zero().
///When parameter_normalization=0, reparameterization will not done, and
///the evaluation at the same parameter has the same values before and after
///of approximate_as_LBRep.
void approximate_as_LBRep(
	MGLBRep& lb,///<Approximated obrep is set, can be this.
	int ordr = 0,	///<new order. When this is MGLBRep, if ordr=0,
				///<ordr=order() will be assumed, else ordr=4 is assumed.
	int parameter_normalization = 0,
	///<Indicates how the parameter normalization be done:
	///<  =0: no parameter normalization.
	///<  =1: normalize to range=(0., 1.);
	///<  =2: normalize to make the average length of the 1st derivative 
	///<    is as equal to 1. as possible.
	bool neglectMulti = false///<Indicates if multiple knots be kept.
		///< true: multiplicity is removed.
		///< false: multiplicity is kept.
)const override;

///Gets new B-Rep by shrinking the original one into a part.

///New one is exactly the same as the original except that it is partial.
///id1 and id2 are id's of this knot_vector, and indicate the parameter range
///as from t[id1] to t[id2]. Here t=this->knot_vector().
///shrinkToKnots() employs the partial knot vector of t and this B-coefficients.
///And so, knot multiplicity of start and end of the new knot vector is not guaranteed.
///It depends on the original one.
void shrinkToKnots(
	int id1,///< start id of this knot vector.
	int id2,///< End id of this knot vector.
	MGLBRep& NewBrep///<shirinked B-Rep　is output, can be this.
)const;

///Gets new B-Rep by computing a part of the original.

///New one is exactly the same as the original except that it is partial.
///If multiple==true(!=0), knot(i)=t1 and knot(n+i)=t2 for i=0,..., k-1
///are guaranteed. Here, n=bdim() and k=order().
///Both t1 and t2 must be inside te range of this.
void shrinkToParameters(
	double t1,///<New parameter range. t1 must be less than t2.
	double t2,///< End of the parameter range. 
	MGLBRep& NewBrep,	///<shirinked B-Rep　is output, can be this.
	int multiple = 0///<Indicates if start and end knot multiplicities
					///<are necessary. =0:unnecessary, !=0:necessary.
)const;

///Gets a new B-Rep by adding knots to an original B-Rep.
void addKnots(const MGKnotArray& knots);///<Knots to add.

///Returns B-Rep Dimension.
int bdim() const { return m_line_bcoef.length(); };

///Return minimum box that includes the partial line.
MGBox box_limitted(const MGInterval& l) const;

///Return minimum box that includes the whole line.
const MGBox& box_unlimit() const;

///Changing this object's space dimension.
void change_dimension(
	int sdim,		///< new space dimension
	int start1=0,	///< Destination order of new object.
	int start2=0 	///< Source order of this object.
);

///Change parameter range, be able to change the direction by providing
///t1 greater than t2.
void change_range(
	double t1,	///<Parameter value for the start of original. 
	double t2	///<Parameter value for the end of original. 
);

///Change order of the B-Rep.
///When new order is greater than the original,
///new B-rep is guaranteed to be the same line as the original. However,
///if new order is less than the original one, new line is not the same
///in general.
void change_order(
	int order	///<New order number. 
);

///Change order of the B-Rep by approximation.
void change_order_by_approximation(
	int ordr	///<New order number. 
);

///Access to (i,j)th element of coef(left-hand side version).
double& coef(int i, int j){return m_line_bcoef(i,j);}

///Access to (i,j)th element of coef (right hand side version).
double coef(int i, int j) const{return m_line_bcoef(i,j);}

///Extract (i,j)element for 0<=j<sdim().
MGVector coef(int i) const{return m_line_bcoef(i);}

///Returns a pointer to the line b-coef data.
const double* coef_data(int i=0, int j=0)const{return m_line_bcoef.data(i,j);}
double* coef_data(int i=0, int j=0){return m_line_bcoef.data(i,j);}

///Compute whole box of the curve.
void compute_box(MGBox& bx) const;

///Connect brep2 to this brep to make one B-Representation.

///This parameter range will not be changed, instead brep2's range
///will be so changed that brep2 has the same 1st derivative magnitude
///as the original this brep's at the connecting point(start or end point of
///this).
///continuity and which can be obtained using the fucntion continuity().
void connect(
	int continuity,	///<continuity. must be>=0.
	int which,	///<which point of this to which of brep2.
				///< =0: start of this and start of brep2.
				///< =1: start of this and end of brep2.
				///< =2: end of this and start of brep2.
				///< =3: end of this and end of brep2.
	const MGLBRep& brep2	///<B-Rep 2.
);

///Compute continuity with brep2.

/// Function's return value is:
/// -1: G(-1) continuity, i.e. two lines are discontinuous.
///  0: G0 continuity, i.e. two lines are connected,
///     but tangents are discontinuous
///  1: G1 continuity, i.e. two lines are connected,
///     and tangents are also continuous.
///  2: G2 continuity, i.e. two lines are connected,
///     and tangents and curvatures are also continuous.
int continuity(
	const MGLBRep& brep2,///<The target 2nd curve.
	int& which,	///<Indicates which point of this is connected
				///< to which of brep2, is meaingfull when continuity>=0.
				///< =0: start of this to start of brep2.
				///< =1: start of this to end of brep2.
				///< =2: end of this to start of brep2.
				///< =3: end of this to end of brep2.
	double& ratio///< Ratio of 1st derivatives of the two line will be returned.
			///< ratio= d2/d1, where d1=1st deriv of this and d2=of brep2
) const;

///Exchange ordering of the coordinates.
///Exchange coordinates (j1) and (j2).
void coordinate_exchange(int j1, int j2);

///Construct new curve object by copying to newed area.
///User must delete this copied object by "delete".
virtual MGLBRep* clone()const;

///Convert this curve to Bezier curve.

///If this is MGLBRep or MGStraight, the shape is exactly the same
///as the original. Otherwise, this is apporoximated by MGLBRep.
///The result MGLBRep is of order 2 if original order is 2 
///and is of order4 otherwise.
///Output bezier can be this.
void convert_to_Bezier(MGLBRep& bezier)const;

///copy as a newed curve.

///The new curve will be MGLBRep or MGRLBRep.
///When original curve was a MGRLBRep, the new curve will be a MGRLBRep.
///Otherwise,  the new curve will be a MGLBRep.
///Returned object must be deleted.
MGLBRep* copy_as_nurbs() const override{return clone();};

///Construct new curve object by changing the original object's space dimension.

///User must delete this copied object by "delete".
MGLBRep* copy_change_dimension(
	int sdim,			///< new space dimension
	int start1=0, 		///< Destination order of new line.
	int start2=0 		///< Source order of this line.
)const;

///Construct a newed curve object.

///Copy is done by limitting the parameter range to prange.
///Returned is newed object and must be deleted.
MGCurve* copy_limitted(const MGInterval& prange) const override;

///Compute curvilinear integral of the 1st two coordinates.
///This integral can be used to compute area sorounded by the curve.
double curvilinear_integral(double t1, double t2) const;

///Divide this curve at the designated knot multiplicity point.
///Function's return value is the number of the curves after divided.
int divide_multi(
	std::vector<UniqueCurve>& crv_list,	//divided curves are appended.
	int multiplicity=-1	///<designates the multiplicity of the knot to divide at.
						///<When multiplicity<=0, order()-1 is assumed.
						///<When multiplicity>=order(), order() is assumed.
)const override;

///Display control polygons using mgVBO::MGDrawPointSeq()
void display_control_polygon(mgSysGL& sgl)const;

///Draw using vbo.
virtual void drawSE(
	mgVBO& vbo,///<The target graphic object.
	double t0,///<Start parameter value of the curve.
	double t1///<End parameter value of the curve.
			///<Draw will be performed from t0 to t1.
)const;

/// Evaluate right continuous n'th derivative data.
/// nderiv=0 means positional data evaluation.
MGVector eval(
	double,			///< Parameter value.
	int nderiv=0,///< Order of Derivative.
	int leftcon=0	///< Left continuous(leftcon=true)
					///< or right continuous(leftcon=false).
)const;

///Compute position, 1st and 2nd derivatives.
/// パラメータ値を与えて位置、一次微分値、二次微分値をもとめる。
void eval_all(
	double,		///< Input parameter value(パラメータ値)
	MGPosition&,///< Position(位置)
	MGVector&,	///< 1st derivative(1次微分値)
	MGVector&	///< 2nd derivative(2次微分値)
)const;

///Evaluate all of i'th derivative data for 0<=i<=nderiv.

///Output will be put on deriv[j+i*sdim()]
///for 0<=i<=nderiv and 0<=j<sdim(), i.e. 
///deriv[j+i*sdim()] is i-th derivative data for 0<=j<sdim(). 
void eval_all(
	double tau,		///< Parameter value to evaluate.
	int nderiv,	///< Order of Derivative.
	double* deriv,	///< Output area of size (nderiv+1)*sdim().
	int leftcon=0	///<Left continuous(leftcon=true) or right continuous(leftcon=false).
)const;

///Evaluate line data at data point seq.(BLELIN)
void eval_line(
	MGENDCOND begin,		///<Begin end condition
	MGENDCOND end,			///<End end conditoion 
	const MGNDDArray& tau,	///<Data points.
	MGBPointSeq& value		///<Values evaluated.
)const;

///Evaluate line data at data point seq.(BLELIN)
void eval_line(
	const MGNDDArray& tau,	///<Data points.
	MGBPointSeq& value		///<Values evaluated.
)const{eval_line(MGENDCOND::MGENDC_NO, MGENDCOND::MGENDC_NO,tau,value);};

///Extrapolate the curve by the chord length.
///The extrapolation is C2 continuous if the order >=4.
void extend(
	int start,		///<Flag of start or end poit of the line.
					///<If start is true extend on the start point.
	double length,	///<chord length to extend. 
	double dk     ///<Coefficient of how curvature should vary at the connecting point.
///<    extrapolation start point. When dk=0, curvature keeps same, i.e.
///<    dK/dS=0. When dk=1, curvature becomes zero at length extrapolated point,
///<    i.e. dK/dS=-K/length at extrapolation start point.
///<    (S=parameter of arc length, K=Curvature at start point)
///<    That is, when dk reaches to 1 from 0, curve changes to flat.
);

///Extrapolate this curve by an (approximate) chord length.
///The extrapolation is C2 continuous.
void extend(
	double length,	///<approximate chord length to extend. 
	bool start=false///<Flag of which point to extend, start or end point of the line.
					///<If start is true extend on the start point.
);

///Extrapolate the curve by the parameter value.
void extend_with_parameter(
	double tau,	///<The parameter value at the end of extended point.
				///<When tau<param_s(), extension will be done at the starting point.
				///<When tau>param_e(), extension will be done at the end point.
	double dk	///<Coefficient of how curvature should vary at the connecting point.
				///<See extend();
);

///Extracts control points.
///Fucntion's return value is 
///true if control points was obtained, false if not.
bool get_control_points(
	MGBPointSeq& cpoints	///<Control points will be output.
)const;

/// Return This object's typeID
long identify_type() const;

///Test if this is a Bezier Curve.

///Functions's return value is MGLBRep* if Bezier, null if not.
///If input ordr>=2, order is also tested if this Bezier's order is the same as input order.
///If input ordr<=1, any ordr>=2 is allowed for Bezier curve.
///Bezier curve is defined as follows. Here t=knot_vector(), k is this LBRep's order,
///n=bdim(), and m=(n-k)/(k-1).
///(1) n=k+(k-1)*m.
///(2) t(0)=t(1)=,...,=t(k-1)=0
///(3) t(i)=t(i+1)=,...,=t(i+k-2)=j+1
///         for i=k, k+(k-1),...,k+j*(k-1) and j=0,...,m-1.
///(4) t(n)=t(n+1)=,...,=t(n+k-1)=m+1
const MGLBRep* is_Bezier(int ordr=0)const;

///Test if this cure is co-planar with the 2nd curve curve2.

///MGPlane expression will be out to plane if this is co-planar.
///Function's return value is true if co-planar.
bool is_coplanar(const MGCurve& curve2, MGPlane& plane)const;

///Test if this cure is planar or not.

///MGPlane expression will be out to plane if this is planar.
///Function's return value is true if planar.
bool is_planar(MGPlane& plane)const;

///Intersection point of spline and curve.
MGCCisects isect(const MGCurve& curve2) const override;
MGCCisects isect(const MGStraight& curve2) const override;
MGCCisects isect(const MGLBRep& curve2) const override;
MGCCisects isect(const MGSurfCurve& curve2) const override;

///Intersection of MGLBRep and Surface
MGCSisects isect(const MGSurface& surf) const;
MGCSisects isect(const MGPlane& surf) const;

///Access to i-th element of knot(left-hand side version).
double& knot(int i){return m_knot_vector(i);};

///Access to i-th element of knot(right hand side version).
double knot(int i) const{return m_knot_vector(i);};

///Returns the raw pointer to the knot vector data.
const double* knot_data() const{return m_knot_vector.data();};
double* knot_data(){return m_knot_vector.data();};

///Returns the knot vector(RHS version).
const MGKnotVector& knot_vector() const{return m_knot_vector;};

///Returns the knot vector(LHS version).
MGKnotVector& knot_vector(){return m_knot_vector;};

///Update this by limiting the parameter range of the curve.
virtual void limit(const MGInterval& i1);

///Returns the B-coef's(RHS version).
const MGBPointSeq& line_bcoef() const{return m_line_bcoef;};

///Returns the B-coef's(LHS version).
MGBPointSeq& line_bcoef(){
	invalidateBox();
	return m_line_bcoef;
};

///Modify the line by moving move_point to to_point.

///fix_point can be
///applied according to move_kind.
/// move_kind=1: Start and end point of the line are fixed. The line is modified
///              linearly so that move_point_param point is the maximum move.
///          =2: The point fix_point[0] is fixed and the other end of
///				the move_point_param side is moved. In this case, maximum move
///              is the end point of the line.
///          =3: fix_point[0]<move_point_param<fix_point[1], and two point
///              fix_point[.] are fixed. The line is modified
///              linearly so that move_point_param point is the maximum move.
///   otherwise: Two fix point fix_point[.] are computed so that the modify
///	            range is the minimum. Other move is same as move_kind=3.
/// Restriction: For the case move_kind=3, actual fix point is wider than
///  specified range. The range is the smallest one possible including
///  fix_point[].
void move(
	int move_kind,			///<Indicates how to move line.
	double move_point_param,///<indicate object point to move by the parameter value.
	const MGPosition& to_point,	///<destination point of the abve source point.
	const double fix_point[2]///<See function explanation.
);

///Change direction of the line.
virtual void negate();

///Obtain parameter value if this curve is negated by "negate()".
double negate_param(double t)const;

///Returns the order.
int order() const{return m_knot_vector.order();};

/// Return ending parameter value.
double param_e() const;

/// Return starting parameter value.
double param_s() const;

///Normalize parameter value t to the nearest knot if their distance is
///within tolerance.
double param_normalize(double t) const;

///Compute part of this curve from parameter t1 to t2.

///Returned is the pointer to newed object, and so should be deleted
///by calling program, or memory leaked.
MGCurve* part(
	double t1,///<from parameter t1.
	double t2,///<To parameter t2.
	int multiple=0	///<Indicates if start and end knot multiplicities
					///<are necessary. =0:unnecessary, !=0:necessary.
)const;

///Return all the foots of the  straight lines that is perpendicular
///to the line.
MGCParam_list perps(
	const MGPosition& point	///< 与ポイント
)const;

///Compute all the perpendicular points of this curve and the second one.

///If f(s) and g(t) are the points of the two curves f and g,
///then obtains points where the following conditions are satisfied:
///  fs*(f-g)=0.    gt*(g-f)=0.
///Here fs and gt are 1st derivatives at s and t of f and g.
///MGPosition P in the MGPosition_list contains this and crv's parameter
///as:     P(0)=this curve's parameter, P(1)=crv2's parameter value.
MGPosition_list perps(const MGCurve& crv2)const;
MGPosition_list perps(const MGStraight& crv2)const;
MGPosition_list perps(const MGRLBRep& crv2)const;
MGPosition_list perps(const MGEllipse& crv2)const;
MGPosition_list perps(const MGLBRep& crv2)const;
MGPosition_list perps(const MGSurfCurve& crv2)const;
MGPosition_list perps(const MGBSumCurve& crv2)const;

///Check if the line B-rep is planar.

///Funtion's return value is;
/// 0: Not planar, nor a point, nor straight line.
/// 1: B-Rep is a point.		2: B-Rep is a straight line.
/// 3: B-Rep is planar.
int planar(
	MGPlane& plane	///<When Brep is not straight line nor a point, plane is returned.
	///<Even when not planar(return value is 0), plane nearest is returned.
	, MGStraight& line	///<When Brep is a line, line is returned.
	, MGPosition& point	///<When Brep is a point, point is returned.
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
int project(
	const MGPlane& plane,	//given surface.
	std::vector<UniqueCurve>& vec_crv_uv,	//uv projection curve will be appended.
	std::vector<UniqueCurve>& vec_crv,	//3d projection curve will be appended.
	const MGVector& vec	//projection vector.
						//if vec = NULL then calculate perpendicular project.
)const;

///Change the B-Rep by decreasing B-Rep dimension by ndec.
///This is an approximation of the origimal B-Rep. Return value is error flag.
int reduce(
	int ndec ///<Number of B-rep dimension to decrease 
);

///Change an original B-Rep to new one with subdivided knot configuration.
///Knots t must be subdivided knots.
void refine(
	const MGKnotVector& t	///<knot vector
);

///removal knot. line_zero tolerance is used.
void remove_knot(){remove_knot(0,0);};

///Remove knot if removed line has the difference less than line_zero().

///The difference is checked only for the space id of coef(.,j+k)
///of j=0, ..., snum-1. When snum=0, snum is set as sdim();
void remove_knot(int j, int snum);

///ノット削除関数(1つのノット)
///戻り値は、削除したノットの数
///When snum!=0, tolerance of totalTol is checked only for coef(.,sid+j),
///where j=0, ..., snum-1. When snum=0, snum is set as sdim();
int remove_knot_one(
	double line0,	///<Tolerance allowed for the knot removal.
					///<When line0=<0., the knot will be removed unconditionally.
	int	nKnot,	///<削除しようとするノットの番号
	double& totalTol,///<誤差合計
	int& num_knot,///<Remained knot number at knot(id) after removed.
	int sid=0,	///<Space dimension start id of this LBRep's B-coef.
	int snum=0	///<Num of space dimension for the totalTol tolerance check.
);

///Returns the space dimension.
int sdim() const{return m_line_bcoef.sdim();};

///Make a sweep surface from crv.

///Returned is a newed MGSurface, must be deleted.
///The sweep surface is defined as:
///This curve(say c(t)) is the rail and the straight line segments from
///C(t)+start_dist*uvec to C(t)+end_dist*uvec are the generatrix.
MGSurface* sweep(
	const MGUnit_vector& uvec,	///<Sweep Direction.
	double start_dist,			///<distance to start edge.
	double end_dist				///<distance to end edge.
)const;

///Return the curve type.
MGCURVE_TYPE type() const{return MGCURVE_TYPE::MGCURVE_SPLINE;};

///Given polar coordinate LBRep, update this to ordinary coordinates system MGLBRep.

///This curve's (x,y) coordinates are polar coordinates system(r,theta), 
///where r is the distance from origin and theta is the angel with x coordinate.
///When return from the function (x,y) are ordinary coordinate system.
///The space dimension of this curve must be >=2;
///If this space dimension is lager than 2, the remaining coordinates are unchanged.
void updatePolarCoordinates2Ordinary();

/// ｌｉｍｉｔをはずす。
MGCurve& unlimit(){return *this;};

///Unlimit parameter range of the curve to the end point direction
MGCurve& unlimit_end(){return *this;};

///Unlimit parameter range of the curve to the start point direction
MGCurve& unlimit_start(){return *this;};

///IGES output function. PD126.
int out_to_IGES(
	MGIgesOfstream& igesfile,
	int SubordinateEntitySwitch=0
)const;

///Debug Function
virtual std::ostream& toString(std::ostream&)const;

///Get the name of the class.
std::string whoami()const{return "LBRep";};

protected:

///Compute intersection point of 1D sub B-Rep of original B-rep.

///Parameter values of intersection point will be returned.
///isect_1D covers this LBRep's C0 ontinuity.
MGCParam_list intersect_1D(						
	double f,			///< Coordinate value
	int coordinate=0	///< Coordinate kind of the data f(from 0).
						///< Id of m_line_bcoef.
)const;	

///@cond
///Obtain so transformed 1D curve expression of this curve that
///f(t)={sum(xi(t)*g[i]) for i=0(x), 1(y), 2(z)}-g[3], where f(t) is the output
///of oneD and xi(t) is i-th coordinate expression of this curve.
///This is used to compute intersections with a plane g[4].
std::unique_ptr<MGCurve> oneD(
	const double g[4]	///<Plane expression(a,b,c,d) where ax+by+cz=d.
)const;
///@endcond

///メンバデータを読み出す関数
/// 戻り値boolは正常に読み出しが出来ればtrue、失敗すればfalseになる
void ReadMembers(MGIfstream& buf);

///メンバデータを書き込む関数
/// 戻り値boolは正常に書き込みが出来ればtrue、失敗すればfalseになる
void WriteMembers(MGOfstream& buf)const;

private:

//Build the MGLBRep from the output of compute_smoothed_p.
void buildLBRepFromSmoothedP(
	double p,
	const MGNDDArray& tau,	//Data point abscissa
	const MGBPointSeq& y,	//Data point ordinates.
	const std::vector<double>& dx,//dx[i]=tau[i+1]-tau[i]] for i=0,.., n-2, is input
	const double* dy,		// dy[i] is the weight  at tau[i] for i=0,..., tau.length()-1.
	const MGBPointSeq& U,	//df2(tau[i])/(6.*p) is input, where df2 is 2nd derivative
					//and p is the output of the compute_smoothed_p.
	const MGBPointSeq& QU	//Q*U is input where Q is the tridiagonal matrix of order n
					//havine the general row [1/dx[i-1], -(1/dx[i-1]+1/dx[i]), dx[i]]
					//for i=0,...,n-1
);

//Build MGLBRep, given PP-Rep(convert PP-Rep to MGLBRep).
//Knot Vector input in m_knot_vector.
//Each knot of the knot vector is break point of pprep.
//The continuities at all the break points must be C(k-2) where
//k is the order of pprep.
void buildLBfromPPRep(const MGPPRep& pp);

///Function for BLUMIX constructor. Given 3D point F,
///compute correct point of F that is closest to F.
MGPosition closest_mix(
	int coordinate1,	///<Missing coordinate kind of this line.
	int coordinate2,	///<Missing coordinate kind of Point P.
	double tau,				///<Parameter value of P of other line,
							///<used only when coordinate2=3.
	const MGPosition& P,	///<Point of 2nd line.
	const MGPosition& F		///<3D point, used to coose closest point to F
							///<when more than one point are found.
)const;

///Draw this by converting straight line segments.
void drawgl(
	mgVBO& vbo,///<The target graphic object.
	double tstart, double tend	///<start and end parameter value of this.
)const;

///Provide divide number of curve span for function intersect.
int intersect_dnum()const;

///compute isects by splitting this curve to sub-curves that do not have
///C0 points in it.
MGCCisects isect_by_split_to_C1(const MGCurve& crv2)const;

///isect for this LBRep that does not include C0 continuity points in it.
///For general intersection computation, use isect.
MGCCisects C1isect(const MGCurve& crv2) const;
MGCCisects C1isect(const MGStraight& crv2) const;

///perps for this LBRep that does not include C0 continuity points in it.
///lb2 may include them. For general perps computation, use perps.
MGPosition_list C1perps(const MGCurve& crv2)const;

///isect_order2 is a private function for isect, computes intersections of
///LBRep(this) of order 2, i.e. polyline B-Rep, with crv or surface.
MGCCisects isect_order2(const MGCurve& crv2) const;
MGCSisects isect_order2(const MGSurface& crv2) const;

///isect with SurfCurve whose m_curve is not a MGTrimmedCurve of MGCompositeCurve.
MGCCisects isect_with_noCompoSC(const MGSurfCurve& curve2)const;

///Compute intersections with MGLBRep curve2 that does not have C0 continuity in it.
MGCCisects isect_withC1LB(const MGLBRep& crv2)const;

///compute isects by splitting this curve to sub-curves that do not have
///C0 points in it.
MGCSisects isect_by_split_to_C1(const MGSurface& surf)const;

///perps_by_split_to_C1 is a private function for perps, computes perpendicular
///points of LBRep(this) with crv.
MGPosition_list perps_by_split_to_C1(const MGCurve& crv2)const;

///Perpendicular points with C1 conitnuity LBRep lbC1.
///MGPosition P in the MGPosition_list contains this and crv's parameter
///as:     P(0)=this curve's parameter, P(1)=crv's parameter value.
MGPosition_list perps_withC1LB(
   const MGLBRep& lbC1
)const;

///Perpendicular points with SurfCurve
///whose m_curve is not a MGTrimmedCurve of MGCompositeCurve.
MGPosition_list perps_with_noCompoSC(const MGSurfCurve& curve2)const;

///perps_order2 is a private function for isect, computes intersection of
///LBRep(this) of order 2, i.e. polyline B-Rep, with crv.
MGPosition_list perps_order2(const MGCurve& crv)const;

///isect for this LBRep that does not include C0 continuity points in it.
///lb2 may include them. For general perpendicular computation, use isect.
MGPosition_list perps_1span(const MGLBRep& lb2)const;

friend class MGCurve;
friend class MGStraight;
friend class MGRLBRep;

};

/** @} */ // end of GEO group
