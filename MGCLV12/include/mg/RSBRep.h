/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#ifndef _MGRSBRep_HH_
#define _MGRSBRep_HH_

#include "mg/Position.h"
#include "mg/SBRep.h"

// MGRSBRep.h
//

// Forward Declaration
class  MGSPointSeq;
class  MGKnotArray;
class  MGMatrix;
class  MGTransf;
class  MGStraight;
class  MGRLBRep;
class  MGPlane;
class  MGCSisects;
class  MGSSisects;
class  MGIfstream;
class  MGOfstream;
/** @file */
/** @addtogroup GEO
 *  @{
 */

/// Defines Surface B-Representation of rational form.

/// This NURBS is of homogeneous form, i.e., B-Coefficients have
/// weight included values. 
/// When usual(non-homogeneous) NURBS form is (xij, yij, zij, wij) ,
/// MGRSBRep form is (xij*wij, yij*wij, zij*wij, wij)
///				 for i=0,..., m-1, and j=0,..., n-1.
class MG_DLL_DECLR MGRSBRep: public MGSurface {

public:

///Scaling.
MG_DLL_DECLR friend MGRSBRep operator* (double scale, const MGRSBRep& sb);


////////Special member functions/////////
MGRSBRep()=default;
~MGRSBRep()=default;
MGRSBRep(const MGRSBRep&)=default;///Copy constructor.
MGRSBRep& operator= (const MGRSBRep&)=default;///Copy assignment.
MGRSBRep(MGRSBRep&&)=default;		///Move constructor.
MGRSBRep& operator= (MGRSBRep&&)=default;///Move assignment.

/// Convert from Non ratoinal form to Rational form.

///When homogeneous==true(non zero), brep is homogeneous form MGSBRep.
///When homogeneous==false(zero), brep is ordinary MGSBRep and
///will be converted to MGRSBRep. That is, weight=1 elements will be
///added to the last space dimension element.
///***** This is the fundamental constructor when homogeneous==1. *****
explicit MGRSBRep(
	const MGSBRep& brep,///<Original SBRep. This can be ordinary SBRep, or 
	///<homogeneous form of MGRSBRep. When homogeneous form,
	///<the last space dimension elements are weights.
	int homogeneous=0	///<true(non zero): homogeneous form,
						///<false(zero):ordinary SBRep.
);

/// Construct a Surface B-Rep by changing space dimension and ordering of coordinates.
MGRSBRep(
	int dim,				///< New space dimension.
	const MGRSBRep& sbrep,	///< Original Surface B-rep.
	int start1=0, 		///< Destination order of new Surface.
	int start2=0 		///< Source order of original Surface.
);

///Approximate an original B-Rep by a new knot configuration.

///The new knot vectors are input from the member data.
///Use setKnotVector() before buildByNewKnotVectorWithKTV() to set the knotvector.
///The new knot config must be inside the range of the original B-Rep
///parameter. However new knots may be coarse or fine.
///Function's return value is Error flag.
///Error is detected only when ut (=2) or vt(=12) is illegal.
///When error!=0, the original old_brep is copied to this.
int buildByNewKnotVectorWithKTV(
	const MGRSBRep& old///<Original B-Rep.
);

///Construct this MGRSBRep, given all the necessary member data.

///SP=MGSPointSeq, KVU and KVV=MGKnotVector.
///MGSPointSeq SP must be homogeneous form(include weights).
///This is move operation conforming. When one of arguments vertex, tu, or tv
///is movable, use this form after void ctor.
template<class SP, class KVU, class KVV>
void buildRSBRepFromMemberData(
	SP&& vertex,///<Control Vertex of Surface B-Rep
	KVU&& tu,	///<knot vector of u-direction
	KVV&& tv	///<knot vector of v-direction
){
	assert(tu.bdim()==vertex.length_u());
	assert(tv.bdim()==vertex.length_v());

	m_surface.m_surface_bcoef=std::forward<SP>(vertex);
	m_surface.m_uknot=std::forward<KVU>(tu);
	m_surface.m_vknot=std::forward<KVV>(tv);
};

///Given bcoef without weight(non-homogeneous form),
///set homogeneous coefficients in m_surface knot vectors.
void buildRSBRepFromMemberDataWithKTV(
	const MGSPointSeq& bcoef,
	//Control Vertex of rational surface B-Rep that does not includes weights.
	const MGSPointSeq& weights	//weights, weights(i,j,0) is for bcoef(i,j,.)
);

///Build MGRSBRep, given rib curves.

///Let v0=start parameter value, v1=terminate parameter value along v, then
///v=v0 const parameter line is curves[0], and v=v1 const parameter line is
///curves[n-1], where n=curves.size(). n must be greater or equal to 2.
///When n==2, the surface is a ruled surface(that is, this->order_u() is 2).
///
///If MGRLBRep's in vecPtrRibRLBReps may have different knot configurations,
///use the global function createSurfaceFromRibs(declared in MGSBRep.h).
void buildByRibRLBRep(
	const std::vector<const MGRLBRep*>& vecPtrRibRLBReps,
	bool direction_adjustment=true//=true, curves[.] direction are adjusted to line
								//to the same direction.
);

///Construct surface of revolution, given planar MGRLBRep and rotation axis sl.

///Parameterization of the surface is:
///	u=const parameter line generates given rlb(when u=0.).
/// v=const parameter line generates a circle whose center is sl.
void buildRevolutionSurface(
	const MGRLBRep& rlb,	///<Planar MGRLBRep to rotate.
	const MGStraight& sl,	///<Rotation axis. This is treated as infinite
			///<one, even if it is not.
	double angle			///<Rotation angle in radian,
	///<-2*pai<=angle<=2*pai,
	///<If angle is positive, circle is anti-clockwise around direction Vector N
	///<of sl. If negative, circle is clockwise around N.
);

///Construct MGRSBRep by sweep NURBS and sweep length.

//rail(say c(u)) is the rail and the straight line segments
//from C(u)+start_dist*uvec to C(u)+end_dist*uvec are the generatrix.
//The surface is expressed as: S(u,v)=c(u)+uvec*v,
//for rail.param_s()<=u<=rail.param_e(), start_dist<=v<=end_dist.
void buildSweep(
	const MGRLBRep& rlbrep,		///<Sweep crv.
	const MGUnit_vector& uvec,	///<Sweep Direction.
	double start_dist,			///<distance to start edge.
	double end_dist			///<distance to end edge.
);

///Assignment.
///When the leaf object of this and srf2 are not equal, this assignment
///does nothing.
MGRSBRep& operator=(const MGGel& gel2);
MGRSBRep& operator=(MGGel&& gel2);

/// 曲面の平行移動を行いオブジェクトを生成する。
///Translation.
MGRSBRep operator+ (const MGVector& ) const;

/// 曲面の逆方向に平行移動を行いオブジェクトを生成する。
///Translation.
MGRSBRep operator- (const MGVector& ) const;

/// 与えられたスケーリングで曲面の変換を行いオブジェクトを生成する。
///Scaling.
MGRSBRep operator* (double) const;

/// 与えられた変換で曲面の変換を行いオブジェクトを生成する。
///Matrix transformation.
MGRSBRep operator* (const MGMatrix& ) const;

/// 与えられた変換で曲面のトランスフォームを行いオブジェクトを生成する。
///General transformation.
MGRSBRep operator* (const MGTransf& ) const;

///Object transformation.
MGRSBRep& operator+=(const MGVector& v);
MGRSBRep& operator-=(const MGVector& v);
MGRSBRep& operator*=(double scale);
MGRSBRep& operator*=(const MGMatrix& mat);
MGRSBRep& operator*=(const MGTransf& tr);

///Comparison of two curves.
bool operator==(const MGRSBRep& gel2)const;
bool operator==(const MGGel& gel2)const;
bool operator<(const MGRSBRep& gel2)const;
bool operator<(const MGGel& gel2)const;
bool operator!=(const MGGel& gel2)const{return !(gel2==(*this));};
bool operator!=(const MGRSBRep& gel2)const{return !(gel2==(*this));};
bool operator==(const MGSBRep& sb)const;

//////////// Member Function ////////////

///Set(update) the knot vector, KV=MGKnotVector.
///This is move operation conforming.
template<class KVU, class KVV>
void setKnotVector(
	KVU&& tu,///<knot vector along u-direction
	KVV&& tv ///<knot vector along v-direction
){
	m_surface.setKnotVector(std::forward<KVU>(tu),std::forward<KVV>(tv));
};

///Gets new B-Rep by adding knots to an original B-Rep.
void addKnots(
	const MGKnotArray& uknots,	///<Knots to add for u-direction
	const MGKnotArray& vknots	///<Knots to add for v-direction.
);

///Returns B-Rep Dimension of u.
int bdim_u() const{return m_surface.bdim_u();}

///Returns B-Rep Dimension of v.
int bdim_v() const{return m_surface.bdim_v();}
	
///Compute minimum box that includes the surface.
MGBox box_limitted(const MGBox& bx) const;///Limited surface be the parameter box.

///Changing this object's space dimension.
void change_dimension(
	int sdim,	///< new space dimension
	int start1=0,///< Destination order of new object.
	int start2=0 ///< Source order of this object.
);

///Change parameter range, be able to change the direction by providing
///t1 greater than t2.
MGRSBRep& change_range(
	int is_u,		///<if true, (t1,t2) are u-value. if not, v.
	double t1,		///<Parameter value for the start of original. 
	double t2		///<Parameter value for the end of original. 
);

///Access to (i,j)th element of coef. Left-hand side version.
double& coef(int i, int j, int k)
{	assert(i<bdim_u() && j<bdim_v() && k<=sdim());
	return m_surface.coef(i,j,k);}

///Access to (i,j)th element of coef. right-hand side version.
double coef(int i, int j, int k)const{return m_surface.coef(i,j,k);}

///Extract (i,j,k) elements for 0<=k<sdim() as a vector.
MGVector coef(int i, int j) const{return m_surface.coef(i,j);}

///Returns a pointer to the surface b-coef data.
const double* coef_data(int i=0, int j=0, int k=0) const
{return m_surface.coef_data(i,j,k);}

///Construct new surface object by copying to newed area.
///User must delete this copied object by "delete".
MGRSBRep* clone() const;

///Return minimum box that includes whole of the surface.
///Compute the box from the scratch.
void compute_box(MGBox& bx) const;

///Construct new surface object by changing
///the original object's space dimension.
///User must delete this copied object by "delete".
MGRSBRep* copy_change_dimension(
	int sdim,			///< new space dimension
	int start1=0, 		///< Destination order of new line.
	int start2=0 		///< Source order of this line.
)const;

///Display control polygons using mgVBO::MGDrawPointSeq(sp)
void display_control_polygon(mgSysGL& sgl)const;

///uまたはv方向に折れ(マルチノット)があるとき面を分割する
///戻り値は、分割数を返却する
int divide_multi_knot(
    std::vector<UniqueSurface>& srfl	///Divided objects are appended.
) const;

///Evaluate right continuous ndu'th and ndv'th derivative data.

///Function's return value is (d(ndu+ndv)f(u,v))/(du**ndu*dv**ndv).
/// ndu=0 and ndv=0 means positional data evaluation.
MGVector eval(
	double u,///< U parameter value of the surface.
	double v,///< V parameter value of the surface.
	int ndu=0,	///< Order of Derivative along u.
	int ndv=0	///< Order of Derivative along v.
) const;

///Evaluate surface data.
MGVector eval(
	const MGPosition& uv///< Parameter value of the surface.
	, int ndu=0			///< Order of derivative along u.
	, int ndv=0			///< Order of derivative along v.
) const;

///Evaluate right continuous surface data.

///Evaluate all positional data and 1st and 2nd derivatives.
void eval_all(
	double u, 		///< U Parameter value of the surface.
	double v,		///< V Parameter value of the surface.
	MGPosition& f,			///< Positional data.
	MGVector&   fu,			///< df(u,v)/du
	MGVector&   fv,			///< df/dv
	MGVector&   fuv,		///< d**2f/(du*dv)
	MGVector&   fuu,		///< d**2f/(du**2)
	MGVector&   fvv			///< d**2f/(dv**2)
) const;

///Evaluate all of i and j'th derivative data for 0<=i<=ndu, 0<=j<=ndv.

/// Output: (d(i+j)f(u,v))/(du**i*dv**j) in deriv[r+j*dim+i*ndv*dim]
///for 0<=r<dim=sdim(), 0<=i<=nderiv and 0<=j<sdim(). 
void eval_all(
	double u, 		///< U Parameter value of the surface.
	double v,		///< V Parameter value of the surface.
	int ndu,	///<Order of Derivative along u.
	int ndv,	///<   along v direction.
	double* deriv	///<Output. (d(i+j)f(u,v))/(du**i*dv**j) in
					///<deriv[r+j*dim+i*(ndv+1)*dim] for 0<=r<dim=sdim().
					///<for 0<=i<=ndu and 0<=j<=ndv.
					///<deriv is an array of deriv[ndu+1][ndv+1][r].
) const;

/// Exchange parameter u and v.
MGSurface& exchange_uv(){m_surface.exchange_uv(); return *this;};

///Modify the original Surface by extrapolating the specified perimeter.

///The extrapolation is C2 continuous if the order >=4.
///The extrapolation is done so that extrapolating length is "length"
///at the position of the parameter value "param" of the perimeter.
MGRSBRep& extend(
	int perimeter,	///<perimeter number of the Surface,
					///< =0:v=min, =1:u=max, =2:v=max, =3:u=min.
	double param,	///< parameter value of above perimeter.
	double length,	///<chord length to extend at the parameter param of the perimeter.
	double dk=0.  ///<Coefficient of how curvature should vary at
///<    extrapolation start point. When dk=0, curvature keeps same, i.e.,
///<    dK/dS=0. When dk=1, curvature becomes zero at length extrapolated point,
///<    i.e. dK/dS=-K/length at extrapolation start point,
///<    (S=parameter of arc length, K=Curvature at start point)
///<    That is, when dk reaches to 1 from 0, curve changes to flat.
);

///Return homogeneous Surface B-Representation of the rational B-Spline.
const MGSBRep& homogeneous() const {return m_surface;};

/// Return This object's typeID
long identify_type() const;

///Test if input parameter value is inside parameter range of the surface.
bool in_range(double u, double v) const{return m_surface.in_range(u,v);};
bool in_range(const MGPosition& uv) const{return m_surface.in_range(uv);};

///Access to i-th element of u knot( left-hand side version).
double& knot_u(int i){return m_surface.knot_u(i);}

///Access to i-th element of u knot(right-hand side version).
double knot_u(int i) const{return m_surface.knot_u(i);}

///Access to i-th element of v knot(left-hand side version).
double& knot_v(int i){return m_surface.knot_v(i);}

///Access to i-th element of v knot(right-hand side version).
double knot_v(int i) const{return m_surface.knot_v(i);}

///Returns a pointer to the u knot vector data.
const double* knot_data_u() const{return m_surface.knot_data_u();}

///Returns a pointer to the v knot vector data.
const double* knot_data_v() const{return m_surface.knot_data_v();}

///Returns the u knot vector.
const MGKnotVector& knot_vector_u() const
{return m_surface.knot_vector_u();}
MGKnotVector& knot_vector_u(){return m_surface.knot_vector_u();}

///Returns the v knot vector.
const MGKnotVector& knot_vector_v() const
{return m_surface.knot_vector_v();}
MGKnotVector& knot_vector_v(){return m_surface.knot_vector_v();}

///Compare two parameter values. If uv1 is less than uv2, return true.
///Comparison is done after prjected to i-th perimeter of the surface.
bool less_than(
	int i,	///<perimeter number.
	const MGPosition& uv1,///<1st paramete.
	const MGPosition& uv2///< 2nd.
)const{return m_surface.less_than(i,uv1,uv2);};

///Obtain partial Surface B-Rep restricted by sub interval of u and v parameter range.
///Both u-range and v-range must be inside the range of the original.
///New one is exactly the same as the original except that it is partial.
void limit(const MGBox& uvrange){
	invalidateBox();
	m_surface.limit(uvrange);
}

///Change direction of the surface.
void negate(			
	int is_u		///< Negate along u-direction if is_u is ture,
					///< else along v-direction.
){m_surface.negate(is_u);}

///Obtain parameter value if this surface is negated by "negate()".
/// Negate along u-direction if is_u is ture,
/// else along v-direction.
MGPosition negate_param(const MGPosition& uv, int is_u=1)const
{ return m_surface.negate_param(uv,is_u);}

///Return non_homogeneous B-Coefficients with weights of
///the rational Surface B-Spline. This MGSPointSeq includes weights.
MGSPointSeq non_homogeneous_bcoef() const;

///Test if this is actually non_rational, i.e. , all of the weights are
///same values.
int non_rational() const;

///Returns the B-Rep order(u-direction).
int order_u() const{return m_surface.knot_vector_u().order();}

///Returns the B-Rep order(v-direction).
int order_v() const{return m_surface.knot_vector_v().order();}

/// Return ending parameter value.
MGPosition param_e() const{return m_surface.param_e();}
double param_e_u() const{return m_surface.param_e_u();}
double param_e_v() const{return m_surface.param_e_v();}

/// Compute parameter curve.
///Returned is newed area pointer, and must be freed by delete.
MGCurve* parameter_curve(
	int is_u			///<Indicates x is u-value if is_u is true.
	, double x			///<Parameter value.
						///<The value is u or v according to is_u.
) const;

/// Compute parameter line.
MGRLBRep parameter_line(
	int is_u			///<Indicates x is u-value if is_u is true.
	, double x			///<Parameter value.
						///<The value is u or v according to is_u.
) const;

/// パラメータ範囲を返す。
///Return parameter range.
MGBox param_range() const;

/// Return starting parameter value.
MGPosition param_s() const{return m_surface.param_s();}
double param_s_u() const{return m_surface.param_s_u();}
double param_s_v() const{return m_surface.param_s_v();}

///Compute part of the surface limitted by the parameter range bx.
///bx(0) is the parameter (us,vs) and bx(1) is (ue,ve).
///That is u range is from us to ue , and so on.
MGRSBRep* part(
	const MGBox& bx,///<Target parameter box.
	int multiple=0	///<Indicates if start and end knot multiplicities
					///<are necessary. =0:unnecessary, !=0:necessary.
)const;

///Retrieve perimeter i of this surface.
/// Compute perimeter Rational line B-Rep.
/// i is perimeter number:
/// =0: v=min line, =1: u=max line, =2: v=max line, =3: u=min line
MGRLBRep perimeter(int i) const;

///Return how many perimeters this surface has.
int perimeter_num() const{return 4;};

///Test if the RSBRep is planar or not.
///Returned is 0(false) if this is not planar, 1(true) if this planar.
int planar(
	MGPlane& plane,		///<Plane that might be closest to this.
						///<Plane is always output even if not planar.
	double& deviation	///<maximum deviation of this from the output plane.
) const;

///Test if part of the surface is planar or not within the tolerance tol.
///The part of the surface is input by the surface parameter range uvbox.
///Returned is 0(false) if this is not planar, 1(true) if planar.
int planar(
	const MGBox& uvbox,///<This surface parameter range.
	double tol,	///<maximum deviation allowed to regard the sub surface as a plane.
	int* divideU=0///<Direction to subdivide will be output, if this was not planar,
				///<=1: u direction, =0: v direction.
)const;

/// 入力パラメータをパラメータ範囲でまるめて返却する。
///Round the input parameter value uv into 
///the parameter range of the surface.
MGPosition range(const MGPosition& uv)const{return m_surface.range(uv);}

///Rebuild this MGRSBRep. Rebuild means:
/// Change the parameterization.
std::unique_ptr<MGSurface> rebuild(
	int how_rebuild=1,
		///< intdicates how rebuild be done.
		///<  =0: no approximation(only parameter change)
		///<  =1: Reconstructed with new knot configuration again as rational spline(MGRSBRep).
		///<  =2: approximated by non-rational spline(MGSBRep) with new knot configuration.
	int parameter_normalization=2,
		///<Indicates how the parameter normalization be done:
		///<  =0: no surface parameter normalization.
		///<  =1: normalize to u_range=(0., 1.), and v_range=(0.,1.);
		///<  =2: normalize to make the average length of the 1st derivative along u and v 
		///<     of the base surface is as equal to 1. as possible.
	double tol=-1.,	///<tolerance allowed for the approximation.
		///<When tol<=0., MGTolerance::line_zero() will be employed.
	int* order=0///<order of the new MGSBRep, >=4 is recomended.
		///<order[0]:u-order, [1]:v-order.
		///<When how_rebuild!=2, order is not used.
		///<When order=0 is input, order[0]=order[1]=4 are assumed.
)const;

///Change the B-Rep by decreasing B-Rep dimension by ndec. This is
///an approximation of the origimal B-Rep. Return value is error flag.
int reduce(
	int is_u,	///<if true, reduce b-rep dimension of u-direction.
	int ndec	///<Number of B-rep dimension to decrease .
){return m_surface.reduce(is_u,ndec);}

///Change an original B-Rep to new one with subdivided knot configuration.
///Knots t must be subdivided knots.
MGRSBRep& refine(
	const MGKnotVector& uknot,	///< new knot of u-direction
	const MGKnotVector& vknot	///< new knot of v-direction
){m_surface.refine(uknot,vknot); return *this;}

///ノット削除関数
///トレランスはline_zeroを使用する。元のノットが細かいものほど削除しやすい
///removal knot. line_zero tolerance is used.
void remove_knot();

///Obtain the partial Surface B-Rep restricted by sub interval of u and v parameter range.

///New one is exactly the same as the original except that it is partial.
///If multiple==true(!=0), knot_u(i)=t1 and knot_u(n+i)=t2 for i=0,..., k-1
///will be guaranteed. Here, n=bdim_u(), k=order_u(),
///t1=uvrange(0).low_point(), and t2=uvrange(0).high_point().
///About knot_v(j), the same.
///Both u-range and v-range must be inside the range of this.
void shrinkToParameters(
	const MGBox& uvrange,	///<u and v parameter range.
	MGRSBRep& newBrep,///<Shrinked surface is output, which can be this.
	int multiple=0 ///<Indicates if start and end knot multiplicities
					///<are necessary. =0:unnecessary, !=0:necessary.
)const;

///Returns the space dimension.
int sdim() const{return m_surface.sdim()-1;}

///Shrink this surface to the part limitted by the parameter range of uvbx.

///New parameter range uvbx2 is so determined that uvbx2 is the smallest
///box tha includes uvbx, and all of the u or v values of uvbx2 is one of 
///the values of u or v knots of the surface knotvector.
///uvbx(0) is the parameter (us,ue) and uvbx(1) is (vs,ve).
///That is u range is from us to ue , and so on.
void shrink_to_knot(
	const MGBox& uvbx,///<The target parameter box.
	int multiple=0	///<Indicates if start and end knot multiplicities
					///<are necessary. =0:unnecessary, !=0:necessary.
){m_surface.shrink_to_knot(uvbx,multiple);};

///Returns the B-coef's.
///Right hand side version.
const MGSPointSeq& surface_bcoef() const{return m_surface.surface_bcoef();}	

///Returns the B-coef's(Left hand side version).
MGSPointSeq& surface_bcoef(){return m_surface.surface_bcoef();}

///Compute surface integral of the 1st two coordinates.

///This integral can be used to compute volume sorounded by the surface.
///double surface_integral(const MGBox&) const;

///Return the surface type.
MGSURFACE_TYPE type() const{return MGSURFACE_TYPE::MGSURFACE_RSPLINE;}

///Unlimit the parameter range. Return the same.
MGSurface& unlimit(){return *this;};

public:

///output to IGES, PD128
int out_to_IGES(
	MGIgesOfstream& igesfile,
	int SubordinateEntitySwitch=0
)const;

///Debug Function
std::ostream& toString(std::ostream&) const;

///Get the name of the class.
std::string whoami()const{return "RSBRep";};

protected:

///Intersection of Surface and a straight line.
MGCSisects isectSl(
	const MGStraight& sl,///<Target straight.
	const MGBox& uvbox=mgNULL_BOX ///<indicates if this surface is restrictied to the parameter
					///<range of uvbox. If uvbox.is_null(), no restriction.
)const override;

///メンバデータを読み込む関数
/// 戻り値boolは正常に読み出しが出来ればtrue、失敗すればfalseになる
void ReadMembers(MGIfstream& buf);
	
///メンバデータを書き込む関数
/// 戻り値boolは正常に書き込みが出来ればtrue、失敗すればfalseになる
void WriteMembers(MGOfstream& buf) const;

private:

//////////////Member Data//////////////
	MGSBRep m_surface;		/// Maximum space dimension id is for weights.

///Obtain coefficient's space dimension.
///This function is used in isect_start etc.
int coef_sdim() const{return m_surface.sdim();};

//Convert th bspline coefficients to homogeneous form.
void convertToHomogeneous();
///The following two function will be used in perps or isect
///to decide how many division of the surface along u or v direction
///should be applied before using perp_guess or isect_guess.
int intersect_dnum_u() const;
int intersect_dnum_v() const;

///"isect1_incr_pline" is a dedicated function of isect_start_incr, will get
/// shortest parameter line necessary to compute intersection.
MGCurve* isect_incr_pline(
	const MGPosition& uv,///<last intersection point.
	int kdt,			///<Input if u=const v-parameter line or not.
						///< true:u=const, false:v=const.
	double du, double dv,///<Incremental parameter length.
	double& u,			///<next u value will be output
	double& v,			///<next v value will be output
	int incr=1		///<Incremental valuse of B-coef's id.
) const;

///Return order of intersection line order of MGLBRep.
int isect_order() const;

///Obtain 1D surface rep. of this surf which can be used for
///isect(const MGPlane& pl). This surf1D is used in isect for
///the argument of isect_startPlane, which will use surf1D to compute isect(pl).
///surf1D=0.(intersection with x=0. plane) is the intersection lines.
std::unique_ptr<MGSBRep> surf1D(const MGPlane& pl)const override;

///u方向に折れ(マルチノット)があるとき面を分割する
///戻り値は、分割数を返却する
int divide_multi_knot_u(
    std::vector<UniqueRSBRep>& srfl) const;  ///Divided MGRSBRep are appended.

///v方向に折れ(マルチノット)があるとき面を分割する
///戻り値は、分割数を返却する
int divide_multi_knot_v(
	std::vector<UniqueRSBRep>& srfl	///Divided objects are appended.
) const;

};

namespace MGCL{
///Compute binominal coefficients.
	
///Let i=degree, then bc(i,j) contains j-th coefficient of the degree i.
///bc(i,j) for 0<=i<=m and 0<=j<=i in bc[(m+1)*i+j].
///bc is an arrary of length (m+1)*(m+1).
///m must be greater than or equal to 1.
extern MG_DLL_DECLR void Binominal(int m, double* bc);

};

/** @} */ // end of GEO group
#endif
