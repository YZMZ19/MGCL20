/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#ifndef _MGBSumSurf_HH_
#define _MGBSumSurf_HH_

class MGIfstream;
class MGOfstream;
#include <iosfwd>
#include "mg/MGCL.h"
#include "mg/CSisects.h"
#include "mg/Surface.h"

// MGBSumSurf.h

/** @addtogroup GEO
 *  @{
 */

///Defines Boolean sum surface.

///Boolean sum surface is defined by three surfaces g1, g2, and g12 whose parameter
///ranges are exactly the same, as f(u,v)=g1+g2-g12. Typically Gordon surface is
///a boolean sum surface. See "Curves and Srufaces for CAGD" by Gerald Farin.
class MG_DLL_DECLR MGBSumSurf: public MGSurface{
 
public:


////////Special member functions/////////
MGBSumSurf():MGSurface(),m_g1(nullptr), m_g2(nullptr), m_g12(nullptr){;};
~MGBSumSurf();
MGBSumSurf(const MGBSumSurf&);///<Copy constructor.
MGBSumSurf& operator= (const MGBSumSurf&);///<Copy assignment.
MGBSumSurf(MGBSumSurf&&);		///<Move constructor.
MGBSumSurf& operator= (MGBSumSurf&&);///<Move assignment.

///Default constructor.

///construct from newed MGSurface. The ownership of g1, g2, and g12 will be
///transfered to this MGBSumSurf.
MGBSumSurf(
	MGSurface* g1,
	MGSurface* g2,
	MGSurface* g12
);

///construct from three MGSurface.
MGBSumSurf(
	const MGSurface& g1,
	const MGSurface& g2,
	const MGSurface& g12
);

									
//////////// Operator overload. 演算子多重定義 ////////////

///Assignment.
///When the leaf object of this and srf2 are not equal, this assignment
///does nothing.
MGBSumSurf& operator=(const MGGel& gel2);
MGBSumSurf& operator=(MGGel&& gel2);

///Object transformation.
MGBSumSurf& operator+=(const MGVector& v);
MGBSumSurf& operator-=(const MGVector& v);
MGBSumSurf& operator*=(double scale);
MGBSumSurf& operator*=(const MGMatrix& mat);
MGBSumSurf& operator*=(const MGTransf& tr);

///Comparison of two curves.
bool operator==(const MGBSumSurf& gel2)const;
bool operator==(const MGGel& gel2)const;
bool operator<(const MGBSumSurf& gel2)const;
bool operator<(const MGGel& gel2)const;
bool operator!=(const MGGel& gel2)const{return !(gel2==(*this));};
bool operator!=(const MGBSumSurf& gel2)const{return !(gel2==(*this));};

//////////// Member Function ////////////

int bdim_u() const;	///<Returns B-Rep Dimension of u.
int bdim_v() const;	///<Returns B-Rep Dimension of v.

/// 入力のパラメータ範囲の曲線部分を囲むボックスを返す。
///Return minimum box that includes limitted surface by uvrange.
MGBox box_limitted(
	const MGBox& uvrange	///< Parameter Range of the curve.
)const;

///Changing this object's space dimension.
void change_dimension(
	int sdim,	///< new space dimension
	int start1=0,///< Destination order of new object.
	int start2=0///< Source order of this object.
);

///Change parameter range, be able to change the direction by providing
///t1 greater than t2.
MGBSumSurf& change_range(
	int is_u,	///<if true, (t1,t2) are u-value. if not, v.
	double t1,	///<Parameter value for the start of original. 
	double t2	///<Parameter value for the end of original. 
);

///Construct new surface object by copying to newed area.
///User must delete this copied object by "delete".
MGBSumSurf* clone() const;

///Return minimum box that includes whole of the surface.
void compute_box(MGBox& bx) const;

///Construct new surface object by changing
///the original object's space dimension.
///User must delete this copied object by "delete".
MGBSumSurf* copy_change_dimension(
	int sdim,	///< new space dimension
	int start1=0,///< Destination order of new line.
	int start2=0 ///< Source order of this line.
)const;

///Evaluate surface data.
///Currently ndu=ndv=0 is assumed.
MGVector eval(
	double u,///< Parameter value of the surface.
	double v ///< must be 0<=u,v<=1.
	, int ndu=0		///< Order of derivative along u.
	, int ndv=0		///< Order of derivative along v.
) const;

///Evaluate surface data.
MGVector eval(
	const MGPosition& uv	///< Parameter value of the surface.
	, int ndu			///< Order of derivative along u.
	, int ndv			///< Order of derivative along v.
) const{return eval(uv[0],uv[1],ndu,ndv);}

/// Exchange parameter u and v.
MGSurface& exchange_uv();

/// Return This object's typeID
long identify_type() const{return MGBSUMSURF_TID;};

bool in_range(double u, double v) const;

///Surface and Surface intersection.

///Return order of intersection line order of MGLBRep.
///The default is 4.
int isect_order() const;

///Access to i-th element of u knot
double knot_u(int i) const;

///Access to i-th element of v knot
double knot_v(int i) const;

///Returns the u knot vector.
const MGKnotVector& knot_vector_u() const;
MGKnotVector& knot_vector_u();

///Returns the v knot vector.
const MGKnotVector& knot_vector_v() const;
MGKnotVector& knot_vector_v();

///Compare two parameter values. If uv1 is less than uv2, return true.
///Comparison is done after projected to i-th perimeter of the surface.
bool less_than(
	int i,	///<perimeter number.
	const MGPosition& uv1,///<1st parameter,
	const MGPosition& uv2///<2nd parameter.
)const{return m_g1->less_than(i,uv1,uv2);};

///Negate direction of surface.
void negate(
	int is_u///< Negate along u-direction if is_u is ture,
			///< else along v-direction.
);

///Returns the order of u.
int order_u() const;

///Returns the order of v.
int order_v() const;

///Output to IGES stream file.
///MGBSumCurve is approximated as MGSBRep and output as PD128.
///int out_to_IGES(MGIgesOfstream& igesfile)const;

///Debug Function
/// Output function.
std::ostream& toString(std::ostream& ostrm) const;

/// Compute parameter curve.
///Returned is newed area pointer, and must be freed by delete.
MGCurve* parameter_curve(
	int is_u	///<Indicates x is u-value if is_u is true.
	, double x	///<Parameter value.
				///<The value is u or v according to is_u.
)const;

/// Return ending parameter value.
double param_e_u() const;
double param_e_v() const;

/// パラメータ範囲を返す。
///Return parameter range.
MGBox param_range() const;

/// Return starting parameter value.
double param_s_u() const;
double param_s_v() const;

///Compute part of the surface limitted by the parameter range bx.
///bx(0) is the parameter (us,vs) and bx(1) is (ue,ve).
///That is u range is from us to ue , and so on.
///Retured is newed object, must be deleted.
MGBSumSurf* part(
	const MGBox& bx,///<Target parameter range.
	int multiple=0	///<Indicates if start and end knot multiplicities
					///<are necessary. =0:unnecessary, !=0:necessary.
)const;

///Return how many perimeters this surface has.
int perimeter_num() const{return 4;};

///Get space dimension.
int sdim()const{return m_g1->sdim();};

///Shrink this surface to the part limitted by the parameter range of uvbx.
///New parameter range uvbx2 is so determined that uvbx2 is the smallest
///box tha includes uvbx, and all of the u or v values of uvbx2 is one of 
///the values of u or v knots of the surface knotvector.
///uvbx(0) is the parameter (us,ue) and uvbx(1) is (vs,ve).
///That is u range is from us to ue , and so on.
void shrink_to_knot(
	const MGBox& uvbx,///<Target parameter range.
	int multiple=0	///<Indicates if start and end knot multiplicities
					///<are necessary. =0:unnecessary, !=0:necessary.
);

/// 面のTypeを返却する。
///Return the surface type.
MGSURFACE_TYPE type() const{return MGSURFACE_TYPE::MGSURFACE_BSUM;}

///メンバデータを読み込む関数
void ReadMembers(MGIfstream& buf);
	
///メンバデータを書き込む関数
void WriteMembers(MGOfstream& buf) const;

///Get the name of the class.
std::string whoami()const{return "BSumSurf";};

//////////// Member Data ////////////

private:

	MGSurface* m_g1;
	MGSurface* m_g2;
	MGSurface* m_g12;

///The following two function will be used in perps or isect
///to decide how many division of the surface along u or v direction
///should be applied before using perp_guess or isect_guess.
int intersect_dnum_u() const;
int intersect_dnum_v() const;

///"isect1_incr_pline" is a dedicated function of isect_start_incr, will get
/// shortest parameter line necessary to compute intersection.
MGCurve* isect_incr_pline(
	const MGPosition& uv,	///<last intersection point.
	int kdt,				///<Input if u=const v-parameter line or not.
							///< true:u=const, false:v=const.
	double du, double dv,///<Incremental parameter length.
	double& u,				///<next u value will be output
	double& v,				///<next v value will be output
	int incr=0			///<Incremental valuse of B-coef's id.
)const;

///Obtain 1D surface rep. of this surf which can be used for
///isect(const MGPlane& pl). This surf1D is used in isect for
///the argument of isect_startPlane, which will use surf1D to compute isect(pl).
///surf1D=0.(intersection with x=0. plane) is the intersection lines.
std::unique_ptr<MGSBRep> surf1D(const MGPlane& pl)const override;

};

/** @} */ // end of GEO group
#endif
