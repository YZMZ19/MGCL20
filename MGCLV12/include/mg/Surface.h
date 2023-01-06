/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/

#ifndef _MGSurface_HH_
#define _MGSurface_HH_
#include <memory>
#include <vector>
#include "mg/Default.h"
#include "mg/Geometry.h"
#include "mg/Position_list.h"
#include "mg/FSurface.h"
#include "mg/SSisects.h"

// MGSurface.h

class MGDrawParam;
class MGBPointSeq;
class MGPosition;
class MGVector;
class MGTransf;
class MGUnit_vector;
class MGCurve;
class MGEllipse;
class MGStraight;
class MGCompositeCurve;
class MGSurfCurve;
class MGLBRep;
class MGSPointSeq;
class MGCSisects;
class MGPlane;
class MGSBRep;
class MGSSisect;
class MGSSisects;
class MGPosition_list;
class MGIfstream;
class MGOfstream;
class MGRLBRep;
class MGCParam_list;
class MGNDDArray;
class MGFace;
class MGShell;
class MGHHisects;
class mgTL2Triangles;

/** @addtogroup GEO
 *  @{
 */

/// MGSurface is an abstract class of 3D surface.

/// Surface is represented using two parameter u and v:f(u,v).
class MG_DLL_DECLR MGSurface: public MGGeometry, public MGFSurface{

public:


////////Special member functions/////////
MGSurface(void)=default;///Void Constructor. 
virtual ~MGSurface()=default;			///Destructor.
MGSurface(const MGSurface& ssi)=default;///Copy constructor.
MGSurface& operator= (const MGSurface& rhs)=default;///Copy assignment.
MGSurface(MGSurface&& ssi)=default;		///Move constructor.
MGSurface& operator= (MGSurface&& rhs)=default;///Move assignment.

///Object transformation.
MGSurface& operator+=(const MGVector& v)=0;
MGSurface& operator-=(const MGVector& v)=0;
MGSurface& operator*=(double scale)=0;
MGSurface& operator*=(const MGMatrix& mat)=0;
MGSurface& operator*=(const MGTransf& tr)=0;

///comparison
virtual bool operator==(const MGGel& gel2)const=0;
virtual bool operator<(const MGGel& gel2)const=0;

///Return space dimension
int sdimFS() const override{	return sdim();};

///Approximate this surface as an MGSBRep.
std::unique_ptr<MGSBRep> approximate_as_SBRep(
	int parameter_normalization=2,
		///<Indicates how the parameter normalization be done:
		///<  =0: no surface parameter normalization.
		///<  =1: normalize to u_range=(0., 1.), and v_range=(0.,1.);
		///<  =2: normalize to make the average length of the 1st derivative along u and v 
		///<      of the base surface is as equal to 1. as possible.
	double tol=-1.,	///<tolerance allowed for the approximation.
		///When tol<=0., MGTolerance::line_zero() will be employed.
	int* order=0///<order of the new MGSBRep, >=4 is recommended.
		///<order[0]:u-order, [1]:v-order.
		///<When order=0 is input, order[0]=order[1]=4 are assumed.
)const;

///Compute average chord length along u(is_u==false) or v(is_u==true).

///Function's return value is the chord length.
double average_chord_length(
	int is_u,///< =0, or 1. Indicates if para is u-value(is_u=1) or v(is_u=0).
	const double para[3],	///<three parameter value of srf, start, middle, and end
							///<of u parameter value if is_u==1, of v if is_u==0.
	const MGNDDArray& tau	///<data points to evaluate on is_u parameter line.
							///<tau[.] is v values if is_u==1, u values if 0.
)const;

///Generate arrow data of the tangent along u and v and the normal.

///Arrow at the parameter value (u,v) of the surface.
///data[0] is the origin of the u-tangent arrow, data[1] is the top of the u-tangent arrow,
///data[2], [3] are two bottoms of u-tangent arrowhead.
///data[0], [4], [5], [6] are the points of v-tangent arrow.
///data[0], [7], [8], [9] are the points of v-tangent arrow.
void arrow(double u,double v, MGPosition data[10])const;
void arrow(const MGPosition& uv, MGPosition data[10])const{
	arrow(uv[0],uv[1],data);
}

///Generate arrow data, given box. The length of the arrows are defined from box.len().
void arrow(const MGBox& box, double u,double v, MGPosition data[10])const;

///Returns B-Rep Dimension of u.
virtual int bdim_u() const {return 1;}

///Returns B-Rep Dimension of v.
virtual int bdim_v() const {return 1;}

///Return minimum box that includes limitted surface by uvrange.

/// ���͂̃p�����[�^�͈͂̋Ȑ��������͂ރ{�b�N�X��Ԃ�.
virtual MGBox box_limitted(
	const MGBox& uvrange	///< Parameter Range of the curve.
) const=0;

///Return box of the parameter space of the surface.
MGBox box_param() const;

///Return box of the parameter space of the FSurface after trimmed one.
const MGBox box_param2() const{return box_param();};

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

///Change parameter range, able to change the direction by providing t1 greater than t2.
virtual MGSurface& change_range(
	int is_u,				///<if true, (t1,t2) are u-value. if not, v.
	double t1,				///<Parameter value for the start of original. 
	double t2				///<Parameter value for the end of original. 
)=0;


///Obtain coefficient's space dimension.

///This function is used in isect_start etc.
virtual int coef_sdim() const{return sdim();};

///Get the clone of this MGFSurface.
MGFSurface* clone_fsurface()const override{return copy_surface();};

///Get the clone of this as a MGFace.

///If this is MGSurface, it is converted to MGFace.
MGFace* clone_as_face()const;

///Compute the closest point parameter value (u,v)of this surface from a point.
virtual MGPosition closest(const MGPosition& point) const;


///Compute the closest point on all the perimeters of the surface.

///The point is returned as the parameter value (u,v) of this surface.
virtual MGPosition closest_on_perimeter(const MGPosition& point)const;
virtual MGPosition closest_on_perimeter(const MGStraight& sl)const;

///Compute closest point from a line to the boundary of the MGFSurface.

///Returned is the parameter value of the FSurface that is closest to point.
virtual MGPosition closest_on_boundary(const MGStraight& sl) const{
	return closest_on_perimeter(sl);
}


///Construct new surface object by copying to newed area.

///User must delete this copied object by "delete".
virtual MGSurface* clone() const override =0;


///compute sample point of the surface to get the approximate plane.

///The approximate plane is of the surface
///within the parameter range (u0,v0) to (u1, v1).
void compute_sample_point(
	double u0,///<Parameter range.
	double u1,///<Parameter range.
	double v0,///<Parameter range.
	double v1,///<Parameter range.
	MGPosition Pn[9],	///<sample points will be output.
	MGPosition& center,	///<center of the sample points will be output.
	MGUnit_vector& normal,///<average normal of Nn[] will be output.
	MGVector* Nn_in=0		///<9 normals of the surface will be output when Nn_in specified.
)const;


///Construct new surface object by changing the original object's space dimension.

///User must delete this copied object by "delete".
virtual MGSurface* copy_change_dimension(
	int sdim,			///< new space dimension
	int start1=0, 		///< Destination order of new line.
	int start2=0 		///< Source order of this line.
)const=0;


///Compute surface curvatures.

///value[0]=K:Gaussian curvature=k1*k2, value[1]=H:Mean curvature=(k1+k2)/2,
///value[2]=k1:minimum curvature, and value[3]=k2=maximum curvature.
///N is the unit normal vector at position (u,v).
void curvatures(
	const MGPosition& uv, double value[4], MGUnit_vector& N) const;
void curvatures(
	double u, double v, double value[4], MGUnit_vector& N) const;

///Compute direction unit vector of the geometry.
MGUnit_vector direction(const MGPosition& param) const;


///Draw 3D curve in world coordinates.

///The object is converted to curve(s) and is drawn.
virtual void drawWire(
	mgVBO& vbo,///<The target graphic object.
	int line_density=1	///<line density to draw a surface in wire mode.
)const{drawWireFS(vbo,line_density);};


///Construct new curve object by copying to newed area.

///User must delete this copied object by "delete".
MGSurface* copy_surface() const;

///display function.
virtual void display_arrows(mgSysGL& sgl)const;

///u�܂���v�����ɐ܂�(�}���`�m�b�g)������Ƃ��ʂ𕪊�����
///�߂�l�́A��������ԋp����
virtual int divide_multi_knot(
	std::vector<UniqueSurface>& srfl///Divided objects are appended.
)const;


///Compute if MGSurfCurve scurve(*this, param_curve) has the same direction to world_curve.

///Assuming that scurve and world_curve are the same curve.
///Function's return value is:
///1: same direction, -1:oppositie direction.
int equal_direction(
	const MGCurve& param_curve,	///<(u,v) parameter representation curve of this.
	const MGCurve& world_curve	///<world representation curve.
)const;

///Evaluate surface data.
virtual MGVector eval(
	double u,	///< U Parameter value of the surface.
	double v		///< V Parameter value of the surface.
	, int ndu=0			///< Order of derivative along u.
	, int ndv=0			///< Order of derivative along v.
) const=0;
MGVector eval(
	const MGPosition& uv,		///<FSurface parameter value(u,v)
	int ndu=0			///< Order of derivative along u.
	, int ndv=0			///< Order of derivative along v.
)const{return eval(uv.ref(0),uv.ref(1),ndu,ndv);}


///Evaluate all the points (ui, vj) into spoint(i,j,.).

///Here ui=utau(i) for 0<=i<utau.length() and vj=vtau(j) for 0<=j<vtau.length().
virtual void eval_spoint(
	const MGNDDArray&	utau,	///<u�����̃f�[�^�|�C���g
	const MGNDDArray&	vtau,	///<v�����̃f�[�^�|�C���g
	MGSPointSeq&		spoint	///<evaluated data will be output to spoint.
)const;

///Evaluate right continuous surface data.

///Evaluate all positional data, 1st and 2nd derivatives.
virtual void eval_all(
	double u,	///< U Parameter value of the surface.
	double v,	///< V Parameter value of the surface.
	MGPosition& f,		///< Positional data.
	MGVector&   fu,		///< df(u,v)/du
	MGVector&   fv,		///< df/dv
	MGVector&   fuv,	///< d**2f/(du*dv)
	MGVector&   fuu,	///< d**2f/(du**2)
	MGVector&   fvv		///< d**2f/(dv**2)
)const;

///Evaluate right continuous surface data.

///Evaluate all positional data, 1st and 2nd derivatives.
virtual void eval_all(
	const MGPosition& uv,	///< Parameter value of the surface.
	MGPosition& f,			///< Positional data.
	MGVector&   fu,			///< df(u,v)/du
	MGVector&   fv,			///< df/dv
	MGVector&   fuv,		///< d**2f/(du*dv)
	MGVector&   fuu,		///< d**2f/(du**2)
	MGVector&   fvv			///< d**2f/(dv**2)
	)const;

///evaluate gap between this surface's perimeter iperi and the given curve curve.
double eval_gap(
	const MGCurve& curve, ///< (I/ ) curve to evaluate.
	int            iperi, ///< (I/ ) 0: vmin, 1: umax, 2: vmax, and 3: umin. 
	MGPosition&    uv     ///< ( /O) the parameter of this that had the largest gap.
)const;

///evaluate gap between this surface's perimeters and the given curve curve.

///evaluation is performed for the perimeter i and curve[i] for 0<=i<=4.
///function's return value is the maximum gap.
double eval_gap(
	const MGCurve* curve[4], ///< (I/ ) curves to evaluate.
	MGPosition&    uv       ///< ( /O) the parameter of this that had the largest gap.
)const;

/// Evaluate n'th derivative data. n=0 means positional data evaluation.
MGVector evaluate(
const MGPosition& t,///< Parameter value, t's space dimension is geometry's manifold dimension.
const int* nderiv	///<Order of derivative of i-th parameter in nderiv[i],
			///<When nderiv=null, nderiv[i]=0 is assumed for all i.
)const;

/// Exchange parameter u and v.
virtual MGSurface& exchange_uv()=0;

///Modify the original Surface by extrapolating the specified perimeter.

///The extrapolation is C2 continuous if the order >=4.
///The extrapolation is done so that extrapolating length is "length"
///at the position of the parameter value "param" of the perimeter.
virtual MGSurface& extend(
	int perimeter,	///<perimeter number of the Surface.
					///< =0:v=min, =1:u=max, =2:v=max, =3:u=min.
	double param,	///< parameter value of above perimeter.
	double length,	///<chord length to extend at the parameter param of the perimeter.
	double dk=0.  ///<Coefficient of how curvature should vary at
///<    extrapolation start point. When dk=0, curvature keeps same, i.e.
///<    dK/dS=0. When dk=1, curvature becomes zero at length extrapolated point,
///<    i.e. dK/dS=-K/length at extrapolation start point
///<    (S=parameter of arc length, K=Curvature at start point)
///<    That is, when dk reaches to 1 from 0, curve changes to flat.
){return *this;};

///Get the MGFSurface pointer if this is MGSurface or MGFace.
const MGFSurface* fsurface()const{return this;};
MGFSurface* fsurface(){return this;};

///Compute 1st and 2nd fundamental quantities of the surface.

///In Q, 1st and 2nd fundamental quantities are returned as:
///Q[0]=E=fu%fu, Q[1]=F=fu%fv, Q[2]=G=fv%fv,
///Q[3]=L=fuu%UN, Q[4]=M=fuv%UN=fvu%UN, Q[5]=N=fvv%UN.
void fundamentals(
	const MGPosition&uv,	///<Surface parameter value (u,v)
	double Q[6],			///<Fundamental quantities will be returned.
	MGUnit_vector& UN		///<Normal vector at uv will be returned.
)const;
void fundamentals(
	double u,	///<Surface parameter u value of (u,v)
	double v,	///<Surface parameter v value of (u,v)
	double Q[6],			///<Fundamental quantities will be returned.
	MGUnit_vector& N		///<Normal vector at uv will be returned.
)const;

///Compute the approximate plane in the parameter range from (u0, v0) to (u1,v1).

///The plane's origin is center point of the plane when the surface is mapped
///onto the plane. The uderiv of the plane is the direction from the point(u0, v0) to (u1,v0).
void get_approximate_plane(
	double u0	///<u range from u0 to u1.
	,double u1,	///<u range from u0 to u1.
	double v0,///<v range from v0 to v1.
	double v1,///<v range from v0 to v1.
	MGPlane& plane,		///<The plane will be output.
	double* width=0,	///<The width and the height of the plane that include all the data
	double* height=0	///<for the surface point to map onto the plane will be output
)const;

///get face pointer if this is MGFace, else null will be returned.
MGFace* get_face_pointer(){return (MGFace*)0;};
const MGFace* get_face_pointer()const{return (const MGFace*)0;};

///Compute common curve part of this surface's perimeter and the crv.

///Function's returned value is the number of common part curve part,
///which is 2 at most.
int getPerimeterCommon(
	const MGCurve& crv,///<crv must be a cotinuous one curve, must not be MGCompositeCurve.
	std::vector<double> pspan[2],///<parameter range of crv and perimeters will be output.
	int peri_num[2]	///<perimeter number of pspan[i] will be output in peri_num[i],
		///<(pspan[i], peri_num[i]) is one pair.
)const;

///get surface pointer. Null will never be returned if this is valid MGFSurface.

///That is, if this is MGFace, base surface will be returned.
MGSurface* get_surface_pointer(){return this;};
const MGSurface* get_surface_pointer()const{return this;};

///Get number of inner boundaries as the output of the function.
virtual int get_number_of_boundaries()const{return 0;};

///Given world curve wcrv on this face, get the parameter space representation pcrv.

///Function's return value is pcrv, which is newed one. Must be deleted.
MGCurve* get_parameterCurve(const MGCurve& wcrv)const;

///Given MGSBRep or MGRSBRep as srf, compute normalize MGKnotVector along u and v.
void get_new_surface_knots(
	int parameter_normalization,
		//Indicates how the parameter normalization be done:
		//=0: no surface parameter normalization.
		//=1: normalize to u_range=(0., 1.), and v_range=(0.,1.);
		//=2: normalize to make the average length of the 1st derivative along u and v 
		//    of the base surface is as equal to 1. as possible.
	MGKnotVector& uknots,//normalized MGKnotVector will be returned.
	MGKnotVector& vknots,
	double* Oldparameter=0//Old parameterization will be returned. [.]=(us, ue, vs,ve).
)const;


///Test if input (u,v) is parameter value on a perimeter of the surface.

///If u or v is on a perimeter, (u,v) will be updated to the perimeter value.
virtual bool on_a_perimeter(
	double& u,	///<Surface parameter u of (u,v)
	double& v,	///<Surface parameter v of (u,v)
	int& perim_num///<if function returns true,the perimete number will be output,
		///<If function returns false, the nearest perimeter number will be output.
)const;

///Test if input x is parameter value on a perimeter of the surface.

///If x is on a perimeter, x will be updated to the perimeter value.
///Function's return value is true if on a perimeter.
bool on_a_perimeter2(
	int is_u,	///<specify if x is u or v value, is_u!=0(true) means u value.
	double& x,	///<Surface parameter (u,v)
	int& perim_num///<if function returns true,the perimete number will be output.
)const;

///Test if this and 2nd object has common area about their box(),taking error into account.
virtual bool has_commonFS(const MGObject& obj2) const{return has_common(obj2);};

///Shade the object in world coordinates.
virtual void shade(
	mgVBO& vbo,
	const MGDrawParam& para,
	MGCL::DRAW_TARGET target= MGCL::SHADING
)const;

///Triangulate this object(MGShell, MGFace, or MGSurface is the target).
virtual void triangulate(
	const MGDrawParam& para,
	MGCL::TL_DATA_KIND dkind,
	std::vector<mgTL2Triangles>& trisVec
)const;

///Compute the approximate plane in the parameter range from (u0, v0) to (u1,v1).

///Also, test if the surface is within surface_tol and angle from the plane.
///The plane's origin is center point of the plane when the surface is mapped�@onto the plane.
///the uderiv is the direction from the point(u0, v0) to (u1,v0).
///Function's return value is true when the surface is within the tolerance surface_tol,
///and the surface normals are within angle from the plane's normal
///plane, width, and height are valid only when function's return value is true.
bool test_and_get_approximate_plane(
	double u0,///<u range from u0.
	double u1,///< to u1.
	double v0,///<v range from v0.
	double v1,///< to v1.
	double surface_tol,	///<tolerance allowed for the deviation from the plane to the surface.
	double angle,		///<angle allowed for the normal of the plane and the normals of the
						///<surface.
	MGPlane& plane,		///<The plane will be output.
	double& width,	///<The width and the height of the plane that include all the data.
	double& height	///<for the surface point to map onto the plane.
)const;

/// Return This object's typeID.
virtual long identify_type()const=0;

///Test if the parameter(u,v) is in this surface's parameter range.
virtual bool in_range(double u, double v)const=0;
bool in_range(const MGPosition& uv)const{return in_range(uv[0], uv[1]);};

///Test if (u,v) is inside the face.

///Function's return value is:
///  0:outside the face.
///  1:unknown.
///  2:inside the face, not on a boundary.
///  <0:(u,v) is on an inner boundary, and abs(return code) is the loop id.
///  4:(u,v) is on the outer boundary.
///  >=10: (u,v) is on a perimeter, (10+perimeter number) will be returned.
int in_range_with_on(const MGPosition& uv)const;

///Test if surface limitted by the parameter range bx is flat and small.

///That is, surface is flat within surftol from the average plane,
///and all of the egedes are small compared with melen2.
bool is_flat_and_small(
	const MGBox& bx,//Paramete range of the surface.
	double surftol,	//Input the maximum deviation allowed from a plane.
	double melen2	//square of maximum edge length allowed of the surface edge.
					//When melen2<=0., no edge length check done.
)const;

///Obtain i-th inner_boundary curves(world coordinates representation) of the FSurface.

///Let the output of inner_boundary(i) be wcurves and
///of inner_boundary_param(i) be pcurves, then wcurves[j] corresponds
///to pcurves[j] one to one. Number of inner_boundary can be obtained
///by the function number_of_inner_boundary().
virtual std::vector<UniqueCurve> inner_boundary(int i)const;

///Obtain i-th inner_boundary curves.

///Obtain i-th inner_boundary curves(world coordinates representation)
///of the FSurface. Let the output of inner_boundary(i) be wcurves and
///of inner_boundary_param(i) be pcurves, then wcurves[j] corresponds
///to pcurves[j] one to one. Number of inner_boundary can be obtained
///by the function number_of_inner_boundary().
virtual std::vector<UniqueCurve> inner_boundary_param(int i)const;

///Default surface-curve intersection function.

///Restriction for this surface and curve:
///	1. this surface and curve must not have C0 continuity in it.
///	2. param_s() and param_e()  of this surface and curve must return
///     real start and end parameter values. That is, must be finite.
MGCSisects intersect(const MGCurve& curve) const;

///Default surface-curve intersection function.
MGCSisects intersect(const MGEllipse& el) const;

///Curve and Surface intersection.

/// Surface �� Curve �̌�_�����߂�B
virtual MGCSisects isect(const MGCurve& curve) const;
virtual MGCSisects isect(const MGStraight& sl) const;

///Surface and Surface intersection.

///Surface �� Surface �̌�������߂�B
virtual MGSSisects isect(const MGSurface& srf2)const;
virtual MGSSisects isect(const MGPlane& srf2)const;
MGSSisects isect(const MGFace& f) const;
MGHHisects isect(const MGShell& shl) const;
MGSSisects isect(const MGFSurface& fsurf) const;

MGSSisects isectFS(const MGFace& fsurf) const override{ return isect(fsurf); }
MGSSisects isectFS(const MGSurface& fsurf) const override{ return isect(fsurf); }
MGCSisects isectFS(const MGCurve& curv) const{ return isect(curv); }

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

///Return order of intersection line order of MGLBRep.
///The default is 4.
virtual int isect_order() const=0;

///@cond
///isect_startH compute one intersection line of two surfaces, this and sf2,
/// given starting intersetion point uvuv((u1,v1) of this and (u2,v2) of sf2).
/// isect_startH halts the computation when intersection
/// reached to a boundary of this or sf2, or reached to one of the points
/// in uvuv_list.
///The function's return value is:
/// =0: Intersection was not obtained.
/// !=0: Intersection was obtained as follows:
///    =1: End point is a point on a perimeter of one of the surfaces.
///    =3: End point is one of boundary points in uvuv_list.
///    =4: End point is the starting point.
///    =7: isect_start halted the computation since intersection was lost
///     during the computation.
int isect_startH(
	const MGPosition& uvuv_startIn, ///<Starting point of the intersection line.
	MGPosition_list& uvuv_list,	///<isect_start will halt when ip reached one of 
		///<the point in uvuv_list. isect_start does not change uvuv_list(actually
		///<uvuv_list is const). uvuv's space dimension is at least 4,
		///<and the first 2 is (u,v) of this and the next 2 is (u,v) of sf2. 
	const MGSurface& sf2,	///<2nd surface.
	MGSSisect& ssi,			///<Surface-surface intersection line will be output.
	MGPosition_list::iterator& uvuv_id
		///<When the end point of ip was one of the points of uvuv_list, that is, 
		///<when the function's return value was 3, uvuv_list's iterator
		///<of the point will be returned,
		///<When the end point was not a point of uvuv_list, end() of uvuv_list
		///<will be returned.
)const;

///Compute intersection points of this face's boundary(outer and inners) with
///face2. If intersection points are found and the boundary is a loop,
///the point's edge pointer(of this) will be stored in a member uvuv of uvuvs.
///uvuv[7] is the edge pointer. If the boundary is not a loop(that is, a perimeter of
///Surfaces), uvuv.sdim()==7 and an edge pointer is not returned.
///When uvuv.sdim()==8, the edge pointer of uvuv[7] is accessed through union mgEdgeP.
///uvuvs[i] is i-th intersection points.
virtual int isect_boundary(
	const MGFSurface& face2,///<2nd surface.
	MGPosition_list& uvuvs,
	///<id1 and id2 are the ids of uvuv where this face's and f2's parameters
	///<are to be stored in a member of uvuvs.
	///<This face's (u,v) is stored in uvuv(id1) and (id1+1).
	///<f2's (u,v) is stored in uvuv(id2) and (id2+1).
	///<id2=0 if id1=2, and id2=2 if id1=0.
	int id1=0///< id of above.
)const{return isect_outcurves(face2,uvuvs,id1);};

///Compute intersection points between the boundary of iid-th inner boundary
///of this face and face2 to compute intersections of face with face2.
///Function's return value is the number of ip's obtained before appending
///into uvuv_list, may not be equal to the enlarged size of uvuv_list.
virtual int isect_incurves(
	const MGFSurface& face2,///<2nd surface.
	int iid,	///<Inner loop id of this face(from 0).
	MGPosition_list& uvuv_list,	///<intersection points will be appended,
		///<One member in the list is of sdim 8,
		///<(4,5,6) is the direction vector, and (7) is Edge pointer of the point.
	int id1			///<id of uvuv(a member of uvuv_list),
		///<uvuv(id1) for this face parameter uvuv(id2) for face2 parameter,
		///<id2=0 if id1=2, and id2=2 if id1=0.
)const{return 0;};

///Compute intersection points of outer boundary curves of this face 
///with face2 to compute intersections.
///Function's return value is the number of ip's obtained(appended)
///into uvuv_list, may not be equal to the enlarged size of uvuv_list.
int isect_outcurves(
	const MGFSurface& face2,///<2nd surface.
	MGPosition_list& uvuv_list,	///<intersection points will be appended,
		///<One member in the list is of sdim 7,
		///<and the last three elements are the ip direction vector.
	int id1			///<id of uvuv(a member of uvuv_list),
		///<uvuv(id1) for this face parameter uvuv(id2) for srf or face2 parameter.
		///<id2=0 if id1=2, and id2=2 if id1=0.
)const;
///@endcond

///Access to i-th element of u knot.
virtual double knot_u(int i) const{return 0.0;}

///Access to i-th element of v knot.
virtual double knot_v(int i) const{return 0.0;}

///Returns the u knot vector.
virtual const MGKnotVector& knot_vector_u() const=0;
virtual MGKnotVector& knot_vector_u()=0;

///Returns the v knot vector.
virtual const MGKnotVector& knot_vector_v() const=0;
virtual MGKnotVector& knot_vector_v()=0;


///Compare two parameter values. If uv1 is less than uv2, return true.

///Comparison is done after projected to i-th perimeter of the surface.
virtual bool less_than(
	int i,	///<perimeter number.
	const MGPosition& uv1,///<1st point parameter.
	const MGPosition& uv2///<2nd point parameter.
)const;

///@cond
///This is a newed MGSurface object.
///If this is a MGSurface, construct a newed MGFace using this newed MGSurface,
///and returns the MGFace*.
///****THIS MUST BE A NEWED OBJECT***********
MGFace* make_face();
///@endcond

///Get manifold dimension.
int manifold_dimension() const{return 2;};


///Negate direction of surface.

/// �ʂ̕����𔽓]����B
///Change direction of the surface.
virtual void negate(){ exchange_uv();};
virtual void negateFS(){negate();};

///Negate direction of surface.
virtual void negate(
	int is_u	///< Negate along u-direction if is_u is ture,
				///< else along v-direction.
)=0;

///Transform the coordinates of boundary of this geometry.

///Transformation is done so that
///new coordinate of boundary is the same coordinate as the new one of
///this geometry after negate() of this geometry is done.
///That is, boundary coordinates are of parameter world of this geometry.
void negate_transform(MGGeometry& boundary)const;

///Compute normal vector(not unit) at uv.
virtual MGVector normal(double u, double v) const;

///Compute normal vector(not unit) at uv.
virtual MGVector normal(const MGPosition& uv) const;

///Get the object point of this MGFSurface.
virtual const MGObject* object_pointer()const{return this;};
virtual MGObject* object_pointer(){return this;};


///Surface offset. 

///���I�t�Z�b�g�֐�
///�I�t�Z�b�g�����́A�m�[�}�������𐳂Ƃ���B�ȗ����a���傫���I�t�Z�b�g�͍s��Ȃ��B
///�߂�l�́A�I�t�Z�b�g�Ȗʃ��X�g���ԋp�����B�G���[�̂Ƃ�����0�̋Ȗʃ��X�g���Ԃ�B
///�g�������X��line_zero()���g�p���Ă���B
///Surface offset. positive offset value means offset normal direction.line_zero() is used.
///the radius of curvature is larger than offset value.
virtual std::vector<UniqueSurface> offset(
	double ofs_value,	///<�I�t�Z�b�g��
	int& error			///<�G���[�R�[�h 0:���� -2:�ȗ����a�ȏ�̃I�t�Z�b�g�s�� -3:�ʐ����R���X�g���N�^�G���[
)const;


///Offset.

///distance is plus value if the direction is toward normal vector of the
///FSurface. Minus if against the normal vector.
///�G���[�R�[�h 0:���� -1:�ȗ����a�ȏ�̃I�t�Z�b�g�s�� -3:�ʐ����R���X�g���N�^�G���[
int offset_fs(
	double distance,
	std::vector<UniqueFSurface>& vecOfsFSurface	//Offset MGFSurfaces are appended.
)const;

///C1 continuous Surface offset. 

///C1�A���Ȗʂ̈��I�t�Z�b�g�֐�
///�I�t�Z�b�g�����́A�m�[�}�������𐳂Ƃ���B�ȗ����a���傫���I�t�Z�b�g�͍s��Ȃ��B
///�߂�l�́A�I�t�Z�b�g�����Ȗʂ̃I�[�g�|�C���^���ԋp�����B�G���[�̂Ƃ��k�����Ԃ�B
///�g�������X��line_zero()���g�p���Ă���B
///C1 continuous Surface offset. positive offset value means offset normal direction.
///the radius of curvature is larger than offset value.line_zero() is used.
virtual std::unique_ptr<MGSurface> offset_c1(
	double ofs_value,	///<�I�t�Z�b�g��
	int& error			///<�G���[�R�[�h 0:���� -1:�ʂɂ��ꂪ����
						///< -2:�ȗ����a�ȏ�̃I�t�Z�b�g�s�� -3:�ʐ����R���X�g���N�^�G���[
)const;

///Test if point P is ont the surface or not.

/// �^����ꂽ�덷���œ_���ʏ�ɂ��邩�ǂ����e�X�g����B
///Test if point P is ont the surface or not. Even if P is not on the
///surface, return parameter of the nearest point of the surface.
virtual bool on(
	const MGPosition& P, ///<A point. �w��_                        
	MGPosition&			 ///<Parameter of the surface will be returned.
)const;


///Test if input (u,v) is on the perimeter perim_num.

///If u or v is on a perimeter, true will be returned.
virtual bool on_the_perimeter(
	int perim_num,	///<a perimete number is input.
	double u,	///<Surface parameter u of (u,v)
	double v	///<Surface parameter v of (u,v)
)const;
	
///Test the uvcurve is on a perimeter.

///If on a perimeter, true will be returned.
virtual bool on_perimeter(
	const MGCurve& uvcurve,	///<curve of surface parameter (u,v)
	int& perim_num///<if function returned true, the perimete number will be output.
)const;

///Returns the order of u.
virtual	int order_u() const{return 1;}	

///Returns the order of v.
virtual	int order_v() const{return 1;}	

/// Output virtual function.
virtual std::ostream& toString(std::ostream& ostrm) const;

/// Output virtual function.
std::ostream& outFS(std::ostream& ostrm) const{return toString(ostrm);};

///Obtain outer_boundary curves(world coordinates representation) of the FSurface.

///Let the output of outer_boundary() be wcurves and of outer_boundary_param()
///be pcurves, then wcurves[i] corresponds to pcurves[i] one to one.
///The output curves can be considered as a continuous counter-clockwise ordered
///boundary of the surface.
virtual std::vector<UniqueCurve> outer_boundary()const;

///Obtain boundary curves(parameter space representation) of the FSurface.

///Let the output of boundary() be wcurves and of boundary_parameter()
///be pcurves, then wcurves[i] corresponds to  pcurves[i] one to one.
virtual std::vector<UniqueCurve> outer_boundary_param()const;

///Return surface parameter value of a point on the surface. 

///If input point is not on the surface, return the nearest point on the
///surface.
virtual MGPosition param(
	const MGPosition &	///<Point.
)const;

///Obtain the parameter of pcurve that represent the same point as wcurve.eval(tw).

///Let wcurve be a world curve rep that lies on this surface, and
///pcurve is parameter (u,v) expression of wcurve. That is,
///wcurve==MGSurfCurve pline(*this,pcurve). Then, param_of_pcurve() will obtain
///Let S() is this surface, fp() is pcurve, and fw() is wcurve.
///Then S(fp(tp))=fw(tb).
double param_of_pcurve(
	double tw,			///<point parameter of wcurve to get the pcurve parameter.
	const MGCurve& wcurve,///<world curve that lies on this surface.
	const MGCurve& pcurve,///<This surface's parameter rep of wcurve.
	const double* guess=0///<guess parameter value to compute tp. When guess=null,
						///<param_of_pcurve will define the guess parameter.
)const;


///Compute parameter value of given point.

/// If input point is not on the geometry, return the nearest point on the
/// geometry.
MGPosition parameter(
	const MGPosition& P	///<Point
)const;

/// Compute parameter curve, that is parameter t=const world curve.
///Where t=u or v of the surface parameter.

///Returned is newed area pointer, and must be freed by delete.
virtual MGCurve* parameter_curve(
	int is_u				///<Indicates x is u-value if is_u is true.
	, double x				///<Parameter value.
							///<The value is u or v according to is_u.
)const=0;

///Obtain parameter curves.

///In the case of MGSurface, parameter curve is only one. However, in the case
///of MGFSurface,  number of parameter curves are more than one.
std::vector<UniqueCurve> parameter_curves(
	int is_u,		///<True(!=0) if x is u-value.(i.e. obtain u=const line)
	double x	///<parameter value. u or v-value accordint to is_u.
)const;

///Obtain parameter space error.
virtual double param_error() const;
virtual double param_error_u() const;
virtual double param_error_v() const;

/// Return ending parameter value.
virtual double param_e_u() const=0;
virtual double param_e_v() const=0;

///Return parameter value of the middle point of the surface.

///The middle point is defined as the parameter (u,v) where
///u=(param_s_u()+param_e_u())/2, and v likewise.
MGPosition param_mid()const{return center_param();};

///Return parameter range.
virtual MGBox param_range() const;

///Return parameter range of the geometry(�p�����[�^�͈͂�Ԃ�).
MGBox parameter_range() const;

/// Return starting parameter value.
virtual double param_s_u() const=0;
virtual double param_s_v() const=0;

///Compute square of parameter span length from (u.min, v.min) to (u.max, v.max).
virtual double param_span() const;

///Compute part of the surface limitted by the parameter range bx.

///bx(0) is the parameter (us,vs) and bx(1) is (ue,ve).
///That is u range is from us to ue , and so on.
///Retured is newed object, must be deleted.
virtual MGSurface* part(
	const MGBox& bx,///<Target box.
	int multiple=0	///<Indicates if start and end knot multiplicities
					///<are necessary. =0:unnecessary, !=0:necessary.
)const=0;

///Retrieve perimeter i of this surface.

/// i must be < perimeter_num().
///When perimeter_num()==0, this function is undefined.
///Retured is newed object, must be deleted.
virtual MGCurve* perimeter_curve(int i) const;

///Retrieve (u,v) parameter rep curve of perimeter i.
/// i must be < perimeter_num().
///When perimeter_num()==0, this function is undefined.
///The direction of the output is the same as perimeter_curve().
UniqueCurve perimeterUVCurve(int i) const;

///Return how many perimeters this surface has.
virtual int perimeter_num() const=0;

/// Construct perimeter i's (u,v) parameter position.
virtual MGPosition perimeter_uv(int i,double t) const;

///Compute a perpendicular point from a point P, given guess parameter value uvguess.

///Function's return value is:
///		true if uv is obtained, false if uv is not obtained.
virtual int perp_guess(
const MGPosition& uv0,///<parameter range of this surface, from uv0 to uv1.
const MGPosition& uv1,///<parameter range of this surface,
				///<When uv0(0)>=uv1(0) or uv0(1)>=uv1(1),
				///<no limit for this parameter range.
const MGPosition& P,		///<Point
const MGPosition& uvguess,	///< guess parameter value of surface.
MGPosition& uv				///< Parameter value will be returned.                    
)const;

///Return the foot of the perpendicular straight line from P.

///Computation is done from the guess parameter value.
///Function's return value is whether point is obtained(true) or not(false).
bool perp_guess(
	const MGPosition& P,		///<Point
	const MGPosition& uvguess,	///< guess parameter value of the shell
	MGPosition& uv				///< Parameter value will be returned.
)const;

///Compute perpendicular points of a curve and a surface.

///Given a guess starting paramter value, compute perpendicular points.
///Function's return value is:
///   perp_guess=true if perpendicular points obtained,
///   perp_guess=false if perpendicular points not obtained,
virtual int perp_guess(
const MGPosition& uv0,///<parameter range of this surface, from uv0 to uv1.
const MGPosition& uv1,///<parameter range of this surface.
		///<When uv0(0)>=uv1(0) or uv0(1)>=uv1(1),
		///<no limit for this parameter range.
const MGCurve& curve,///<curve.
double t0,///<parameter range of curve, from t0,
double t1,///<parameter range of curve, to t1.
		///<When t0>=t1, no limit for curve2 parameter range.
const MGPosition& tuvg,///<Guess parameter value of curve and this surface.
MGPosition& tuv	///<perpendicular points' parameter values will be output.
///<tuv(0): curve's parameter, (tuv(1),tuv(2)):this surface's parameter.
)const;

///Compute perpendicular points of a curve and a surface.
virtual int perp_guess(
const MGPosition& uv0,///<parameter range of this surface, from uv0 to uv1.
const MGPosition& uv1,///<When uv0(0)>=uv1(0) or uv0(1)>=uv1(1),
		///<no limit for this parameter range.
const MGCompositeCurve& crv,///<MGCompositeCurve.
double t0,///<parameter range of curve, from t0 to t1.
double t1,///<When t0>=t1, no limit for curve2 parameter range.
const MGPosition& tuvg,	///<Guess parameter value of curve and this surface.
MGPosition& tuv	///<perpendicular points' parameter values will be output,
///<tuv(0): curve's parameter, (tuv(1),tuv(2)):this surface's parameter.
)const;

///Compute perpendicular points of a curve and a surface.

///Given guess a starting paramter value, compute perpendicular points.
///given guess starting paramter values.
///Function's return value is:
///   perp_guess=true if perpendicular points obtained,
///   perp_guess=false if perpendicular points not obtained,
virtual bool perp_guess(
	const MGCurve& curve,	///<curve.
	const MGPosition& uvguess,	///<Guess parameter value of the FSurface.
	double tguess,	///<Guess parameter value of the curve.
	MGPosition& uv,	///<perpendicular point's parameter values of the shell will be output.
	double& t		///<Curve's point's parameter value will be output.		
)const;

///Return the foot of the perpendicular straight line from P.

///Return the foot that is nearest to point uvguess.
/// �^�_�ɂ����Ƃ��߂��A�^�_����ʂɐ����Ȗʏ�̓_�����߂�B
///Function's return value is whether point is obtained(1) or not(0).
virtual int perp_point(
	const MGPosition &P,	///< �w��_                       
	MGPosition& uv,			///<Parameter value of the plane will be output
							///<�p�����[�^�l
	const MGPosition* uvguess=0	///< guess parameter value of the surface
)const;

///Return all foots of perpendicular straight lines from P.
virtual MGPosition_list perps(
	const MGPosition& P				///< Point of a space(�w��_)
)const;

///Compute the parameter value of the closest point from the straight to this object.

///sl is the eye projection line whose direction is from yon to hither, and if
///sl had multiple intersection points, The closest point to the eye will be selected.
MGPosition pick_closest(const MGStraight& sl)const;

///Round the input parameter value uv.

/// ���̓p�����[�^���p�����[�^�͈͂ł܂�߂ĕԋp����B
///Round the input parameter value uv into the parameter range of the surface.
virtual MGPosition range(const MGPosition&) const;

///removal redundant knots within the tolerance line_zero().

///�m�b�g�폜�֐�(B�\���Ȑ��̂�)
///�g�������X��line_zero���g�p����B���̃m�b�g���ׂ������̂قǍ폜���₷��
virtual void remove_knot();

///Rotate this surface around the straight line.

///Rotate the surface around the straight line whose direction is vec and
///that passes through origin.
virtual MGSurface& rotate_self(
	const MGVector& vec,
	double angle,
	const MGPosition& origin=mgORIGIN
);

///Return the surface type.
virtual int sdim() const=0;

///Shrink this surface to the part limitted by the parameter range of uvbx.

///New parameter range uvbx2 is so determined that uvbx2 is the smallest
///box tha includes uvbx, and all of the u or v values of uvbx2 is one of 
///the values of u or v knots of the surface knotvector.
///uvbx(0) is the parameter (us,ue) and uvbx(1) is (vs,ve).
///That is u range is from us to ue , and so on.
virtual void shrink_to_knot(
	const MGBox& uvbx,///<Target parameter box.
	int multiple=0	///<Indicates if start and end knot multiplicities
					///<are necessary. =0:unnecessary, !=0:necessary.
){;};

///split this fsurface at the parameter param.
virtual void split(
	double param,///<parameter value of this fsurface. if is_u is true, param is u-value,
				///<else v-value.
	bool is_u,	///<indicates if param is u or v of the surface parameter (u,v).
	std::vector<UniqueFSurface>& surfaces///<splitted surfaces will be output.
)const;

///Construct a tangent plane LBRep of a surface's perimeter.
///Output MGLBRep�@tp's knot vector is the same as srf's perimeter's 
///if the perimeter is MGLBRep. Let ft(t) is pcrv and g(t) be the output tp ,
///then g(t)=normal at srf(u,v) where (u,v)=f(t).
virtual void TPatPerimeter(
	int perimeterNum,
	MGLBRep& tp///<Obtained tp is output.
)const;

///Return the surface type.
virtual MGSURFACE_TYPE type() const=0;

///Compute unit normal vector at uv.
MGUnit_vector unit_normal(const MGPosition& uv) const;

///Compute unit normal vector at uv.
MGUnit_vector unit_normal(double u,double v) const;

///Get the name of the class.
virtual std::string whoami()const{return "Surface";};

protected:

///Test if the surface is flat or not.
	
///Test is done within the parameter value rectangle of uvbox.
///Function's return value is:
///	true: if the surface is flat
///  false: if the surface is not falt.
///When this is not falt, the direction that indicates which direction the surface
///should be divided will be output.
///***** the flatness is tested only approximately. This is for exclusive use of
///planar().
virtual bool flat(
	const MGBox& uvbox,///<Target parameter range.
	double tol,		///<Tolerance allowed to regart flat
					///<(Allowed distance from a plane).
	int& direction,	///<   1: u-direction is more non flat.
					///<   0: v-direction is more non flat.
	MGPosition& P,	///<Position of the flat plane will be output.
	MGUnit_vector& N///<Normal of the flat plane will be output.
)const;

///Default intersection program of MGSurface.

///It is assumed that both this and srf2 are not a plane.
MGSSisects intersect(const MGSurface& srf2) const;

///Default intersection program of MGSurface with a plane.
virtual MGSSisects intersectPl(const MGPlane& srf2) const;

///Compute intersection points of an inner parameter line of this surface and sf2.

///The intersection point is used to compute surface to surface intersection lines.
///Function's return value is at most one intersection point un uvuv_list.
///One member of uvuv_list is (u1,v1,u2,v2), where (u1,v1) is a parameter of
///this surface and (u2,v2) is a parameter of surf.
MGPosition_list intersectInner(
	const MGSurface& sf2		///<The second surface.
)const;

///@cond
///The following two function will be used in perps or isect
///to decide how many division of the surface along u or v direction
///should be applied before using perp_guess or isect_guess.
virtual int intersect_dnum_u() const=0;
virtual int intersect_dnum_v() const=0;

///isect_area_length() returns initial area length for the intersection
///line.
virtual int isect_area_length() const{return (bdim_u()+bdim_v())*2;};

///isect_direction() is used by isect_startPt() to define which constant
///parameter line should be used to compute intersection, and what
///incremental value be used for the parameter.
///Function's return value is direction to get next intersection(with dt).
///When =1: u=const direction, =0: v=const, =-1: cannot get intersection.
virtual int isect_direction(
	const MGFSurface& sf2,	///<Second surface for the intersection.
	int m1,		///<id of uvuvS that indicates this surface's parameter
		///<position in uvuvS. (uvuvS(m1), uvuvS(m1+1))=(u,v) of this surface.
	MGPosition& uvuvS,///<start parameter (u,v) pair of this surface and sf2.
	double& du,	///<Incremental value of the parameter kind of kdt will be output.
	double& dv, ///<Right dt will be output according to the function's output =0,1.
	double acuRatio=1.	///acuracy ratio.
)const{return MGFSurface::isect_direction(sf2,m1,uvuvS,du,dv,acuRatio);};

///isect_div_id_max is maximum id of array of sect_div defined in
///isect_dt_coef. That is, isect_div_id_max+1 is the length of the array
///sect_div.
int isect_div_id_max()const;

///"isect1_incr_pline" is a dedicated function of isect_start_incr, will get
/// shortest parameter line necessary to compute intersection.
virtual MGCurve* isect_incr_pline(
	const MGPosition& uv,	///<last intersection point.
	int kdt,				///<Input if u=const v-parameter line or not,
							///< true:u=const, false:v=const.
	double du,	///<Incremental parameter length of u.
	double dv,	///<Incremental parameter length of v.
	double& u,				///<next u value will be output.
	double& v,				///<next v value will be output.
	int incr=0			///<Incremental valuse of B-coef's id.
)const=0;

///"isect_inner_dt" is a dedicated function of isect_startPt,
/// comutes adequate incremental parameter value(du,dv) and parameter line kind
///kdt(u=const or v=const).
virtual void isect_inner_dt(
	int n,		///<num of i.p. obtained so far(not include uvnow).
	const MGPosition& uvnow,///<intersection point obtained last(of this).
	double& du,	///<incremental u length from previous to uvnow is input.
				///<New du or dv will be output according to kdt's return value.
	double& dv,	///<incremental v length.
	int& kdt,	///<Parameter kind used so far is input, will be output as:
				///<=1:parameter line kind(u=const), =0: v=const,
				///<=-1:should halt computation since incremental value is zero.
	double acuRatio=1.	///<Accurate ratio.
)const{MGFSurface::isect_inner_dt(n,uvnow,du,dv,kdt,acuRatio);};

///<Compute intersections with MGLBRep lb that does not have C0 continuity in it.
virtual MGCSisects isect_withC1LB(const MGLBRep& lb)const;

///isect with SurfCurve whose m_curve is not a MGTrimmedCurve of MGCompositeCurve.
virtual MGCSisects isect_with_noCompoSC(const MGSurfCurve& scrv)const;
///@endcond

///Intersection of Surface and a straight line.
virtual MGCSisects isectSl(
	const MGStraight& sl,///<Target straight.
	const MGBox& uvbox=mgNULL_BOX ///<indicates if this surface is restrictied to the parameter
					///<range of uvbox. If uvbox.is_null(), no restriction.
)const;

///Obtain 1D surface rep. of this surf.

///This surf1D is used in isect for
///the argument of isect_startPlane, which will use surf1D to compute isect(pl).
///surf1D=0.(intersection with x=0. plane) is the intersection lines.
virtual std::unique_ptr<MGSBRep> surf1D(const MGPlane& pl)const=0;

///�����o�f�[�^��ǂݍ��ފ֐�
virtual void ReadMembers(MGIfstream& buf);

///�����o�f�[�^���������ފ֐�
virtual void WriteMembers(MGOfstream& buf) const;

private:

///Define curve division number when a curve crv be projected onto this MGFSurface.
///The result is used in prj2GetParamRange().
int get_proj_divnum(const MGCurve& crv)const;

///Compute perpendicular points of a curve and a surface,
///given guess starting paramter values.
///Function's return value is:
///   perp_guess=true if perpendicular points obtained,
///   perp_guess=false if perpendicular points not obtained,
int perp_guess_general(
	const MGPosition& uv0,
	const MGPosition& uv1,///<parameter range of this surface.
		///<When uv0(0)>=uv1(0) or uv0(1)>=uv1(1),
		///<no limit for this parameter range.
	const MGCurve& curve,///<curve.
	double t0, double t1,///<parameter range of curve,
		///<When t0>=t1, no limit for curve2 parameter range.
	const MGPosition& tuvg,///<Guess parameter value of curve and this surface.
	MGPosition& tuv	///<perpendicular points' parameter values
		///<will be output, tuv(0): curve's parameter, (tuv(1),tuv(2)):this surface's parameter.
) const;

///C1�A���Ȗʂ̈��I�t�Z�b�g�֐��Ŏg�p����m�b�g�x�N�g���ƃf�[�^�|�C���g�����߂�֐�
void offset_calc_knot_vec(
	MGKnotVector& knot_vec_u,			///<���܂���u�m�b�g�x�N�g��
	MGKnotVector& knot_vec_v,			///<���܂���v�m�b�g�x�N�g��
	MGNDDArray& data_point_u,			///<���܂���u�f�[�^�|�C���g
	MGNDDArray& data_point_v	///<���܂���v�f�[�^�|�C���g
) const;

///�ȗ����a�ƃI�t�Z�b�g�l�̔�����s��
///error�̂Ƃ�false���ԋp�����
int offset_check_curva(
	double ofs_value///<�I�t�Z�b�g��
)const;	

///�ȗ����a�ƃI�t�Z�b�g�l�̔�����s��(u=const�̃p�����[�^�Ȑ����g�p����)
///error�̂Ƃ�false���ԋp�����
int offset_check_curva_one(
	double ofs_value	///<�I�t�Z�b�g��
)const;

///�I�t�Z�b�g����T���v���|�C���g��1�p�b�`���Ƃ̕����������߂�
///�S�Ẵp�b�`���̕������ōő�̒l��Ԃ�
virtual int offset_div_num() const;

///�I�t�Z�b�g����T���v���|�C���g��1�p�b�`���Ƃ̕����������߂�
///������n = sqrt(1 / tol) * sqrt((M1 + 2 *M2 + M3) / 8)�͈ȏ�̎��ŋ��܂�
///Compute number of division for offset in range.
int offset_div_num_one(
	const MGBox& param_range) const;

//Get nearest (u,v) that is not degenerated point.
//That is, nearest point |fu*fv|>MGTolerance::mach_zero().
void nearest_non_degenerated(double& u,double& v) const;

friend class MGStraight;
friend class MGLBRep;
friend class MGSurfCurve;
friend class MGSBRep;
friend class MGPlane;
friend class MGCylinder;
friend class MGFSurface;
friend class MGFace;
friend class MGBSumSurf;

///@cond
///Check if the intersection line lineb's start and end tangent vectors are accurate
///enough. If they do not have enough accuracy, isect_start_tan returns
///which end did not have the accuracy.
/// 1:start, 2:end, 3:start and end. If both ends had enough accuracy, returns 0.
MG_DLL_DECLR friend int isect_start_tan(
	const MGFSurface& sf1,///<Surface 1.
	const MGFSurface& sf2,///<Surface 2.
	const MGLBRep& lineb,///<the intersection line.
	MGVector* Tse[2]	///<If an end had not the accuracy, accurate tangent
						///<will be output. Tse[0]:start, Tse[1]:end.
);

///Update lineb so as to have the tangent tan for start or end according to ngtan.
MG_DLL_DECLR friend void isect_start_adjustSE(
	int ngtan,	///<Return value of isect_start_tan, indicates which end be
				///<adjusted. =1: start tangent, =2:end tangent, =3:both tangent.
	MGNDDArray& tau,	///<data point abcissa.
	MGBPointSeq& point,	///<data point ordinate.
	MGLBRep& lineb,	///<line b-rep obtained so far. tangent adjusted new B-rep
					///<will be output.
	MGVector* tan[2]///<accurate tangent data obtained by isect_start_tan.
);
///@endcond

};

///Update lineb so as to have the tangent tan for start or end according to ngtan.
MG_DLL_DECLR void isect_start_adjustSE(
	int ngtan,	///<Return value of isect_start_tan, indicates which end be
				///<adjusted. =1: start tangent, =2:end tangent, =3:both tangent.
	MGNDDArray& tau,	///<data point abcissa.
	MGBPointSeq& point,	///<data point ordinate.
	MGLBRep& lineb,	///<line b-rep obtained so far. tangent adjusted new B-rep
					///<will be output.
	MGVector* tan[2]///<accurate tangent data obtained by isect_start_tan.
);

namespace MGCL{

///���u�Ȑ��񂩂�ʂ��쐬����

///���u�Ȑ��̑S�Ẵm�b�g�������X�v���C���̎�MGSBRep��MGRSBRep���ԋp�����
///����ȊO�̏ꍇ�́A�Ȑ���LBRep�ōč\�����Ėʂ��쐬����̂�MGSBRep���ԋp�����
///�쐬����ʂ̃m�b�g�x�N�g���̓��u�Ȑ��̌�����u,���u�������v�Ƃ���
///curves[i] must have the same direction.
///Let v0=start parameter value, v1=terminate parameter value along v, then
///v=v0 const parameter line is curves[0], and v=v1 const parameter line is
///curves[n-1], where n=curves.size(). n must be greater or equal to 2.
///When n==2, the surface is a ruled surface(that is, order_u() is 2).
MG_DLL_DECLR std::unique_ptr<MGSurface> createSurfaceFromRibs(
	const std::vector<const MGCurve*>& curves,///<���u�Ȑ���
	bool direction_adjustment=true///<  =true, curves[.] direction are adjusted to line
								///<to the same direction.
);

/// Creates a ruled surface.
MG_DLL_DECLR std::unique_ptr<MGSurface> create_ruled_surface(
	const MGCurve& cross1,    ///< a curve as Edge No.1.
	const MGCurve& cross2,     ///< another curve as Edge No.2.
	bool direction_adjustment=true///<  =true, curves[.] direction are adjusted to line
								///<to the same direction.
);


///Creates a surface of revolution.

///Parameterization of the surface is:
///	u=const parameter line generates given curve(when u=0.).
///  v=const parameter line generates a circle whose center is axis.
MG_DLL_DECLR std::unique_ptr<MGRSBRep> create_revolved_surface(
	const MGCurve& curve,     ///< generatrix curve
	const MGStraight& axis,   ///< revolution axis
	double angle = mgDBLPAI   ///< revolution angle
);


///Test if pline has the same direction to world_curve.

///Assuming that pline=MGSurfCurve(MGSurface srf, parameter_curve of the srf) and
///world_curve are the same curve.
///Function's return value is:
///1: same direction, -1:oppositie direction.
int SurfCurve_equal_direction(
	const MGCurve& pline,	///<MGSurfCurve(MGSurface srf, parameter_curve of the srf)
	const MGCurve& world_curve	///<world representation curve.
);

///Evaluate which direction is longer, u or v, from the 9 sample points.

///Points are at the parameter range square of the surface.
///Function's return value is square of the length of the longer direction.
double get_length(
	MGPosition Pn[9],	//9 sample points.
					//Pn[.]={(umin, vmin), (umid, vmin), (umax,vmin),
					//       (umin, vmid), (umid, vmid), (umax,vmid),
					//       (umin, vmax), (umid, vmax), (umax,vmax),
	bool& direction	//which direction was longer is returned.
					//True: when u direction is longer than v direction.
);

}

/** @} */ // end of GEO group
#endif
