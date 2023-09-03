/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#ifndef _MGFace_HH_
#define _MGFace_HH_

#include <vector>
#include "mg/Default.h"
#include "mg/drawParam.h"
#include "mg/BPointSeq.h"
#include "mg/Unit_vector.h"
#include "mg/Position.h"
#include "mg/Curve.h"
#include "mg/FSurface.h"
#include "mg/Surface.h"

#include "topo/Cell.h"
#include "topo/Loop.h"
#include "topo/LEPoint.h"
#include "topo/FOuterCurve.h"

class MGSSisects;
class MGCSisects;
class MGPosition_list;
class MGLBRep;
class MGCell;
class MGLoop;
class MGLEPoint;
class MGShell;
class MGFace;
class MGHHisects;

//
//Define MGFace Class.

/** @addtogroup TOPO
 *  @{
 */


//Used for tessellatiion(shade).
//Represents  n polylines of n edges of a loop.
//SHLL_COM_EDGES[j] is j-th edge's polyline of a loop.
//Generally a face has multiple SHLL_COM_EDGES.
typedef std::vector<const MGLBRep*> SHLL_COM_EDGES;

///MGFace is a trimmed surface.

///MGFace is an instance of MGCell, can be a constituent of MGShell.
///Many useful functions are provided in MGFSurface. See MGFSurface.
class MG_DLL_DECLR MGFace: public MGCell, public MGFSurface{
private:
	///vector of boundaries who bound this cell.
	std::vector<UniqueLoop> m_boundaries;

	mutable MGBox m_box_param;///<Box of parameter space of the face.

public:

///Iterator definition.
typedef std::vector<UniqueLoop>::iterator iterator;
typedef std::vector<UniqueLoop>::const_iterator const_iterator;
typedef std::vector<UniqueLoop>::reverse_iterator reverse_iterator;
typedef std::vector<UniqueLoop>::const_reverse_iterator const_reverse_iterator;

///Faceのスケーリングを行い，Faceを作成する。
///Scaling of the Face by a double.
MG_DLL_DECLR friend MGFace operator* (double s, const MGFace& face);

///////// Constructor /////////

///Null face.
MGFace() = default;
~MGFace();

///Copy constructor.
MGFace(const MGFace& face);
MGFace(MGFace&& face);

///Constructor.
MGFace(const MGFace& face, bool boundaryIsNecessary);

///Fundamental constructor.
///Construct a face from geometry of manifold dimension 2
///and the boundaries.
///The constructor takes the ownership of geo and MGLoop in boundaries.
///boundaries must be loops.
MGFace(
	MGSurface* geo,
	std::vector<UniqueLoop>&& boundaries
);
MGFace(
	MGSurface* geo,
	const std::vector<UniqueLoop>& boundaries
);

///Face of whole surface of no boundary.
MGFace(const MGSurface& surf)
:MGCell(surf),m_box_param(surf.param_range())
{;}

///Conversion constructor from MGFSurface to MGFace.
MGFace(const MGFSurface& surf);

///Face of whole surface of no boundary.
///Ownership of surf is transfered to the face.
///(that is surf must be a newed object.)
explicit MGFace(MGSurface* surf)
:MGCell(surf),m_box_param(surf->param_range())
{;}

///Construct a face by copying boundaries(only parameter rep of the boundary)
///from argument boundaries.
///Second form is to input a newed surface. The constructor takes the ownership
///of the surf and boundaries.
MGFace(const MGSurface& surf,
	const std::vector<MGLoop*>& boundaries);
MGFace(UniqueSurface&& surf,
	std::vector<UniqueLoop>&& boundaries);

///////// operator overload/////////

///Assignment.
///When the leaf object of this and cell2 are not equal, this assignment
///does nothing.
MGFace& operator=(const MGGel& gel2);
MGFace& operator=(const MGFace& gel2);
MGFace& operator=(MGFace&& f2);

/// Faceに平行移動を行ないオブジェクトを生成する。
///Translation of the Face
MGFace operator+ (const MGVector& v) const;

/// Faceに逆方向の平行移動を行ないオブジェクトを生成する。
///Translation of the Face
MGFace operator- (const MGVector& v) const;

///Faceのスケーリングを行い，Faceを作成する。
///Scaling of the Face by a double.
MGFace operator* (double s) const;

/// 与えられた変換でFaceの変換を行い，Faceを作成する。
///Transformation of the Face by a matrix.
MGFace operator* (const MGMatrix& mat) const;

/// 与えられた変換によってトランスフォームをおこないFaceを生成する。
///Transformation of the Face by a MGTransf.
MGFace operator* (const MGTransf& tr) const;

///Complexのスケーリングを行い，Complexを作成する。
///Scaling of the Complex by a double.
MGFace operator/ (double s) const{return (*this)*(1./s);};

///Comparison of two faces.
bool operator==(const MGFace& gel2)const;
std::partial_ordering operator<=>(const MGFace& gel2)const;

//gel2 must be the same class as this.
bool equal_test(const MGGel& gel2)const override;

//gel2 must be the same class as this.
std::partial_ordering ordering_test(const MGGel& gel2)const override;

///Return space dimension
int sdimFS() const override{	return MGCell::sdim();};

///PD144=MGFace. Output to PD144(Trimmed surface).
int out_to_IGES(
	MGIgesOfstream& igesfile,
	int SubordinateEntitySwitch=0
)const;

///Debug Function
std::ostream& toString(std::ostream& ostrm) const;

/// Output virtual function.
std::ostream& outFS(std::ostream& ostrm) const{return toString(ostrm);};

/////////Member function/////////

///Add a new loop to this face as aboundary.
///When the old loops that are outside the nloop will be removed from this.
///nloop can be inner or outer.
///nloop must be a newed MGLoop, and the ownership is transfered to this.
void add_boundary(MGLoop* nloop);

///Generate arrow data of the tangent along u and v and the normal
///at the parameter value (u,v) of the surface.
///data[0] is the origin of the u-tangent arrow, data[1] is the top of the u-tangent arrow,
///data[2], [3] are two bottoms of u-tangent arrowhead.
///data[0], [4], [5], [6] are the points of v-tangent arrow.
///data[0], [7], [8], [9] are the points of v-tangent arrow.
void arrow(double u,double v, MGPosition data[10])const;
void arrow(const MGPosition& uv, MGPosition data[10])const{
	arrow(uv[0],uv[1],data);
}

///Append new one boundary to boundary vectors.
///Returned is the number of boudaries after appending.
///bound must be a newed MGLoop, and the ownership is transfered to this.
///*** append_boundary does not check validity with other loops
///(e.x. already existed loops will be outside the new boudanry bound).
///If the validity check is necessary, use add_boudanry().
int append_boundary(MGLoop* bound);

///Prepend new one boundary to boundary vectors.
///Returned is the number of boudaries after prepending.
virtual int prepend_boundary(MGLoop* bound);

///Obtain binder of this if this is a binder(null when not).
const MGBCell* binderCell()const override{return nullptr;};
MGBCell* binderCell()override{return nullptr;};

///Obtain i-th boundary pointer.
UniqueLoop& boundary(int i) { return m_boundaries[i]; };
const UniqueLoop& boundary(int i) const { return m_boundaries[i]; };

///Obtain iterator of m_boundaries.
const_iterator boundaryIterator(const MGLoop* bnd) const;
iterator boundaryIterator(MGLoop* bnd);

///Obtain boundaries of this cell.
const std::vector<UniqueLoop>& boundaries() const{ return m_boundaries; };

///Obtain all the boundary curves(world coordinates representation)
///of the face.
///That is, all of the outer boundaries and all of the inner boundaries.
std::vector<UniqueCurve> face_boundaries()const;

///If this had boundary binders, free them. As the result this
///will have no neighbours.
void free_neighbours() override;

///Return box of the parameter space of the face.
///After trimmed one.
const MGBox& box_param() const;

///Return box of the parameter space of the FSurface.
///After trimmed one.
const MGBox box_param2() const{return box_param();};

///Build a loop of this face, given a closed curve crv on this face. Although crv
///is generally a MGCompositeCurve, this may be not the case. Returned MGLoop is not
///added into this face as a boundary. User must add it after the direction is adjusted.
///That is, the output loop can be an outer or inner loop.
std::unique_ptr<MGLoop> build_loop(
	const MGCurve& crv	///<curve of world coordinates.
		///<Generally this is not on face and always is projectd onto the face.
)const;

///Obtain the center parameter value of this cell.
MGPosition center_param() const;

///Make a clone of the cell.
///clone() does not copy the binder cell of this.
MGFace* clone() const override;

///Make a clone of the cell that does not have boundaries.
///Does not copy the binder cell of this.
MGFace* cloneWithoutBoundary() const override;

///Get the clone of this MGFSurface.
MGFSurface* clone_fsurface()const override{return clone();};

///Get the clone of this as a MGFace.
///If this is MGSurface, it is converted to MGFace.
MGFace* clone_as_face()const{return clone();};

///Compute closest point from a point.
///Returned is the parameter value of the face that is closest to point.
MGPosition closest(const MGPosition& point) const;
	
///Compute closest point from a line to the boundary of the MGFSurface.
///Returned is the parameter value of the FSurface that is closest to point.
MGPosition closest_on_boundary(const MGStraight& sl) const;

///compute box of the cell in bx.
///Currently this does not compute correct box, compute m_extent box.
void compute_box(MGBox& bx) const;

///////display member function.

///Display direction arrows on the surface.
void display_arrows(mgSysGL& sgl)const;

///Display control polygon of the surface.
void display_control_polygon(mgSysGL& sgl)const;

///Draw 3D curve in world coordinates.
///The object is converted to curve(s) and is drawn.
void drawWire(
	mgVBO& vbo,///<Target graphic object.
	int line_density=1	///<line density to draw a surface in wire mode.
)const{drawWireFS(vbo,line_density);};

///Draw 3D point(vertex) in world coordinates.
///The object is converted to point(s) and is drawn.
///This is valid only for topology objects or MGPoint.
void drawVertex(mgVBO& vbo)const;

///Shade the object in world coordinates.
void shade(
	mgVBO& vbo,///<Target graphic object.
	const MGDrawParam& para,///<Parameter to draw.
	MGCL::DRAW_TARGET target=MGCL::SHADING///<Shading target.
)const;

///Triangulate this object(MGShell, MGFace, or MGSurface is the target).
virtual void triangulate(
	const MGDrawParam& para,
	MGCL::TL_DATA_KIND dkind,
	std::vector<mgTL2Triangles>& trisVec
)const;

///Test if directions of parameter curve and world curve of the face boundary
///is equal or not. This function can be used to test the pair of
///the output of outer_boundary() and outer_boundary_param(), or the pair of
///inner_boundary() and inner_boundary_param().
///Return is:
///true if equal direction, false if opposite direction.
bool equal_direction(
	const std::vector<UniqueCurve>& wcurves,
		///<output of outer_boundary() or inner_boundary().
	const std::vector<UniqueCurve>& pcurves,
		///<output of outer_boundary_param() or inner_boundary_param().
	int i
		///<id of the curve in wcurves and pcurves to test the direction.
)const;

///Erase i-th boundary.
///erase_boundary remove from this cell's bounary and destruct the boundary.
void erase_boundary(iterator i);
void erase_boundary(int i);

///erase_boundary removes from this cell's bounary and destruct the boundary.
void erase_boundary(MGLoop* bnd);

///Evaluate.
///Input parameter value is not checked if it is in_range() or not.
///Even if it is not in_range(), surface evaluation will be executed.
MGVector eval(
	double u,	///<Face parameter value u of (u,v).
	double v,	///< v of (u.v).
	 int ndu=0///<Order of derivative along u.
	 , int ndv=0///< along v.
) const;

///Evaluate at uv.
MGVector eval(
	const MGPosition& uv,	///<Face parameter value(u,v)
	int ndu=0,///<Order of derivative along u.
	int ndv=0 ///<Order of derivative along v.
) const;

///Extract all the loops of this face.
void extract_loops(std::vector<const MGLoop*>& loops)const;

///Extract sub  face that is bounded by networks loops.
///Extracted sub face is the smallest closed part of this face bounded by
///the networks that includes the parameter position uv(u,v).
void extract_sub_face(
	const std::vector<UniqueLoop>& networks,///<(u,v) representation networks.
	const MGPosition& uv,///<Parameter value to indicate which part of the face to extract.
	std::unique_ptr<MGFace>& face///<Result extracted face will be output.
)const;

///Get the MGFSurface pointer if this is MGSurface or MGFace.
const MGFSurface* fsurface()const{return this;};
MGFSurface* fsurface(){return this;};

///Get inner_aboundary loops included in the input box.
std::vector<const MGLoop*> get_inner_boundary_loops(const MGBox& uvbox) const;

///Judge if the display list for vmode is made or not.
//bool displayList_is_made(MGCL::VIEWMODE vmode)const;

///get face pointer if this is MGFace, else null will be returned.
MGFace* get_face_pointer(){return this;};
const MGFace* get_face_pointer()const{return this;};

///get surface pointer. Null will never be returned if this is valid MGFSurface.
///That is, if this is MGFace, base surface will be returned.
MGSurface* get_surface_pointer(){return surface();};
const MGSurface* get_surface_pointer()const{return surface();};

///Get number of inner boundaries as the output of the function.
int get_number_of_boundaries()const{return number_of_boundaries();};

///Test if this and 2nd object has common area about their box(),
///taking error into account.
bool has_commonFS(const MGObject& obj2)const{return has_common(obj2);};

///Test if this face has boundary loops or not in the specified box.
///If this has one, return true.
bool hasLoop(const MGBox& uvbox)const;

///Test if this face has an inactive loop.
///If this has one, return true.
bool hasInactiveLoop()const;

///Test if this face has the outer boundary loop instead of perimeter boundary
///loops. If this has the outer boundary loop and has not perimeter boundary loops,
///return true.
bool hasOuterBoundaryLoop()const;

///Test if this face has perimeter boundary loops or not.
///If this has one, return true.
bool hasPerimeterBoundaryLoop() const;

///Obtain i-th inner_boundary curves(world coordinates representation)
///of the face. Let the output of inner_boundary(i) be wcurves and
///of inner_boundary_param(i) be pcurves, then wcurves[j] corresponds
///to pcurves[j] one to one. Number of inner_boundary can be obtained
///by the function number_of_inner_boundary().
std::vector<UniqueCurve> inner_boundary(int i)const override;

///Obtain i-th inner_boundary curves(world coordinates representation)
///of the face. Let the output of inner_boundary(i) be wcurves and
///of inner_boundary_param(i) be pcurves, then wcurves[j] corresponds
///to pcurves[j] one to one. Number of inner_boundary can be obtained
///by the function number_of_inner_boundary().
std::vector<UniqueCurve> inner_boundary_param(int i)const override;

///Return Object's type ID (TID)
long identify_type()const;

///Test if parameter value (u,v) is in the range of the face parameter.
bool in_range(double u, double v)const;
bool in_range(const MGPosition& uv)const;

///Test if (u,v) is inside the face.

///Function's return value is:
///  0:outside the face.
///  1:unknown.
///  2:inside the face, not on a boundary.
///  <0:(u,v) is on an inner boundary, and abs(return code) is the loop id.
///  4:(u,v) is on the outer boundary.
///  >=10: (u,v) is on a perimeter, (10+perimeter number) will be returned.
int in_range_with_on(const MGPosition& uv)const;

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

///Intersection.
MGCSisects isect(const MGCurve& curv) const;
MGSSisects isect(const MGFSurface& fsurf) const;
MGSSisects isect(const MGFace& fsurf) const;
MGSSisects isect(const MGSurface& fsurf) const;
MGHHisects isect(const MGShell& shell2) const;

MGSSisects isectFS(const MGFace& fsurf) const override{ return isect(fsurf); }
MGSSisects isectFS(const MGSurface& fsurf) const override{ return isect(fsurf); }
MGCSisects isectFS(const MGCurve& curv) const{ return isect(curv); }

///Access to i-th element of u knot
double knot_u(int i) const{return surface()->knot_u(i);};

///Access to i-th element of v knot
double knot_v(int i) const{return surface()->knot_v(i);};

///Returns the u knot vector.
const MGKnotVector& knot_vector_u()const{return surface()->knot_vector_u();};
MGKnotVector& knot_vector_u(){return surface()->knot_vector_u();};

///Returns the v knot vector.
const MGKnotVector& knot_vector_v()const{return surface()->knot_vector_v();};
MGKnotVector& knot_vector_v(){return surface()->knot_vector_v();};

///Obtain i-th boundary loop of the face.
UniqueLoop& loop(int i){ return m_boundaries[i]; };
const UniqueLoop& loop(int i)const{ return m_boundaries[i]; };

///This is a newed MGFace or MGSurface object.
///If this is a MGFace, returns this pointer.
///If this is a MGSurface, construct a newed MGFace using this newed MGSurface,
///and returns the MGFace*.
MGFace* make_face(){return this;};

///Make outer boundary if not existed.
void make_outer_boundary();

///Get manifold dimension.
int manifold_dimension() const{return 2;};

///Negate the face.
void negate();
void negateFS(){negate();};

///Compute normal vector(not unit) at uv.
MGVector normal(const MGPosition& uv) const
{	return surface()->normal(uv.ref(0), uv.ref(1));};

///Compute normal vector(not unit) at (u,v).
MGVector normal(double u,double v) const
{	return surface()->normal(u, v);};

///Test if no outer boundary except the surface perimeters.
///That is, test if the following two conditions are satisfied:
///         1. no perimeter boundaries.
///         2. no outer boundary.
bool no_outer_boundaries()const;

///Return number of boundaries.
int number_of_boundaries() const{ return int(m_boundaries.size()); };

///Get number of inner boundaries as the output of the function.
int number_of_inner_boundaries()const
{	int dummy; return number_of_inner_boundaries(dummy);}

///Get number of inner boundary loops.
///Returned i is the id of the first inner boundary loop if inner boundaries
///exist.
int number_of_inner_boundaries(int& i)const;

///Compute number of active loops.
int number_of_loops()const;

///Get number of perimeter boundary loop.
int number_of_perimeter_boundaries()const;

///Return MGObject pointer if this MGGel is an MGObject, else return null.
MGObject* object_pointer(){ return this; };
const MGObject* object_pointer()const{return this;};

///Offset.
///distance is plus value if the direction is toward normal vector of the
///face. Minus if against the normal vector.
///エラーコード 0:成功 -1:曲率半径以上のオフセット不可 -3:面生成コンストラクタエラー
MGFace offset(double distance, int& error)const;

///Offset.
///distance is plus value if the direction is toward normal vector of the
///face. Minus if against the normal vector.
///エラーコード 0:成功 -1:曲率半径以上のオフセット不可 -3:面生成コンストラクタエラー
int offset(
	double distance,
	std::vector<UniqueFace>& vecOfsFace	//Offset faces are returned(not appended).
)const;

///Offset.
///distance is plus value if the direction is toward normal vector of the
///FSurface. Minus if against the normal vector.
///エラーコード 0:成功 -1:曲率半径以上のオフセット不可 -3:面生成コンストラクタエラー
int offset_fs(double distance, std::vector<UniqueFSurface>& vecOfsFSurface)const;

///Test if a point P is on the face.
///Returned is true if the point P is on the face.
///false(0) if P was not on the face.
bool on(
	const MGPosition& P,///<Target postion data.
	MGPosition& uv	///<Parameter value of the face is returrned.
					///<Even if P is not on the face, nearest point
					///<parameter value will be returned.
)const;

///Test if input (u,v) is parameter value on a perimeter of the base surface.
///If u or v is on a perimeter, they will be updated to the perimeter value.
bool on_a_perimeter(
	double& u,	///<Surface parameter u of (u,v).
	double& v,	///<Surface parameter v of (u,v).
	int& perim_num	///<if function returns true,
					///<the perimete number is output.
)const;

///Obtain outer_boundary curves(world coordinates representation) of the FSurface.
///Let the output of outer_boundary() be wcurves and of outer_boundary_param()
///be pcurves, then wcurves[i] corresponds to pcurves[i] one to one.
///The output curves can be considered as a continuous counter-clockwise ordered
///boundary of the surface.
std::vector<UniqueCurve> outer_boundary()const;

///Obtain boundary curves(parameter space representation) of the face.
///Let the output of boundary() be wcurves and of boundary_parameter()
///be pcurves, then wcurves[i] corresponds to  pcurves[i] one to one.
std::vector<UniqueCurve> outer_boundary_param()const;

///Obtain parameter value of the face whose world coordinates are P.
MGPosition param(const MGPosition& P)const;

///Obtain parameter curves.
///In the case of surface, parameter curve is only one. However, in the case
///of face,  number of parameter curves are more than one.
std::vector<UniqueCurve> parameter_curves(
	int is_u,		///<True(!=0) if x is u-value.(i.e. obtain u=const line)
	double x	///<parameter value. u or v-value accordint to is_u.
)const;

/// パラメータ範囲を返す。
///Return parameter range.
MGBox param_range() const{return box_param();};

/// Return ending parameter value.
double param_e_u()const;
double param_e_v()const;

/// Return starting parameter value of the base surface.
double param_s_u()const;
double param_s_v()const;

///Obtain parent shell that this face belongs to.
MGShell* parent_shell();
const MGShell* parent_shell()const;

///Obtain perimeter boundadary loop's curve representation.
///Returned are curves of perimeter boundaries, do not contain perimeter
///of the surface.
std::vector<UniqueCurve> PBloop_curves() const;

///Return the foot of the perpendicular straight line from P.
///Computation is done from the guess parameter value.
///Function's return value is whether point is obtained(true) or not(false).
bool perp_guess(
	const MGPosition& P,		///<Point
	const MGPosition& uvguess,	///< guess parameter value of the shell
	MGPosition& uv				///< Parameter value will be returned.
) const;

///Compute perpendicular points of a curve and the face,
///given guess starting paramter values.
///Function's return value is:
///   perp_guess=true if perpendicular points obtained,
///   perp_guess=false if perpendicular points not obtained,
bool perp_guess(
	const MGCurve& curve,	///<curve.
	const MGPosition& uvguess,	///<Guess parameter value of the face.
	double tguess,			///<Guess parameter value of the curve.
	MGPosition& uv,			///<perpendicular point's parameter values of the shell
	double& t				///<will be output.
) const;

///指定点から最も近い、垂線の足とパラメータ値を返す。
///Return the foot of the perpendicular straight line from p that is 
///nearest to point p.
/// Function's return value is whether point is obtained(1) or not(0)
int perp_point (
	const MGPosition& p,		///< 指定点(point)
	MGPosition& uv,		///<Parameter value of the surface will be returned.
	const MGPosition* uvguess=0	///< guess parameter value of surface
)const;

///Compute perpendicular points on the face from a point P((x,y,z)).
///MGPosition uv in the MGPosition_list is:
///uv(0): u parameter, and uv(1): v parameter of the face.
///Generally number of uv are more than one.
MGPosition_list perps(const MGPosition& P) const;

///Compute the parameter value of the closest point from the straight to
///this object.
///sl is the eye projection line whose direction is from yon to hither, and if
///sl had multiple intersection points, The closest point to the eye will be selected.
MGPosition pick_closest(const MGStraight& sl)const;

///Round the input parameter (u,v) of the face to the nearest point of
///the face parameter range.
MGPosition range(const MGPosition& uv) const;

///Rebuild this face. Rebuild means:
/// 1) Rebuild the base surface.
/// 2) Rebuild the parameter edges.
/// 3) Rebuild the binder edges from the rebuilt parameter edges.
/// This rebuild does not change the number of edges. All of the old edges
/// does exist after the rebuild.
///Base surface rebuild is done only when the base surface is MGSBRep, or MGRSBRep.
std::unique_ptr<MGFace> rebuild(
	int how_rebuild=1,
		///<intdicates how rebuild be done.
		///< =0: no base surface approximation(only parameter change)
		///< =1: If base surface is MGRSBRep, reconstruct it with new knot configuration
		///<     again as rational spline(MGRSBRep). If MGSBRep, reconstruct with new knot
		///<     configuration as MGSBRep.
		///< =2: approximated by non-rational spline(MGSBRep) with new knot configuration,
		///<     even base surface is MGRSBRep. If the base surface is MGSBRep, same as =1.
	int parameter_normalization=2,
		///<Indicates how the parameter normalization be done:
		///< =0: no surface parameter normalization.
		///< =1: normalize to u_range=(0., 1.), and v_range=(0.,1.);
		///< =2: normalize to make the average length of the 1st deriv along u and v 
		///<     of the base surface is as equal to 1. as possible.
		///< =3: Specify parameter range in range[4].
	double tol=-1.,///<tolerance allowed for the approximation.
		///<When tol<=0., MGTolerance::line_zero() will be employed.
	int* order=0, ///<order of the new MGSBRep, >=4 is recomended.
		///<order[0]:u-order, [1]:v-order.
		///<When how_rebuild!=2, order is not used.
		///<When order=0 is input, order[0]=order[1]=4 are assumed.
	double* range=0///<valid only when parameter_normalization=3,
		///< and range[]={umin, umax, vmin, vmax}.
		///<When parameter_normalization=3 and range=0, parameter_normalization=2
		///<is assumed.
)const;

///Remove inactive loops from this face.
void remove_inactive_loops();

///Rotate only the boundary of this face, but do not rotate the base surface.
///This is designed for the face of scalePolar().
void rotateBoundary(const MGMatrix& mat);

///Execute polar-scaling to this MGFace.
///Parameter curve's (u,v) of the boundary loops are updated, and
///the binder edges extents are removed if exist.
///The base surface is unchanged.
///The updated result parameter curve is always MGLBRep.
///
///Rotation is performed from the angle range (angleBase,angle1) to
///(angleBase,angle2).
///That is, when angle1=angle2, no change is done.
///When angle2 is angleBase, all the data will lie on the straight of from origin to
///(cos(angleBase), sin(angleBase)).
///angle1-angleBase must be >MGTolerance::angle_zero().
void scalePolar(
	double angleBase,	///<base angle.
	double angle1,		
	double angle2
);

///Shrink the base surface of this face to the part limitted by the parameter range of uvbx.
///New parameter range uvbx2 is so determined that uvbx2 is the smallest
///box tha includes uvbx, and all of the u or v values of uvbx2 is one of 
///the values of u or v knots of the surface knotvector.
///uvbx(0) is the parameter (us,ue) and uvbx(1) is (vs,ve).
///That is u range is from us to ue , and so on.
void shrink_base_surface_to_knot(
	int multiple=0	///<Indicates if start and end knot multiplicities
					///<are necessary. =0:unnecessary, !=0:necessary.
);

///Sort boundary occurreces in m_boundaries.
///Sorting is done according to operator< of MGLoop.
///parameter space box will be set.
void sort_boundaries();

///split this fsurface at the parameter param.
void split(
	double param,///<parameter value of this fsurface. if is_u is true, param is u-value,
				///<else v-value.
	bool is_u,	///<indicates if param is u or v of the surface parameter (u,v).
	std::vector<UniqueFSurface>& surfaces///<splitted surfaces will be output.
)const override;

///Split the face giving networks loops. Splitting is done by finding the smallest
///closed areas out of networks.
void split(
	const std::vector<UniqueLoop>& networks,///<Network to split.
	std::vector<UniqueFace>& faces///<Result trimmed face(s) will be appended.
)const;

///Get base surface pointer.
MGSurface* surface();
const MGSurface* surface() const;

///Trim the face by the projection of a curve along a vector direction.

///If mgNULL_VEC is specified as direction, surface normal projection
///will be employed.
///crv has a direction. That is, when this face is divided into two faces,
///left part face of the crv is the face selected.
///
///When the projected curve on the face is connected to already existed
///boundary curve, no new boundary is generated, and inserted(connected)
///to the old boundary. However new projection is floating, that is, not
///connected to any old boundaries, or this boundary is the first boundary,
///new boundary is generated.
///Function's return value is error code:
///0= normal return
///(not error, this includes the case of inactive loop generation)
///2= tried to generate outer boundary loop inside perimeter boudary.
///3= tried to generate inner boundary loop that incudes active loop inside.
///4= tried to generate perimeter boudary loop that inactivates perimeter
///   boundary loops existed.
///5= tried to generate a loop outside the face.
int trim_projection(
	const MGCurve& crv,		///<curve of world coordinates.
		///<Generally this is not on face and always is projectd to the face.
	const MGVector& direction=mgNULL_VEC///<Projection directin vector.
);

///Trim given parameter curve of the face.

///crv has a direction. That is, when this face is divided into two faces,
///left part face of the crv is the face selected.
///
///When the pcrv is connected to already existed boundary curve,
///no new boundary is generated, and inserted(connected) to the old boundary.
///However, new projection is floating, that is, not
///connected to any old boundaries, or this boundary is the first boundary,
///new boundary is generated.
///Function's return value is error code:
///0= normal return
///(not error, this includes the case of inactive loop generation).
///1= input pcrv includes a part that is outside surface parameter range.
///2= tried to generate outer boundary loop inside perimeter boudary.
///3= tried to generate inner boundary loop that incudes active loop inside.
///4= tried to generate perimeter boudary loop that inactivates perimeter
///   boundary loops existed.
///5= tried to generate a loop outside the face.
int trim(
	const MGCurve& pcrv	///<parameter(u,v) space curve of the face.
);

///Trim given a loop new_loop that does not have the parent face.

///new_loop must be parrameter representaion of this face and
///must not have intersections with the loops of this face except the end points
///of new_loop.
///
///Function's return value is error code:
///0= normal return
///(not error, this includes the case of inactive loop generation).
///2= tried to generate outer boundary loop inside perimeter boudary.
///3= tried to generate inner boundary loop that incudes active loop inside.
///4= tried to generate perimeter boudary loop that inactivates perimeter
///   boundary loops existed.
///5= tried to generate a loop outside the face.
int trim(MGLoop&& new_loop_in);

///Trim the face giving networks loops. Trimming is done by removing the smallest
///closed area out of networks that includes the parameter position uv(u,v).
void trim(
	const std::vector<UniqueLoop>& networks,///<network to trim the face.
	const MGPosition& uv,///< position parameter data to indicate which part of the face to trim.
	std::vector<UniqueFace>& faces///<Result trimmed face(s) will be appended.
)const;

///Compute unit normal vector at uv.
MGUnit_vector unit_normal(const MGPosition& uv) const
{	return surface()->unit_normal(uv.ref(0), uv.ref(1));};

///Compute unit normal vector at (u,v).
MGUnit_vector unit_normal(double u,double v) const
{	return surface()->unit_normal(u, v);};

///Get the name of the class.
std::string whoami()const{return "Face";};

protected:

///Write Object's Member Data
void WriteMembers(MGOfstream& buf) const;

///Read Object's member data.
void ReadMembers(MGIfstream& buf);

///Free specified boundary(bound) from a member of parent cell's boundaries.
///Return MGComplex* bound if freed normally.
///If bound was not a member of the boundaries, return 0.
///Only free, does not destruct the boundary.
virtual MGComplex* free_boundary(const MGComplex* bound);

///Get all the MGPCell* of the all the boundaries of this.
std::vector<const MGPCell*> getBoundaryPcells()const override;

private:

///Copy all boundaries of cellin into this(but does not copy own binder cell relation),
///and the copied boundarys' MGPCell* association  of the original(key) and the new one into cmap.
//1st(key) is original MGPCell*(MGEdge) of boundaries(m_boundaries(MGLoop's)).
//2nd is copyied new.
void copy_all_boundaries(
	const MGCell& cellin,
	std::map<const MGPCell*, MGPCell*>* cmap = 0
) override;

///Test in_range if this is a face, if not, do nothing.
///This is to accelerate the test of in_range in isect_guess().
///See isect_guess(MGCurve).
bool in_range_face(const MGPosition& uv)const{return in_range(uv);};

///Obtain coefficient's space dimension.
///This function is used in isect_start etc.
int coef_sdim() const{return surface()->coef_sdim();};

///Compute parameter range box.
void compute_box_param() const;

///Negate the boundary.
void negate_boundary()override;

///Define curve division number when a curve crv be projected onto this MGFSurface.
///The result is used in prj2GetParamRange().
int get_proj_divnum(const MGCurve& crv)const;

///set box as null(to set the box as initial)
void invalidateBox()const override;

///Test if (u,v) is inside inner boundary.

///inside means not on the
///boundary and not included inside the face.
///If true is returned, the id of m_boundaies is returned.
///Function's return value is:
///  0:outside of all the inner loops(not on the loop)
///  1:unknown
///  2:inside an inner loop(not on the loop), and the loop id is returned in id.
/// otherwise:on the loop(int(MGEdge* of parameter edge))will be returned, 
///			 and the loop id is returned in id.
size_t inside_inner_boundary(const MGPosition& uv, int& id)const;

///Test if (u,v) is inside the outer boundary.

///Inside the outer boundary means that inside outer_boudary_param() or not.
///***Caution***
///(1)This must not be used for faces that do not have perimeter or outer boundary
///loop.
///(2)inside_outer_boundary does not check about the inner loops.
///
///Function's return value is:
///  0:outside the outer boundary(not on a loop)
///  1:unknown
///  2:inside the outer boundary(not on a loop)
/// otherwise:on the outer boundary loop
size_t inside_outer_boundary(
	const MGPosition& uv
)const;

///Intersection with MGFSurface.
MGSSisects intersect(const MGFSurface& face2)const;

///Compute intersection points of an inner parameter line of this face and f2.
///The intersection point is used to compute surface to surface intersection lines.
///Function's return value is at most one intersection point in uvuv_list.
///One member of uvuv_list is (u1,v1,u2,v2), where (u1,v1) is a parameter of
///this surface and (u2,v2) is a parameter of surf.
MGPosition_list intersectInner(
	const MGFace& f2		///<The second face.
) const;

///Compute intersection points of an inner parameter line of this face and sf2.
///The intersection point is used to compute surface to surface intersection lines.
///Function's return value is at most one intersection point in uvuv_list.
///One member of uvuv_list is (u1,v1,u2,v2), where (u1,v1) is a parameter of
///this surface and (u2,v2) is a parameter of surf.
MGPosition_list intersectInner(
	const MGSurface& sf2		///<The second surface.
) const;

///isect_area_length() returns initial area length for the intersection
///line.
int isect_area_length() const{return surface()->isect_area_length();};

///Compute intersection points of this face's boundary(outer and inners) with
///face2. If intersection points are found and the boundary is a loop,
///the point's edge pointer(of this) will be stored in a member uvuv of uvuvs.
///uvuv[7] is the edge pointer. If the boundary is not a loop(that is, a perimeter of
///Surfaces), uvuv.sdim()==7 and an edge pointer is not returned.
///When uvuv.sdim()==8, the edge pointer of uvuv[7] is accessed through union mgEdgeP.
///uvuvs[i] is i-th intersection points.
int isect_boundary(
	const MGFSurface& face2,
	MGPosition_list& uvuvs,
	///<id1 and id2 are the ids of uvuv where this face's and f2's parameters
	///<are to be stored in a member of uvuvs.
	///<This face's (u,v) is stored in uvuv(id1) and (id1+1).
	///<f2's (u,v) is stored in uvuv(id2) and (id2+1).
	///<id2=0 if id1=2, and id2=2 if id1=0.
	int id1=0
)const;

///"isect1_incr_pline" is a dedicated function of isect_start_incr, will get
/// shortest parameter line necessary to compute intersection.
MGCurve* isect_incr_pline(
	const MGPosition& uv,	///<last intersection point.
	int kdt,				///<Input if u=const v-parameter line or not.
							///< true:u=const, false:v=const.
	double du, double dv,	///<Incremental parameter length.
	double& u,				///<next u value will be output
	double& v,				///<next v value will be output
	int incr=0			///<Incremental valuse of B-coef's id.
)const;

///Compute intersection points between the boundary of iid-th inner boundary
///of this face and face2 to compute intersections of face with face2.
///Function's return value is the number of ip's obtained before appending
///into uvuv_list, may not be equal to the enlarged size of uvuv_list.
int isect_incurves(
	const MGFSurface& face2,
	int iid,	///<Inner loop id of this face(from 0)
	MGPosition_list& uvuv_list,	///<intersection points will be appended.
		///<One member in the list is of sdim 8,
		///<(4,5,6) is the direction vector, and (7) is Edge pointer of the point.
	int id1			///<id of uvuv(a member of uvuv_list).
		///<uvuv(id1) for this face parameter uvuv(id2) for face2 parameter.
		///<id2=0 if id1=2, and id2=2 if id1=0.
)const;

///Compute intersection points of outer boundary curves of this face 
///with face2 to compute intersections.
///Function's return value is the number of ip's obtained(appended)
///into uvuv_list, may not be equal to the enlarged size of uvuv_list.
int isect_outcurves(
	const MGFSurface& face2,
	MGPosition_list& uvuv_list,	///<intersection points will be appended.
		///<One member in the list is of sdim 7,
		///<and the last three elements are the ip direction vector.
	int id1			///<id of uvuv(a member of uvuv_list).
		///<uvuv(id1) for this face parameter uvuv(id2) for srf or face2 parameter.
		///<id2=0 if id1=2, and id2=2 if id1=0.
)const;

///Compute intersection points between loop lp of this face and face2
///to compute intersections of face with face2.
///Function's return value is the number of ip's obtained before appending
///into uvuv_list, may not be equal to the enlarged size of uvuv_list.
int isect_lpcurves(
	const MGFSurface& face2,		///<srf!=null, face2=null.
	const MGLoop& lp,				///<Loop id of this face.
	MGPosition_list& uvuv_list,	///<intersection points will be appended.
		///<One member in the list is of sdim 8,
		///<(4,5,6) is the direction vector, and (7) is Edge pointer of the point.
	int id1			///<id of uvuv(a member of uvuv_list).
		///<uvuv(id1) for this face parameter uvuv(id2) for face2 parameter.
		///<id2=0 if id1=2, and id2=2 if id1=0.
)const;

///Compute intersection points between loop lpid of this face and face2
///to compute intersections of face with face2.
///Function's return value is the number of ip's obtained before appending
///into uvuv_list, may not be equal to the enlarged size of uvuv_list.
int isect_lpcurves(
	const MGFSurface& face2,
	int lpid,				///<Loop id of this face.
	MGPosition_list& uvuv_list,	///<intersection points will be appended.
		///<One member in the list is of sdim 8,
		///<(4,5,6) is the direction vector, and (7) is Edge pointer of the point.
	int id1			///<id of uvuv(a member of uvuv_list).
		///<uvuv(id1) for this face parameter uvuv(id2) for face2 parameter.
		///<id2=0 if id1=2, and id2=2 if id1=0.
) const;

///Compute intersection lines, given end points of the i.l.
///isectEP does not execute has_common() with srf2. Users of isectEP are recommended
///to execute has_common().
MGSSisects isectEP(
	MGPosition_list& uvuv_list,	///<End points list of the intersection.
		///<On return, uvuv_list.size() will be 0.
	const MGFSurface& fsrf2		///<The second surface.
) const;

///Obtain parameter u(kcod=1), or v(kcod=0) of the intersection point of
///v=x(const) line(kcod=1), or u=x(const) line (kcod=0) with
///all the face boundaries.
std::vector<double> isect1D_with_boundaries(
	double x,							///<coordinate value of kcod.
	int kcod					///<Coordinate kind, =0:u, =1: v.
)const;

///Obtain outer boundary curve expression as the combination of
///loop pointers and perimeter id's.
std::vector<MGFOuterCurve> outer_curve() const;

///Obtain perimeter i's parameter range.
///Let rvec be std::vector<MGInterval> of the fucntion's output,
///then rvec[j] is j-th parameter range of i-th perimeter. 
std::vector<MGInterval> perimeter_param_range(int i) const;

///Dedicated function of range.
///Will check if point (u,v) is inside inner boundary or not.
///If inside an inner boundary, obtain the closest point to the boundary.
MGPosition range_check_inner_boundary(const MGPosition& uv) const;

///Remove parameter uv from uvs that is outside face parameter range.
void remove_outside_param(MGPosition_list& uvs)const;

///Shade the object in world coordinates.
void shade(
	mgVBO& vbo,///<Target graphic object.
	const MGDrawParam& para,///<Parameter to draw.
	std::vector<SHLL_COM_EDGES>* polylines,///<Shared edge lines when this is a part of shell
	MGCL::DRAW_TARGET target=MGCL::SHADING///<Shading target.
)const;

///Proprietry routine for make_outer_boundary(), will append newed straight line
///edge to the end of lp if param (uv1 and uv2) are far away enough compared with error
///err.
void sl_edge_append(
	MGLoop*& lp,	///<lp that the edge should be append.
					///<If lp is null, new lp will be generated.
	int id,			///<perimeter num of the staight line edge.
	const MGPosition& uv2,	///<face parameter of the end point of sl.
	const MGPosition& uv1,	///<face parameter of the start point of sl.
	double err_sqr		///<Error square allowed to regard as same points.
);

///ボックス枠に囲まれる交点を持つUV曲線を生成する
void getTrimCrv(
	double uerror, double verror,///<u and v parameter error.
	const MGBox& box,	///<parameter (u,v) box of the surface.
	std::vector<UniqueCurve>& vecCrv	///<paremter curves will be output
) const;

friend class MGShell;
friend class MGSurface;

};

///@cond

///mgEdgeP is used in output of  isect_lpcurves() to store Edge pointer in
///uvuv_list.
union MG_DLL_DECLR mgEdgeP{const MGEdge* pointer; double doubleV;};

///@endcond

///Obtain the closest point from point uv to vector of curves.
///MGClosest_to_curves does not change wc_zero, and so calling program of
///MGClosest_to_curves should change it if necessary.
MG_DLL_DECLR MGPosition MGClosest_to_curves(
	const MGPosition& uv,				///<Point.
	const std::vector<UniqueCurve>& curves	///<vector of curves.
);

/** @} */ // end of TOPO group
#endif
