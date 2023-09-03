/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#ifndef _MGPlane_HH_
#define _MGPlane_HH_
#include <memory>
#include "mg/drawParam.h"
#include "mg/Position.h"
#include "mg/Unit_vector.h"
#include "mg/Surface.h"
#include "mg/SBRep.h"
#include "mgGL/VBO.h"

// MGPlane.h
// Header for class MGPlane
class MGBPointSeq;
class MGCCisect;
class MGCSisect;
class MGTransf;
class MGStraight;
class MGSBRep;
class MGIfstream;
class MGOfstream;
class MGFace;
/** @file */

/** @addtogroup GEO
 *  @{
 */

///MGPlane is infinite plane in 3D space.

///Using orthonormal two vector m_uderiv and m_vderiv in 3D space,
///plane function f(u,v)= m_root_point + u*m_uderiv + v*m_vderiv,
///where u and v are two parameter of surface representation. 
///m_d(the distance from the origin) and m_normal(plane's unit normal) are
///also hold to save the computation.
/// MGPlane�N���X�͂R������Ԃɂ����镽�ʂ�\���N���X�ł���B
/// MGPlane�N���X�ł͈ȉ��̂悤�ȃp�����[�^�\�����g�p���܂��B
/// Point(u,v) = m_root_point + u * m_uderiv + v * m_vderiv
class MG_DLL_DECLR MGPlane :public MGSurface{

	MGUnit_vector	m_normal;	///<Normal of the plane ���ʂ̖@���x�N�g��.
	double		m_d;
	///<Distance from the origin(0,0,0) ���ʂ̉A�֐��\��
	///< (m_d=ax+by+cz+....where m_normal=(a,b,c,....))
	///< ���Ȃ킿m_d�͌��_�ƕ��ʂ̋���.
	mutable std::unique_ptr<MGKnotVector> m_uknotV;
	mutable std::unique_ptr<MGKnotVector> m_vknotV;
	///<When knot_vector_u,v() is invoked, the knot vector will be set,
	///<These two variables must be initialize to 0(null).
	
protected:
	MGPosition	m_root_point;	///<A point on the plane,
	///<���ʂ̃p�����[�^�\���̊�_.
	MGVector	m_uderiv;	///<U direction vector,
	///<���ʂ̃p�����[�^�\���̂�����.
	MGVector	m_vderiv;	///<V direction vector,
	///< ���ʂ̃p�����[�^�\���̂�����.

public:

///Vector translation.
MG_DLL_DECLR friend MGPlane operator+ (const MGVector& v, const MGPlane& pl);

///Scaling.
MG_DLL_DECLR friend MGPlane operator* (double scale, const MGPlane& pl);

////////Special member functions/////////
MGPlane();
virtual ~MGPlane()=default;
MGPlane(const MGPlane&);///Copy constructor.
MGPlane& operator= (const MGPlane&);///Copy assignment.
MGPlane(MGPlane&&)=default;		///Move constructor.
MGPlane& operator= (MGPlane&&)=default;///Move assignment.

/// Construct a plane by changing this space dimension or ordering the coordinates.
MGPlane(
	int dim,			///< New space dimension.
	const MGPlane& plane,///< Original Plane.
	int start1=0, 		///< Destination order of new Surface.
	int start2=0		///< Source order of original Surface.
); 

/// Construct a plane from the coefficients of the plane equation.
///     a*x+b*y+c*z=d.
/// Coefficients a,b,c,d are provided by double array g[4].
MGPlane(
	const double g[4],	///<coefficients g[0]=a, b=g[1], c=g[2], d=g[3].
	const double* root_point=0///<When root_point!=0, root_point[.] are (x,y,z) values of the root point.
);

///Plane from normal of the plane and the distance from the origin(0,0,0).
MGPlane(
	const MGUnit_vector& normal,///<Normal of the plane.
	double d			///<distance from origin of the plane.
						///<When normal=(a,b,c,....), d=a*x+b*y+c*z+.... .
);

///Plane from a point on the plane and the normal.
MGPlane(
	const MGUnit_vector& normal,
	const MGPosition& p
);

///Plane from a point and a straight line on the plane.
MGPlane(
	const MGStraight& st,
	const MGPosition& point
);

///Plane from a point on the plane, u and v direction vector of the plane.

///The plane's (uderiv, vderiv) are transfomred to a orthnormal system.
///***** This is the fundamental constructor.*****
MGPlane(
	const MGVector& uderiv,
	const MGVector& vderiv,
	const MGPosition &origin
);

/// Construct a plane by interpolating two planes.

/// If two planes intersect, interpolation is rotation around the
/// intersection line.
///If two planes are parallel, interpolation is parallel move.
MGPlane(
	const MGPlane& plane1,///<Target plane1 1.
	const MGPlane& plane2,///< Plane 2.
	double t	///<Input ratio.
				///< When t=0, the plane will be plane1,
				///< When t=1, the plane will be plane2.
);

///Construct a plane from three points on the plane.
MGPlane(
	const MGPosition& P1,
	const MGPosition& P2,
	const MGPosition& P3
);


///Assignment.
///When the leaf object of this and srf2 are not equal, this assignment
///does nothing.
MGPlane& operator=(const MGGel& gel2);
MGPlane& operator=(MGGel&& gel2);

///Transformation object construction
MGPlane operator+ (const MGVector& v) const;
MGPlane operator- (const MGVector& v) const;
MGPlane operator* (double scale) const;
MGPlane operator* (const MGMatrix& mat) const;
MGPlane operator* (const MGTransf& tr) const;

///Object transformation.
MGPlane& operator+=(const MGVector& v);
MGPlane& operator-=(const MGVector& v);
MGPlane& operator*=(double scale);
MGPlane& operator*=(const MGMatrix& mat);
MGPlane& operator*=(const MGTransf& tr);

///Comparison of two curves.
bool operator==(const MGPlane& gel2)const;
std::partial_ordering operator<=>(const MGPlane& gel2)const;

//gel2 must be the same class as this.
bool equal_test(const MGGel& gel2)const override;

//gel2 must be the same class as this.
std::partial_ordering ordering_test(const MGGel& gel2)const override;

///Output to IGES stream file(=PD190).
int out_to_IGES(
	MGIgesOfstream& igesfile,
	int SubordinateEntitySwitch=0
)const;

///String output function.
///Output to ostream �����o�f�[�^��W���o�͂ɏo�͂���B
std::ostream& toString(std::ostream &) const;

//////////////Member function �����o�֐�////////////

///Gets parameters(a,b,c,d) of the plane expression a*x+b*y+c*z=d.
void abcd(double g[4]) const;///< g[.]=(a,b,c,d)

///Box that includes limitted plane by box.
MGBox box_limitted(
	const MGBox& uvrange	///< Parameter Range of the surface.
) const;

///Obtain ceter coordinate of the geometry.
virtual MGPosition center() const{return m_root_point;};

///Changing this object's space dimension.
void change_dimension(
	int sdim,			///< new space dimension
	int start1=0, 		///< Destination order of new object.
	int start2=0 		///< Source order of this object.
);

///Change parameter range, be able to change the direction by providing
///t1 greater than t2.
/////********** MGPlane does not accept change_range, does nothing.****///
MGPlane& change_range(
	int is_u,				///<if true, (t1,t2) are u-value. if not, v.
	double t1,				///<Parameter value for the start of original. 
	double t2				///<Parameter value for the end of original. 
){return *this;};

///Change root point.
void change_root_point(
	const MGPosition& new_point
);

///Compute the closest point parameter value (u,v) of this surface
///from a point.
MGPosition closest(const MGPosition& point) const;

///Construct new surface object by copying to newed area.
///User must delete this copied object by "delete".
MGPlane* clone() const;

///Return minimum box that includes whole of the surface.
///Compute the box from the scratch.
void compute_box(MGBox& bx) const;

///Construct new surface object by changing
///the original object's space dimension.
///User must delete this copied object by "delete".
MGPlane* copy_change_dimension(
	int sdim,			///< new space dimension
	int start1=0, 		///< Destination order of new line.
	int start2=0 		///< Source order of this line.
)const;

void display_arrows(mgSysGL& sgl)const;

///Return the distace between plane and the origin(0,0,0).
double distance() const{return m_d;};

///Return the distace between plane and the point.
double distance(const MGPosition& point) const;

///Draw 3D curve in world coordinates.
///The object is converted to curve(s) and is drawn.
virtual void drawWire(
	mgVBO& vbo,///<The target graphic object.
	int line_density=1	///<line density to draw a surface in wire mode.
)const{drawWirePlane(vbo,line_density, MGCL::WIRE);};

///Draw 3D curve in world coordinates.
///The object is converted to curve(s) and is drawn.
void drawWirePlane(
	mgVBO& vbo,///<The target graphic object.
	int line_density=1,	///<line density to draw a surface in wire mode.
	MGCL::DRAW_TARGET target= MGCL::WIRE///<target vbo element.
)const;

///Shade the object in world coordinates.
virtual void shade(
	mgVBO& vbo,///<The target graphic object.
	const MGDrawParam& para,///<Parameter to draw.
	MGCL::DRAW_TARGET target= MGCL::SHADING///<target vbo element.
)const{drawWirePlane(vbo,1,target);};

///Evaluate surface data.
MGVector eval(
	double u, 	///< U Parameter value of the surface.
	double v	///< V Parameter value of the surface.
	, int ndu=0	///< Order of derivative along u.
	, int ndv=0	///< Order of derivative along v.
)const;

///Evaluate surface data.
MGVector eval(
	const MGPosition& uv	///< Parameter value of the surface.
	, int ndu=0			///< Order of derivative along u.
	, int ndv=0			///< Order of derivative along v.
)const{return eval(uv.ref(0),uv.ref(1),ndu,ndv);}

///Evaluate right continuous surface data.
///Evaluate all positional data, 1st and 2nd derivatives.
void eval_all(
	double u, 	///< U Parameter value of the surface.
	double v,	///< V Parameter value of the surface.
	MGPosition& f,		///< Positional data.
	MGVector&   fu,		///< df(u,v)/du
	MGVector&   fv,		///< df/dv
	MGVector&   fuv,	///< d**2f/(du*dv)
	MGVector&   fuu,	///< d**2f/(du**2)
	MGVector&   fvv		///< d**2f/(dv**2)
)const;

/// Exchange parameter u and v.
MGSurface& exchange_uv();

/// Return This object's typeID
long identify_type() const;

bool in_range(double u, double v) const{ return 1;};
bool in_range(const MGPosition& uv) const{ return 1;};

///The following two function will be used in perps or isect
///to decide how many division of the surface along u or v direction
///should be applied before using perp_guess or isect_guess.
int intersect_dnum_u() const{ return 2;};
int intersect_dnum_v() const{ return 2;};

///Compute intesection of Plane and Curve.
MGCSisects isect(const MGCurve& curve)const;
MGCSisects isect(const MGEllipse& curve)const;

///Surface �� Surface �̌�������߂�B
///Surface and Surface intersection.
MGSSisects isect(const MGSurface& srf2)const;
MGSSisects isect(const MGFSurface& srf2)const;
MGSSisects isect(const MGFace& f2)const;
MGSSisects isect(const MGPlane& srf2)const{ return intersectPl(srf2); };

MGHHisects isect(const MGShell& shell)const;

///isect_startHPL compute one intersection line of two surfaces, this and sf2,
/// given starting intersetion point uvuv((u1,v1) of this and (u2,v2) of sf2).
/// isect_startHPL halts the computation when intersection
/// reached to a boundary of this or sf2, or reached to one of the points
/// in uvuv_list.
///The function's return value is:
/// =0: Intersection was not obtained.
/// !=0: Intersection was obtained as follows:
///    =1: End point is a point on a perimeter of one of the surfaces.
///    =3: End point is one of boundary points in uvuv_list.
int isect_startHPL(
	const MGPosition& uvuv_startIn, ///<Starting point of the intersection line.
	MGPosition_list& uvuv_list,	///<isect_startHPL will halt when ip reached one of 
		///<the point in uvuv_list. isect_startHPL does not change uvuv_list(actually
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
) const;

///Return knot value of (infinite-minus, infinite-plus)
double knot_u(int)const ;
double knot_v(int)const ;

///Returns the u knot vector.
const MGKnotVector& knot_vector_u() const;
MGKnotVector& knot_vector_u();

///Returns the v knot vector.
const MGKnotVector& knot_vector_v() const;
MGKnotVector& knot_vector_v();

///Make a display list of this gel.
virtual void make_display_list(
	MGCL::VIEWMODE vmode=MGCL::DONTCARE
)const;

///Negate the normal of the plane,���ʂ𔽓]����B�m�[�}�����t�����ɂ���.
void negate(
	int is_u///< Negate along u-direction if is_u is ture, else along v-direction.
);

///Obtain parameter value if this surface is negated by "negate()".
/// Negate along u-direction if is_u is ture,
/// else along v-direction.
MGPosition negate_param(const MGPosition& uv, int is_u=1)const;

///Return the normal of the plane, ���ʂ̖@����ԋp����.
MGVector normal(double u, double v) const{return m_normal;}
MGVector normal(const MGPosition& uv) const{return m_normal;}
const MGUnit_vector& normal() const{return m_normal;}

///Update this plane so that m_uderiv, m_vderiv, m_normal
/// construct a orthonormal system.
void normalize();

///Surface offset. positive offset value is offset normal direction.
///the radius of curvature is larger than offset value.line_zero() is used.
///C1�A���Ȗʂ̈��I�t�Z�b�g�֐�
///�I�t�Z�b�g�����́A�m�[�}�������𐳂Ƃ���B�g�������X��line_zero()���g�p���Ă���B
///�߂�l�́A�I�t�Z�b�g�ʂ̃I�[�g�|�C���^�[���ԋp�����B
std::unique_ptr<MGSurface> offset_c1(
	double ofs_value,	///<�I�t�Z�b�g��.
	int& error			///<�G���[�R�[�h 0:���� -1:�ʂɂ��ꂪ����
						///< -2:�ȗ����a�ȏ�̃I�t�Z�b�g�s�� -3:�ʐ����R���X�g���N�^�G���[.
) const;

///Test if a point is on the plane. If on the plane, return true.
/// �w��_���ʏ�ɂ��邩���ׂ�B�i�ʏ�Ȃ��true�j.
bool on(const MGPosition& point) const;
bool on(
	const MGPosition& point,	///<A point to test �w��_.
	MGPosition& puv	///<Parameter value of point on the plane will be returned either
					///<point is on the plane or not.
) const;

///Test if a straight line is on the plane. Return true if on.
/// ���������ʏ�ɂ��邩���ׂ�B�i���ʏ�Ȃ��true�j
bool on(const MGStraight&) const;

///Test if input (u,v) is parameter value on a perimeter of the surface.
///If u or v is on a perimeter, (u,v) will be updated to the perimeter value.
bool on_a_perimeter(
	double& u,///<Surface parameter
	double& v,///<Surface parameter (u,v).
	int& perim_num	///<if function returns true,
						///<the perimete rnumber is output.
)const{return false;};

///Test if input (u,v) is on the perimeter perim_num.
///If u or v is on a perimeter, true will be returned.
bool on_the_perimeter(
	int perim_num,	///<a perimete number is input.
	double u,	///<Surface parameter
	double v	///<Surface parameter (u,v).
)const{return false;};

///Return plane parameter value of a point on the plane. 
///If input point is not on the plane, returned is
///the nearest point parameter of the plane.
MGPosition param(
	const MGPosition&	///< �w��_
) const;

///Obtain parameter space error.
double param_error() const;
double param_error_u() const;
double param_error_v() const;

/// Return ending parameter value.
virtual double param_e_u() const{return mgInfiniteVal;}
virtual double param_e_v() const{return mgInfiniteVal;}

/// �p�����[�^�͈͂�Ԃ��B
///Return parameter range of the plane(Infinite box).
virtual MGBox param_range() const;

/// Return starting parameter value.
virtual double param_s_u() const{return -mgInfiniteVal;}
virtual double param_s_v() const{return -mgInfiniteVal;}

/// Compute parameter curve.
///Returned is newed area pointer, and must be freed by delete.
virtual MGCurve* parameter_curve(
	int is_u		///<Indicates x is u-value if is_u is true.
	, double x		///<Parameter value.
					///<The value is u or v according to is_u.
)const;

///Compute part of the surface limitted by the parameter range bx.
///bx(0) is the parameter (us,ue) and bx(1) is (vs,ve).
///That is u range is from us to ue , and so on.
MGSurface* part(
	const MGBox& uvbox,///<This plane parameter range.
	int multiple=0	///<Indicates if start and end knot multiplicities
					///<are necessary. =0:unnecessary, !=0:necessary.
)const;

/// i must be < perimeter_num().
///When perimeter_num()==0, this function is undefined.
virtual MGCurve* perimeter_curve(int i) const{return 0;};

///Return how many perimeters this surface has.
virtual int perimeter_num() const{return 0;};

/// Construct perimeter (u,v) parameter position.
/// i is perimeter number:
/// =0: v=min line, =1: u=max line, =2: v=max line, =3: u=min line
/// t is perimeter parameter line's parameter value of u or v.
virtual MGPosition perimeter_uv(int i,double t) const;

/// �^����ꂽ�_�ɂ����Ƃ��߂��ʏ�̓_��ԋp����B�p����
/// �[�^�l���ԋp����B
///Return the nearest point of the plane from P.
///Function's return value is always true.
int perp_point(
	const MGPosition& P,///< �^����ꂽ�_.
	MGPosition& uv,		///<Parameter value of the plane will be output.
	const MGPosition* uvguess=NULL	///< guess.
) const;

///Return all(actually one) foots of perpendicular straight lines from P.
MGPosition_list perps(
	const MGPosition& P	///< Point of a space(�w��_).
) const;

///Test if the surface is planar or not.
///Returned is 0(false) if this is not planar, 1(true) if this planar.
int planar(
	MGPlane& plane,	///<Plane that might be closest to this.
					///<Plane is always output even if not planar.
	double& deviation	///<maximum deviation of this from the output plane.
) const;

///Test if part of the surface is planar or not within the tolerance tol.
///The part of the surface is input by the surface parameter range uvbox.
///Returned is 0(false) if this is not planar, 1(true) if planar.
///For plane, planar always returns true.
int planar(
	const MGBox& uvbox,///<This surface parameter range.
	double tol,	///<maximum deviation allowed to regard the sub surface as a plane.
	int* divideU=0///<Direction to subdivide will be output, if this was not planar,
				///<=1: u direction, =0: v direction.
)const;

///�^����ꂽ�Ȑ������g�ɖʒ��܂��̓x�N�g�����e���ċȐ����X�g�����߂�B
///���e�Ȑ��͖ʏ�̃p�����[�^�Ȑ���3�����Ȑ��Ƃ��Ă��ꂼ�ꏇ�ԂɁA
///vec_crv_uv, vec_crv�Ɋi�[�����B
///uv�Ȑ��̃g�������X��rc_zero()���A3�����Ȑ���line_zero()�����ꂼ��g�p���Ă���B
///get perpendicular or vector projection curve list.
///uv projection curves are put into vec_crv_uv(rc_zero() is used),
///3d projection curves are put into vec_crv(line_zero() is used) respectively.
///�����F
///		const MGCurve&			crv,		(I/ )	given curve.
///		std::vector<UniqueCurve>&	vec_crv_uv,		( /O)	uv projection curve.
///		std::vector<UniqueCurve>&	vec_crv,		( /O)	3d projection curve.
///		const MGVector&			vec=mgNULL_VEC	(I/ )	projection vector.
///�߂�l�F
///		���e�Ȑ��̐�:		���e�Ȑ������܂���
///		0:			���e�Ȑ������܂�Ȃ�����
///		-1:			���������G���[
///		-2:			���������G���[�i�������Ȃ������j
///�ǋL�F����vec���^�����Ȃ�(null�j�Ƃ��A�ʒ����e����B
///Obtain the projected curve of a curve onto the surface.
///The direction of the projection is along the vector vec if the vec is not
///NULL, and normal to the surface if the vec is NULL.
///Output of 'project' is two kind of curves:
///one is general world coordinate curves('vec_crv'), and the other is (u,v) curves of
///the parameter space of the surfaces(vec_crv_uv).
///vec_crv_uv.size() is equal to vec_crv.size(). Let the size be n, then
/// (vec_crv_uv[i], vec_crv[i]) is one pair for 0<=i<n.
///Function's return value is:
/// >=0: number of curves obtained, <0 : Some error detected.
int project(
	const MGStraight& sl,	///<given curve.
	std::vector<UniqueCurve>& vec_crv_uv,	///<uv projection curve will be appended.
	std::vector<UniqueCurve>& vec_crv,	///<3d projection curve will be appended.
	const MGVector& vec= mgNULL_VEC	///<projection vector,
						///<if vec = NULL then calculate perpendicular project.
)const;

///�^����ꂽ�Ȑ������g�ɖʒ��܂��̓x�N�g�����e���ċȐ����X�g�����߂�B
///���e�Ȑ���3�����Ȑ��Ƃ���vec_crv�Ɋi�[�����B
///3�����Ȑ���tolerance��line_zero()���g�p���Ă���B
///get perpendicular or vector projection curve list.
///3d projection curves are put into vec_crv(line_zero() is used).
///�����F
///		const MGCurve&			crv,		(I/ )	given curve.
///		std::vector<UniqueCurve>&	vec_crv,		( /O)	3d projection curve.
///		const MGVector&			vec=mgNULL_VEC	(I/ )	projection vector.
///�߂�l�F
///		���e�Ȑ��̐�:		���e�Ȑ������܂���
///		0:			���e�Ȑ������܂�Ȃ�����
///		-1:			���������G���[
///		-2:			���������G���[�i�������Ȃ������j
///�ǋL�F����vec���^�����Ȃ�(null�j�Ƃ��A�ʒ����e����B
///Obtain the projected curve of a curve onto the surface.
///The direction of the projection is along the vector vec if the vec is not
///NULL, and normal to the surface if the vec is NULL.
///Output of 'project' is general world coordinate curves('vec_crv')
///Function's return value is:
/// >=0: number of curves obtained, <0 : Some error detected.
int project(
	const MGCurve& crv,				///<given curve.
	std::vector<UniqueCurve>& vec_crv,	///<3d projection curve will be appended.
	const MGVector& vec = mgNULL_VEC///<projection vector,
							///<if vec = NULL then calculate perpendicular project.
) const;

/// ���̓p�����[�^���p�����[�^�͈͂ł܂�߂ĕԋp����.
///Round the input uv into parameter range of the plane, 
///return the same value as input.
MGPosition range(const MGPosition& uv) const{return uv;}

/// �^����ꂽ���ʂƂ̊֌W��Ԃ�.
///Relation of two planes.
MGPSRELATION relation(
	const MGPlane&,	///<Second plane
	MGStraight&		///<Intersection line of the two planes
					///<will be output if intersect.
) const;

/// �^����ꂽ�����Ƃ̊֌W��Ԃ��B
///Relation between plane and straight line.
MGPSRELATION relation(
	const MGStraight&,	///<Straight line
	MGCSisect&			///<Intersection point if intersect.
) const;

/// ���ʂ̃p�����[�^�\���̋N�_��ԋp����B
///Return root point of the plane.
const MGPosition& root_point() const{return m_root_point;}

///Return the space dimension.
int sdim() const;

///Obtain boundary and main parameter lines of the FSurface.
///skeleton includes boundary() and inner parameter lines.
///density indicates how many inner parameter lines are necessary
///for both u and v directions.
virtual std::vector<UniqueCurve> skeleton(int density=1)const;

///Obtain all the parameter curves at knots of u and v knot vector.
std::vector<UniqueCurve> skeleton_at_knots()const;

///split this fsurface at the parameter param.
void split(
	double param,///<parameter value of this fsurface. if is_u is true, param is u-value,
				///<else v-value.
	bool is_u,	///<indicates if param is u or v of the surface parameter (u,v).
	std::vector<std::unique_ptr<MGFSurface>>& surfaces///<splitted surfaces will be output.
)const;

///Construct a tangent plane LBRep of a surface's perimeter.

///Output MGLBRep�@tp's knot vector is the same as srf's perimeter's 
///if the perimeter is MGLBRep. Let ft(t) is pcrv and g(t) be the output tp ,
///then g(t)=normal at srf(u,v) where (u,v)=f(t).
void TPatPerimeter(
	int perimeterNum,
	MGLBRep& tp///<Obtained tp is output.
)const override{tp=MGLBRep();};

/// �Ȗʃ^�C�v��ԋp����B
///Return surface type of the plane.
virtual MGSURFACE_TYPE type() const{return MGSURFACE_TYPE::MGSURFACE_PLANE;}

/// ���ʂ̃p�����[�^�\���̂�������ԋp����B
///Return u-direction vector of the plane.
const MGVector& u_deriv() const{return m_uderiv;}

/// ���ʂ̃p�����[�^�\���̂�������ԋp����B
///Return v-direction vector of the plane.
const MGVector& v_deriv() const{return m_vderiv;}

/// �_�𕽖ʂɓ��e�����_�̕��ʂ̃p�����[�^�\��(u,v)�����߂�B
///Return uv parameter of the point projected from point p to the plane.
MGPosition uv(const MGPosition& p) const;

/// Vector�𕽖ʂɓ��e����Vector�̕��ʂ̃p�����[�^�\��(u,v)�����߂�B
///Return uv parameter of the point projected from point p to the plane.
///p is the end of the vector v originated from the root_point().
MGVector uv(const MGVector& v) const;

///Get the name of the class.
std::string whoami()const{return "Plane";};

protected:

///Test if the surface is flat or not within the parameter value rectangle of uvbox.
///Function's return value is:
///	true: if the surface is flat
///  false: if the surface is not falt.
///When this is not falt, the direction that indicates which direction the surface
///should be divided will be output.
///***** the flatness is tested only approximately. This is for exclusive use of
///planar().
bool flat(
	const MGBox& uvbox,///<Target parameter box.
	double tol,		///<Tolerance allowed to regard flat
					///<(Allowed distance from a plane).
	int& direction,	///<   1: u-direction is more non flat,
					///<   0: v-direction is more non flat.
	MGPosition& P,	///<Position of the flat plane will be output.
	MGUnit_vector& N///<Normal of the flat plane will be output.
)const;

///Test if the surface is flat or not within the parameter value rectangle of uvbox.
///Function's return value is:
///	true: if the surface is flat
///  false: if the surface is not falt.
///When this is not falt, the direction that indicates which direction the surface
///should be divided will be output.
///This is the same as flat except that this does not have the arguments P, N.
///***** the flatness is tested only approximately. This is for exclusive use of
///tessellation.
bool flat_tess(
	double u0,///<u range from u0.
	double u1,///<  to u1
	double v0,///<v range from v0.
	double v1,///<  to v1
	double tol,		///<Tolerance allowed to regart flat
					///<(Allowed distance from a plane).
	bool& direction	///<   1: u-direction is more non flat,
					///<   0: v-direction is more non flat.
)const;

///Intersection of Surface and a straight line.
MGCSisects isectSl(
	const MGStraight& sl,///< Target straight.
	const MGBox& uvbox=mgNULL_BOX ///<indicates if this surface is restrictied to the parameter
					///<range of uvbox. If uvbox.is_null(), no restriction.
)const override;

///�����o�f�[�^��ǂݍ��ފ֐�
void ReadMembers(MGIfstream& buf);
	
///�����o�f�[�^���������ފ֐�
void WriteMembers(MGOfstream& buf) const;

private:

///get the display arrow length vector.
void get_uv_display_vector(
	MGVector& u,
	MGVector& v
)const;

///Compute the intersection line of this and the plane pl.
MGSSisects intersectPl(const MGPlane& pl) const override;

///isect_area_length() returns initial area length for the intersection
///line.
int isect_area_length() const{return 10;};

/// "isect_guess" computes one intersection point of surface and a curve,
/// given initail guess parameter values of surface and curve.
///Function's return value is 1(true) if i.p. obtained, and 0(false) if not.
int isect_guess(
	const MGCurve& crv,		///<Curve.
	const MGPosition& uvi,	///<Input initial guess parameter value
						///< of the i.p. of the surface,
	double ti,			///<Input initial guess parameter value of the line.
	MGPosition& uv,		///< Output parameter value obtained. 
	double& t			///< Output parameter value obtained. 
)const;

/// "isect_guess_straight" computes one intersection point of surface and
///a straight line, given initail guess parameter values of the surface and 
///the straight line.
///Function's return value is 1(true) if i.p. obtained, and 0(false) if not.
int isect_guess_straight(
	const MGStraight& sl,	///<Straight line.
	double ti,			///<Initial guess parameter value of the straight.
	const MGPosition& uvi,	///<Input initial guess parameter value,
						///< of the i.p. of the surface. 
	double& t,			///<Straight parameter obtained.
	MGPosition& uvout		///<Surface parameter value obtained(u,v). 
) const;

///"isect1_incr_pline" is a dedicated function of isect_start_incr, will get
/// shortest parameter line necessary to compute intersection.
MGCurve* isect_incr_pline(
	const MGPosition& uv,	///<last intersection point.
	int kdt,				///<Input if u=const v-parameter line or not,
							///< true:u=const, false:v=const.
	double du, double dv,///<Incremental parameter length.
	double& u,				///<next u value will be output
	double& v,				///<next v value will be output
	int incr=1		///<Incremental valuse of B-coef's id.
) const;

///"isect_inner_dt" is a dedicated function of isect_startPt,
/// comutes adequate incremental parameter value(du,dv) and parameter line kind
///kdt(u=const or v=const).
void isect_inner_dt(
	int n,				///<num of i.p. obtained so far(not include uvnow).
	const MGPosition& uvnow,///<intersection point obtained last(of this).
	double& du, double& dv,	///<incremental length from previous to uvnow is input.
				///<New du or dv will be output according to kdt's return value.
	int& kdt,	///<Parameter kind used so far is input, will be output as:
				///<=1:parameter line kind(u=const), =0: v=const,
				///<=-1:should halt computation since incremental value is zero.
	double acuRatio=1.	///<Accurate ratio.
) const;

///Return order of intersection line order of MGLBRep.
///The default is 4.
int isect_order() const{return 4;}

///Compute intersections with MGLBRep lb that does not have C0 continuity in it.
MGCSisects isect_withC1LB(const MGLBRep& lb)const;

///isect with SurfCurve whose m_curve is not a MGTrimmedCurve of MGCompositeCurve.
MGCSisects isect_with_noCompoSC(const MGSurfCurve& curve2)const;

///�I�t�Z�b�g����T���v���|�C���g��1�p�b�`���Ƃ̕����������߂�
///�S�Ẵp�b�`���̕������ōő�̒l��Ԃ�
int offset_div_num() const;

///Obtain 1D surface rep. of this surf which can be used for
///isect(const MGPlane& pl). This surf1D is used in isect for
///the argument of isect_startPlane, which will use surf1D to compute isect(pl).
///surf1D=0.(intersection with x=0. plane) is the intersection lines.
std::unique_ptr<MGSBRep> surf1D(const MGPlane& pl)const override{
	assert(false); return std::unique_ptr<MGSBRep>(nullptr);
};

friend class MGStraight;

};

/** @} */ // end of GEO group
#endif
