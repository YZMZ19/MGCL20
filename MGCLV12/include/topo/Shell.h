/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#ifndef _MGShell_HH_
#define _MGShell_HH_

#include <vector>
#include "mg/MGCL.h"
#include "mg/Default.h"
#include "mg/drawParam.h"
#include "mg/Position_list.h"
#include "mg/FPoint.h"
#include "topo/Complex.h"
#include "topo/CFisects.h"
#include "topo/Face.h"
#include "topo/HHisect.h"
#include "topo/HHisects.h"

class MGFSurface;
class MGIgesOfstream;
/** @file */
/** @addtogroup TOPO
 *  @{
 */

///MGShell is a composition of MGFace's(trimmed surface).

///See topology structure.
class MG_DLL_DECLR MGShell: public MGComplex{

public:

/////////Constructor/////////

///Default constructor.
MGShell()=default;
~MGShell();

///Copy and move constructor.
MGShell(const MGShell& shell) = default;
MGShell(MGShell&& shell) = default;

///Construct a shell of one face(face copy version).
MGShell(const MGFace& face);

///Construct a shell of one face.

///face is a pointer to newed face, and the ownership
///will be transfered to MGShell.
explicit MGShell(MGFace* face);

///Construct a shell of one face.

///Copy version of face.
MGShell(const MGFSurface& face);

///Construct a shell of one face.

///face is a pointer to newed face, and the ownership
///will be transfered to MGShell.
explicit MGShell(MGFSurface* face);

///Construct from boundary complex(i.e. MGLoop).
///This constructor takes the ownership of MGCell* in boundary.
MGShell(
	std::list<MGCell*> faces
			///<Boundary data of the super class MGComplex. List of faces.
);

///Assignment.
///When the leaf object of this and bnd2 are not equal, this assignment
///does nothing.
MGShell& operator=(const MGGel& gel2);
MGShell& operator=(const MGShell& shell)=default;
MGShell& operator=(MGShell&& shell) = default;

///Object transformation.
MGShell& operator+=(const MGVector& v){MGComplex::operator+=(v);return *this;};
MGShell& operator-=(const MGVector& v){MGComplex::operator-=(v);return *this;};
MGShell& operator*=(double scale){MGComplex::operator*=(scale);return *this;};
MGShell& operator*=(const MGMatrix& mat){MGComplex::operator*=(mat);return *this;};
MGShell& operator*=(const MGTransf& tr){MGComplex::operator*=(tr);return *this;};

/// Complexに平行移動を行ないオブジェクトを生成する。
///Translation of the Complex
MGShell operator+ (const MGVector& v) const;

/// Complexに逆方向の平行移動を行ないオブジェクトを生成する。
///Translation of the Complex
MGShell operator- (const MGVector& v) const{ return operator+(-v); };

///Complexのスケーリングを行い，Complexを作成する。
///Scaling of the Complex by a double.
MGShell operator* (double s) const;

/// 与えられた変換でComplexの変換を行い，Complexを作成する。
///Transformation of the Complex by a matrix.
MGShell operator* (const MGMatrix& mat) const;

/// 与えられた変換によってトランスフォームをおこないComplexを生成する。
///Transformation of the Complex by a MGTransf.
MGShell operator* (const MGTransf& tr) const;

///Complexのスケーリングを行い，Complexを作成する。
///Scaling of the Complex by a double.
MGShell operator/ (double s) const{ return (*this)*(1./s); };

///comparison
bool operator<(const MGShell& gel2)const;
bool operator<(const MGGel& gel2)const;

///String output.
std::ostream& toString(std::ostream& ostrm) const;

///IGES output function
///Function's return value is the directory entry id created.
int out_to_IGES(
	MGIgesOfstream& igesfile,
	int SubordinateEntitySwitch=0
)const;

///Test if this is an active boundary.
bool active() const{return false;};

///Append a face this is independent or connected to the already member faces of this
///shell. append_face() does not make any face join process.
///f must be a newed object, and the ownership will be transfered to this shell.
void append_face(MGFace* f){append_pcell(f);};

///Make a clone.
///Returned is pointer of newed object, must be deleted.
///When parent is specified, clone's parent is set to the parent.
MGShell* cloneWithParent(MGCell& parent) const;
MGShell* clone() const;

///Test if this is closed boundary.
bool closed() const{return false;};

///Compute the closest point from a point to this shell.
MGFPoint closest(const MGPosition& point) const;

///Judge if the display list for vmode is made or not.
//virtual bool displayList_is_made(MGCL::VIEWMODE vmode)const;

///////display member function.
virtual void display_arrows(mgSysGL& sgl)const;
virtual void display_control_polygon(mgSysGL& sgl)const;

///Get the face pointer from its iterator in MGComplex.
MGFace* face(iterator i);
const MGFace* face(const_iterator i)const;

///Get the face pointer from its pcell id in MGComplex.
MGFace* face(int i);
const MGFace* face(int i)const;

/// Return This object's typeID
long identify_type() const;

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

///Intersection of a shell and a curve.
MGCFisects isect(const MGCurve& curve) const;

///Intersection of two shells.
MGHHisects isect(const MGShell& shell2)const;

///Intersection of a shell and a face.
///This shell's face is face1 in HHisect and face2 is face.
///The name of isectFS is to get uniqueness to isect(const MGObject&).
MGHHisects isectFS(const MGFSurface& face) const;

///Make 2 types of display list of this gel(wire and shading).
void make_display_list(
	MGCL::VIEWMODE vmode=MGCL::DONTCARE
)const;

///Get manifold dimension.
int manifold_dimension() const{return 2;};

///Merge a face at free edges of this shell.

///Function's return value is
///   false: merge was not done because no common edges were found.
///   true: merge of the shell and the face was done.
///      Case of true includes that merge was done only partialy.
///      This case occurs when after some common edges are merged,
///      some common edges are found to have contradictionary direction.
///      Merge is done only for the first common edges found to have
///      same direction.
///When the function's return value is false, the face will not be added to
///the shell. Merge negates the input face direction when this shell's direction
///is not the same as the input face's one along the common edge(s).
///The second form is a pointer version. The face must be newed object
///and will be destructed after merge when function's return value is true
///(merge was processed).
///This shell may be dummy shell. In this case, the face is added to the 1st
///face in the shell.
bool merge_at_common_edge(const MGFace& face);
bool merge_at_common_edge(MGFace* face);
bool merge_at_common_edge(const MGFSurface& face);
bool merge_at_common_edge(MGFSurface* face);

///Get the number of faces included in this shell.
int number_of_faces()const{return number_of_pcells();};

///Test if a point is on the shell or not.

///Function's return value is true if the point is on the shell, and false if not.
///The point parameter of the shell is returned in fp if true is returned.
///If false is returned, the closest point of the shell will be returned in fp.
bool on(
	const MGPosition& point,	///<Target point.
	MGFPoint& fp				///<Shell's point parameter value.
) const;

///Obtain perpendicular points of a shell from a point.
std::vector<MGFPoint> perps(const MGPosition& point) const;

///Obtain the projected curves of a curve onto the shell.
///The direction of the projection is along the vector vec if the vec is not
///NULL, and normal to the shell if the vec is NULL.
///Output of 'project' is two kind of curves:
///one is general world coordinate curves(iline() of the MGHHisect members of lines),
///and the other is (u,v) curves of the parameter space of the surfaces
///(uvline1() of the MGHHisect members of lines ).
///*** uvline2() of the MGHHisect members of lines is a deque of length zero
///(not used).
///Function's return value is the number of curves obtained.
///When <0 is returned, some internal error occured.
int project(
	const MGCurve& crv,		///<given world coordinate curve to project.
	MGHHisects& lines,
			///<World coordinates (x,y,z) lines and (u,v) lines of the projected
			///<curves will be returned in lines.
	const MGVector& vec=mgNULL_VEC	///<projection direction.
			///<if vec = NULL, then projection that is normal to the shell.
)const;

///Shade the object in world coordinates.
void shade(
	mgVBO& vbo,
	const MGDrawParam& para,
	MGCL::DRAW_TARGET target=MGCL::SHADING
)const;

///Triangulate this object(MGShell, MGFace, or MGSurface is the target).
virtual void triangulate(
	const MGDrawParam& para,
	MGCL::TL_DATA_KIND dkind,
	std::vector<mgTL2Triangles>& trisVec
)const;

///Obtain boundary and main parameter lines of the FSurface.
///skeleton includes boundary() and inner parameter lines.
///density indicates how many inner parameter lines are necessary
///for both u and v directions.
std::vector<UniqueCurve> skeleton(int density=1)const;

///Obtain all the parameter curves at knots of u and v knot vector.
std::vector<UniqueCurve> skeleton_at_knots()const;

///Get the name of the class.
virtual std::string whoami()const{return "Shell";};

private:

///Get the partner point uvuv_id of the uvuv_list[j].
///If found, return true, if not, return false.
///Shell*(face or surface) version.
///The other one must not be shell, a surface or a face.
bool isect_partner(
	const MGPosition& uvuv,///<The point to get the partner.
	MGPosition_list* uvuv_list,
	int& j,		///<Partner edge's face id(of uvuv_list) will be returned.
	MGPosition_list::iterator& uvuvItr
		///<Partner point's iterator of uvuv_list[j] will be returned.
)const;

///Get the partner point uvuv_id of the uvuv_list[j].
///If found, return true, if not, return false.
///Shell*shell intersection version.
bool isect_partner(
	const MGPosition& uvuv,///<The point to get the partner.
	MGPosition_list* uvuv_list,
	int& i,
	int& j,		///<Partner edge's face id(of uvuv_list) will be returned.
	MGPosition_list::iterator& uvuvItr
		///<Partner point's iterator of uvuv_list[i+nf1*j] will be returned.
)const;

};

///Set up common edges world coordinate line data in polylines.
///polylines[i] is std::vector<SHLL_COM_EDGES> for face(i) for i=0,...,number_of_faces()-1.
void MG_DLL_DECLR set_up_shell_shade(
	const MGShell& shell,///<Target shell
   	MGDrawParam& para,///< input tessellation parameter.
				///< When para.maximum_edge_length_tess()<=0. is input, maximum_edge_length is
				///<computed and set in para.
	std::vector<UniqueLBRep>& boudaries,///<Container of the boundary data.
		///<Polylined boundaries will be held in boundaries.
	std::vector< std::vector<SHLL_COM_EDGES> >& polylines///<polylines[i] holds boundary
		///<data of shell.face(i)'s loops polylined boundaries data.
);

/** @} */ // end of TOPO group
#endif
