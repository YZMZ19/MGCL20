/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#ifndef _MGObject_HH_
#define _MGObject_HH_

#include "mg/MGCL.h"
#include "mg/Box.h"
#include "mg/Position.h"
#include "mg/drawParam.h"
#include "mg/isects.h"
#include "mg/AttribedGel.h"
#include "mgGL/Appearance.h"

//
//Define MGObject Class.
class MGIfstream;
class MGOfstream;
class MGVector;
class MGMatrix;
class MGTransf;
class MGObject;
class MGGeometry;
class MGPoint;
class MGStraight;
class MGCurve;
class MGFSurface;
class MGSurface;
class MGFace;
class MGShell;
class MGAppearance;
class mgVBO;

/** @defgroup MGObjectRelated Object Related class
 *  MGObject is top abstract class for MGPoint, MGCurve, and MGSurface.
 *  @{
 */

///Is an abstract class which represents a whole geometry and a topology.
class MG_DLL_DECLR MGObject:public MGAttribedGel{
private:
	std::unique_ptr<MGAppearance> m_appearance;///<MGAppearance pointer. .

protected:
	mutable MGBox m_box;

public:

////////Special member functions/////////
MGObject():m_appearance(nullptr){;};	///void constructor.
virtual ~MGObject();///Destructor.
MGObject(const MGObject&);///Copy constructor.
MGObject& operator= (const MGObject& obj2){return set_object(obj2);};;///Copy assignment.
MGObject(MGObject&&);		///Move constructor.
MGObject& operator= (MGObject&&);///Move assignment.

///Object transformation.
virtual MGObject& operator+=(const MGVector& v)=0;
virtual MGObject& operator-=(const MGVector& v)=0;
virtual MGObject& operator*=(double scale)=0;
virtual MGObject& operator*=(const MGMatrix& mat)=0;
virtual MGObject& operator*=(const MGTransf& tr)=0;

////////////Member Function////////////

///Return minimum box that includes whole of the geometry.
const MGBox& box() const;
bool box_is_null()const;

///Output virtual function.
virtual std::ostream& toString(std::ostream&) const;

///Get the MGAppearance pointer of this object. If not defined, null will be
///returned.
///See ensure_appearance().
MGAppearance* appearance(){return m_appearance.get();};
const MGAppearance* appearance()const{return m_appearance.get();};

///Construct new object by copying to newed area.
///User must delete this copied object by "delete".
virtual MGObject* clone() const override =0;

///Compute box of the geometry.
///Compute the box from the scratch.
virtual void compute_box(MGBox& bx) const = 0;

///display member function.
virtual void display_arrows(mgSysGL& sgl)const{ ; };
virtual void display_break_points(mgSysGL& sgl)const{ ; };
virtual void display_control_polygon(mgSysGL& sgl)const{ ; };
virtual void display_curvatures(
	mgSysGL& sgl,	///<Sgl to make the pictures in.
	int		density,///<densitiy of the graph.
	bool	use_radius,///<true:radius display, false:curvature display.
	double	scale = 1.	///<scaling of the graph.
)const{	;};

///Draw the object in wire mode, in the world coordinates.
///The object is converted to curve(s) and is drawn.
virtual void drawWire(
	mgVBO& vbo,///<The target graphic object.
	int line_density=1	///<line density to draw a surface in wire mode.
)const=0;

///Compute the intersections of two objects.
///Intersections are obtained from two objects, which are known using
///the MGisects::object1() and object2().
///****NOTE****
///When two objects' manifold dimension are the same, object1 is this object
///at the invocation of MGObject::intersection(), and object2 is the argument
///object.
///However, their manifold dimension are not the same, object1 is always
///the lower dimension's object and object2 is the higer dimension's object.
virtual MGisects isect(const MGObject& obj2)const=0;
MGisects isect(const MGFSurface& f)const;
virtual MGisects isect(const MGPoint& obj2)const{return MGisects();};

///Shade the object in world coordinates.
virtual void shade(
	mgVBO& vbo,
	const MGDrawParam& para,
	MGCL::DRAW_TARGET target= MGCL::SHADING
)const{drawWire(vbo);};

///make this group has appearance and get the MGAppearance pointer.
///See appearance().
MGAppearance* ensure_appearance();

///Make a display list of this gel.
virtual void make_display_list(
	MGCL::VIEWMODE vmode=MGCL::DONTCARE
)const;

///Test if this and 2nd object has common area about their box(),
///taking error into account.
bool has_common(const MGObject& obj2) const;

///Compute the parameter value of the closest point from the straight to
///this object.
///sl is the eye projection line whose direction is from yon to hither, and if
///sl had multiple intersection points, The closest point to the eye will be selected.
virtual MGPosition pick_closest(const MGStraight& sl)const{
	return MGPosition();
}

///Remove the MGAppearance of this MGAttribedGel.
std::unique_ptr<MGAppearance> remove_appearance() override;

///Get the MGFSurface pointer if this is MGSurface or MGFace.
virtual const MGFSurface* fsurface()const{return (const MGFSurface*)0;};
virtual MGFSurface* fsurface(){return (MGFSurface*)0;};

//set the copy of appr2 to this MGAttribedgel.
void set_appearance(const MGAppearance& appr2)override;
void set_appearance(MGAppearance* appr2)override;

///Transform the gel by the argument.
virtual void transform(const MGVector& v){(*this)+=v;};///translation
virtual void transform(double scale){(*this)*=scale;};///scaling.
virtual void transform(const MGMatrix& mat){(*this)*=mat;};///matrix transformation.
virtual void transform(const MGTransf& tr){(*this)*=tr;};///general transformation.

protected:

///Mark this as updated.
virtual void invalidateBox()const;

///Read all member data.
virtual void ReadMembers(MGIfstream& buf);

///Write all member data
virtual void WriteMembers(MGOfstream& buf)const;

///Assignment.
///When the leaf object of this and gel2 are not equal, this assignment
///does nothing.
MGObject& set_object(const MGObject& gel2);

private:


friend class MGIfstream;
friend class MGOfstream;
friend MGCell;

};

/** @} */ // end of MGObjectRelated group
#endif
