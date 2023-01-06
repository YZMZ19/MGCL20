/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#ifndef _MGCell_HH_
#define _MGCell_HH_

#include <vector>
#include <map>
#include "mg/Box.h"
#include "mg/Position.h"
#include "mg/Geometry.h"
#include "mgGL/VBO.h"

//
//Defines MGCell Class.

class MGBox;
class MGComplex;
class MGGeometry;
class MGPCell;

/** @addtogroup TOPO
 *  @{
 */

///Cell is a cell without boundaries(No Boundaries).

///MGCell ia an abstract class and the super class of MGEdge and MGFace.
///There are two types of cells. One is parameter cell(pcell) and
///the other is binder cell(bcell). They are exclusive, that is, if
///a cell is a parameter cell, the cell cannot be binder cell and
///vice versa.
///MGCell cell is a constituent of a complex and MGCell's are stored in MGComplex.
///Binder cell is a binder of parameter cells. Plural cells are connected through a binder.
class MG_DLL_DECLR MGCell:public MGObject{

protected:
	/// World of this cell.
	/// When m_parent_complex is null, or m_parent_complex's parent boundary is null
	/// this is ordinary world cell, not parameter world.
	mutable MGComplex* m_parent_complex;

	///Geometry of this cell, which is a curve or surface.
	std::unique_ptr<MGGeometry> m_extent;	///<Geometry

	mutable double m_perror;///<Error allowed for the parameter space of the cell.

public:

///////Constructor/////////

///Special constructors.
MGCell();
virtual ~MGCell();

MGCell(const MGCell& rhs);
MGCell(MGCell&& rhs);//Move constructor.

MGCell& operator=(const MGCell& rhs);
MGCell& operator=(MGCell&& rhs)=default;//Move assignment.

///Cell of whole geometry(no boundary).
///The second form that input MGGeometry* takes the ownership of the geo
///into the MGCell, must not delete the object and the object must be
///newed one.
MGCell(const MGGeometry& geo);
explicit MGCell(
	MGGeometry* geo///<Geometry of this MGCell. Must be a newed object and the
		///<ownership is transfered to this.
);

///////operator overload//////

///Object transformation.
virtual MGCell& operator+=(const MGVector& v);
virtual MGCell& operator-=(const MGVector& v){ return operator+=(-v); };
virtual MGCell& operator*=(double scale);
virtual MGCell& operator*=(const MGMatrix& mat);
virtual MGCell& operator*=(const MGTransf& tr);

/////////////Member Function///////////////

///Return space dimension
int sdim() const;

virtual const MGBCell* binderCell()const=0;
virtual MGBCell* binderCell()=0;

///Obtain the center of this cell.
MGPosition center() const;

///Obtain the center parameter value of this cell.
virtual MGPosition center_param() const=0;

///Make a clone of the cell.
///clone() does not copy the binder cell of this.
virtual MGCell* clone() const override =0;

///Make a clone of the cell that does not have boundaries.
///Does not copy the binder cell of this.
virtual MGCell* cloneWithoutBoundary() const = 0;

///Obtain the direction of the cell.
virtual MGUnit_vector direction() const;

///Draw 3D point(vertex) in world coordinates.
///The object is converted to point(s) and is drawn.
///This is valid only for topology objects or MGPoint.
virtual void drawVertex(mgVBO& vbo)const=0;

///Get extent geometry, may be null if this does not have extent.
const std::unique_ptr<MGGeometry>& extent() const{ return m_extent;};
std::unique_ptr<MGGeometry>& extent() {return m_extent;};

///If this had boundary binders, free them. As the result this
///will have no neighbours.
virtual void free_neighbours() {;};

///Free(but does not delete) the extent geometry.
///Freed extent is returned as the function's return value.
MGGeometry* free_extent();

///Free from membership of the parent complex.
///free_from_parent() does not maintain the box of the complex this cell
///belonged to. And so, users of free_from_parent() must do it.
///free_from_parent() frees this cell's boudary bindness, that is,
///if this had neibours(otherwords, if this boundaries had binders)
///they are freed.
MGComplex* free_from_parent();

///Return Object's type ID (TID)
virtual long identify_type()const=0;

///Return true if this is a binder cell
virtual bool is_bcell()const{ return false; };

///Obtain manifold dimension.
virtual int manifold_dimension() const=0;

///Negate the direction of the cell.
virtual void negate();

///Return parameter space error of the cell.
virtual double parameter_error()const;

///Obtain parent complex.
const MGComplex* parent_complex() const{return m_parent_complex;};
MGComplex* parent_complex() {return m_parent_complex;};

///Set extent of this cell.
virtual void set_extent(std::unique_ptr<MGGeometry>&& extent);
virtual void set_extent_as_null();

/// Output virtual function.
virtual std::ostream& toString(std::ostream&) const;

protected:

///check if boundary's binder transformation is necessary or not.
bool bn_binder_tr_necessary()const{	return !parent_complex();}

///Copy all boundaries of cellin into this(but does not copy own binder cell relation),
//and register cellin's boundary MGPCell* association into cmap.
//1st(key) is original MGPCell* of boundaries of this,
//2nd is copyied new.
virtual void copy_all_boundaries(
	const MGCell& cellin,
	std::map<const MGPCell*, MGPCell*>* cmap = 0
)=0;

///Free specified boundary(bound) from a member of parent cell's boundaries.
///Return MGBoundary if freed normally.
///If bound was not a member of the boundaries, return 0.
///Only free, does not destruct the boundary.
virtual MGComplex* free_boundary(const MGComplex* bound)=0;

///Get all the MGPCell* of the all the boundaries of this.
virtual std::vector<const MGPCell*> getBoundaryPcells()const=0;

///Cell comparison.
bool compare(const MGCell& cell2)const;

///Read Object's member data.
virtual void ReadMembers(MGIfstream& buf);

///set member datas.
///virtual void set_members(
///	MGGeometry* geo);

///Write Object's Member Data
virtual void WriteMembers(MGOfstream& buf) const;

private:

///Generate a copied MGCell of this by newing.
///This is a proprietry routine of MGComplex copy.
///Copy is performed by registering all the original boundary data's MGPCell
///and new MGPCell association in cmap, (but does not copy own binder cell relation).
///My parent is not copied to cloned cell.
UniqueCell cloneWithMap(std::map<const MGPCell*, MGPCell*>& cmap)const;

///Negate the boundary.
virtual void negate_boundary()=0;

friend class MGPCell;
friend class MGComplex;
friend class MGEdge;

};

/** @} */ // end of TOPO group
#endif
