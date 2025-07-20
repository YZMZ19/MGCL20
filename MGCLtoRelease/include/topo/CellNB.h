/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#ifndef _MGCellNB_HH_
#define _MGCellNB_HH_

#include <vector>
#include <map>
#include "mg/Box.h"
#include "mg/Position.h"
#include "mg/Geometry.h"
#include "topo/Topology.h"

//
//Defines MGCellNB Class.

class MGBox;
class MGComplex;
class MGGeometry;
class MGPCell;
class MGBVertex;

/** @addtogroup TOPO
 *  @{
 */

///CellNB is a cell without boundaries(No Boundaries).

///MGCellNB ia an abstract class and the super class of MGEdge and MGFace.
///There are two types of cells. One is parameter cell(pcell) and
///the other is binder cell(bcell). They are exclusive, that is, if
///a cell is a parameter cell, the cell cannot be binder cell and
///vice versa.
///MGCellNB cell is a constituent of a complex and MGCellNB's are stored in MGComplex.
///Binder cell is a binder of parameter cells. Plural cells are connected through a binder.
class MG_DLL_DECLR MGCellNB:public MGTopology{

protected:
	/// World of this cell.
	/// When m_parent_complex is null, or m_parent_complex's parent boundary is null
	/// this is ordinary world cell, not parameter world.
	mutable MGComplex* m_parent_complex;

	///Geometry of this cell, which is a curve or surface.
	std::unique_ptr<MGGeometry> m_extent;	///<Geometry

	mutable double m_perror;///<Error allowed for the parameter space of the cell.

	mutable MGBox m_box;	///<Box of this cell.
		///<Initially this is null,
		///<and will be computed by compute_box() when necessary.

public:

///////Constructor/////////

///Special constructors.
MGCellNB();
MGCellNB(const MGCellNB& rhs);
MGCellNB(MGCellNB&& rhs);//Move constructor.
virtual ~MGCellNB() = default;
MGCellNB& operator=(const MGCellNB& rhs);
MGCellNB& operator=(MGCellNB&& rhs);//Move assignment.

///CellNB of whole geometry(no boundary).
///The second form that input MGGeometry* takes the ownership of the geo
///into the MGCellNB, must not delete the object and the object must be
///newed one.
MGCellNB(const MGGeometry& geo);
explicit MGCellNB(
	MGGeometry* geo///<Geometry of this MGCellNB. Must be a newed object and the
		///<ownership is transfered to this.
);

///////operator overload//////

///Object transformation.
virtual MGCellNB& operator+=(const MGVector& v);
virtual MGCellNB& operator-=(const MGVector& v){ return operator+=(-v); };
virtual MGCellNB& operator*=(double scale);
virtual MGCellNB& operator*=(const MGMatrix& mat);
virtual MGCellNB& operator*=(const MGTransf& tr);

/////////////Member Function///////////////

///Return space dimension
int sdim() const;

///Obtain the box of the cell.
const MGBox& box() const;

virtual const MGBCell* binderCell()const=0;
virtual MGBCell* binderCell()=0;

///Obtain the center of this cell.
MGPosition center() const;

///Obtain the center parameter value of this cell.
virtual MGPosition center_param() const=0;

///Make a clone of the cell.
///clone() does not copy the binder cell of this.
virtual MGCellNB* clone() const=0;

///Obtain the direction of the cell.
virtual MGUnit_vector direction() const;

///Draw 3D point(vertex) in world coordinates.
///The object is converted to point(s) and is drawn.
///This is valid only for topology objects or MGPoint.
virtual void drawVertex(mgVBO& vbo)const=0;

///Get extent geometry, may be null if this does not have extent.
const std::unique_ptr<MGGeometry>& extent() const{ return m_extent;};
std::unique_ptr<MGGeometry>& extent() {return m_extent;};

///Free(but does not delete) the extent geometry.
///Freed extent is returned as the function's return value.
MGGeometry* free_extent();

///Free from membership of the parent complex.
///free_from_parent() does not maintain the box of the complex this cell
///belonged to. And so, users of free_from_parent() must do it.
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
virtual void set_extent_as_null(){ m_extent.reset(); };

///Obtain star cells.
const MGCellNB* star() const;
MGCellNB* star();

/// Output virtual function.
virtual std::ostream& toString(std::ostream&) const;

///Get the name of the class.
virtual std::string whoami()const{return "CellNB";};

protected:

///check if boundary's binder transformation is necessary or not.
bool bn_binder_tr_necessary()const{	return !parent_complex();}

///compute box of the cell in m_box.
virtual void compute_box() const;

///Copy all boundaries of cellin into this(but does not copy own binder cell relation),
//and register cellin's boundary MGPCell* association into cmap.
//1st(key) is original MGPCell* of boundaries of m_pcells,
//2nd is copyied new.
virtual void copy_all_boundaries(
	const MGCellNB& cellin, std::map<const MGPCell*, MGPCell*>* cmap = 0
)=0;

///Free specified boundary(bound) from a member of parent cell's boundaries.
///Return MGBoundary if freed normally.
///If bound was not a member of the boundaries, return 0.
///Only free, does not destruct the boundary.
virtual MGComplex* free_boundary(const MGComplex* bound)=0;

///Get all the MGPCell* of the all the boundaries of this.
virtual std::vector<const MGPCell*> getBoundaryPcells()const=0;

///Cell comparison.
bool is_less_than(const MGCellNB& cell2)const;

///set box as null(to set the box as initial)
virtual void set_box_as_null()const{ m_box.set_null(); };

///Read Object's member data.
virtual void ReadMembers(MGIfstream& buf);

///set member datas.
///virtual void set_members(
///	MGGeometry* geo);

///Write Object's Member Data
virtual void WriteMembers(MGOfstream& buf) const;

private:

///Generate a copied MGCellNB of this by newing.
///This is a proprietry routine of MGComplex copy.
///Copy is performed by registering all the original boundary data's MGPCell
///and new MGPCell association in cmap, (but does not copy own binder cell relation).
UniqueCellNB cloneWithMap(std::map<const MGPCell*, MGPCell*>& cmap)const;

///Make sure that this has an extent expression.
///When this did not have an extent, make the extent from the partner
///member's parameter expression and the star cell.
///This must be a binder cell that has partner members that are
///boundaries. When this is not the case or this had an extent already,
///it does nothing.
virtual void make_extent() const=0;

///Negate the boundary.
virtual void negate_boundary()=0;

friend class MGPCell;
friend class MGComplex;
friend class MGEdge;

};

/** @} */ // end of TOPO group
#endif
