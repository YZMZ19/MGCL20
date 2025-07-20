/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#ifndef _MGComplex_HH_
#define _MGComplex_HH_

#include <list>
#include <map>
#include "mg/Box.h"
#include "mg/Object.h"
#include "topo/Cell.h"

//
//Define MGComplex Class.

class MGPosition;
class MGCell;
class MGPCell;

/** @addtogroup TOPO
 *  @{
 */

///MGComplex is a container of parameter cells and binder cells.
class MG_DLL_DECLR MGComplex: public MGObject{
public:
typedef std::list<UniqueCell> container_type;
typedef container_type::iterator iterator;
typedef container_type::const_iterator const_iterator;

private:
	mutable MGCell* m_parent_cell;///<Cell that has this boundary as a boundary.
	container_type m_pcells;///< list of pcell elements.

public:

///Void constructor.
MGComplex();
virtual ~MGComplex();

///Copy constructor. Copy as a boundary complex of parent.
///When parent is not specified, this is ordinary world complex.
/// Not a boundary complex.
MGComplex(const MGComplex& complex);
MGComplex(MGComplex&& complex);

///Construct of one cell.
explicit MGComplex(const MGCell& cell);

///This form takes the ownership of the cell, must be newed object.
explicit MGComplex(MGCell* cell);

/////////operator overload/////////

///Assignment.
///When the leaf object of this and topo2 are not equal, this assignment
///does nothing.
virtual MGComplex& operator=(const MGGel& gel2);
virtual MGComplex& operator=(const MGComplex& comp2);
virtual MGComplex& operator=(MGComplex&& comp2)noexcept;

///Object transformation.
virtual MGComplex& operator+=(const MGVector& v);
virtual MGComplex& operator-=(const MGVector& v){ return operator+=(-v); };
virtual MGComplex& operator*=(double scale);
virtual MGComplex& operator*=(const MGMatrix& mat);
virtual MGComplex& operator*=(const MGTransf& tr);

/////////Member Function/////////

///Obtain binders of all the pcells of the boundary.
///i-th binder of the fucntion's returned value is the binder
///of i-th pcell of the boundary, and may be null.
std::vector<MGBCell*> binders() const;

///Compute barycenter of all the vertex(binder cell of 0D manifold
///dimension).
MGPosition center() const;

///Construct new object by copying to newed area.
///User must delete this copied object by "delete".
virtual MGComplex* clone() const override =0;

///Compute the box from the scratch.
void compute_box(MGBox& bx) const;

///Draw the object in wire mode, in the world coordinates.
///The object is converted to curve(s) and is drawn.
virtual void drawWire(
	mgVBO& vbo,///<The target graphic object.
	int line_density = 1	///<line density to draw a surface in wire mode.
)const override;

///Erase first pcell element, including binder cells
///that will be free when the first pcell is erased.
void erase_first_pcell(){m_pcells.pop_front();};	///erase first pcell.

///Erase last pcell element, including binder cells
///that will be free when the last pcell is erased.
void erase_last_pcell(){m_pcells.pop_back();};	///erase last pcell.

///Erase the pcell element, including binder cells
///that will be free when the last pcell is erased.
void erase(iterator pcelli);

///Get fisrt pcell pointer.
const MGCell* first_pcell()const;
MGCell* first_pcell();

//Get boundary binders of all the boundaries of member cells(m_pcells).
//Shared binders are stored once in cvec.
//cvec contais MGBVertex if this is MGLoop, and
//MGEdge(as MGBCell) if MGShell.
void get_all_boundary_binders(std::vector<MGBCell*>& cvec) const;

///Return Object's type ID (TID)
virtual long identify_type()const;

///Test if this complex includes the MGCell cell as a contituent.
///Returns true if cell is included in this complex.
bool includes(const MGCell* cell)const;

///Get last pcell pointer.
const MGCell* last_pcell()const;
MGCell* last_pcell();

///Get manifold dimension.
virtual int manifold_dimension() const=0;

///Reverse the direction of the boundary.
///(Coordinate transformation is not performed.)
virtual void negate();

//Negate the boundary according to the parent cell negation.
//That is,
//1. Transform the coordinates of the bondary cell.
//(This transfromation depends on how the parent cell is transformed
//when negate() is invoked. So, the member cells of this boundary
//are transformed by negate_transoform of the parent cell.)
//2. Reverse the direction of the parameter cells(negate each cell).
//3. Reverse the ordering of the parameter cells.
//(*****Does not negate the binders*****)
void negate_as_boundary(const MGCell* parent);

///count number of pcells of the complex
int number_of_pcells() const{return int(m_pcells.size());};

/// Output virtual function.
virtual std::ostream& toString(std::ostream&) const;

///Obtain first pcell iterator.
const_iterator pcell_begin() const{ return m_pcells.begin();};
iterator pcell_begin() { return m_pcells.begin();};

///Obtain end pcell iterator(next of the last pcell).
const_iterator pcell_end() const{ return m_pcells.end();};
iterator pcell_end() { return m_pcells.end();};

///Obtain i-the pcell. MGCell version.
const MGCell* pcelli(int i) const;
MGCell* pcelli(int i);

///Obtain i-th pcell iterator.
const_iterator pcellIterator(int i) const;
iterator pcellIterator(int i);

///Cehck if pcell exist.
bool pcell_exist()const;

///Obtain pcells that constitute the boundary.
///Let pcellvec[.] be pcells' return value and bindervec[.] be
///binders' return value. Then pcellvec[i] corresponds to bindervec[i].
///bindervec[i] is binder cell of i-th pcell element of the boundary.
container_type& pcells(){ return m_pcells; };
const container_type& pcells() const{ return m_pcells; };

///Compute the parameter value of the closest point from the straight to
///this object.
///sl is the eye projection line whose direction is from yon to hither, and if
///sl had multiple intersection points, The closest point to the eye will be selected.
MGPosition pick_closest(const MGStraight& sl)const;

///Set parent cell.
///Returned is the conventional parent cell attached to
///before execution of this set_parent.
MGCell* set_parent(MGCell& new_parent) const;

///Get star cell pointer if this complex is a boundary of a cell.
///Else, null will be returned.
const MGCell* star() const{ return m_parent_cell; };
MGCell* star(){ return m_parent_cell; };

protected:

///Fundamental constructor.
///Construct from a list of pcells.
///This constructor takes the ownership of all pcells in pcells.
///Fundamental constructor of MGComplex does not take box as input since
///this constructor has to do binder append treatment that is attached to
///pcells.
explicit MGComplex(std::list<MGCell*>& pcells);

///Insert a Cell(binder or parameter) before the position loc.
///loc is the iterator of m_pcells or m_bcells according to pcell.
iterator add_cell(
	MGCell* cell,	///<Cell to insert
	iterator loc	///<Iterator that indicates insert position. Insert cell befoer 
					///<the iterator loc.
);

///Append a PCell to the end of pcell sequence.
iterator append_pcell(MGCell* cell);

//Copy m_pcells of comp into this. Binders of m_pcells are not copied, but deleted.
//Each MGPCell in m_pcells association of the original(1st and key) and
//newed one are registered in cmap when m_pcells is a vector of MGPCell's subclass.
//cmap contains MGEdge associations if MGLoop.
//When this is MGShell, cmap contains nothing.
void copy_all_elements(const MGComplex& comp, std::map<const MGPCell*, MGPCell*>* cmap = 0);

///Copy boundary.
///This boundary data is cleared and bnd's boundary is copied into this.
virtual void copy_boundary(const MGComplex& bnd);

///Prepend a PCell to the end of pcell sequence.
iterator prepend_pcell(MGCell* cell);

///Get the name of the class.
virtual std::string whoami()const{return "Complex";};

///Read Object's member data.
virtual void ReadMembers(MGIfstream& buf);

///Write Object's Member Data
virtual void WriteMembers(MGOfstream& buf) const;

friend class MGCell;
friend class MGPCell;
friend class MGCell;
friend class MGEdge;
friend class MGFace;
};

/** @} */ // end of TOPO group
#endif
