/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#include "StdAfx.h"
#include "mg/Geometry.h"
#include "mg/Point.h"
#include "mg/Position.h"
#include "mg/Curve.h"
#include "mg/Straight.h"
#include "topo/BCell.h"
#include "topo/PVertex.h"
#include "topo/Complex.h"
#include "topo/Cell.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//
//Implement MGComplex Class.

//Constructor

//Void constructor (under a boundary).
MGComplex::MGComplex():m_parent_cell(nullptr){;}

//Construct of one cell.
//The second form takes the ownership of the cell, must be newed object.
MGComplex::MGComplex(const MGCell& cell):m_parent_cell(nullptr){
	add_cell(cell.clone(), m_pcells.begin());
}
MGComplex::MGComplex(MGCell* cell):m_parent_cell(nullptr){
	add_cell(cell, m_pcells.begin());
}

//Fundamental constructor.
//Construct from a list of pcells.
//This constructor takes the ownership of all pcells in pcells.
MGComplex::MGComplex(std::list<MGCell*>& pcells)
	:m_parent_cell(nullptr){
	for(auto i:pcells){
		add_cell(i, m_pcells.end());
	}
}

//Copy constructor. Binders of member cells(m_pcells) are deleted.
MGComplex::MGComplex(
	const MGComplex& complex)	//Original complex.
:MGObject(complex), m_parent_cell(nullptr){
	copy_all_elements(complex);
}
MGComplex::MGComplex(
	MGComplex&& complex	//Original complex.
):MGObject(std::move(complex)), m_parent_cell(nullptr)
, m_pcells(std::move(complex.m_pcells)){
}
MGComplex::~MGComplex() {
	m_parent_cell = nullptr;
	for (UniqueCell& ci : m_pcells)
		ci.reset();
}

///////operator overload///////

// Complexに平行移動を行ない自身のComplexとする。
//Translation of the Complex
MGComplex& MGComplex::operator+= (const MGVector& v){
	iterator cs, ce, csave;
	ce=m_pcells.end();
	for(cs=m_pcells.begin(); cs!=ce; cs++)
		**cs +=v;

	std::vector<MGBCell*> bvec;
	get_all_boundary_binders(bvec);
	for(auto& bvi:bvec)
		bvi->binder_tr(v);
	if(!box_is_null())
		m_box +=v;
	return *this;
}

//Complexのスケーリングを行い自身のComplexとする。
//Scaling of the Complex by a double.
MGComplex& MGComplex::operator*= (double s){
	iterator cs, ce, csave;
	ce=m_pcells.end();
	for(cs=m_pcells.begin(); cs!=ce; cs++)
		**cs *=s;
	invalidateBox();
	return *this;
}

// 与えられた変換でComplexの変換を行い自身のComplexとする。
//Transformation of the Complex by a matrix.
MGComplex& MGComplex::operator*= (const MGMatrix& mat){
	iterator cs, ce, csave;
	ce=m_pcells.end();
	for(cs=m_pcells.begin(); cs!=ce; cs++) **cs *=mat;

	invalidateBox();
	return *this;
}

// 与えられた変換によってトランスフォームをおこない自身のComplexにする。
//Transformation of the Complex by a MGTransf.
MGComplex& MGComplex::operator*= (const MGTransf& tr){
	iterator cs, ce, csave;
	ce=m_pcells.end();
	for(cs=m_pcells.begin(); cs!=ce; cs++) **cs *=tr;

	invalidateBox();
	return *this;
}

MGComplex& MGComplex::operator=(const MGComplex& comp2){
	if(m_parent_cell){
		m_parent_cell->free_boundary(this);
		m_parent_cell = 0;
	}
	m_pcells.clear();
	MGObject::operator =(comp2);
	copy_all_elements(comp2);
	return *this;
}
MGComplex& MGComplex::operator=(MGComplex&& comp2)noexcept {
	if(m_parent_cell){
		m_parent_cell->free_boundary(this);
		m_parent_cell = 0;
	}
	MGObject::operator =(std::move(comp2));
	m_pcells=std::move(comp2.m_pcells);
	return *this;
}
MGComplex& MGComplex::operator=(const MGGel& gel2){
	const MGComplex* gel2_is_this=dynamic_cast<const MGComplex*>(&gel2);
	if(gel2_is_this)
		operator=(*gel2_is_this);
	return *this;
}

//Set parent cell.
//Returned is the conventional parent cell attached to
//before execution of this set_parent.
MGCell* MGComplex::set_parent(MGCell& new_parent)const{
	MGCell* parentOld = m_parent_cell;
	if(parentOld){
		if(parentOld==&new_parent)
			return parentOld;
		new_parent.free_boundary(this);
	}
	m_parent_cell = &new_parent;
	return parentOld;
}

//Cell comparison.
bool MGComplex::operator<(const MGComplex& cell2)const{
	if(this==&cell2)
		return false;

	const MGCell* c1=star();
	const MGCell* c2=cell2.star();
	return number_of_pcells()<cell2.number_of_pcells();
}
bool MGComplex::operator<(const MGGel& gel2)const{
	const MGComplex* comlx2=dynamic_cast<const MGComplex*>(&gel2);
	if(comlx2)
		return operator<(*comlx2);
	return identify_type() < gel2.identify_type();
}

//Member Function

//Insert a Cell(binder or parameter) before the position loc.
//loc is the iterator of m_pcells .
MGComplex::iterator MGComplex::add_cell(
	MGCell* cell,	//Cell to insert
	MGComplex::iterator loc
					//Iterator that indicates insert position. Insert cell befoer 
					//the iterator loc.
){
	MGComplex::iterator newloc=loc;
	MGComplex* compold=cell->parent_complex();
	if(compold!=this){
		//1. Free from the old parent(if exist).
		if(compold){
			cell->free_from_parent();
			compold->invalidateBox();			
		}

		//2. Add own cell.
		newloc= m_pcells.insert(loc, UniqueCell(cell));
		cell->m_parent_complex=this;
		invalidateBox();
	}
	return newloc;
}

//Append a PCell to the end of pcell sequence.
MGComplex::iterator MGComplex::append_pcell(MGCell* cell){
	assert(!cell || (cell&&cell->extent()));//Extent must be attached.
	return add_cell(cell, m_pcells.end());
}

//Obtain binders of all the pcells of the boundary.
//i-th binder of the fucntion's returned value is the binder
//of i-th pcell of the boundary, and may be null.
std::vector<MGBCell*> MGComplex::binders() const{
	int n=number_of_pcells();
	std::vector<MGBCell*> bdrs(n);
	const_iterator cell=pcell_begin();
	for(int i=0; i<n; i++){
		const UniqueCell& pcell=*cell++;
		if(pcell)
			bdrs[i]=pcell->binderCell();
		else
			bdrs[i]=0;
	}
	return bdrs;
}

//Compute barycenter of all the vertex(binder cell of 0D manifold
//dimension).
MGPosition MGComplex::center() const{
	MGPosition cntr;
	int num=0;
	for(auto& pcell:pcells()){
		cntr += pcell->center();
		num++;
	}
	if(num) cntr/=double(num);
	return cntr;
}

//Get boundary biders of all the boundaries of member cells(m_pcells).
//Shared binder pointers are stored once in cvec.
//cvec contais MGBVertex if this is MGLoop, and
//MGEdge(as MGBCell) if MGShell.
void MGComplex::get_all_boundary_binders(std::vector<MGBCell*>& cvec) const{
	for(auto& pcelli : m_pcells){	//Traverse over all the MGCell in m_pcells.
		std::vector<const MGPCell*> bndryPcells=pcelli->getBoundaryPcells();
		for(auto& bndryPcelli:bndryPcells){
			SharedBCell& bcelli = bndryPcelli->binder();
			if(!bcelli)
				continue;
			if(bndryPcelli==bcelli->first_partner_member())
				cvec.push_back(bcelli.get());
		}
	}
}

//Copy m_pcells of comp into this. Binders of m_pcells are also copied to maintain the binderness.
//Each MGPCell in m_pcells association of the original(1st and key) and
//newed one are registered in cmap when m_pcells is a vector of MGPCell's subclass.
//cmap contains MGEdge associations if MGLoop.
//When this is MGShell cmap contains nothing.
void MGComplex::copy_all_elements(
	const MGComplex& comp,
	std::map<const MGPCell*, MGPCell*>* cmap
){
	//1. Clone all the MGCell's(without binders) of m_pcells.

    // 1st(key) is original MGPCell* of boundaries of m_pcells, 2nd is copyied new.
	std::map<const MGPCell*, MGPCell*> cmapBndry;
	for(auto& pcelli : comp.m_pcells){
		//pcelli's boundary's MGPCell* into cmapBndry.
		UniqueCell pcellNew(pcelli->cloneWithMap(cmapBndry));
		if(cmap){
			//If cmap given, register m_pcells association of original and new in cmap.
			MGPCell* pcelOld = dynamic_cast<MGPCell*>(pcelli.get());
			if(pcelOld)
				cmap->insert(std::make_pair(pcelOld, dynamic_cast<MGPCell*>(pcellNew.get())));
		}
		append_pcell(pcellNew.release());
	}

	//2. Clone the binders of the boundaries of m_pcells.
	std::vector<MGBCell*> bcells;
	comp.get_all_boundary_binders(bcells);
	//bcells are binders of the boundaries of comp.m_pcells,
	//not binders of the member pcells.
	//E.g. if this is MGLoop, m_pcells are edges, and the boundarise are MGPVertex.
	//So bcells contains MGBVetex.
	//If this is MGShell, m_pcells are faces, and the boundaries are MGLoop.
	//So bcells contains MGEdge.

	for(auto& bcelli : bcells){
		const MGBCell& bcell = *bcelli;
		size_t n=bcell.number_of_partner_members();
		if (!n)
			continue;
		std::vector<const MGPCell*> newPcells(n);
		size_t i = 0;
		for(auto& pcelli:bcell.partner_members()){
			newPcells[i++]= cmapBndry[pcelli];
		}
		SharedBCell newBcelli = bcell.cloneWithPartnerMembers(newPcells);
		for(auto& newPcelli : newPcells){
			newPcelli->m_binder = newBcelli;
		}
	}
	invalidateBox();
}

//Copy boundary.
//This boundary data is cleared and bnd's boundary is copied into this.
void MGComplex::copy_boundary(const MGComplex& bnd){
	MGComplex::operator=(bnd);
}

///Erase the pcell element, including binder cells
///that will be free when the last pcell is erased.
void MGComplex::erase(iterator pcelli){ m_pcells.erase(pcelli); };

//Get fisrt pcell pointer.
const MGCell* MGComplex::first_pcell() const{
	return m_pcells.front().get();
}
MGCell* MGComplex::first_pcell(){
	return m_pcells.front().get();
}

//Test if this complex includes the MGCell cell as a contituent.
//Returns true if cell is included in this complex.
bool MGComplex::includes(const MGCell* cell)const{
	const_iterator ps=m_pcells.begin(), pe=m_pcells.end();
	for(; ps!=pe; ps++){	
		if((*ps).get()==cell)
			return true;
	}
	return false;
}

//Get last pcell pointer.
const MGCell* MGComplex::last_pcell() const{
	return m_pcells.back().get();
}
MGCell* MGComplex::last_pcell(){
	return m_pcells.back().get();
}

//Reverse the direction of the boundary.
//(Coordinate transformation is not performed.)
void MGComplex::negate(){
	iterator i, cs, ce;
	i = cs = pcell_begin(); ce = pcell_end();
	//Negate each parameter cell.
	for(; i!=ce; i++)
		(*i)->negate();

	//Reverse the ordering of the parameter cells.
	std::reverse(cs, ce);
}

//Negate the boundary according to the parent cell negation.
//That is,
//1. Transform the coordinates of the bondary cell.
//(This transfromation depends on how the parent cell is transformed
//when negate() is invoked. So, the member cells of this boundary
//are transformed by negate_transoform of the parent cell.)
//2. Reverse the direction of the parameter cells(negate each cell).
//3. Reverse the ordering of the parameter cells.
//(*****Does not negate the binders*****)
void MGComplex::negate_as_boundary(const MGCell* parent){
	const MGCell* parent_cell = parent;
	if(!parent){
		assert(star()); //This must be a boudary of a cell when parent not specified.
		parent_cell = star();
	}
	const UniqueGeometry& parent_geo = parent_cell->extent();

	iterator cs, ce;
	cs = pcell_begin(); ce = pcell_end();
	for(; cs!=ce; cs++){
		//1. Transform the coordinates of each cell.
		parent_geo->negate_transform(*((*cs)->extent()));
		//2. Negate each parameter cell.
		(*cs)->negate();
	}
	//3. Reverse the ordering of the parameter cells.
	std::reverse(pcell_begin(), ce);
	invalidateBox();
}

// Output virtual function.
std::ostream& MGComplex::toString(std::ostream& ostrm) const{
	ostrm << "parent=";
	if (m_parent_cell) {
		ostrm << (const MGGel*)m_parent_cell;
	}else {
		ostrm << "Null";
	}
	MGObject::toString(ostrm);

	size_t n = m_pcells.size();
	std::string pcellName = n ? m_pcells.front()->whoami():"PCell";
	ostrm<<", "<<std::endl<<"Number of "<<pcellName<<"s = "<<n<<"::";
	int i=0;
	ostrm<<std::endl;
	for(auto& pcelli:m_pcells)
		ostrm<<(*pcelli).whoami()<<i++<<"="<<*pcelli<<std::endl;

	return ostrm;
}

//Cehck if pcell exist.
bool MGComplex::pcell_exist()const{
	return m_pcells.begin()!=m_pcells.end();
}

//Obtain i-the pcell. MGCell version.
const MGCell* MGComplex::pcelli(int i) const{
	const_iterator pcellp=pcell_begin();
	std::advance(pcellp,i);
	return pcellp->get();
}
MGCell* MGComplex::pcelli(int i){
	iterator pcellp=pcell_begin();
	std::advance(pcellp,i);
	return pcellp->get();
}

//Obtain i-th pcell iterator.
MGComplex::const_iterator MGComplex::pcellIterator(int i) const{
	const_iterator pcellp=pcell_begin();
	std::advance(pcellp,i);
	return pcellp;
}
MGComplex::iterator MGComplex::pcellIterator(int i){
	iterator pcellp=pcell_begin();
	std::advance(pcellp,i);
	return pcellp;
}

//Compute the parameter value of the closest point from the straight to
//this object.
//sl is the eye projection line whose direction is from yon to hither, and if
//sl had multiple intersection points, The closest point to the eye will be selected.
MGPosition MGComplex::pick_closest(const MGStraight& sl)const{
	const_iterator ps=m_pcells.begin(), pe=m_pcells.end();
	MGPosition prm;
	if(ps==pe) return prm;
	prm=(**ps).pick_closest(sl);
	double t=sl.closest((**ps).extent()->evaluate(prm));
	for(ps++; ps!=pe; ps++){
		MGPosition prm2=(**ps).pick_closest(sl);
		double t2=sl.closest((**ps).extent()->evaluate(prm2));		
		if(t2>t){
			t=t2; prm=prm2;
		}
	}
	return prm;
}

//Prepend a PCell to the end of pcell sequence.
MGComplex::iterator MGComplex::prepend_pcell(MGCell* cell){
	assert(!cell || (cell&&cell->extent()));//Extent must be attached.
	return add_cell(cell, m_pcells.begin());
}

//Compute the box from the scratch.
void MGComplex::compute_box(MGBox& bx) const{
	bx.set_null();
	const_iterator ps=m_pcells.begin(), pe=m_pcells.end();
	for(;ps!=pe;ps++)
		bx |=(*ps)->box();
}
