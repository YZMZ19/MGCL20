/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#include "StdAfx.h"
#include "mg/Tolerance.h"
#include "mg/Unit_vector.h"
#include "mg/Geometry.h"
#include "mg/Point.h"
#include "topo/PCell.h"
#include "topo/Cell.h"
#include "topo/Complex.h"
#include "mgGL/Appearance.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//
//Implements MGCell Class.
//Cell is a cell without boundaries(No Boundaries).

//There are two types of cells. One is parameter cell(pcell) and
//the other is binder cell(bcell). They are exclusive, that is, if
//a cell is a parameter cell, the cell cannot be binder cell and
//vice versa.
//Parameter cell is a constituent of a complex.
//Binder cell is a binder of parameter cells. Plural cells are connected
//through a binder.
//MGCell ia an abstract class.

//////Constructor///////

//Void constructor. Constructor of pcell.
MGCell::MGCell():m_parent_complex(nullptr), m_perror(-1.){;}

//Copy constructor. Result cell is not a member of any complex.
//Partners of cell will not be copied.
MGCell::MGCell(const MGCell& cell)
:MGObject(cell), m_parent_complex(nullptr), m_perror(cell.m_perror){
	if(cell.m_extent)
		m_extent.reset(cell.m_extent->clone());
}

//Move constructor. Result cell is not a member of any complex.
//Partners of cell will not be copied.
MGCell::MGCell(MGCell&& cell)
:MGObject(std::move(cell)), m_parent_complex(nullptr), m_perror(cell.m_perror),
 m_extent(std::move(cell.m_extent)){
}

//Cell of whole geometry(no boundary), under parent.
//The second form that input MGGeometry* takes the ownership of the geo
//into the MGCell, must not delete the object and the object must be
//newed one.
MGCell::MGCell(const MGGeometry& geo)
:m_extent(geo.clone()),m_parent_complex(nullptr), m_perror(-1.){
	copy_appearance(geo);
	m_box = geo.box();
}
MGCell::MGCell(MGGeometry* geo)
:m_extent(geo),m_parent_complex(nullptr), m_perror(-1.){
	if(geo){
		copy_appearance(*geo);
		geo->remove_appearance();
		m_box = geo->box();
	}
}

MGCell::~MGCell() {
	free_from_parent();
}

//////operator overload/////
//does not change parent complex.
MGCell& MGCell::operator=(const MGCell& cell2){
	MGObject::operator=(cell2);
	if(cell2.m_extent)
		m_extent.reset(cell2.m_extent->clone());
	else
		m_extent.reset();
	m_box = cell2.m_box;
	m_perror = cell2.m_perror;
	return *this;
}

///Generate a copied MGCell of this by newing.
///This is a proprietry routine of MGComplex copy.
///Copy is performed by registering all the original boundary data's MGPCell
///and new MGPCell association in cmap, (but does not copy own binder cell relation).
///My parent is not copied to cloned cell.
UniqueCell MGCell::cloneWithMap(
	//1st(key) is original and 2nd is copyied new.
	std::map<const MGPCell*, MGPCell*>& cmap
)const{
	UniqueCell cell(cloneWithoutBoundary());
	cell->copy_all_boundaries(*this, &cmap);
	return cell;
}

// Cellに平行移動を行ない自身のCellとする。
//Translation of the cell
MGCell& MGCell::operator+= (const MGVector& v){
	if(m_extent)
		*m_extent +=v;
	m_box+=v;
	return *this;
}

//Cellのスケーリングを行い自身のCellとする。
//Scaling of the cell by a double.
MGCell& MGCell::operator*= (double s){
	if(m_extent)
		*m_extent *=s;
	invalidateBox();
	return *this;
}

// 与えられた変換でCellの変換を行い自身のCellとする。
//Transformation of the cell by a matrix.
MGCell& MGCell::operator*= (const MGMatrix& mat){
	if(m_extent)
		*m_extent *=mat;
	invalidateBox();
	return *this;
}

// 与えられた変換によってトランスフォームをおこない自身のCellにする。
//Transformation of the cell by a MGTransf.
MGCell& MGCell::operator*= (const MGTransf& tr){
	if(m_extent)
		*m_extent *=tr;
	invalidateBox();
	return *this;
}

///Return space dimension
int MGCell::sdim() const{
	return m_extent ? m_extent->sdim():0; 
}

//Cell comparison.
bool MGCell::compare(const MGCell& cell2)const{
	long thisID=identify_type(),  ID2=cell2.identify_type();
	if (thisID != ID2)//If not the same type cells.
		return thisID < ID2;

	if(this==&cell2)
		return false;

	const MGComplex* comp1=parent_complex();
	if(!comp1)
		return true;
	const MGComplex* comp2=cell2.parent_complex();
	if(!comp2)
		return false;

	if(comp1!=comp2)
		return (*comp1)<(*comp2);

	//Now this and cell2 have the same parent complex.
	//Comparison is done by the appearance order of this or cell2 in the complex.
	MGComplex::const_iterator i=comp1->pcell_begin(), ie=comp1->pcell_end();
	for(auto& pcelli:comp1->pcells()){
		const MGCell* celli = pcelli.get();
		if(celli==this)
			return true;
		if(celli==&cell2)
			return false;
	}
	return true;//This statement will never be executed.
}

///////////Member Function/////////////

//Obtain the center of this cell.
MGPosition MGCell::center() const{
	MGPosition cntr;
	if(m_extent)
		cntr=m_extent->evaluate(center_param());
	return cntr;
}

//Obtain the direction of the cell.
MGUnit_vector MGCell::direction() const{
	MGPosition param=center_param();
	return extent()->direction(param);
}

//Free(but does not delete) the extent geometry.
//Freed extent is returned as the function's return value.
MGGeometry* MGCell::free_extent(){
	return m_extent.release();
}

///Free from membership of the parent complex.
///free_from_parent() does not maintain the box of the complex this cell
///belonged to. And so, users of free_from_parent() must do it.
///free_from_parent() frees this cell's boudary bindness, that is,
///if this had neibours(otherwords, if this boundaries had binders)
///they are freed.
MGComplex* MGCell::free_from_parent(){
	free_neighbours();
	MGComplex* parent=m_parent_complex;
	if(!parent) return parent;

	//free this cell from parent complex.
	MGComplex::iterator itri, itrend;
	itri=parent->m_pcells.begin(); itrend=parent->m_pcells.end();
	for(; itri!=itrend; itri++){
		if(this==(*itri).get()) {
			itri->release();
			parent->m_pcells.erase(itri);
			m_parent_complex=nullptr;
			return parent;
		}
	}	
	return parent;
}

//Negate the direction of the cell.
void MGCell::negate(){
	if(m_extent){
		//1. Negate each boudary.
		negate_boundary();

		//2. Negate own extent.
		m_extent->negate();
		invalidateBox();
	}
}

//Return parameter space error of the cell.
double MGCell::parameter_error()const{
	if(m_perror<=0.){
		if(m_extent)
			m_perror = m_extent->parameter_error();
		else 
			return MGTolerance::wc_zero();
	}
	return m_perror;
}

//Set extent of this cell.
void MGCell::set_extent(std::unique_ptr<MGGeometry>&& extent){
	m_extent = std::move(extent); invalidateBox();
	MGComplex* parent = parent_complex();
	if(parent)
		parent->m_box.set_null();
}
void MGCell::set_extent_as_null(){
	m_extent.reset(); 
	invalidateBox();
}

// Output virtual function.
std::ostream& MGCell::toString(std::ostream& ostrm) const{
	MGObject::toString(ostrm);
	if(m_parent_complex){
		std::string parentName = m_parent_complex->whoami();
		ostrm<<", parent"<<parentName<<"="<<(const MGGel*)m_parent_complex<<std::endl;
	}
	ostrm<<", perror="<<m_perror;
	if(m_extent){
		std::string myName = whoami();
		ostrm<<","<<std::endl<<myName<<"Extent="<<*m_extent;
	}
	return ostrm;
}
