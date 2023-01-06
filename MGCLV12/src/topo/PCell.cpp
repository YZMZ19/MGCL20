/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#include "StdAfx.h"
#include "mg/Box.h"
#include "mg/Gel.h"
#include "topo/PCell.h"
#include "topo/BCell.h"
#include "topo/Cell.h"
#include "topo/Complex.h"
#include "topo/Edge.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//
//Define MGPCell Class.
//A MGPCell includes only a binder cell pointer if exists.

//There are two types of cells. One is parameter cell(pcell) and
//the other is binder cell(bcell). They are exclusive, that is, if
//a cell is a parameter cell, the cell cannot be binder cell and
//vice versa.
//Binder cell is a binder of parameter cells. Plural cells are connected
//through a binder.
//MGPCell ia an abstract class.


//////////Virtual Destructor//////////

MGPCell::~MGPCell(){
	resetBinder();
}

///When the leaf objects of this and gel2 are not equal, this assignment
///does nothing.
MGGel& MGPCell::operator=(const MGGel& gel2){
	MGPCell* pcel = dynamic_cast<MGPCell*>(this);
	if(pcel)
		(*this) = *pcel;
	return *this;
}

MGPCell& MGPCell::operator=(const MGPCell& rhs){
	if(this==&rhs)
		return *this;

	resetBinder();
	return *this;
}
MGPCell& MGPCell::operator=(MGPCell&& rhs){//Move assignment.
	if(this==&rhs)
		return *this;

	resetBinder();
	return *this;
}

///Bind this cell to cell2.
///This cell is a pcell of a boundary of a higher manifold dimension's cell A,
///and cell2 is also is a pcell of a boundary of another cell B.
///That is, this cell is a part of a boundary of cell A,
///and cell2 is a part of a boundary of cell B.
///If cell A's manifold dimension is n, then this cell's manifold dimension is n-1.
///B's manifold dimension is same as A's. and cell2's manifold dimension is n-1.
void MGPCell::bind(MGPCell& cell2){
	assert(manifold_dimension()==cell2.manifold_dimension());

	MGCell* star_cell1=star();
	MGCell* star_cell2=cell2.star();

	MGComplex* complex=star_cell1->parent_complex();
	MGComplex* complex2=star_cell2->parent_complex();
	assert(!complex || !complex2 || complex==complex2);
	//If both has complexes, they must be the same.

	if(complex&&!complex2){
		complex->append_pcell(star_cell2);
	}else if(!complex && complex2){
		complex2->append_pcell(star_cell1);
	}
	merge_2bcells(*this, cell2);
}

///Merge two binder cells of pcell1 and pcell2.
///When both had binders, binder cell of pcell2 will be destructed.
void merge_2bcells(MGPCell& pcell1, MGPCell& pcell2){
	assert(pcell1.manifold_dimension()==pcell2.manifold_dimension());
	SharedBCell& bcell1=pcell1.binder();
	SharedBCell& bcell2=pcell2.binder();
	if(bcell1){
		if(bcell1==bcell2)
			return;
		if(bcell2){//If both had binder, merge them.
			//Merge extent.
			if(!bcell1->extentBC() && bcell2->extentBC())
				bcell1->set_extentBC(UniqueGeometry(bcell2->release_extentBC()));

			//Merge own binder.
			size_t npatnr=bcell2->number_of_partner_members();
			for(size_t i=0; i<npatnr; i++){
				//Change binder from bcell2 to this.
				const MGPCell* pcell2Partner=bcell2->partner_member((int)i);
				pcell2Partner->resetBinder(bcell1);
			//We do not update bcell2's old MGBCell partner members
			//since the MGBCell will be destructed.
			}
		}else{//If bcell1 had binder and bcell2 had not.
			pcell2.resetBinder(bcell1);
		}
	}else if(bcell2){//If bcell2 had binder and bcell1 had not.
		pcell1.resetBinder(bcell2);
	}else{//Both had not binders.
		pcell2.resetBinder(pcell1.make_binder());
	}
	assert(bcell1->number_of_partner_members()==bcell1.use_count());
	assert(bcell1.get()==bcell2.get());
}

///Make a binder associated with the extent geometry rep.
///If this already had the binder with the extent,
///make_binder_with_curve() only returns the reference.
///If this had the binder without extent, make_binder_with_extent() will
///generate the extent.
std::shared_ptr<MGBCell>& MGPCell::make_binder_with_extent()const{
	make_binder();
	if(!m_binder->extentBC()){
		//If binder did not have curve representation, make it.
		m_binder->set_extentBC(make_binder_extent());
	}
	return m_binder;
}

//Return number of partners.
//Partners do not inclue own cell.
int MGPCell::number_of_partners() const{
	int n=0;
	if(m_binder)
		n=int(m_binder->number_of_partner_members())-1;
	//decrease 1 since Partners do not inclue own cell.
	return n;
}

///Obtain partner cells.
///Partners represent same world's(same cell's parameter) coordinates.
///The partners do not include this pcell except when star cell is
///connected to the star cell itself(closed only by the star cell).
///Let prtnrs[.] is the function's output, then prtners[0] is
///the netxt partner of this, and prtnrs[last] is the previous one of this
///in the partners array of the binder.
std::vector<const MGPCell*> MGPCell::partners() const{
	const std::shared_ptr<MGBCell>& bcell=binder();
	if(bcell){//If have binder, get bcell's partner.
		int m=number_of_partners();
		std::vector<const MGPCell*> prtnrs(m);
		int i,j=0; const MGPCell* prtnr;
		//Find this in bcell->m_partners.
		for(i=0; i<=m; i++){
			prtnr=bcell->partner_member(i);
			if(prtnr==this)
				break;
		}
		int mp1=m+1;
		for(j=0; j<m; j++){
			int id=++i%mp1;
			prtnrs[j]=bcell->partner_member(id);
		}
		return prtnrs;
	}else
		return std::vector<const MGPCell*>();
}

///Reset the binder.
void MGPCell::resetBinder(){
	SharedBCell& bdr = binder();
	if(bdr){
		bdr->free_partner_member(this);
		bdr.reset();
	}
}

///Set the binder.
void MGPCell::resetBinder(SharedBCell& newBinder)const{
	if(m_binder)
		m_binder->free_partner_member(this);
	m_binder = newBinder;
	if(m_binder)
		m_binder->add_partner_member(*this);
}

///Rest the binder by newBinder. newBinder's ownership is transfered to this as
///SharedBCell.
void MGPCell::resetBinder(MGBCell* newBinder)const{
	auto bndr = SharedBCell(newBinder);
	resetBinder(bndr);
}

//Make this binder the same as pcell2's binder.
//pcell2 must have binder already.
void MGPCell::shareBinder(const MGPCell & pcell2){
	assert(pcell2.m_binder);
	m_binder = pcell2.m_binder;
	m_binder->add_partner_member(*this);
	assert(m_binder.use_count()==m_binder->number_of_partner_members());
}

//Make this binder the same as bcell's partner members
//if bcell already had a partner member.
//If bcell did not have partner members, reset bcell.
void MGPCell::shareBinderIfHadPartner(MGBCell* bcell){
	const MGPCell* prtnr = bcell->first_partner_member();
	prtnr ?	shareBinder(*prtnr) : resetBinder(bcell);
}

// Output virtual function.
std::ostream& MGPCell::toString(std::ostream& ostrm) const{
	ostrm << ",";
	if(m_binder){
		const MGBCell* bcell = m_binder.get();
		ostrm<<bcell->whoami();
		ostrm << "=" << (MGGel*)bcell;
	}else
		ostrm << "binder=Null";
	return ostrm;
}
