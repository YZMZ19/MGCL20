/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#include "StdAfx.h"
#include "mgGL/Appearance.h"
#include "topo/BVertex.h"
#include "topo/Edge.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//
//Implements MGBCell Class.
//Cell is a cell without boundaries(No Boundaries).
//Boundaries are defined in the subclass, MGEdge, or MGFace.

//There are two types of cells. One is a parameter cell(pcell) and
//the other is a binder cell(bcell). They are exclusive, that is, if
//a cell is a parameter cell, the cell cannot be binder cell and
//vice versa.
//Binder cell is a binder of parameter cells. Plural cells are connected
//through a binder.
//MGBCell is an abstract class.

MGBCell::MGBCell(const MGBCell& cell){ ; }//We release all of the partner members.
MGBCell& MGBCell::operator=(const MGBCell& gel2){
	if(m_partners) m_partners->clear();
	return *this;
}
MGBCell::MGBCell(const MGPCell* pcell)
	:m_partners(std::make_unique< PCellVec>(1, pcell)){ ; }
MGBCell::MGBCell(std::vector<const MGPCell*>& pcells)
	:m_partners(std::make_unique< PCellVec>(pcells)){;}//Construct of partner member of pcells.
MGBCell::MGBCell(MGBCell&& rhs){ ; }//Move constructor.
MGBCell& MGBCell::operator=(MGBCell&& rhs){ m_partners.reset(); return *this; }//Move assignment.

///Assignment.

///When the leaf objects of this and gel2 are not equal, this assignment
///does nothing.
MGGel& MGBCell::operator=(const MGGel& gel2){
	MGBCell* bcel = dynamic_cast<MGBCell*>(this);
	if(bcel)
		(*this) = *bcel;
	return *this;
}

//Add partner to this binder cell.
void MGBCell::add_partner_member(const MGPCell& partner){
	if (m_partners) {
		auto comp = [](const MGPCell* cell1, const MGPCell* cell2) {
			return cell1->is_less_than(*cell2); };
		std::vector<const MGPCell*>::iterator
			i = std::lower_bound(m_partners->begin(), m_partners->end(), &partner, comp);
		m_partners->insert(i, &partner);
	}
	else
		m_partners = std::make_unique< PCellVec>(1,&partner);
	
}

//Obtain the 1st partner member.
//When this had no partners, null is returned.
const MGPCell* MGBCell::first_partner_member()const{
	const MGPCell* pcel = nullptr;
	if(number_of_partner_members())
		pcel= partner_member(0);
	return pcel;
}

//Free specified partner(cellin).
void MGBCell::free_partner_member(const MGPCell* cellin) const{
	if (m_partners) {
		const_iterator itr = m_partners->begin();
		const_iterator itre = m_partners->end();
		for (; itr != itre; itr++) {
			if (cellin == *itr) {
				m_partners->erase(itr);
				break;
			}
		}
	}
}

// Output virtual function.
std::ostream& MGBCell::toString(std::ostream& ostrm) const{
	size_t n = number_of_partner_members();
	ostrm << ", "<<n<<" membr";
	if(n){
		const MGPCell* pcel = partner_member(0);
		std::string pcname;
		if(pcel)
			pcname = pcel->whoami();
		ostrm << pcname;
		MGBCell::const_iterator ps = m_partners->begin(), pe = m_partners->end();
		if (ps != pe) {
			const MGPCell* pcel0 = *ps;
			ostrm << "=[";
			ostrm << (MGGel*)pcel0;
			for (ps++; ps != pe; ps++)
				ostrm << "," << (MGGel*)(*ps);
			ostrm << "]";
		}
	}
	return ostrm;
}
