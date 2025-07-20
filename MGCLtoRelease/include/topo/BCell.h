/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#ifndef _MGBCell_HH_
#define _MGBCell_HH_

#include <vector>
#include "mg/Position.h"
#include "topo/PCell.h"

//
//Defines MGBCell Class.

class MGBox;
class MGGeometry;
using PCellVec = std::vector<const MGPCell*>;

/** @addtogroup TOPO
 *  @{
 */

///MGBCell is a binder cell that connects MGPCell partners.

///MGBCell ia an abstract class and the known subclass are MGBVertex and MGEdge.
///There are two types of cells. One is parameter cell(MGPCell) and
///the other is binder cell(MGBCell). They are exclusive, that is, if
///a cell is a parameter cell, the cell cannot be binder cell and
///vice versa.
///Binder cell is a binder of parameter cells. Plural parameter cells are bound through a binder.
class MG_DLL_DECLR MGBCell:virtual public MGGel{

protected:
///The vector of partnerÅ@member parameter cells(MGPCell) who share this bcell.
mutable std::unique_ptr<PCellVec> m_partners;

public:

///Iterator definition.
typedef std::vector<MGPCell*>::iterator iterator;
typedef std::vector<const MGPCell*>::const_iterator const_iterator;
typedef std::vector<MGPCell*>::reverse_iterator reverse_iterator;
typedef std::vector<const MGPCell*>::const_reverse_iterator const_reverse_iterator;

///////Constructor/////////

///Void constructor. Constructor of bcell.
MGBCell(const MGBCell& cell);
MGBCell& operator=(const MGBCell& gel2);
explicit MGBCell(const MGPCell* pcell);//Construct as one partner member of pcell.
MGBCell(std::vector<const MGPCell*>& pcells);//Construct of partner member of pcells.
MGBCell(MGBCell&& rhs);//Move constructor.
MGBCell& operator=(MGBCell&& rhs);//Move assignment.
virtual ~MGBCell()=default;

///Assignment.

///When the leaf objects of this and gel2 are not equal, this assignment
///does nothing.
virtual MGGel& operator=(const MGGel& gel2);

///Transform the binder.
virtual void binder_tr(const MGVector& v) = 0;//Translation.
virtual void binder_tr(double s) = 0;//Scaling.
virtual void binder_tr(const MGMatrix& mat) = 0;//Matrix transform.
virtual void binder_tr(const MGTransf& tr) = 0;//Trans transform.

///Get extent geometry, may be null if this does not have extent.
virtual const MGGeometry* extentBC() const=0;
virtual MGGeometry* extentBC()=0;

//Obtain the 1st partner member.
//When this had no partners, null is returned.
const MGPCell* first_partner_member()const;

///Make sure that this has an extent expression.
///When this did not have an extent, make the extent from the partner
///member's parameter expression and the star cell.
///This must be a binder cell that has partner members that are
///boundaries. When this is not the case or this had an extent already,
///it does nothing.
virtual void make_extent() const=0;

///Obtain the i-th member partner. This must be a binder cell.
virtual const MGPCell* partner_member(int i)const{return (*m_partners)[i];};

///Obtain member partners. This must be a binder cell.
///before use of partner_members number_of_partner_members()>=1 must be assured.
const std::vector<const MGPCell*>& partner_members()const{return *m_partners;};
std::vector<const MGPCell*>& partner_members(){return *m_partners;};

///Return nummber of partners stored in m_partners.
size_t number_of_partner_members() const{return m_partners?m_partners->size():0;};

///Release the extent of this binder cell.
virtual MGGeometry* release_extentBC()=0;

///Set extent of this binder cell.
virtual void set_extentBC(std::unique_ptr<MGGeometry>&& extent)=0;
virtual void set_extent_as_nullBC() = 0;

/// Output virtual function.
virtual std::ostream& toString(std::ostream&) const;

protected:

MGBCell()=default;

///Add partner to this m_partners vector.
///Do not update partner's m_binder.
void add_partner_member(const MGPCell& partner);

///Clone this BCell, building the new partner membership of partners.
virtual SharedBCell cloneWithPartnerMembers(std::vector<const MGPCell*>& newPartners)const = 0;

///Read Object's member data.
virtual void ReadMembers(MGIfstream& buf) override;

///Write Object's Member Data
virtual void WriteMembers(MGOfstream& buf) const override;

private:

///Free specified partner(cellin).
void free_partner_member(const MGPCell* cellin) const;

///Merge two binder cells of pcell1 and pcell2.
///Boht had binder, binder cell of pcell2 will be destructed.
friend void merge_2bcells(MGPCell& pcell1, MGPCell& pcell2);

friend class MGPCell;
friend class MGEdge;
friend class MGComplex;
friend class MGOfstream;
friend class MGIfstream;

};

/** @} */ // end of TOPO group
#endif
