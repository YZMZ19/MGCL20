/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#ifndef _MGBVertex_HH_
#define _MGBVertex_HH_

#include "mg/Default.h"
#include "mg/Position.h"
#include "mg/Point.h"
#include "topo/BCell.h"

class MGBox;
class MGOfstream;
class MGIfstream;
class MGPVertex;

//
//Define MGBVertex Class.

/** @addtogroup TOPO
 *  @{
 */

///MGBVertex is 0 manifold dimension binder cell, is an point.

///MGBVertex is a binder cell of MGPVertex, and has manifold dimension 0.
///MGBVertex is not used as a parameter cell(boundary of an edge).
///Since MGBVertex's manifold dimension is 0, MGBVertex does not have
///boundaries.
class MG_DLL_DECLR MGBVertex:public MGBCell{
	friend MGPVertex;

private:
	mutable std::unique_ptr<MGPoint> m_point;

public:
///////// Constructor /////////

///Void constructor.

///Special constructors.
MGBVertex() = default;
MGBVertex(const MGBVertex& bv);
MGBVertex(MGBVertex&& rhs)=default;//Move constructor.
~MGBVertex() = default;
MGBVertex& operator=(const MGBVertex& gel2);
MGBVertex& operator=(MGBVertex&& rhs) = default;//Move assignment.

//Constructor that holds pv as a partner member.
MGBVertex(const MGPVertex& pv);

///Comparison of two objects.

///Comparison.
bool operator==(const MGBVertex& gel2)const;
std::partial_ordering operator<=>(const MGBVertex& gel2)const;

//gel2 must be the same class as this.
bool equal_test(const MGGel& gel2)const override;

//gel2 must be the same class as this.
std::partial_ordering ordering_test(const MGGel& gel2)const override;

/////////Member Function/////////

///Transform the binder.
void binder_tr(const MGVector& v) override;//Translation.
void binder_tr(double s) override;//Scaling.
void binder_tr(const MGMatrix& mat) override;//Matrix transform.
void binder_tr(const MGTransf& tr) override;//Trans transform.

///Generate copied gel of this gel.
///Returned is a newed object. User must delete the object.
MGBVertex* clone()const override{return new MGBVertex(*this);};

///Draw point(vertex) in this coordinates.
void drawVertex(mgVBO& vbo)const;

///Get extent geometry, may be null if this does not have extent.
const MGPoint* extentBC() const override;
MGPoint* extentBC()  override;

///Return Object's type ID (TID)
long identify_type()const override;

///Make sure that this has an extent expression.
///When this did not have an extent, make the extent from the partner
///member's parameter expression and the star cell.
///This must be a binder cell that has partner members that are
///boundaries. When this is not the case or this had an extent already,
///it does nothing.
void make_extent() const override;

///Get manifold dimension.
int manifold_dimension()const override{return 0;};

///Obtain the i-th member partner PVertex.
const MGPVertex* partner_member_vertex(int i)const;

/// Output function.
std::ostream& toString(std::ostream&) const;

///Compute the parameter value of the closest point from the straight to
///this object.
///sl is the eye projection line whose direction is from yon to hither, and if
///sl had multiple intersection points, The closest point to the eye will be selected.
///This will be never invoked.
//MGPosition pick_closest(const MGStraight& sl)const;

///Return extent data(i.e. MGPoint*), which may be null if this has no extent.
const MGPoint* point() const;

///Get the position data of this vertex.
///Returns MGPoint data if this has extent. Otherwise obtain from partner's
///star edge by evaluating the edge's data.
MGPosition position() const;

///Release the extent of this binder cell.
MGGeometry* release_extentBC() override;

///Set extent of this binder cell.
void set_point(const MGPosition& P);
void set_extentBC(std::unique_ptr<MGGeometry>&& extent)override;
void set_extent_as_nullBC() override;

protected:

//Constructor that holds newPartners as partner members.
MGBVertex(std::vector<const MGPCell*>& newPartners);

///Get the name of the class.
std::string whoami()const override{ return "BV"; };

///Read Object's member data.
void ReadMembers(MGIfstream& buf);

///Write Object's Member Data
void WriteMembers(MGOfstream& buf) const;

private:

///Clone this BCell, building the new partner membership of partners.
SharedBCell cloneWithPartnerMembers(std::vector<const MGPCell*>& newPartners)const override;

/// Serialization.
MG_DLL_DECLR friend MGOfstream& operator<< (MGOfstream& buf, const MGBVertex& bv);
MG_DLL_DECLR friend MGIfstream& operator>> (MGIfstream& buf, MGBVertex& bv);

};

/** @} */ // end of TOPO group
#endif
