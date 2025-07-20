/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#ifndef _MGGelPosition_HH_
#define _MGGelPosition_HH_

#include <iosfwd>
#include <vector>
#include "mg/MGCL.h"
#include "mg/Group.h"
#include "mg/AttribedGel.h"
class MGShell;

//
//Define MGGelPosition Class.

/** @file */
/** @addtogroup GelRelated
 *  @{
 */

///MGGelPosition is a class to locate where a target gel is in a group hierarchy.

/// Generally, A group includes other groups,
/// and the included groups include other groups. In this way, the groups make
/// a group hierachy. MGGelPosition represents this hierarcy.
///
/// Let n=m_Ghierarchy.size().
/// When n==0, m_group is the MGGroup that includes the object m_object .
/// When n >0, m_group includes m_Ghierarchy[0] (which is MGShell, MGGroup, or MGAttrib).
/// group m_Ghierarchy[i-1] includes m_Ghierarchy[i] for i=1,...,n-2.
/// m_Ghierarchy[n-1] includes m_object(if m_object!=0).
/// m_object is the leaf MGObject pointer if exist.
/// 
/// There are 3 kinds of target gel of MGGelPosition.
/// 1) The target gel is MGObject that is (directly) included in MGGroup.
///    The MGObject may be MGShell, independent MGFace that does not consititute MGShell),
///    or other general MGObject.
///    In this case, all of the members of m_Ghierarchy are MGGroup and
///    m_object is the target MGObject( is not null).
///    This case is tested by leaf_isObject().
/// 2) The target is MGFace that constitutes a MGShell, or MGShell.
///    In this case, m_Ghierarchy[i] for i=0, ..., n-2 are MGGroup,
///    m_Ghierarchy[n-1] is MGShell*, and m_object is MGFace*.
///    This case is tested by is_shell_face().
/// 3) The target is not MGObject but MGGel(currently MGGroup, or MGAttrib).
///    In this case, m_Ghierarchy[i] for i=0, ..., n-2 are MGGroup,
///    m_Ghierarchy[n-1] is the target MGGel(MGGroup*, or MGAttrib*),
///    and m_object is nullptr.
///    This case is tested by leaf_isGroup(), or leaf_isAttrib().
/// 
/// When n>0 and m_object=0, this means that leaf_isGroup() or leaf_isAttrib() is true and
/// the positioning target is MGGroup, or MGAttrib.
/// Although m_Ghierarchy[i] for i=0,...,n-2 are always MGGroup, m_Ghierarchy[n-1] may be
/// MGShell that includes MGFace, MGAttrib.
/// 
class MG_DLL_DECLR MGGelPosition{

public:

///String stream function.
MG_DLL_DECLR friend std::ostream& operator<< (std::ostream&, const MGGelPosition&);


////////Special member functions/////////
MGGelPosition();
virtual ~MGGelPosition()=default;
MGGelPosition(const MGGelPosition&)=default;//Copy constructor.
MGGelPosition& operator= (const MGGelPosition&)=default;///Copy assignment.
MGGelPosition(MGGelPosition&&)=default;		///Move constructor.
MGGelPosition& operator= (MGGelPosition&&)=default;///Move assignment.


///Constructor of no hierarched group(m_Ghierarchy.size()==0).
///When obj is not null, obj must belong to group.
explicit MGGelPosition(MGGroup* group, MGObject* obj=0);

///Equal operator
bool operator== (const MGGelPosition& gelp2) const;
bool operator!= (const MGGelPosition& gelp2) const{return !(*this==gelp2);}

///Comparison operator.
bool operator<(const MGGelPosition& gp2)const;
bool operator>(const MGGelPosition& gp2)const{return gp2<(*this);};
bool operator<=(const MGGelPosition& gp2)const{return !((*this)>gp2);};
bool operator>=(const MGGelPosition& gp2)const{return !(gp2>(*this));};

////////////Member Function////////////

///Append lower level group or shell data.
//gel must be a MGShell or MGGroup.
void append_lower_gel(MGGel* gel);

/// <summary>
/// Build MGGelPosition of MGGroup hierarchy whose top group is groupRoot.
/// Although gel may not belong to groupRoot, must belong to a member group
/// (or a member group of the member group, ....) of groupRoot.
/// </summary>
void buildGroupHierarchy(MGGroup* groupRoot, MGObject* obj);

///Return the MGGel i;
const MGGel* gel(int i)const{return m_Ghierarchy[i];};
MGGel* gel(int i){return m_Ghierarchy[i];};

///Generate a newed clone object.
virtual MGGelPosition* clone()const;

///perform add operation of this gel position.
///(push back the object of m_object to the group that includes m_object)
/// When is_shell_face(), 
void do_add();

///perform remove operation of this gel position.
///(Release the object of m_object from the group that includes m_object,
///but does not delete the gel).
void do_remove();

///Get the group pointer that includes leaf_gel();
const MGGroup* bottom_group()const;
MGGroup* bottom_group();

///Get the top group pointer.
const MGGroup* top_group()const{return m_group;};
MGGroup* top_group(){return m_group;};

//Get the shell pointer when is_shell_face() is true.
//When is_shell_face() is false, behavior is undefined.
MGShell* get_shell_of_shell_face()const;

///Test if this is null.
bool is_null() const;

///Test if this is MGGelpotion that point shell and the member face.
///That is, m_object is MGFace and top_object() is MGShell.
bool is_shell_face()const;

//Test if this is one of the types of types.
bool is_type(const MGAbstractGels& types)const;

//Test if this leaf is MGAttrib.
bool leaf_isAttrib()const;

/// <summary>
/// Test if this leaf is a MGAttribedGel.
/// </summary>
bool leaf_isAttribedGel()const;

//Test if this leaf is MGGroup.
bool leaf_isGroup()const;

/// <summary>
/// Test if this leaf is MGObject*, and not MGFace* that constitute MGShell.
/// </summary>
bool leaf_isObject()const;

/// Get the leaf object or group pointer of this.
/// Returned is MGObject or MGGroup.
/// When is_shell_face() is true, MGShell is not returned, but MGFace.
const MGAttribedGel* leafAttribedGel()const;
MGAttribedGel* leafAttribedGel();

/// Get the leaf gel pointer of this.
/// When is_shell_face() is true, MGShell is not returned, but MGFace.
const MGGel* leafGel()const;
MGGel* leafGel();

/// Get the leaf object pointer of this.
/// When is_shell_face() is true, MGShell is not returned, but MGFace.
const MGObject* leaf_object()const{return m_object;};
MGObject* leaf_object(){return m_object;};

///Set the leaf MGAttrib data.
void set_leafAttrib(MGAttrib* agel);

///Set the leaf MGGroup data.
void set_leafGroup(MGGroup* grp);

///Set the leaf object data.
void set_leaf_object(MGObject* obj);

///Set the group data.
void set_top_group(MGGroup* group){m_group=group;};

///Set this as null.
void set_null();

///Test if this is symmetric to gel2.
///Symmetric means:
///Both leaf objects are MGObject and they have the same manifold dimension.
bool symmetric(const MGGelPosition& gp2)const;

/// Get the top object pointer of this.
/// When is_shell_face() is true, MGShell is returned.
/// When is_shell_face() is false, top_object()== leaf_object().
const MGObject* top_object()const;
MGObject* top_object();

protected:
	MGGroup* m_group;///<The top group pointer which includes
		/// the object m_object if m_Ghierarchy.size()==0.
		/// If m_Ghierarchy.size()>0, m_group includes m_Ghierarchy[0].
		/// Generally m_group is MGGroup of .mgl file.

	std::vector<MGGel*> m_Ghierarchy;
		/// The group hierarchy that belong to m_group.
		///  Let n=m_Ghierarchy.size(),
		/// then group m_Ghierarchy[i-1] includes m_Ghierarchy[i] for n=0,...,n-2.
		/// m_Ghierarchy[n-1] includes m_object(if m_object!=0);
		/// Although m_Ghierarchy[i] for i=0,...,n-2 are always MGGroup,
		/// m_Ghierarchy[n-1] may be MGShell that includes a MGFace.
		/// In this case, m_object is the MGFace.

	MGObject* m_object;	///<The leaf MGObject pointer if the target gel is an MGObject
		/// (MGShell, MGFace that does or does not consititute MGShell, or a other MGObject).

private:
	//Get the object pointer of this.
	const MGObject* top_object_sub()const;

	//Get the lowerest level of group pointer of this.
	const MGGroup* bottom_group_sub()const;
};

/** @} */ // end of GelRelated group
#endif
