/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/

#include "StdAfx.h"
#include "mg/GelPosition.h"
#include "topo/Shell.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//
//Implements MGGelPosition Class.
//MGGelPosition is a class which expresses which group a gel belongs to.
//


///Void constructor.
MGGelPosition::MGGelPosition()
:m_group(nullptr), m_object(nullptr){
}

///Constructor of no hierarched group(m_Ghierarchy.size()==0).
MGGelPosition::MGGelPosition(MGGroup* group, MGObject* obj)
:m_group(group), m_object(obj){
}

//Equal operator
bool MGGelPosition::operator==(const MGGelPosition& gelp2) const{
	if(m_group!=gelp2.m_group)
		return false;
	if(m_object!=gelp2.m_object)
		return false;
	size_t n1=m_Ghierarchy.size(), n2=gelp2.m_Ghierarchy.size();
	if(n1!=n2)
		return false;
	for(size_t i=0; i<n1; i++)
		if(m_Ghierarchy[i]!=gelp2.m_Ghierarchy[i])
			return false;
	return true;
}

bool MGGelPosition::operator<(const MGGelPosition& gelp2)const{
	const MGGroup* g1=top_group();
	const MGGroup* g2=gelp2.top_group();
	if(g1!=g2)
		return g1<g2;

	size_t n1=m_Ghierarchy.size(), n2=gelp2.m_Ghierarchy.size();
	if(n2>n1)
		return false;
	if(n1>n2)
		return true;

	for(size_t i=0; i<n1; i++){
		if(m_Ghierarchy[i]!=gelp2.m_Ghierarchy[i])
			return m_Ghierarchy[i]<gelp2.m_Ghierarchy[i];
	}
	if(m_object && gelp2.m_object)
		return m_object->ordering_test(*gelp2.m_object)<0;
	else
		return m_object <gelp2.m_object;
}

//String output function.
std::ostream& operator<<(std::ostream& out, const MGGelPosition& gelp){
	out<<"MGGelPosition::top_group="<<gelp.m_group
		<<",obj="<<gelp.m_object<<std::endl;
	size_t n=gelp.m_Ghierarchy.size();
	if(n){
		out<<"num of mid hierary="<<n<<":";
		for(size_t i=0; i<n; i++){
			out<<i<<":"<<gelp.m_Ghierarchy[i];
			if(i<(n-1))
				out<<", ";
		}
		out<<std::endl;
	}
	return out;
}

///Append lower level group or shell data.
void MGGelPosition::append_lower_gel(MGGel* gel){
	assert(dynamic_cast<MGGroup*>(gel) || dynamic_cast<MGShell*>(gel)
	|| dynamic_cast<MGAttrib*>(gel));
	m_Ghierarchy.push_back(gel);
}

bool buildGroupHierarchyOneHierarcy(
	const MGGroup& group,//scanning target group.
	MGObject* obj,
	std::vector<MGGroup*>& groupHierarchy, //Hierarchy so far.
	MGGelPosition& gp
) {
	for (auto& member : group.m_list) {
		MGGel* geli = member.get();

		MGObject* memberObj = dynamic_cast<MGObject*>(geli);
		if (memberObj && memberObj == obj) {
			for (auto& groupH : groupHierarchy)
				gp.append_lower_gel(groupH);
			return true;
		}

		MGGroup* memberGroup = dynamic_cast<MGGroup*>(geli);
		if (memberGroup) {
			groupHierarchy.push_back(memberGroup);
			if(buildGroupHierarchyOneHierarcy(*memberGroup,obj, groupHierarchy, gp))
				return true;
			groupHierarchy.pop_back();
		}
	}
	return false;
}

/// <summary>
/// Build MGGelPosition of MGGroup hierarchy whose top group is groupRoot.
/// Although gel may not belong to groupRoot, must belong to a member group
/// (or a member group of the member group, ....) of groupRoot.
/// </summary>
void MGGelPosition::buildGroupHierarchy(MGGroup* groupRoot, MGObject* obj){
	m_group = groupRoot;
	m_Ghierarchy.clear();
	m_object = obj;

	std::vector<MGGroup*> groupHierarchy;//Group hierarchy except groupRoot.
	if (!buildGroupHierarchyOneHierarcy(*groupRoot, obj, groupHierarchy, *this))
		set_null();
}

//Get the lowerest level of group pointer of this.
const MGGroup* MGGelPosition::bottom_group_sub()const{
	const MGGroup* grp=m_group;
	int nhierk=(int)m_Ghierarchy.size();
	if(nhierk--){
		const MGGel* gel=m_Ghierarchy[nhierk--];
		const MGGroup* grp2=dynamic_cast<const MGGroup*>(gel);
		if(grp2)
			grp=grp2;
		else{
			if(nhierk>=0){
				gel=m_Ghierarchy[nhierk];
				grp= dynamic_cast<const MGGroup*>(gel);;
			}
		}
	}
	return grp;
}

///Get the group pointer.
const MGGroup* MGGelPosition::bottom_group()const{
	return bottom_group_sub();
}
MGGroup* MGGelPosition::bottom_group(){
	const MGGroup* grp=bottom_group_sub();
	return const_cast<MGGroup*>(grp);
}

//Generate a newed clone object.
MGGelPosition* MGGelPosition::clone()const{
	return new MGGelPosition(*this);
}

//Get the shell pointer when is_shell_face() is true.
//When is_shell_face() is false, behavior is undefined.
MGShell* MGGelPosition::get_shell_of_shell_face()const {
	assert(is_shell_face());
	return dynamic_cast<MGShell*>(m_Ghierarchy.back());
}

///Test if this is null.
bool MGGelPosition::is_null()const{
	if(m_Ghierarchy.size())
		return false;
	return m_object == nullptr;
}

///Test if this is MGGelpotion that point shell and the member face.
///That is, m_object is MGFace and top_object() is MGShell.
bool MGGelPosition::is_shell_face()const{
	if(!m_object)
		return false;
	if(top_object_sub()==m_object)
		return false;//Since m_object=MGFace and m_Ghierarchy.back()=MGShell.

	assert(dynamic_cast<const MGFace*>(m_object)&& 
		dynamic_cast<const MGShell*>(m_Ghierarchy.back())) ;
	return true;
}

//Test if this is one of the types of types.
bool MGGelPosition::is_type(const MGAbstractGels& types)const{
	const MGGel* gel = leafGel();
	bool isType= gel ? gel->type_is(types):false;

	if(!isType && is_shell_face())
		isType= top_object()->type_is(types);

	return isType;
}

//Test if this leaf is MGGroup.
bool MGGelPosition::leaf_isAttrib()const{
	if(m_object)
		return false;
	size_t nhierk=m_Ghierarchy.size();
	if(!nhierk)
		return false;
	return dynamic_cast<const MGAttrib*>(m_Ghierarchy.back())!=nullptr;
}

bool MGGelPosition::leaf_isAttribedGel() const {
	return leafAttribedGel()!=nullptr;
}

//Test if this leaf is MGGroup.
bool MGGelPosition::leaf_isGroup()const {
	if (m_object)
		return false;
	size_t nhierk = m_Ghierarchy.size();
	if (!nhierk)
		return false;
	return dynamic_cast<const MGGroup*>(m_Ghierarchy.back()) != nullptr;
}

bool MGGelPosition::leaf_isObject() const{
	if (m_object == nullptr || is_shell_face())
		return false;
	return true;
}

///Get the leaf object or group pointer of this.
///Returned is MGObject that is not is_shell_face(), or MGGroup.
const MGAttribedGel* MGGelPosition::leafAttribedGel()const {
	const MGAttribedGel* agel = nullptr;
	if (m_object)
		agel= m_object;
	else if (m_Ghierarchy.size()) 
		agel=dynamic_cast<const MGGroup*>(m_Ghierarchy.back());
	return agel;
}
MGAttribedGel* MGGelPosition::leafAttribedGel() {
	const MGGelPosition& gp = *this;
	return const_cast<MGAttribedGel*>(gp.leafAttribedGel());
}

///Get the leaf gel pointer of this.
///Returned is MGGel excluding Shell.
const MGGel* MGGelPosition::leafGel()const {
	const MGGel* gel = nullptr;
	if (m_object)
		gel = m_object;
	else if (m_Ghierarchy.size())
		gel = dynamic_cast<const MGGroup*>(m_Ghierarchy.back());
	return gel;
}
MGGel* MGGelPosition::leafGel() {
	const MGGelPosition& gp = *this;
	return const_cast<MGAttribedGel*>(gp.leafAttribedGel());
}


///Set the leaf object data.
void MGGelPosition::set_leafAttrib(MGAttrib* agel){
	m_object =nullptr;
	append_lower_gel(agel);
}

///Set the leaf object data.
void MGGelPosition::set_leafGroup(MGGroup* grp) {
	m_object = nullptr;
	append_lower_gel(grp);
}

///Set the leaf object data.
void MGGelPosition::set_leaf_object(MGObject* obj){
	m_object=obj;
}

//Get the object pointer of this.
const MGObject* MGGelPosition::top_object_sub()const{
	const MGObject* obj=0;
	size_t nhierk = m_Ghierarchy.size();
	if(nhierk)
		obj=dynamic_cast<const MGObject*>(m_Ghierarchy[nhierk-1]);
	if(obj)
		return obj;//This must be MGShell.
	return m_object;
}
const MGObject* MGGelPosition::top_object()const{
	return top_object_sub();
}

MGObject* MGGelPosition::top_object(){
	const MGObject* obj=top_object_sub();
	return const_cast<MGObject*>(obj);
}

///Clear the content.
void MGGelPosition::set_null(){
	m_group=0;
	m_Ghierarchy.clear();
	m_object=0;
}

//Test if this is symmetric to gel2.
//Symmetric means:
///Both leaf objects are MGObject and they have the same manifold dimension.
bool MGGelPosition::symmetric(const MGGelPosition& gelp2)const{
	const MGObject* obj1=leaf_object();
	if(!obj1)
		return false;
	const MGObject* obj2=gelp2.leaf_object();
	if(!obj2)
		return false;

	return (obj1->manifold_dimension()==obj2->manifold_dimension());
}

//Perform add operation of this gel position.
//(insert the object of m_object in the lowerest level of group).
//This is valid only for top_object() is not shell
//(unable to add face into shell).
void MGGelPosition::do_add(){
	MGGroup* group_up;
	MGGel* gel_to_add;
	MGGroup* grp=bottom_group();

	if(is_shell_face()){
	//When is_shell_face(), the target to add is MGShell, and not MGFace.
		group_up=grp;
		gel_to_add=top_object();//This is MGShell.
	}else if(leaf_isGroup()){
		gel_to_add=grp;
		size_t nhierak=m_Ghierarchy.size();
		if(nhierak>=2)
			group_up= dynamic_cast<MGGroup*>(m_Ghierarchy[nhierak-2]);
		else
			group_up=m_group;
	}else{
		group_up=grp;
		gel_to_add= m_object;
	}
	group_up->append(gel_to_add);
}

//Perform remove operation of this gel position.
//(Release the gel of m_gel from the group of m_group, but does not delete the gel).
void MGGelPosition::do_remove(){
	MGGroup* group_up;
	MGGel* gel_to_remove;
	MGGroup* grp=bottom_group();

	if(is_shell_face()){
		//When is_shell_face(), the target to remove is MGShell, and not MGFace.
		group_up=grp;
		gel_to_remove=top_object();
	}else if(leaf_isGroup()){
		gel_to_remove=grp;
		size_t nhierak=m_Ghierarchy.size();
		if(nhierak>=2)
			group_up=dynamic_cast<MGGroup*>(m_Ghierarchy[nhierak-2]);
		else
			group_up=m_group;
	}else{
		group_up=grp;
		gel_to_remove=m_object;
	}

	MGGroup::reverse_iterator i=group_up->rbegin(), iend=group_up->rend(), inext;
	while(i!=iend){
		inext=i; inext++;
		if(i->get()==gel_to_remove){
			i->release();
			group_up->erase(inext.base());//release the found gel
			break;
		}
		i=inext;
	}
}
