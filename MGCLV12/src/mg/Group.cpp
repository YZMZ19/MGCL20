/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#include "StdAfx.h"
#include "mg/Box.h"
#include "mg/AttribedGel.h"
#include "mg/Ifstream.h"
#include "mg/Ofstream.h"
#include "mg/Group.h"
#include "mg/Attrib.h"
#include "mg/GelFactory.h"
#include "mg/LBRep.h"
#include "mgGL/Context.h"
#include "mgGL/Name.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//
//Implements MGGroup Class.
//MGGroup is a class which constains MGGel elements.
//MGGroup provides functions: Serialization of MGGel elements.
//


//Construct MGGroup from the file made by MGOfstream.
//error contains error flag as:
//=0: open succeeded.
//=1: file not found, or could not be opened.
//=2: file found, but, the format is not MGCL format.
//When error!=0, MGGroup contains no MGGel's.
MGGroup::MGGroup(const TCHAR* file, int& error){
	MGIfstream infile;
	error=infile.open(file);
	if(error) return;
	load(infile);
}

bool MGGroup::operator<(const MGGroup& gel2)const{return size()<gel2.size();};
bool MGGroup::operator<(const MGGel& gel2)const{
	const MGGroup* gel2_is_this=dynamic_cast<const MGGroup*>(&gel2);
	if(gel2_is_this)
		return operator<(*gel2_is_this);
	return identify_type()<gel2.identify_type();
}

//////////Member Function//////////

//Get the MGAppearance pointer in this group. If not defined, null will be
//returned.
MGAppearance* MGGroup::appearance(){
	MGAppearance* app=nullptr;
	if(size()){
		MGGel* gel=front().get();
		app=dynamic_cast<MGAppearance*>(gel);
		if(!app){
			MGContext* cnt=dynamic_cast<MGContext*>(gel);
			if(cnt)
				app=cnt->appearance();
		}
	}
	return app;
}

const MGAppearance* MGGroup::appearance()const{
	MGGroup* grp=const_cast<MGGroup*>(this);
	return grp->appearance();
}

//Get the box of the group.
//If no objects were included in this group, null box will be returned.
MGBox MGGroup::box()const{
	MGBox bx;
	const_iterator i=begin(), ie=end();	
	for(; i!=ie; i++){
		const MGObject* obj=dynamic_cast<const MGObject*>(i->get());
		if(obj){
			const MGBox& bxobj=obj->box();
			if(bxobj.finite())
				bx|=bxobj;
			continue;
		}
		const MGGroup* grp=dynamic_cast<const MGGroup*>(i->get());
		if(grp){
			bx|=grp->box();
			continue;
		}
	}
	return bx;
}

//Generate copied gel of this gel.
//Returned is a newed object. User must delete the object.
MGGroup* MGGroup::clone()const{
	MGGroup* grpNew=new MGGroup;
	const_iterator i=begin(), ie=end();	
	for(; i!=ie; i++){
		grpNew->append((**i).clone());
	}
	return grpNew;
}

//Get the MGContext pointer stored in this group. If not defined, null will be
//returned.
MGContext* MGGroup::context(){
	if(!size())
		return nullptr;
	return dynamic_cast<MGContext*>(front().get());
}
const MGContext* MGGroup::context()const{
	if(!size())
		return nullptr;
	return dynamic_cast<const MGContext*>(front().get());
}

//////display member function.
void MGGroup::display_arrows(mgSysGL& sgl)const{
	for(auto& gel:m_list){
		auto obj = dynamic_cast<const MGObject*>(gel.get());
		if(obj)
			obj->display_arrows(sgl);
	}
}
void MGGroup::display_break_points(mgSysGL& sgl)const{
	for(auto& gel:m_list){
		auto obj = dynamic_cast<const MGObject*>(gel.get());
		if(obj)
			obj->display_break_points(sgl);
	}
}
void MGGroup::display_control_polygon(mgSysGL& sgl)const{
	for(auto& gel:m_list){
		auto obj = dynamic_cast<const MGObject*>(gel.get());
		if(obj)
			obj->display_control_polygon(sgl);
	}
}
void MGGroup::display_curvatures(
	mgSysGL& sgl,
	int		density,//densitiy of the graph.
	bool	use_radius,//true:radius display, false:curvature display.
	double	scale	//scaling of the graph.
)const{
	for(auto& gel:m_list){
		auto obj = dynamic_cast<const MGObject*>(gel.get());
		if(obj)
			obj->display_curvatures(sgl,density,use_radius,scale);
	}
}

//make this group has appearance and get the MGAppearance pointer.
MGAppearance* MGGroup::ensure_appearance(){
	MGAppearance* app=appearance();
	if(!app){
		app=new MGAppearance();
		push_appearance(app);
	}
	return app;
}

//set the copy of appr2 to this MGAttribedgel.
void MGGroup::set_appearance(const MGAppearance& appr2){
	push_appearance(appr2.clone());
}
void MGGroup::set_appearance(MGAppearance* appr2){
	push_appearance(appr2);
}
//Find the position of the gel in the gel list and the group pointer
//which includes the gel. Searching will be done into the member group gel
//of this list.
MGGroup::const_iterator MGGroup::find(
	const MGGel* gel, const MGGroup*& grp
)const{
	std::vector<const MGGroup*> grps;
	const_iterator i=begin(), ie=end();	
	for(; i!=ie; i++){
		MGGel* gel2=i->get();
		if(gel2==gel){ grp=this; return i;}
		const MGGroup* grp=dynamic_cast<const MGGroup*>(gel2);
		if(grp) grps.push_back(grp);
	}
	size_t ngrp(grps.size());
	for(size_t j=0; j<ngrp; j++){
		const MGGroup* gelj=grps[j];
		const_iterator pos=gelj->find(gel,grp);
		if(grp) return pos;
	}
	grp=0;
	return end();
}
MGGroup::iterator MGGroup::find(
	MGGel* gel, MGGroup*& grp
){
	std::vector<MGGroup*> grps;
	iterator i=begin(), ie=end();	
	for(; i!=ie; i++){
		MGGel* gel2=i->get();
		if(gel2==gel){ grp=this; return i;}
		MGGroup* grp=dynamic_cast<MGGroup*>(gel2);
		if(grp) grps.push_back(grp);
	}
	size_t ngrp(grps.size());
	for(size_t j=0; j<ngrp; j++){
		MGGroup* gelj=grps[j];
		iterator pos=gelj->find(gel,grp);
		if(grp) return pos;
	}
	grp=0;
	return end();
}

MGGroup::iterator MGGroup::search_by_id(MGGEL_TID tid){
	iterator i=begin(), ie=end();
	for(;i!=ie; i++){
		if((*i)->identify_type()==tid) return i;
	}
	return ie;
}
MGGroup::const_iterator MGGroup::search_by_id(MGGEL_TID tid)const{
	const_iterator i=begin(), ie=end();
	for(;i!=ie; i++){
		if((*i)->identify_type()==tid) return i;
	}
	return ie;
}

///Find the 1st gel(MGGroup or MGObject) member in this MGGroup
///that has the given name.
///**********Searching IS performed into the member group recursively.
void MGGroup::nameSearch(
	const MGName& name,
	MGGroup*& groupIncluding,///<If found the group including the gel of the name is set,
					///else, null will be returned.
	MGAttribedGel*& gelFound///<found gel will be returned if found, else null.
){
	search_by_nameSub(name,groupIncluding,gelFound);
}

///Find the 1st gel(MGGroup or MGObject) member in this MGGroup
///that has the given name.
///**********Searching IS performed into the member group recursively.
void MGGroup::nameSearch(
	const MGName& name,
	const MGGroup*& groupIncluding,///<If found, the group including the gel of the name is set,
					///else, null will be returned.
	const MGAttribedGel*& gelFound///<found gel will be returned if found, else null.
)const{
	MGGroup* grp2;
	MGAttribedGel* gel2;
	MGGroup* group=const_cast<MGGroup*>(this);
	group->search_by_nameSub(name,grp2,gel2);
	groupIncluding=grp2;
	gelFound=gel2;
}

///Find the 1st gel(MGGroup or MGObject) that has the given name in this MGGroup.
///**********Searching IS performed into the member group recursively.
void MGGroup::search_by_nameSub(
	const MGName& name,
	MGGroup*& groupIncluding,///<If found, the group including the gel of the name is set,
					///else, null will be returned.
	MGAttribedGel*& gelFound///<found gel will be returned if found, else null.
){
	groupIncluding=0;
	gelFound=0;

	iterator i=begin(), ie=end();
	for(;i!=ie; i++){
		MGGel* geli=i->get();
		MGAttribedGel* agel=dynamic_cast<MGAttribedGel*>(geli);
		if(!agel)
			continue;
		const MGName* namei=agel->get_name();
		if(namei){
			if((*namei) == name){
				groupIncluding=this;
				gelFound=agel;
				return;
			}
		}
		MGGroup* grpi=dynamic_cast<MGGroup*>(agel);
		if(grpi){
			grpi->search_by_nameSub(name,groupIncluding,gelFound);
			if(gelFound)
				return;
		}
	}
}

//Test if this gel includes an object.
//Function's return value is the 1st object found in the gel list of this
//if this includes an object. Otherwise null will be returned.
const MGObject* MGGroup::includes_object()const{
	const_reverse_iterator i=rbegin(), ie=rend();	
	for(; i!=ie; i++){
		MGGel* geli=i->get();
		const MGObject* obj=dynamic_cast<const MGObject*>(geli);
		if(obj)
			return obj;
		const MGGroup* grp=dynamic_cast<const MGGroup*>(geli);
		if(grp){
			obj=grp->includes_object();
			if(obj)
				return obj;
		}
	}
	return 0;
}
MGObject* MGGroup::includes_object(){
	return
		const_cast<MGObject*>(const_cast<const MGGroup*>(this)->includes_object());
}

//Count up how many MGGroup members are included in this group.
//Function's return value is the number of member group.
//Only the members of this group are counted.
int MGGroup::getMemberGroupNumber()const{
	const_iterator i=begin(), ie=end();
	int ngroup=0;
	for(; i!=ie; i++){
		MGGel* geli=i->get();
		const MGGroup* grp=dynamic_cast<const MGGroup*>(geli);
		if(grp)
			ngroup++;
	}
	 return ngroup;
}

//Get i-th MGGroup pointer of this member group.
///If i-th MGGroup except MGAppearance is not found, null is returned.
const MGGroup* MGGroup::get_i_th_MemberGroup(int i)const{
	const_iterator igrp=begin(), igrpend=end();
	for(int j=0; igrp!=igrpend; igrp++){
		MGGel* geli=igrp->get();
		const MGGroup* grp=dynamic_cast<const MGGroup*>(geli);
		if(grp){
			if(j==i)
				return grp;
			j++;
		}
	}
	return 0;
}
MGGroup* MGGroup::get_i_th_MemberGroup(int i){
	const MGGroup* grpC=this;
	const MGGroup* foundG=grpC->get_i_th_MemberGroup(i);
	return const_cast<MGGroup*>(foundG);
}

//Make an MGOfstream file which contains this group.
//The file generated by make_file can be retrieved by the constructor
//MGGroup(const char* file, int error);
//Function's return value is:
//=0: the file is successfully made.
//=1: file could not be opened.
int MGGroup::make_file(const TCHAR* file){
	MGOfstream outfile;
	int error=outfile.open(file);
	if(error)
		return error;
	outfile<<(*this);
	return 0;
}

// Output virtual function.
std::ostream& MGGroup::toString(std::ostream& ostrm) const{
	ostrm<<"MGGroup="<<this<<", number of gels = "<<size()<<std::endl;
	const_iterator i=begin(), ie=end();	
	for(int j=0; i!=ie; i++, j++){
		ostrm<<"gel"<<j<<":"<<(**i)<<std::endl;
	}
	return ostrm;
}

//Get the number of objects included in thie group.
int MGGroup::num_of_objects() const{
	int num=0;
	const_iterator i=begin(), ie=end();	
	for(; i!=ie; i++){
		const MGGel* geli=i->get();
		const MGObject* obj=dynamic_cast<const MGObject*>(geli);
		if(obj){ num++; continue;}
		const MGGroup* grp=dynamic_cast<const MGGroup*>(geli);
		if(grp){num+=grp->num_of_objects();}
	}
	return num;
}

void MGGroup::push_appearance(MGAppearance* appr){
	if(size()){
		MGGel* gel0=front().get();
		MGContext* ctx=dynamic_cast<MGContext*>(gel0);
		if(ctx){
			ctx->set_appearance(appr);
			return;
		}
		MGAppearance* app2=dynamic_cast<MGAppearance*>(gel0);
		if(app2)
			pop_front();
	}
	push_front(MYELM(appr));
}

void MGGroup::push_context(MGContext* cntx){
	if(size()){
		MYELM& gel0=front();
		MGContext* ctx2=dynamic_cast<MGContext*>(gel0.get());
		if(ctx2){
			pop_front();
		}else{
			MGAppearance* app2=dynamic_cast<MGAppearance*>(gel0.get());
			if(app2)
				erase(begin());
		}
	}
	push_front(MYELM(cntx));
}

// push element x at the end.
void MGGroup::append(MGGel* x){
	MGAppearance* app=dynamic_cast<MGAppearance*>(x);
	if(app){
		push_appearance(app);
		return;
	}
	MGContext* ctx=dynamic_cast<MGContext*>(x);
	if(ctx){
		push_context(ctx);
		return;
	}
	push_back(std::unique_ptr<MGGel>(x));
}

// push element x at the first.
void MGGroup::prepend(MGGel* x){
	MGAppearance* app=dynamic_cast<MGAppearance*>(x);
	if(app){
		push_appearance(app);
		return;
	}
	MGContext* ctx=dynamic_cast<MGContext*>(x);
	if(ctx){
		push_context(ctx);
		return;
	}

	if(size()){
		iterator first=begin();
		app=dynamic_cast<MGAppearance*>(first->get());
		if(app){
			insert(++first,MYELM(x));
			return;
		}
		ctx=dynamic_cast<MGContext*>(first->get());
		if(ctx){
			insert(++first,MYELM(x));
			return;
		}
	}
	push_front(MYELM(x));
};

//Remove the MGAppearance of this MGAttribedGel.
std::unique_ptr<MGAppearance> MGGroup::remove_appearance(){
	std::unique_ptr<MGAppearance> apr;
	if(size()){
		MYELM& gel0=front();
		MGAppearance* app=dynamic_cast<MGAppearance*>(gel0.get());
		if(app){
			gel0.release(); pop_front();
			apr.reset(app);
		}else{
			MGContext* cnt=dynamic_cast<MGContext*>(gel0.get());
			if(cnt)
				apr=cnt->remove_appearance();
		}
	}
	return apr;
}

//Read all member data.
void MGGroup::ReadMembers(MGIfstream& buf){
	int n;
	buf>>n;
	for(int i = 0; i<n; i++){
		UniqueGel gel(buf.ReadPointer());
		push_back(std::move(gel));
	}
}
//Write all member data
void MGGroup::WriteMembers(MGOfstream& buf)const{
	int nPlaneImage=0;
	const_iterator i=begin(), ie=end();	
	for(; i!=ie; i++){
		if((**i).identify_type()==MGPLANEIMAGE_TID)
			nPlaneImage++;
	}

	int n=int(size())-nPlaneImage;
	buf<<n;
	i=begin();	
	for(; i!=ie; i++){
		if((**i).identify_type()==MGPLANEIMAGE_TID)
			continue;
		buf.WritePointer(i->get());
	}
}

//Transform the gel by the argument.
void MGGroup::transform(const MGVector& v)//translation
{
	for(auto& geli:m_list){
		auto obj = dynamic_cast<MGObject*>(geli.get());
		if(obj)
			obj->transform(v);
	}
}
void MGGroup::transform(double scale)//scaling.
{
	for(auto& geli:m_list){
		auto obj = dynamic_cast<MGObject*>(geli.get());
		if(obj)
			obj->transform(scale);
	}
}
void MGGroup::transform(const MGMatrix& mat)//matrix transformation.
{
	for(auto& geli:m_list){
		auto obj = dynamic_cast<MGObject*>(geli.get());
		if(obj)
			obj->transform(mat);
	}
}
void MGGroup::transform(const MGTransf& tr)//general transformation.
{
	for(auto& geli:m_list){
		auto obj = dynamic_cast<MGObject*>(geli.get());
		if(obj)
			obj->transform(tr);
	}

}

/// グループ最下層のデータを取得する
void MGGroup::getSmallGroupData(std::vector<MGObject*>& outputObjects){
	const_iterator i=begin(), ie=end();	
	for(; i!=ie; i++){
		MGGel* geli=i->get();
		MGGroup* grp=dynamic_cast<MGGroup*>(geli);
		if (!grp){
			continue;
		}

		bool haveChildGroup = false;
		const_iterator i2=grp->begin(), ie2=grp->end();
		// 子グループのループ
		for(; i2!=ie2; i2++){
			MGGel* geli2=i2->get();
			MGGroup* childGrp=dynamic_cast<MGGroup*>(geli2);

			if(childGrp){
				haveChildGroup = true;
			}else{
				MGObject* obj=dynamic_cast<MGObject*>(geli2);
				if(obj){
					outputObjects.push_back(obj);
				}
			}
		}

		if(haveChildGroup){
			// 子グループが存在する場合は、再帰的に処理する
			grp->getSmallGroupData(outputObjects);
		}	
	}
}

/// グループに含まれるすべてのobjectを取得する
void MGGroup::getAllObjects(std::vector<MGObject*>& allObjects){
	const_iterator i=begin(), ie=end();	
	for(; i!=ie; i++){
		MGObject* obj=dynamic_cast<MGObject*>(i->get());
		if (obj){
			allObjects.push_back(obj);
		}
	}
}

AUTO_GEL_REGISTER(MGGroup, MGGROUP_TID);
