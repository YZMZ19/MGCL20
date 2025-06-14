/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#ifndef _MGGroup_HH_
#define _MGGroup_HH_

#include "mg/AttribedGel.h"

//
//Define MGGroup Class.
//
class MGIfstream;
class MGOfstream;
class MGBox;
class MGObject;
class MGName;
class MGAttrib;
class MGGLAttrib;
class MGAppearance;
class MGContext;
class MGVector;
class MGMatrix;
class MGTransf;
class mgVBO;

/** @addtogroup GelRelated
 *  @{
 */

///MGGroup is a class which constains MGGel elements.

///MGGroup provides functions:
///(1) a container of MGGel(MGAtrribedGel(MGObject or MGGroup), MGAttrib) elements.
///(2) To attach an appearance to the context of the group.
///NOTE: Do NOT use push_back() or push_front() for MGAppearance and MGContex.
///Instead use append() or push_appearance, push_context.
class MG_DLL_DECLR MGGroup: public MGAttribedGel{

public:
///	別名定義

	using MYELM=std::unique_ptr<MGGel>;
	using MYLIST=std::list<MYELM>;

	typedef MYLIST::reference              reference;
	typedef MYLIST::const_reference        const_reference;
	typedef MYLIST::iterator               iterator;
	typedef MYLIST::const_iterator         const_iterator;
	typedef MYLIST::size_type              size_type;
	typedef MYLIST::difference_type        difference_type;
	typedef MYLIST::value_type             value_type;
	typedef MYLIST::reverse_iterator       reverse_iterator;
	typedef MYLIST::const_reverse_iterator const_reverse_iterator;

	MYLIST m_list;

	/////MYLIST's member function/////
	//We cannot use inheritance of std::list<MYELM> to make DLL.
	iterator begin() { return m_list.begin(); };
	const_iterator begin()const { return m_list.begin(); };
	iterator end() { return m_list.end(); };
	const_iterator end()const { return m_list.end(); };
	reverse_iterator rbegin() noexcept { return m_list.rbegin(); };
	reverse_iterator rend() noexcept { return m_list.rend(); };
	const_reverse_iterator rbegin() const noexcept { return m_list.rbegin(); };
	const_reverse_iterator rend() const noexcept { return m_list.rend(); };
	void clear() { m_list.clear(); };
	bool empty() const { return m_list.empty(); };
	template <class... Args>
	void emplace_back(Args&&... args) { m_list.emplace_back(std::forward<Args>(args)...); };
	iterator erase(iterator x) { return m_list.erase(x); };
	iterator erase(iterator first, iterator last) { return m_list.erase(first, last); };
	const MYELM& front() const { return m_list.front(); };
	MYELM& front() { return m_list.front(); };
	const MYELM& back() const { return m_list.back(); };
	MYELM& back() { return m_list.back(); };
	iterator insert(iterator it, MYELM&& x) { return m_list.insert(it, std::move(x)); };
	template <class InputIterator>
	iterator insert(const_iterator it, InputIterator first, InputIterator last) {
		return m_list.insert(it, first, last);
	};
	void pop_back() { m_list.pop_back(); };
	void pop_front() { m_list.pop_front(); };
	void push_back(MYELM&& x) { m_list.push_back(std::move(x)); };//Only move operation is allowed.
	void push_front(MYELM&& x) { m_list.push_front(std::move(x)); };//Only move operation is allowed.
	size_t size() const { return m_list.size(); };
	void sort() { m_list.sort(); };

////////Special member functions/////////
MGGroup()=default;///void constructor.
virtual ~MGGroup()=default;
MGGroup(const MGGroup& rhs)=delete;
MGGroup& operator=(const MGGroup& rhs)=delete;
MGGroup(MGGroup&& rhs)=default;
MGGroup& operator=(MGGroup&& rhs)=default;

///Construct MGGroup from the file made by MGOfstream.
///error contains error flag as:
///=0: open succeeded.
///=1: file not found, or could not be opened.
///=2: file found, but, the format is not MGCL format.
///When error!=0, MGGroup contains no MGGel's.
MGGroup(const TCHAR* file, int& error);

///Comparison.
bool equal_test(const MGGel& gel2)const override;//gel2 must be MGGroup.
std::partial_ordering ordering_test(const MGGel& gel2)const override;//gel2 must be MGGroup.
bool operator==(const MGGroup& g2)const;
std::partial_ordering operator<=>(const MGGroup& gel2)const;

////////////Member Function////////////

///Get the MGAppearance pointer in this group. If not defined, null will be
///returned.
MGAppearance* appearance();
const MGAppearance* appearance()const;

/// push element x at the end.
///The x's ownership is transfered to this MGGroup.
void append(MGGel* x);
void append(MYELM&& x){ append(x.release()); };

/// push element x at the first.
///x must be a newed object, and the ownership will be transfered to thisMGGroup.
void prepend(MGGel* x);

///The x's ownership is transfered to this MGGroup.
void push_appearance(MGAppearance* appr);
void push_context(MGContext* cntx);

///Get the box of the group.
///If no objects were included in this group, null box will be returned.
MGBox box()const;

///Generate copied gel of this gel.
///Returned is a newed object. User must delete the object.
virtual MGGroup* clone()const override;

///Get the MGContext pointer stored in this group. If not defined, null will be
///returned.
MGContext* context();
const MGContext* context()const;

///////display member function.
virtual void display_arrows(mgSysGL& sgl)const;
virtual void display_break_points(mgSysGL& sgl)const;
virtual void display_control_polygon(mgSysGL& sgl)const;
virtual void display_curvatures(
	mgSysGL& sgl,///<Sgl to make pictures in.
	int		density,///<densitiy of the graph.
	bool	use_radius,///<true:radius display, false:curvature display.
	double	scale=1.	///<scaling of the graph.
)const;

///make this group has appearance and get the MGAppearance pointer.
MGAppearance* ensure_appearance();

///Find the position of the gel in the gel list and the group pointer
///which includes the gel. Searching will be done into the member group gel
///of this list.
const_iterator find(
	const MGGel* gel,///<The target gel.
	const MGGroup*& grp	///<The group that includes gel will be returned.
						///<When not found grp=null will be returned
)const;

///Find the position of the gel in the gel list and the group pointer
///which includes the gel. Searching will be done into the member group gel
///of this list.
iterator find(
	MGGel* gel,///<The target gel.
	MGGroup*& grp///<The group that includes gel will be returned,
			///<When not found grp=null will be returned.
);

///Count up how many MGGroup members are included in this group.
///Function's return value is the number of member group.
///Only the members of this group are counted.
int getMemberGroupNumber()const;

///Get i-th MGGroup pointer of this memeber group.
///If i-th MGGroup except MGAppearance is not found, null is returned.
const MGGroup* get_i_th_MemberGroup(int i)const;
MGGroup* get_i_th_MemberGroup(int i);

///Find the position of the 1st member gel of type tid in this MGGroup.
///**********Searching is NOT performed into the member group.
iterator search_by_id(MGGEL_TID tid);
const_iterator search_by_id(MGGEL_TID tid)const;

///Find the 1st gel(MGGroup or MGObject) member in this MGGroup
///that has the given name.
///**********Searching IS performed into the member group recursively.
void nameSearch(
	const MGName& name,///<The target name.
	MGGroup*& groupIncluding,///<If found the group including the gel of the name is set,
					///else, null will be returned.
	MGAttribedGel*& gelFound///<found gel will be returned if found, else null.
);

///Find the 1st gel(MGGroup or MGObject) member in this MGGroup
///that has the given name.
///**********Searching IS performed into the member group recursively.
void nameSearch(
	const MGName& name,///<The target name.
	const MGGroup*& groupIncluding,///<If found, the group including the gel of the name is set,
					///else, null will be returned.
	const MGAttribedGel*& gelFound///<found gel will be returned if found, else null.
)const;

/// Return This object's typeID
virtual long identify_type() const{return MGGROUP_TID;};

///Test if this gel includes an object.
///Function's return value is the 1st object found in the gel list of this
///if this includes an object. Otherwise null will be returned.
const MGObject* includes_object()const;
MGObject* includes_object();

//Load MGGroup data from istrm, which is unloaded by MGGroup::unload.
void load(MGIfstream& istrm);

//Unload MGGroup data into ostrm.
void unload(MGOfstream& ostrm)const;

///Make an MGOfstream file which contains this group.
///The file generated by make_file can be retrieved by the constructor
///MGGroup(const char* file, int error);
///Function's return value is:
///=0: the file is successfully made.
///=1: file could not be opened.
int make_file(const TCHAR* file);

///Make a display list of this gel.
void make_display_list(MGCL::VIEWMODE vmode=MGCL::DONTCARE)const;

///Get manifold dimension.
///MGGroup returns right one, MGGroup return 2, and others return -1.
virtual int manifold_dimension() const{return 2;};

///Delete the mgVBO of the i-th element.
void delete_displayList(const_iterator x)const;

///Delete display list of the sequence [first, last).
void delete_displayList(
	const_iterator first, const_iterator last
)const;

///Delete the mgVBO of gels_to_delete
void delete_displayList(
	const std::vector<const MGGel*>& gels_to_delete
)const;

///Process of draw or render attributes.
void drawAttrib(
	mgVBO& vbo,///<The target graphic object.
	bool no_color=false	///if true, color attribute will be neglected.
)const;

///Get the number of objects included in thie group.
int num_of_objects() const;

///IGES output function.
virtual int out_to_IGES(
	MGIgesOfstream& igesfile,
	int SubordinateEntitySwitch=0
)const;

/// Output virtual function.
virtual std::ostream& toString(std::ostream&) const;

///Remove the MGAppearance of this MGAttribedGel.
std::unique_ptr<MGAppearance> remove_appearance() override;

///Execute polar-scaling to all the MGCurve and MGFace of this group.
///curve's (x,y) are updated. No other coordinates are unchanged.
///The updated result curve is always MGLBRep.
///For MGFace, the boundaries are polar-scaled.
///
///Rotation is performed from the angle range (angleBase,angle1) to
///(angleBase,angle2).
///That is, when angle1=angle2, no change is done.
///When angle2 is angleBase, all the data will lie on the straight of from origin to
///(cos(angleBase), sin(angleBase)).
///angle1-angleBase must be >MGTolerance::angle_zero().
///IF a member gel is not MGCurve nor MGFace, it is unchanged.
void scalePolar(
	double angleBase,	///<base angle.
	double angle1,		
	double angle2
);

//set the copy of appr2 to this MGAttribedgel.
void set_appearance(const MGAppearance& appr2)override;
void set_appearance(MGAppearance* appr2)override;//ownership transfer mode.

/// グループ最下層のデータを取得する??????
void getSmallGroupData(
	std::vector<MGObject*>& outputObjects	///グループ最下層のデータを入れる
);

/// グループに含まれるすべてのobjectを取得する
void getAllObjects(
	std::vector<MGObject*>& allObjects	///オブジェクトデータ
);

///Transform the gel by the argument.

///translation
virtual void transform(const MGVector& v);

///scaling.
virtual void transform(double scale);

///matrix transformation.
virtual void transform(const MGMatrix& mat);

///general transformation.
virtual void transform(const MGTransf& tr);

///Get the name of the class.
virtual std::string whoami()const{return "Group";};

///Read all member data.
virtual void ReadMembers(MGIfstream& buf);

///Write all member data
virtual void WriteMembers(MGOfstream& buf)const;

private:

///Find the 1st gel(MGGroup or MGObject) member in this MGGroup
///that has the given name.
///**********Searching IS performed into the member group recursively.
void search_by_nameSub(
	const MGName& name,
	MGGroup*& groupIncluding,///<If found, the group including the gel of the name is set,
					///else, null will be returned.
	MGAttribedGel*& gelFound///<found gel will be returned if found, else null.
);


friend class MGIfstream;
friend class MGOfstream;
};

/** @} */ // end of GelRelated group
#endif
