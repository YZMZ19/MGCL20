/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#ifndef _MGPosition_list_HH_
#define _MGPosition_list_HH_
/** @file */
/** @addtogroup IsectContainer
 *  @{
 */

#include <list>
#include "mg/Position.h"

//Forward class declaration.
class MGBox;
class MGCurve;
class MGFSurface;
class MGIfstream;
class MGOfstream;
class MGSSisect;

///MGPosition_list provides a list of Positions.

///MGPosition_list is used to store curve or surface parameter as
///two objects intersection data.
class MG_DLL_DECLR MGPosition_list{
//We cannot use inheritance of std::list<MYELM> to make DLL.

public:

	using MYELM = MGPosition;
	using MYLIST = std::list<MYELM>;

	typedef MYLIST::reference              reference;
	typedef MYLIST::const_reference        const_reference;
	typedef MYLIST::iterator               iterator;
	typedef MYLIST::const_iterator         const_iterator;
	typedef MYLIST::size_type              size_type;
	typedef MYLIST::difference_type        difference_type;
	typedef MYLIST::value_type             value_type;
	typedef MYLIST::allocator_type         allocator_type;
	typedef allocator_type::pointer       pointer;
	typedef allocator_type::const_pointer const_pointer;
	typedef MYLIST::reverse_iterator       reverse_iterator;
	typedef MYLIST::const_reverse_iterator const_reverse_iterator;

	MYLIST m_list;

	/////MYLIST's member function/////
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
	iterator insert(const_iterator it, const MYELM& x) { return m_list.insert(it, x); };
	iterator insert(iterator it, MYELM&& x) { return m_list.insert(it, std::move(x)); };
	template <class InputIterator>
	iterator insert(const_iterator it, InputIterator first, InputIterator last) {
		return m_list.insert(it, first, last);
	};
	void pop_back() { m_list.pop_back(); };
	void pop_front() { m_list.pop_front(); };
	void push_back(const MYELM& x) { m_list.push_back(x); };
	void push_front(const MYELM& x) { m_list.push_front(x); };
	void push_back(MYELM&& x) { m_list.push_back(std::move(x)); };
	void push_front(MYELM&& x) { m_list.push_front(std::move(x)); };
	size_t size() const { return m_list.size(); };
	void sort() { m_list.sort(); };

//////////// Constructor ////////////

/// Void constructor(Constructor of length 0).
MGPosition_list(){;};

/// Constructor of length 1.
MGPosition_list(const MGPosition& P): m_list(1,P){;};

///Construct MGPosition_list by replacing each element P of list by Q,
///where, Q=MGPosition(P.sdim(),P,start1,start2).
MGPosition_list(
	const MGPosition_list& list,///<Original MGPosition_list 
	int start1,	///<Start position to store of elements Q of
				///<new MGPosition_list.
	int start2///<Start position to retrieve of elements P of list.
);

///Debug Function
MG_DLL_DECLR friend std::ostream& operator << (std::ostream&, const MGPosition_list& );

//////////// Member Function. ////////////

/// Adds the parameter to the end of the list.
///Function's return value is true if appended, false if not.
bool append(const MGPosition& pos);

///Add parameter pair p of crv1 and 2 to the end of list if p is not
///included in the list already.
///Function's return value is true if appended, false if not.
bool append(
	const MGCurve& crv1,
	const MGCurve& crv2,
	const MGPosition& p
);

///Add parameter all the data in the list to the end of list if the member p is not
///included in the list already.
void append(
	const MGCurve& crv1,
	const MGCurve& crv2,
	const MGPosition_list& list
);

///Add parameter uv of surface to the end of list if uv is not
///included in the list already.
///Function's return value is true if appended, false if not.
bool append(
	const MGFSurface& srf,
	const MGPosition& uv
);

///Add parameter all the uv's of list to the end of this list if the
///uv is not included in this list already.
void append(
	const MGFSurface& srf,
	const MGPosition_list& list
);

///Add parameter pair tuv of crv and surface to the end of list if tuv is not
///included in the list already.
bool append(
	const MGCurve& crv,
	const MGFSurface& srf,
	const MGPosition& tuv
);

///Add parameter pair uvuv of surface srf1 and sur2 to the end of list
/// if uvuv is not included in the list already.
///Function's return value is true if appended, false if not.
bool append(
	const MGFSurface& srf1,
	const MGFSurface& srf2,
	const MGPosition& uvuv
){ return add(true,srf1,srf2,uvuv);}

///Add parameter pair uvuv of surface srf1 and sur2 to the end of list
void append(
	const MGFSurface& srf1,
	const MGFSurface& srf2,
	const MGPosition_list& list
);

/// Returns the number of items that are in the list.
int entries() const{return int(size());};

/// Returns(but does not remove) first element in the list.
/// If list is empty, behavior is undefined.
const MGPosition& first() const{return front();};
MGPosition& first(){return front();};

///Test if one of the points in the list is included in
///the box. If so, return id of the point. The id is of the first point
///found in the list.
///If no points was not included in the list, end() wil be returned in id.
///Function's return value is true if a point is in the box, and false
///if no points are in the box.
///n is the number of space dimension of points in this list to check.
///If n>0 is specified, the first n coordinates of the position data in
///the list are checked if included in the box.
///If n==0 is specified, all the coordinate data are checked.
bool in(const MGBox& box, const_iterator& id, int n=0) const;
bool in(const MGBox& box, iterator& id, int n=0);

///Inserts position at the index i.
///This index must be between zero and the number of items in the list,
/// or behavior is undefined.
void insertAt(iterator i, const MGPosition& pos){insert(i, pos);};

///Return true (1) if there are no items in the list,
/// false(0) otherwise.
bool isEmpty() const{return empty();};

///Add parameter pair uvuv of surface srf1 and sur2 to the beginning of list
///if uvuv is not included in the list already.
bool prepend(
	const MGFSurface& srf1,
	const MGFSurface& srf2,
	const MGPosition& uvuv
){ return add(false,srf1,srf2,uvuv);};

///Remove position in the list that is the same point as P.
///When distace of the two point is within error, they are regarded as same.
///Function's return value is num of points deleted.
int remove(
	double error,		///<square of error allowed to regard same.
	const MGPosition& P,///<Position to remove.
	int n=0			///<Number of space dimension to check,
						///<If n==0, all the coordinates are checked.
);

///Remove parameter pair uv of surface srf1 from the list
/// if uv is included in the list.
int remove(
	const MGFSurface& srf,
	const MGPosition& uv
);

///Remove parameter pair uvuv of surface srf1 and sur2 from the list
/// if uvuv is included in the list.
///uvuv's space dimension is 4. And uvuv(0) and uvuv(1) are srf1's
/// parameter (u,v), uvuv(2) and uvuv(3) are srf2's parameter (u,v).
///Function's return value is true if removal was done, and false if
///no same point was found and removal was not performed.
int remove(
	const MGFSurface& srf1,
	const MGFSurface& srf2,
	const MGPosition& uvuv
);

///Remove same elements in this list(parameter pair uvuv of surface srf1 and sur2)
///if elements are within the tolerance.
///Space dimension of the elements(uvuv) is 4. And uvuv(0) and uvuv(1) are srf1's
/// parameter (u,v), uvuv(2) and uvuv(3) are srf2's parameter (u,v).
///Function's return value is number of removed elements.
int remove_uvuv(
	const MGFSurface& srf1,
	const MGFSurface& srf2
);

///Remove the position and return the position. If i is not valid, 
/// behavior is undefined.
MGPosition removeAt(iterator i);

///Remove the first position in the list and return the position.
///If i is not valid, behavior is undefined.
MGPosition removeFirst();

///Remove the last position in the list and return the position.
///If i is not valid, behavior is undefined.
MGPosition removeLast();

///Remove uvuv that is on the ssi.
///Function's return value is the number of points deleted.
int removeOn(
	const MGFSurface& srf1,
	const MGFSurface& srf2,
	const MGSSisect& ssi
);

///reverse the ordering of the elements in the list.
void reverse_order();

///Sort the positions in the list according to the surface parameter space
///ordering. Positions (uvuv(id),uvuv(id+1)) in the list is a parameter
///position of a face.
void sort_uv_space(int id);

private:

///Dump Functions
int dump_size() const;

///Dump Function
int dump(MGOfstream& ) const;

///Restore Function
int restore(MGIfstream& );

///Add parameter pair uvuv of surface srf1 and sur2 to the end(append==true)
///or the beginning of list if uvuv is not included in the list already.
bool add(
	bool append,			
	const MGFSurface& srf1,
	const MGFSurface& srf2,
	const MGPosition& uvuv
);

};

/** @} */ // end of IsectContainer group
#endif
