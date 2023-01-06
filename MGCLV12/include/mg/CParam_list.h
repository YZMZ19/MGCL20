/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#ifndef _MGCParam_list_HH_
#define _MGCParam_list_HH_
/** @addtogroup IsectContainer
 *  @{
 */
#include <list>
#include <algorithm>
#include "mg/MGCL.h"

//Forward class declaration.
class MGCurve;

/// MGParam_Vector provides a list to store parameters of a curve.
class MG_DLL_DECLR MGCParam_list{
//We cannot use inheritance of std::list<MYELM> to make DLL.

public:
	using MYELM = double;
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
	template <class... Args>
	void emplace_back(Args&&... args) { m_list.emplace_back(std::forward<Args>(args)...); };
	bool empty() const { return m_list.empty(); };
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

private:
const MGCurve *m_curve;	///< Curve.
double m_error;			///<Error in parameter space of the curve.

public:

///String stream Function
MG_DLL_DECLR friend std::ostream& operator << (std::ostream&, const MGCParam_list& );

////////Special member functions/////////
explicit MGCParam_list(const MGCurve *curve=0);
virtual ~MGCParam_list()=default;
MGCParam_list(const MGCParam_list& rhs)=default;
MGCParam_list& operator=(const MGCParam_list& rhs)=default;
MGCParam_list(MGCParam_list&& rhs)=default;
MGCParam_list& operator=(MGCParam_list&& rhs)=default;

/// Adds the parameter to the end of the list.
void append(double param);

/// Adds the parameter list to the end of the list.
void append(const MGCParam_list& lst);

///Returns the pointer to the curve.
const MGCurve* curve() const {return m_curve;}

/// Returns the number of items that are in the list.
int entries() const{return int(size());};

/// Returns(but does not remove) first element in the list.
/// If list is empty, behavior is undefined.
double first() const{return front();};

///Return true (1) if there are no items in the list,
/// false(0) otherwise.
bool isEmpty() const{return empty();};

/// Return(but does not remove) last element in the list.
/// If list is empty, behavior is undefined.
double last() const{return back();};

///Get Lower Bound(this must be sorted)
///if no lower bound return false
bool lower_bound(double param, double& lowerBound)const{
	const_iterator citer = std::lower_bound(begin(), end(), param);
	if(citer == end())return false;
	lowerBound = *citer;
	return true;
}

/// Adds the parameter to the beginning of the list.
void prepend(double param){push_front(param);};

///Remove the parameter and return the parameter. If i is not valid, 
/// behavior is undefined.
double removeAt(iterator i);

///Remove the first parameter int the list and return the parameter.
///If i is not valid, behavior is undefined.
double removeFirst();

///Remove the last parameter in the list and return the parameter.
///If i is not valid, behavior is undefined.
double removeLast();
};

/** @} */ // end of IsectContainer group
#endif
