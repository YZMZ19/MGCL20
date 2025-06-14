/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#ifndef _MGLCisect_vector_HH_
#define _MGLCisect_vector_HH_
/** @file */
/** @addtogroup IsectContainer
 *  @{
 */

#include <vector>
#include "topo/LCisect.h"

//Forward class declaration.
class MGCurve;
class MGLoop;
class MGInterval;

/// MGLCisect_vector defines linked list of MGLCisect.

/// Used to represent Intersection points of a loop and a curve.
class MG_DLL_DECLR MGLCisect_vector{

public:

using MYELM = MGLCisect;
using MYVEC = std::vector<MYELM>;

typedef MYVEC::reference              reference;
typedef MYVEC::const_reference        const_reference;
typedef MYVEC::iterator               iterator;
typedef MYVEC::const_iterator         const_iterator;
typedef MYVEC::size_type              size_type;
typedef MYVEC::difference_type        difference_type;
typedef MYVEC::value_type             value_type;
typedef MYVEC::reverse_iterator       reverse_iterator;
typedef MYVEC::const_reverse_iterator const_reverse_iterator;

MYVEC m_vec;

/////MYVEC's member function/////
//We cannot use inheritance of std::vector<MYELM> to make DLL.
iterator begin() { return m_vec.begin(); };
iterator end() { return m_vec.end(); };
const_iterator begin()const { return m_vec.begin(); };
const_iterator end()const { return m_vec.end(); };
reverse_iterator rbegin() noexcept { return m_vec.rbegin(); };
reverse_iterator rend() noexcept { return m_vec.rend(); };
const_reverse_iterator rbegin() const noexcept { return m_vec.rbegin(); };
const_reverse_iterator rend() const noexcept { return m_vec.rend(); };
void clear() { m_vec.clear(); };
bool empty() const { return m_vec.empty(); };
template <class... Args>
void emplace_back(Args&&... args) { m_vec.emplace_back(std::forward<Args>(args)...); };
iterator erase(iterator x) { return m_vec.erase(x); };
iterator erase(iterator first, iterator last) { return m_vec.erase(first, last); };
const MYELM& front() const { return m_vec.front(); };
MYELM& front() { return m_vec.front(); };
const MYELM& back() const { return m_vec.back(); };
MYELM& back() { return m_vec.back(); };
iterator insert(const_iterator it, const MYELM& x) { return m_vec.insert(it, x); };
iterator insert(iterator it, MYELM&& x) { return m_vec.insert(it, std::move(x)); };
template <class InputIterator>
iterator insert(const_iterator position,
	InputIterator first, InputIterator last) {
	return m_vec.insert(position, first, last);
};
void pop_back() { m_vec.pop_back(); };
void push_back(const MYELM& x) { m_vec.push_back(x); };
void push_back(MYELM&& x) { m_vec.push_back(std::move(x)); };
void reserve(size_type n) { m_vec.reserve(n); };
size_t size() const { return m_vec.size(); };
MYELM& operator[](size_t i) { return m_vec.operator[](i); };
const MYELM& operator[](size_t i)const { return m_vec.operator[](i); };

public:

///String stream Function
MG_DLL_DECLR friend std::ostream& operator<< (std::ostream&, const MGLCisect_vector& );

/// Constructor
MGLCisect_vector();
MGLCisect_vector(const MGLoop& loop);		///Loop

/// Adds the MGLCisect to the end of the list.
void append(const MGLCisect& lcisect);

void append(
	const MGLEPoint& lp,		///<loop's parameter with edge id.
	double t,				///<Curve's parameter value.
	const MGPosition& uv	///<Face's parameter value(u,v) data.
);

/// Adds the MGLCisect_vector to the end of the list.
void append(const MGLCisect_vector& list);

/// Return the number of items that are in the list.
int entries() const{return (int)size();};

/// Return  the first element in the list.
/// If list is empty, behavior is undefined.
const MGLCisect& first() const{return front();};

///Insert MGLCisect at the position i.
void insertAt(iterator i, const MGLCisect& lcisect){insert(i, lcisect);};

/// Return(but does not remove) last element in the list.
/// If list is empty, behavior is undefined.
const MGLCisect& last() const{return back();};

///Return the pointer to loop.
const MGLoop* loop() const {return m_loop;};

///Update MGLEPoint in this LCisect_vector.
///This is to update MGLEPoints in this obtained before MGLoop::make_vertex
///and the loop is updated by MGLoop::make_vertex.
///Generally speaking, when make_vertex is invoked after MGLCisect_vector is obtaind,
///the MGLEPoints in the MGLCisect_vector do not contain correct values, since
///a new edge is inserted int the MGComplex's cell vector.
void update_lepoint(
	const MGLEPoint& lep
);

private:
	const MGLoop* m_loop;
	double m_error_square;		///<Error square allowed to compute isect and 
								///<end point coincidence.

};

/** @} */ // end of IsectContainer group
#endif
