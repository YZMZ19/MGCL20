/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#ifndef _MGKnotArray_HH_
#define _MGKnotArray_HH_
/** @file */
/** @addtogroup BASE
 *  @{
 */

#include <vector>
#include "mg/Knot.h"

// MGKnotArray.h
//

//Forward Declaration
class MGIfstream;
class MGOfstream;

/// Defines Array of Knots.

///This is not knot vector for B-Rep.
///See MGKnotVector.
class MG_DLL_DECLR MGKnotArray{
//We cannot use inheritance of std::vector<MGKnot> to make DLL.

public:

using MYVEC = std::vector<MGKnot>;
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
const MGKnot& front() const { return m_vec.front(); };
MGKnot& front() { return m_vec.front(); };
const MGKnot& back() const { return m_vec.back(); };
MGKnot& back() { return m_vec.back(); };
iterator insert(const_iterator it, const MGKnot& x){return m_vec.insert(it, x); };
iterator insert(iterator it, MGKnot&& x) { return m_vec.insert(it, std::move(x)); };
template <class InputIterator>
iterator insert(const_iterator position,
	InputIterator first, InputIterator last) {
	return m_vec.insert(position, first, last);
};
void pop_back() { m_vec.pop_back(); };
void push_back(const MGKnot& x) { m_vec.push_back(x); };
void push_back(MGKnot&& x) { m_vec.push_back(std::move(x)); };
void reserve(size_type n) { m_vec.reserve(n); };
size_t size() const { return m_vec.size(); };
MGKnot& operator[](size_t i) { return m_vec.operator[](i); };
const MGKnot& operator[](size_t i)const { return m_vec.operator[](i); };

///String stream Function.
MG_DLL_DECLR friend std::ostream& operator<< (std::ostream&, const MGKnotArray&);

////////////Constructor////////////

/// Dummy constructor.
MGKnotArray(){;};

///From Knot.
explicit MGKnotArray(const MGKnot&);

///From knot and the multiplicity.
MGKnotArray(double knot, int mult);

///Dump Functions.

	///Calculate dump size
	int dump_size() const;

	///Dump Function
	int dump(MGOfstream& ) const;

	///Restore Function
	int restore(MGIfstream& );

private:

};

/** @} */ // end of BASE group
#endif
