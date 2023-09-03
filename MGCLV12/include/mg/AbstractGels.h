/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#ifndef _MGAbstractGells_HH_
#define _MGAbstractGells_HH_

#include <vector>
#include <iosfwd>
#include "mg/MGCL.h"
#include "mg/types.h"

/** @file */
/** @addtogroup GelRelated
 *  @{
 */

//
//Define MGAbstractGels Class.

///Is a container of MGAbstractGel, to specify what kind of gels are required.

///MGAbstractGels is a class which constains MGAbstractGel elements as a vector,
///provides OR conditions on the specification of gels.
class MG_DLL_DECLR MGAbstractGels{
//We cannot use inheritance of std::vector<MYELM> to make DLL.

	using MYELM = MGAbstractGel;
	using MYVEC= std::vector<MGAbstractGel> ;
public:

/// Alias definition.
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

///  MYVEC's member function definitions.
iterator begin() { return m_vec.begin(); };
const_iterator begin()const { return m_vec.begin(); };
iterator end() { return m_vec.end(); };
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
size_t size() const { return m_vec.size(); };
MYELM& operator[](size_t i) { return m_vec.operator[](i); };
const MYELM& operator[](size_t i)const { return m_vec.operator[](i); };

///String output function.
MG_DLL_DECLR friend std::ostream& operator<< (std::ostream&, const MGAbstractGels&);

///Construct MGAbstractGels of a MGAbstractGel. This is a conversion contructor.
MGAbstractGels(const MGAbstractGel& agell);

///push elements in agells at the end of this.
void push_back(const MGAbstractGels& agells){
	insert(end(), agells.begin(), agells.end());
}

};

/** @} */ // end of GelRelated group
#endif
