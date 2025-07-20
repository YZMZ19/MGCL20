/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#ifndef _MGGelPositions_HH_
#define _MGGelPositions_HH_

class MGGroup;
#include <algorithm>
#include <vector>
#include "mg/GelPosition.h"

class MGPickObjects;

/** @file */
/** @addtogroup GelRelated
 *  @{
 */

///MGGelPosition Container Class.

///MGGelPositions is a class which constains MGGelPosition elements as a vector.
class MG_DLL_DECLR MGGelPositions{
//We cannot use inheritance of std::vector<MGGelPosition> to make DLL.

public:

using MYVEC=std::vector<MGGelPosition>;
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

///// MGGelPositions's member function /////
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
const MGGelPosition& front() const { return m_vec.front(); };
MGGelPosition& front() { return m_vec.front(); };
const MGGelPosition& back() const { return m_vec.back(); };
MGGelPosition& back() { return m_vec.back(); };
iterator insert(const_iterator it, const MGGelPosition& x) { return m_vec.insert(it, x); };
iterator insert(iterator it, MGGelPosition&& x) { return m_vec.insert(it, std::move(x)); };
template <class InputIterator>
iterator insert(const_iterator position,
	InputIterator first, InputIterator last) {
	return m_vec.insert(position, first, last);
};
void pop_back() { m_vec.pop_back(); };
void push_back(const MGGelPosition& x) { m_vec.push_back(x); };
void push_back(MGGelPosition&& x) { m_vec.push_back(std::move(x)); };
void reserve(size_type n) { m_vec.reserve(n); };
size_t size() const { return m_vec.size(); };
MGGelPosition& operator[](size_t i) { return m_vec.operator[](i); };
const MGGelPosition& operator[](size_t i)const { return m_vec.operator[](i); };

///String output function.
MG_DLL_DECLR friend std::ostream& operator<< (std::ostream&, const MGGelPositions&);


MGGelPositions()=default;

///Construct MGGelPositions of a MGGelPosition.
///When gelp.is_null(), the gelp will not be set.
MGGelPositions(const MGGelPosition& gelp);

///Conversion constructor of MGPickObjects.
MGGelPositions(const MGPickObjects& gelp);

/////////////////operator overloaded////////////////

///Set operation.
MGGelPositions& operator+=(const MGGelPositions& gelps){append(gelps);return *this;};
MGGelPositions& operator+=(const MGGelPosition& gelp){push_back(gelp);return *this;};
MGGelPositions& operator-=(const MGGelPositions& gelps){remove(gelps);return *this;};
MGGelPositions& operator-=(const MGGelPosition& gelp){remove(gelp);return *this;};
MGGelPositions& operator-=(const MGAbstractGels& types){remove(types);return *this;};
MGGelPositions& operator&=(const MGGelPositions& gelps){reset_with_common(gelps);return *this;};

///append MGGelPosition to this if the target is MGObject* and pred is true
///for [first, last).
///&**first must be MGGel*.
template <class InputIterator, class Predicate>
void append_if(
	MGGroup* grp,
	InputIterator first, InputIterator last,
	Predicate pred,//Predicate that is invoked when MGObject* as pred(MGObject*).
	MGGelPosition* gelp=nullptr //input null when grp is top group.
	    //When gelp is input, grp is set in gelp as an lower hierarchy group of top group.
){
	MGGelPosition gelp0(grp);
	if(gelp)
		gelp->append_lower_gel(grp);
	else
		gelp=&gelp0;
	for(; first != last; ++first){
		MGGel* gel=&**first;
		MGGroup* grplow= dynamic_cast<MGGroup*>(gel);
		if(grplow){//When gel is a MGGroup.
			append_if(grplow, grplow->begin(), grplow->end(), pred, gelp);
		}else{
			MGObject* obj=dynamic_cast<MGObject*>(gel);
			if(obj && pred(obj)){
				gelp->set_leaf_object(obj);
				push_back(std::move(*gelp));
			}
		}
	}
}

///push elements in gelps at the end. All of the gel pointers are
///transfered to this. On return, gelps will have no gel pointer in it.
void append(const MGGelPositions& gelps){
	insert(end(), gelps.begin(), gelps.end());
}

///Find the input MGGelposition.
iterator find(const MGGelPosition& gelp){return std::find(begin(),end(),gelp);};
const_iterator find(const MGGelPosition& gelp)const{return std::find(begin(),end(),gelp);};

//Test if there is a MGPickObject that includes input objin
//in this MGPickObjects' member. If input objin is MGShell,
//and a member is_shell_face(), test is performed to the shell.
//Returns true if objin is included in this MGPickObjects.
iterator includes(const MGObject* objin);

///Remove gelp if found in this.
void remove(const MGGelPosition& gelp);

///Remove objects of type from this pickobjects.
void remove(const MGAbstractGels& types);

///Remove gelps from this pickobjects.
void remove(const MGGelPositions& gelps);

///replace this with the common objects of this and pobjs2.
void reset_with_common(const MGGelPositions& pobjs2);

///replace this with symmetric_differecne of this and pobj, that is;
///(1) remove the same MGPickObject from this and pobjss.
///(2) append the result pobjs2 to this.
///On return, pobjs2 will have null sequence.
void reset_with_symmetric_difference(MGGelPositions& pobjs2);

///replace the i-th elemnet to pobj.
void reset(int i, const MGGelPosition& pobj);

///Select objects of input type from this.
///Function's return value is pickobjects selected.
///This is unchanged.
MGGelPositions select(const MGAbstractGels& types)const;

///Select the 1st MGCurve from this.
///Function's return value is MGPickObject of MGCurve 1st encountered in this
///MGPickObject sequence. If this did not includes any MGCurve,
///null MGPickOjbect will be returned.
///This is unchanged.
MGGelPosition select_1st_curve()const;

///Select all the MGCurve from this.
///MGPickObject of MGCurve encountered in this MGPickObject sequence will be appended
///in curves.
///This is unchanged.
void select_curves(MGGelPositions& curves)const;

///Select the 1st MGFSurface from this.
///Function's return value is MGPickObject of MGFSurface 1st encountered in this
///MGPickObject sequence. If this did not includes any MGFSurface,
///null MGPickObject will be returned.
///This will be unchanged.
MGGelPosition select_1st_fsurface()const;

///Select all the MGFSurface from this.
///MGGelPositions of MGFSurface encountered in this MGPickObject sequence will be appended
///in surfaces.
///This will be unchanged.
void select_fsurfaces(MGGelPositions& surfaces)const;

///Test if this is symmetric to gels2.
///Symmetric means:
///(1) number of gels included is the same.
///(2) all of the gels are MGObject and they have the same manifold dimension.
bool symmetric(const MGGelPositions& gels2)const;

private:

};

/** @} */ // end of GelRelated group
#endif
