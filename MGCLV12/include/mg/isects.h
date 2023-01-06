/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#ifndef _MGisects_HH_
#define _MGisects_HH_
/** @file */
/** @addtogroup IsectContainer
 *  @{
 */

#include <memory>
#include <list>
#include "mg/isect.h"

//Forward class declaration.
class MGisect;
class MGObject;
class MGCCisects;
class MGCSisects;
class MGSSisects;
class MGCFisects;
class MGHHisect;
class MGHHisects;

///MGisects defines a list of MGisect.

///MGisects is used to represent an array of intersection lines of 
///two objects.

///Intersections are obtained from two objects, which are known using
///the member functions object1() and object2().
///****NOTE****
///When two objects' manifold dimension are the same, object1 is this object
///at the invocation of MGObject::intersection(), and object2 is the argument
///object.
///However, their manifold dimension are not the same, object1 is always
///the lower dimension's object and object2 is the higer dimension's object.
class MG_DLL_DECLR MGisects{
public:

using MYELM = std::unique_ptr<MGisect>;
using MYLIST=std::list<MYELM>;

protected:
MYLIST m_list;

private:
const MGObject* m_object1;///< Object 1.
const MGObject* m_object2;///< Object 2.

public:

using iterator= MYLIST::iterator;
using const_iterator= MYLIST::const_iterator;
using reverse_iterator= MYLIST::reverse_iterator;
using const_reverse_iterator = MYLIST::const_reverse_iterator;
using reference= MYLIST::reference;
using const_reference= MYLIST::const_reference;
using size_type= MYLIST::size_type;

///String stream Function
MG_DLL_DECLR friend std::ostream& operator << (std::ostream& ostrm, const MGisects& is);

template <class T, class... Args>
void emplace_back(Args&&... args) {
	m_list.emplace_back(new T(std::forward<Args>(args)...));
};

template<class T>
std::unique_ptr<T> releaseFront(){
	auto isec = std::unique_ptr<T>(static_cast<T*>(front().release()));
	pop_front();
	return isec;
};
template<class T>
std::unique_ptr<T> releaseBack(){
	auto isec = std::unique_ptr<T>(static_cast<T*>(back().release()));
	pop_back();
	return isec;
};
template<class T>
std::unique_ptr<T> release(iterator i){
	auto isec = std::unique_ptr<T>(static_cast<T*>(i->release()));
	erase(i);
	return isec;
};

////////Special member functions/////////
~MGisects()=default;//Destructor.
MGisects(const MGisects& rhs)=delete;//Copy constructor.
MGisects& operator=(const MGisects& rhs)=delete;//Copy assignment.
MGisects(MGisects&& rhs)=default;//Move constructor.
MGisects& operator=(MGisects&& rhs)=default;//Move assignment.

///Constructor(of size 0)
MGisects(
	const MGObject* obj1=nullptr,
	const MGObject* obj2=nullptr
);

///Construct from MGCFisects.
//MGisects(const MGCFisects& cfis);

/////MYLIST's member function/////
//We cannot use inheritance of std::list<MYELM> to make DLL.

MGisects::iterator begin() { return m_list.begin(); }
MGisects::const_iterator begin()const { return m_list.begin(); }
MGisects::iterator end() { return m_list.end(); }
MGisects::const_iterator end()const { return m_list.end(); }
MGisects::reverse_iterator rbegin() noexcept { return m_list.rbegin(); };
MGisects::const_reverse_iterator rbegin() const noexcept { return m_list.rbegin(); }

void clear() { m_list.clear(); }
bool empty() const { return m_list.empty(); }
MGisects::iterator erase(iterator x) { return m_list.erase(x); }
MGisects::iterator erase(iterator first, iterator last) { return m_list.erase(first, last); }
const MGisects::MYELM& front() const { return m_list.front(); }
MGisects::MYELM& front() { return m_list.front(); }
const MGisects::MYELM& back() const { return m_list.back(); }
MGisects::MYELM& back() { return m_list.back(); }
MGisects::iterator insert(
	MGisects::iterator it, MYELM&& x) {
	return m_list.insert(it, std::move(x));
}
void pop_back() { m_list.pop_back(); }
void pop_front() { m_list.pop_front(); }
void push_back(MGisects::MYELM&& x) {
	m_list.push_back(std::move(x));
};//Only move operation is allowed.
size_t size() const { return m_list.size(); }

void sort() { m_list.sort(); };
const MGObject* object1()const{ return m_object1; };
const MGObject* object2()const{ return m_object2; };


///Exchange first and second order of MGisect.
void exchange12();

///Adds one MGisect* to the end of the list.
///Transfers the ownership of the isect to this list.
void emplace_back(MGisect* isect) {
	m_list.emplace_back(isect);
}

///append all the member of isects to the end of the list.
///Transfers the ownership of the isect in isects to this list.
void push_back(MGisects&& isects);

///Adds one MGisect* to the fron of the list.
///Transfers the ownership of the isect to this list.
void emplace_front(MGisect* isect){
	m_list.emplace_front(isect);
}

///append all the member of isects to the fron of the list.
///Transfers the ownership of the isect in isects to this list.
void push_front(MGisects&& isects);

/// Output virtual function.
virtual std::ostream& toString(std::ostream& ostrm)const;

};

///Cast the is to each type of isect,
///e.g. auto& csi=isectCast<MGCSisect>(i);
template<class T>
T& isectCast(MGisects::iterator is){ return *static_cast<T*>(is->get()); };

///Cast the is to each type of isect,
///e.g. auto& csi=isectCast<MGCSisect>(i);
template<class T>
const T& isectCast(MGisects::const_iterator is){ return *static_cast<const T*>(is->get()); };

/** @} */ // end of IsectContainer group
#endif