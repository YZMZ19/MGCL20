/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#pragma once 

#ifndef _MGPickObjects_HH_
#define _MGPickObjects_HH_

#include <vector>
#include "mg/AbstractGels.h"
#include "mg/PickObject.h"

class MGCurve;
class MGFSurface;
class MGOfstream;
class MGIfstream;
class MGGelPositions;
class MGPickObject;

/////////////////////////////////////////////////

/** @addtogroup MGObjectRelated
 *  @{
 */

///a container class for MGPickObject.
class MG_DLL_DECLR MGPickObjects{

public:
	typedef std::vector<UniquePickObject> MYVEC;

	/// types definitions.
	typedef MYVEC::reference              reference;
	typedef MYVEC::const_reference        const_reference;
	typedef MYVEC::iterator               iterator;
	typedef MYVEC::const_iterator         const_iterator;
	typedef MYVEC::size_type              size_type;
	typedef MYVEC::reverse_iterator       reverse_iterator;
	typedef MYVEC::const_reverse_iterator const_reverse_iterator;

/// Constructors

////////Special member functions/////////
MGPickObjects(void) = default;///Void Constructor. 
virtual ~MGPickObjects() = default;			///Destructor.
MGPickObjects(const MGPickObjects& rhs);///Copy constructor.
MGPickObjects& operator= (const MGPickObjects& rhs);///Copy assignment.
MGPickObjects(MGPickObjects&& rhs) = default;		///Move constructor.
MGPickObjects& operator= (MGPickObjects&& rhs) = default;///Move assignment.

///Construct MGPickObjects of one pobj.
explicit MGPickObjects(const MGPickObject& pobj);

///Index Operator overload.
const MGPickObject& operator[](size_t i)const{return *m_PickObjects[i];};
MGPickObject& operator[](size_t i){return *m_PickObjects[i];};

///Set operation.
MGPickObjects& operator+=(const MGPickObjects& gelps){push_back(gelps);return *this;};
MGPickObjects& operator+=(const MGPickObject& gelp);
MGPickObjects& operator-=(const MGPickObjects& gelps){remove(gelps);return *this;};
MGPickObjects& operator-=(const MGPickObject& gelp){remove(gelp);return *this;};
MGPickObjects& operator-=(const MGAbstractGels& types){remove(types);return *this;};
MGPickObjects& operator&=(const MGPickObjects& gelps){reset_with_common(gelps);return *this;};

///append the current objects(MGGelPositions).
void append_object(const MGGelPositions& gelps);

///Replace this sequence with [first,last).
void assign(const_iterator first, const_iterator last);

const MGPickObject& front()const{return *m_PickObjects.front();};
MGPickObject& front(){return *m_PickObjects.front();};

const MGPickObject& back()const{return *m_PickObjects.back();};
MGPickObject& back(){return *m_PickObjects.back();};

iterator begin(){return m_PickObjects.begin();};
const_iterator begin()const{return m_PickObjects.begin();};

iterator end(){return m_PickObjects.end();};
const_iterator end()const{return m_PickObjects.end();};

reverse_iterator rbegin(){return m_PickObjects.rbegin();};
const_reverse_iterator rbegin()const{return m_PickObjects.rbegin();};

reverse_iterator rend(){return m_PickObjects.rend();};
const_reverse_iterator rend()const{return m_PickObjects.rend();};

void clear(){m_PickObjects.clear();};
bool empty()const{return m_PickObjects.empty();};

void pop_back(){m_PickObjects.pop_back();};

///find the same pobj in this objects.
iterator find(const MGPickObject& pobj);
const_iterator find(const MGPickObject& pobj)const;

//Get an 1st object to tessellate from this pick objects.
//Function's return value is the object pointer(MGFSurface).
const MGFSurface* get_object_to_tessellate()const;

//Test if there is a MGPickObject that includes input objin
//in this MGPickObjects' member. If input objin is MGShell,
//and a member is_shell_face(), test is performed to the shell.
//Returns is the iterator of found obj.
iterator includes(const MGObject* objin);

//Test if input MGPickObject is includes in this MGPickObjects' member.
//If input objin is MGShell, and a member is_shell_face(),
//test is performed to the shell.
//Returns is the iterator of found obj.
iterator includes(const MGPickObject& pobj);

/// erase sequence [first, last).
void erase(iterator first, iterator last);

/// erase sequence i.
iterator erase(iterator i);

/// erase i-th element.
void erase(int i){erase(begin()+i);};

/// erase after the elments after the front().
///Resutl has length 1 sequence.
void erase_except_front();

///add one pobj.
void push_back(const MGPickObject& pobj);
void push_back(const MGPickObjects& pobjs);
void push_back(MGPickObjects&& pobjs);
void emplace_back(const MGGelPosition& gelp);
void emplace_back(MGGroup* group, MGObject* obj = nullptr);

///Remove pobj if found in this.
void remove(const MGPickObject& pobj);
void remove(const MGPickObjects& pobjs);

///Remove objects of types from this pickobjects.
void remove(const MGAbstractGels& types);

///Remove gelps from this pickobjects.
void remove(const MGGelPositions& gelps);

///replace this with the common objects of this and pobjs2.
void reset_with_common(const MGPickObjects& pobjs2);

///replace this with symmetric_differecne of this and pobj, that is;
///(1) remove the same MGPickObject from this and pobjs2.
///(2) append the result pobjs2 to this.
void reset_with_symmetric_difference(const MGPickObjects& pobjs2);

///reserve the size n, which are all null.
void reserve(size_t n);

//Erase all the objects that is not of specified types from this.
void reset_objects(const MGAbstractGels& types);

///resize the length of the sequence.
void resize(size_t n){m_PickObjects.resize(n);};

///resize the length of the sequence.
void reset(size_t i, const MGPickObject& pobj);

///Select objects of input type from this.
///Function's return value is pickobjects selected.
///This will be unchanged.
MGPickObjects select(const MGAbstractGels& types)const;

///Select the 1st MGCurve from this.
///Function's return value is MGPickObject of MGCurve 1st encountered in this
///MGPickObject sequence. If this did not includes any MGCurve,
///null MGPickOjbect will be returned.
///This will be unchanged.
MGPickObject select_1st_curve()const;

///Select all the MGCurve from this.
///MGPickObject of MGCurve encountered in this MGPickObject sequence will be appended
///in curves.
///This will be unchanged.
void select_curves(MGPickObjects& curves)const;

///Select the 1st MGFSurface from this.
///Function's return value is MGPickObject of MGFSurface 1st encountered in this
///MGPickObject sequence. If this did not includes any MGFSurface,
///null MGPickObject will be returned.
///This will be unchanged.
MGPickObject select_1st_fsurface()const;

///Select all the MGFSurface from this.
///MGPickObjects of MGFSurface encountered in this MGPickObject sequence will be appended
///in surfaces.
///This will be unchanged.
void select_fsurfaces(MGPickObjects& surfaces)const;

///Obtain the pobj number defined.
size_t size()const{return m_PickObjects.size();};

MYVEC& object_vector(){return m_PickObjects;};
const MYVEC& object_vector()const{return m_PickObjects;};

///Set no display for this vector of MGPickObject.
void setNoDisplay()const;

///Set no display for this vector of MGPickObject.
void setDisplay()const;

protected:
	MYVEC m_PickObjects;///vector of MGPickObject.

};

/** @} */ // end of MGObjectRelated group
#endif // _MGPickObjects_HH_
