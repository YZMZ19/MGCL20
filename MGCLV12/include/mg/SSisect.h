/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#ifndef _MGSSisect_HH_
#define _MGSSisect_HH_
/** @addtogroup IsectContainer
 *  @{
 */

#include "mg/MGCL.h"
#include "mg/isect.h"
#include "mg/Curve.h"

// MGSSisect.h
// Header for MGSSisect

//Forward Declaration
class MGCurve;

///MGSSisect represents one intersection line of two surfaces.

///A list of intersection lines are expressed by MGSSisects.
///The behavior of MGSSisect is like a auto_ptr. Copy or assignment
///of MGSSisect means transfer of the ownership of all the included curve
///to copied or assigned MGSSisect and original MGSSisect does not have the
///ownership of the curves any more. User should be aware of it.
/// Surface と Surface の交線を一つのみ表現する。交線の集合は別に表現される。
class MG_DLL_DECLR MGSSisect:public MGisect{
private:	
	std::unique_ptr<MGCurve> m_iline;	///<(x,y,z)coordinate representaion of the line.
					///<交線の座標値よる表現
	std::unique_ptr<MGCurve> m_param1;	///<(u,v) representaion of the line of the first
					///<surface. 交線の第１surface のパラメータ(u,v)による表現
	std::unique_ptr<MGCurve> m_param2;	///<(u,v) representaion of the line of the second
					///<surface. 交線の第2surface のパラメータ(u,v)による表現
	MGSSRELATION m_rel;	///<Two surfaces relationship at the intersection line.
					///<交線における両Surfaceの関係。

public:

////////Special member functions/////////
MGSSisect():m_rel(MGSSRELATION::MGSSREL_UNKNOWN){;};	///void constructor.
~MGSSisect()=default;			///Destructor.
MGSSisect(const MGSSisect& ssi)=delete;///Copy constructor.
MGSSisect& operator= (const MGSSisect& rhs)=delete;///Copy assignment.
MGSSisect(MGSSisect&& ssi);		///Move constructor.
MGSSisect& operator= (MGSSisect&& rhs);///Move assignment.

///Construct providing all the raw data.
///The ownership of iline, param1, and param2 are all transfered to MGSSisect.
///All of these objects must be newed ones.
MGSSisect(
	MGCurve* iline,	///<Pointer of newed object.
	MGCurve* param1,///<Pointer of newed object.
	MGCurve* param2,///<Pointer of newed object.
	const MGSSRELATION r1= MGSSRELATION::MGSSREL_UNKNOWN///<Relation of the two.
):m_iline(iline), m_param1(param1), m_param2(param2), m_rel(r1){;};

///Construct providing all the raw data.
///Copy version. Copy of the three curves will take place.
MGSSisect(
	const MGCurve& iline,
	const MGCurve& param1,
	const MGCurve& param2,
	const MGSSRELATION r1= MGSSRELATION::MGSSREL_UNKNOWN
);

///Comparison operator.
bool operator< (const MGSSisect& ssi2)const;
bool operator> (const MGSSisect& ssi2)const{return ssi2<(*this);};
bool operator<= (const MGSSisect& ssi2)const{return !(ssi2<(*this));};
bool operator>= (const MGSSisect& ssi2)const{return !((*this)<ssi2);};
bool operator== (const MGSSisect& ssi2)const;
bool operator!= (const MGSSisect& ssi2)const{return !operator==(ssi2);};

///Ordering functions.
bool operator< (const MGisect& is)const;
bool operator< (const MGCCisect& is)const{return false;};
bool operator< (const MGCSisect& is)const{return false;};
//bool operator< (const MGCFisect& is)const{return false;};
//bool operator< (const MGFFisect& is)const{return true;};
bool operator== (const MGisect& is)const;

///Exchange 1st and 2nd order of the parameter line representation.
///When 1st object's manifold dimension is less than 2nd one's,
///this does nothing. Valid only when their manifold dimensions are equal.
void exchange12()override;

///Test if two ssi's world curve have common parts (in line_zero()).
///Fucntion's return value is 
///		1:have common part.
///		0>=:no common part(except a point).
int has_common(const MGSSisect& ssi2)const;

///Test if this SSI is null.
bool is_null()const{return m_iline.get()==nullptr;};

///Return the object of the intersection(world coordinates representation).
const MGObject& isect()const{return *m_iline;};

///Return the 1st object's parameter value curve of the intersection.
///*****This function is valid only when manifold_dimension()==1.
const MGCurve* isect1_param1()const{return m_param1.get();};

///Return the 2nd object's parameter value curve of the intersection.
///*****This function is valid only when manifold_dimension()==1.
const MGCurve* isect1_param2()const{return m_param2.get();};

///Return (x,y,z) coordinate representation intersection line.
MGCurve& line() const{return *m_iline;}

///Return the manifold dimension of the intersection.

///That is,
///0: when the intersection is a point,
///1: when                  is a curve,
///2: when                  is a surface.
int manifold_dimension()const{return 1;};

///negate the direction of the intersection line.
void negate();

/// Output function.
std::ostream& toString(std::ostream& ostrm)const;

///Return (u,v) parameter representation intersection line
///of the 1st surface.
MGCurve& param1() const{return *m_param1;}

///Return (u,v) parameter representation intersection line
///of the 2nd surface.
MGCurve& param2() const{return *m_param2;}

///Return the relationship at the intersection line.
MGSSRELATION rel() const{return m_rel;}

///Release each curve pointer from this.

///After the use of release_xxxx(), MGSSisect does not have the ownership of
///the each curve.
MGCurve* release_line(){return m_iline.release(); };
MGCurve* release_param1(){return m_param1.release(); };
MGCurve* release_param2(){return m_param2.release(); };

void set_null();
	friend class MGSSisects;
};

/** @} */ // end of IsectContainer group
#endif
