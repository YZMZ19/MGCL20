/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#ifndef _MGKnot_HH_
#define _MGKnot_HH_
/** @file */
/** @addtogroup BASE
 *  @{
 */
#include "mg/MGCL.h"

// MGKnot.h
//

class MGIfstream;
class MGOfstream;

///Defines knot value and its multiplicity.
class MG_DLL_DECLR MGKnot{

public:

///String stream Function.
MG_DLL_DECLR friend std::ostream& operator<< (std::ostream&, const MGKnot& );

////////////Constructor////////////
MGKnot():m_multiplicity(0){;};	///Default Constructor

///*****This is the fundamental constructor.*****
MGKnot(
	double t,			///<Knot value
	int multiplicity	///<Multiplicity
);

////////////Operator overload.////////////

///Cast operator to double.
operator double()const{ return m_value; };

///Comparison operator.
bool operator< (const MGKnot& kt2)const{return m_value<kt2.m_value;};
bool operator> (const MGKnot& kt2)const{return kt2<(*this);};
bool operator<= (const MGKnot& kt2)const{return !(kt2<(*this));};
bool operator>= (const MGKnot& kt2)const{return !((*this)<kt2);};
bool operator== (const MGKnot& kt2)const;
bool operator!= (const MGKnot& kt2)const{return !operator==(kt2);};

////////////Member Function////////////

///Get the knot value.
double value() const {return m_value;}

///Get the multiplicity.
int multiplicity() const {return m_multiplicity;}

///Dump Functions.
///Calculate dump size
int dump_size() const;

///Dump Function
int dump(MGOfstream& ) const;

///Restore Function
int restore(MGIfstream& );

private:
	double	m_value;			///<Knot value.
	int		m_multiplicity;		///<multiplicity of the value.

};

/** @} */ // end of BASE group
#endif
