/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#ifndef _MGEReal_HH_
#define _MGEReal_HH_
/** @file */
/** @addtogroup BASE
 *  @{
 */
#include "mg/MGCL.h"

// MGEReal.h
//

class MGIfstream;
class MGOfstream;

/// MGEReal is extended real number to define infinity.

///Includes minus infinite and plus infinite, adding to ordinary real number.
class MG_DLL_DECLR MGEReal {

public:

MG_DLL_DECLR friend MGEReal operator+ (double, const MGEReal&);
MG_DLL_DECLR friend MGEReal operator- (double, const MGEReal&);
MG_DLL_DECLR friend MGEReal operator* (double, const MGEReal&);
MG_DLL_DECLR friend MGEReal operator/ (double, const MGEReal&);

///String stream Function.
MG_DLL_DECLR friend std::ostream& operator<< (std::ostream&, const MGEReal& );

////////////Constructor////////////

///Default Constructor
MGEReal(double val=0.0):m_value(val){;};

///infinite=-1 means minus_infinite,
///         +1 means plus_infinite.
MGEReal(MGINFINITE_TYPE infinite);

////////////Operator overload.////////////

///To cast to double.
operator double()const{ return m_value; };

///Operator definitions.
MGEReal operator+ (double) const;
MGEReal operator+ (const MGEReal&) const;

MGEReal& operator+= (double);
MGEReal& operator+= (const MGEReal&);

///Unary minus.
MGEReal operator- () const;

MGEReal operator- (double) const;
MGEReal operator- (const MGEReal&) const;

MGEReal& operator-= (double);
MGEReal& operator-= (const MGEReal&);

MGEReal operator* (double) const;
MGEReal operator* (const MGEReal&) const;

MGEReal& operator*= (double);
MGEReal& operator*= (const MGEReal&);

MGEReal operator/ (double) const;
MGEReal operator/ (const MGEReal&) const;

MGEReal& operator/= (double);
MGEReal& operator/= (const MGEReal&);

///Comparison
bool operator== (const MGEReal& r2) const;
auto operator<=> (const MGEReal& r2) const->std::partial_ordering;
bool operator== (const double t) const;
auto operator<=> (const double t) const->std::partial_ordering;

////////////Member Function////////////

///return -1 if minus_infinite(), 1 if plus_infinite(), else 0.
int infinite_coef()const;

///Test if this is equal to t regarding to base.
bool equal_base(double t, double base)const;
bool equal_base(const MGEReal& t,double base)const;

///Test if this is finite.
bool finite()const{return infinite_coef()==0;};

///Test if this is infinite.
bool infinite()const{return infinite_coef()!=0;};

///Invert the sign of this.
void invert(){m_value*=-1.;};

///Test if this is minus infinite.
bool minus_infinite()const{return (m_value<=(-mgInfiniteVal));};

///Test if this is plus infinite.
bool plus_infinite()const{return (mgInfiniteVal<=m_value);};

///Update this to double value val.
void set_real(double val){m_value=val;};

///Update this to be plus infite.
void set_plus_infinite(){m_value=mgInfiniteVal+1.;};

///Update this to minus infinite.
void set_minus_infinite(){m_value=-mgInfiniteVal-1.;};

///Update this to be zero.
void set_zero(){m_value=0.;};

///Get the double data of this.
double value() const {return m_value;};

///Dump Functions.
///Calculate dump size
int dump_size() const;

///Dump Function
int dump(MGOfstream& ) const;

///Restore Function
int restore(MGIfstream& );
	
////////////Member Data////////////

private:
	double	m_value;///<when m_value<=-mgInfiniteVal, this is minus infinite.
					///<when m_value>=mgInfiniteVal, this is plus infinite.

friend class MGInterval;

};

/** @} */ // end of BASE group
#endif
