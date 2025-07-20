/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#ifndef _MGLEPoint_HH_
#define _MGLEPoint_HH_
/** @file */
/** @addtogroup IsectContainer
 *  @{
 */

#include "mg/Unit_vector.h"
#include "topo/Complex.h"

class MGLCisect;
class MGLPoint;
class MGLoop;
class MGEdge;

//
//Define MGLEPoint Class.

///Is to represent a Loop's point.

///This expreses (i, t), where i is the edge's iterator in a loop,
///and t is the parameter value of the curve of the edge.
class MG_DLL_DECLR MGLEPoint{

public:

///String stream Function
MG_DLL_DECLR friend std::ostream& operator<< (std::ostream&, const MGLEPoint& );

/////////Constructor/////////
MGLEPoint(){;};
MGLEPoint(
	MGComplex::const_iterator i,	///<Loop's edge iterator.
	double t)						///<Parameter value of i-th pcell curve.
	:m_i(i), m_t(t){;};

///Construct from the loop, edge number, and the edge's parameter value.
MGLEPoint(
	const MGLoop& lp,///<The target loop.
	int i,			///<Edge number.
	double t		///<Parameter value of the edge i.
);

///Conversion constructor from MGLCisect.
MGLEPoint(const MGLCisect& lci);

///convert loop1's MGLEPoint lep to loop2's MGLEPoint. loop2 must be a copy of loop1.
MGLEPoint(
	const MGLoop& loop1,///<Loop 1.
	const MGLEPoint& lep,///<Parameter value of i-th pcell curve of loop1.
	const MGLoop& loop2 ///<Loop2.
);

/////////Operator oveload/////////

///Comparison operator.
bool operator< (const MGLEPoint& lp)const;
bool operator> (const MGLEPoint& lp)const;
bool operator<= (const MGLEPoint& lp)const;
bool operator>= (const MGLEPoint& lp)const;
bool operator== (const MGLEPoint& lp)const;
bool operator!= (const MGLEPoint& lp)const{return !operator==(lp);};

/////////Member function/////////

///Obtain the edge pointer.
const MGEdge* edge() const;
MGEdge* edge_to_update() const;

///return loop's edge number.
int edge_num()const;

///Evaluation of the loop at the LEPoint.
///When nderi=0, get the positional data(a parameter (u,v) of the surface)
///at the point.
MGVector eval(int nderi=0)const;

///Evaluation of the star curves of the loop at the point t.
///When nderi=0, get a position of the surface at the boundary point.
///The star curve is SurfCurve(face's surface, loop's curve).
///(The star curve has the same world coordinate with the binder curve's, but
///their direction may be opposite. The star curve has always the same direction
///as the loop.)
MGVector eval_star(int nderi=0)const;

///Test if two LEPoints are of the same edge.
bool equal_edge(const MGLEPoint& le2) const;

///Test if two LEPoints are of the same position.
bool equal_position(const MGLEPoint& le2) const;

///Test if this is the same position as P.
bool equal_position(const MGPosition& P) const;

///Compute a vector at the MGLEPoint point that goes inside the face
///and perpendicular to the boundary loop.
MGUnit_vector inner_vector()const;

///test if this is the end point of the loop.
bool is_end_point()const;

///test if this is the start point of the loop.
bool is_start_point()const;

///test if this is the end point of the loop.
bool is_Edge_end_point()const;

///test if this is the start point of the loop.
bool is_Edge_start_point()const;

///Return the iterator of the edge.
MGComplex::const_iterator iterator()const{return m_i;};

///Get loop pointer.
const MGLoop* loop()const;

///Return edge's curve parameter value.
double param()const{return m_t;};

///Set iteraror of the edge and parameter value of the edge.
void set(
	MGComplex::const_iterator i,///<Iterator of the edge.
	double t	///<parameter of the edge.
);							

///Update only parameter value of the edge.
void set(
	double t///<parameter of the edge.
){m_t=t;};

///Get the end point of the next edge in the loop sequence.
MGLEPoint end_of_pre_edge()const;

///Get the start point of the next edge in the loop sequence.
MGLEPoint start_of_next_edge()const;

private:
	MGComplex::const_iterator m_i;	///<edge iterator in the loop.
	double m_t=0.;						///<edge's curve parameter value.
	

};

/** @} */ // end of IsectContainer group
#endif
