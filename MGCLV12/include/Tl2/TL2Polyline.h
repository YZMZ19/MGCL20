#ifndef _mgTL2Polyline_HH_
#define _mgTL2Polyline_HH_

#include "mg/LBRep.h"
#include "topo/Complex.h"

/****************************************************************/
/*   Copyright (c) 2021 by System fugen G.K.                    */
/*                       All rights reserved.                   */
/****************************************************************/

class MGPosition;
class mgTL2parameter;
class MGEdge;
class MGLEPoint;
class mgTL2LPlines;
class mgTLEdgePoint;
class mgTL2LPPoint;

/** @addtogroup UseTessellation
 *  @{
 */

///mgTL2Polyline holds a parameter line(polyline) of a surface.

///mgTL2Polyline is a proprietry class for Face tessellation.
///mgTL2Polyline holds (u,v) MGLBRep of order 2(polyline).
class mgTL2Polyline: public MGLBRep{

public:

	enum polyline_type{
		WHOLE_INNER=0,	//whole points are free inner points.
		START_BOUNDARY,	//only start point is connected to a boundary line.
		END_BOUNDARY,	//only end point is connected to a boundary line.
		START_END_BOUNDARY,//Only start and end points are connected to a boundary.
						//Inner points except start and end are not on a boundary.
		WHOLE_BOUNDARY	//All of the points are on a boundary.
	};

private:
	polyline_type m_type=WHOLE_INNER;//indicated how this line is connectied to boundaries.
	const mgTL2parameter* m_tlparam=nullptr;//Tessellation parameter.
	short m_start[3], m_end[3];//Ids of m_Bpolylines of mgTL2Face, indicates
			//this Polyline's start or end point. That is,
			//e.g. let m_start[.] be (i,j,k), then k-th point of the polyline
			//m_Bpolylines[i][j] is the start point of this polyline.
	mutable MGPosition m_uv_start;//Start point uv data if is not null.
	mutable MGVector m_dire_start, m_dire_end;//Start and end points direction
			//in the world coordinate if is not null.

public:

//////////// constructor ///////////////
	
mgTL2Polyline(){ ; };

//Construct null mgTL2Polyline.
mgTL2Polyline(const mgTL2parameter& para)
:MGLBRep(),m_tlparam(&para),m_type(WHOLE_INNER){;};

//Construct dummy mgTL2Polyline of lengh n.
mgTL2Polyline(
	const mgTL2parameter& para,
	int n //nubmer of points
);

//Copy constructor.
mgTL2Polyline(const mgTL2Polyline& pline2);
mgTL2Polyline(mgTL2Polyline&& pline2)noexcept;
mgTL2Polyline& operator=(mgTL2Polyline&&)noexcept;

//Construct mgTL2Polyline from (u,v) curve representation of the surface.
//The type of the boundary is WHOLE_INNER.
mgTL2Polyline(
	const mgTL2parameter& para,
	const MGCurve& uvline
);

//Construct mgTL2Polyline of straight line of point number n from uvS to uvE.
//The type of the boundary is WHOLE_INNER, and so necessary to set_endID/set_startID,
//if necessary.
mgTL2Polyline(
	const mgTL2parameter& para,
	int n, //nubmer of points to 
	const MGPosition& uvE,
	const MGPosition& uvS
);

//Construct mgTL2Polyline from the straight from pvS to pvE.
//The straight is guaranteed to keep the deviation from the surface
//within tessellation curve error and within the maximum edge length.
//The type of the boundary is set according the boundary infromation of pvS and pvE.
mgTL2Polyline(
	const MGLEPoint& pvE,
	const MGLEPoint& pvS
);

//Construct mgTL2Polyline from the straight from pvS to pvE.
//The type of the boundary is set according the boundary infromation of pvS and pvE.
mgTL2Polyline(
	const mgTL2LPlines& lines,
	const mgTLEdgePoint& pvE,
	const mgTLEdgePoint& pvS
);

polyline_type boundaryType()const{return m_type;};

//////////// Member Function ///////////////

void setParam(const mgTL2parameter& p){m_tlparam=&p;};

//Compute the angle of this curves's world rep and u=const iso-parameter line world rep
//at the parameter t. The parameter t is normalized value from 0.(start point)
//to 1.(end point).
double angle2Uline(double t)const;

//Compute the angle of this curves's world rep and u=const iso-parameter line world rep
//at the point id i.
double angle2Uline_at_i(int i)const;

///Construct new curve object by copying to newed area.
///User must delete this copied object by "delete".
mgTL2Polyline* clone()const;

//Get the direction of the line in the world coordinate at the normalized parameter t.
//The parameter t is normalized value from 0.(start point) to 1.(end point).
MGVector direWorld(double t)const;

//Get the start or end point direction of the line in the world coordinate.
const MGVector& direWorld_start()const;
const MGVector& direWorld_end()const;

//Evaluate the direction at t(not normalized ordinary parameter value) in its
//world representation.
MGVector direWorld_with_non_normalized(double t)const;

//Get id of m_Bpolylines of mgTL2Face from point id of this polyline.
//Function's return value is
//true if the point of idVertex is a boundary point, false if not.
bool get_id_from_VertexID(int idVertex, short id[3])const;

//Get id of m_Bpolylines of mgTL2Face from the parameter t.
//Function's return value is
//true if the point of t is a boundary point, false if not.
bool get_id(double t, short id[3])const;

//Update this by limiting the parameter range of the curve.
///Limitting is done at the knot parameter for both start and end.
void limit(const MGInterval& i1);

///Change direction of the line.
void negate();

//Get the number of points of this polyline.
int number_of_points()const{return bdim();};

///Debug Function
std::ostream& toString(std::ostream& ostrm)const;

const mgTL2parameter& TL2param()const{return *m_tlparam;};

//Get i-th point surface parameter (u,v) of this polyline
MGPosition uv(int i)const;

//Get start point surface parameter (u,v) of this polyline
const MGPosition& uv_start()const;

//Get i-th point(x,y,z) of this polyline
MGPosition xyz(int i, bool need_normal)const;

void setBoundaryID(
	const mgTL2LPPoint& P, //End point.
	bool isEnd = true
);
void set_type(polyline_type type){ m_type=type;};
void set_startID(const short start[3]);
void set_endID(const short end[3]);

void get_startID(short start[3]);
void get_endID(short end[3]);

//Update the point number. The parameter range is also updated according to the newN.
mgTL2Polyline& updatePointNumber(int newN);

};

//Get curve representation(MGLBRep of uniform B-Spline of order 2).
const mgTL2Polyline* TL2Polyline(const MGEdge* edg);
const mgTL2Polyline* TL2Polyline(MGComplex::const_iterator ei);

/** @} */ // end of UseTessellation group
#endif
