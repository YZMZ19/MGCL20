/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#ifndef _MGDrawParam_HH_
#define _MGDrawParam_HH_

#include <iosfwd>
#include "mg/MGCL.h"

class MGObject;
class MGContext;
class MGOfstream;
class MGIfstream;
class mgTLInputParam;

///@cond
#define MAXDIVNUM 64 ///<Define how many lines a curve of screen width be devided.
			///<When a curve is approximated by lines, this is used.
#define SURFACE_TOL_BY_CURVE 1.5 ///<Default surface tolerance factor for the tessellation.
			///<SURFACE_TOL_BY_CURVE*cruve_tolerace is the surface tolerance.
///@endcond

/** @file */
/** @addtogroup BASE
 *  @{
 */


///\brief Defines parameters to draw MGObject, maily to approximate by lines and facets.

///1. To approximate the objects to draw in wire mode.
///   m_span_length_wire and m_line_desity_wire_face are used for the wire mode.
///2. To approoximate the objects(faces or surfaces) to draw in shading mode.
///   m_curve_tolerance, m_surface_tolerance, and m_maximum_edge_length are 
///   used for the shading mode.
///
///(1) m_span_length_wire: Used for the approximation of a curve by a polyline.
///                        m_span_length_wire is the maximum length of the one
///                        line segment of the polyline.
///(2) m_line_desity_wire_face: Used for the approximation of a surface by polylines.
///                        All of the curves(boundary or inner lines) are approximated
///                        by m_span_length_wire. And m_line_desity_wire_face indicates
///                        how many wires be drawn for the inner lines of
///                        the wire representation of the face.
///(3) m_curve_tolerance:  The maximum deviation from a curve allowed for the tessellaiton
///                        of a surface. The length of the one line segment of this approximation
///                        is limitted by m_maximum_edge_length.
///(4) m_surface_tolerance: The maximum deviation allowed of a traiangle from surfaces
///                        forn the tessellation.
///(5) m_maximum_edge_length: The maximum length of a triangle of the tessellation.
///
/// The minimum parameter to construct MGDrawParam object is m_line_desity_wire_face.
/// All of the other parameters can be derived from objects to draw and m_line_desity_wire_face.
class MG_DLL_DECLR MGDrawParam{

public:

////////////Constructor////////////

MGDrawParam(
	double curve_tolerance=-1.,///<Maximum deviation allowed to approximate curves.
	double surface_tolerance=-1.,///<Maximum deviation allowed to approximate surfaces.
	double maximum_edge_length=-1.,	///<Maximum edge length for the tessellation.
	double span_length_wire=20,///<length of a line segment to approximate curves.
	int line_desity_wire_face=1///<Indicates how many wires be drawn for the inner lines of
		///<the wire representation of the face.
);

///Construct from a object(box) data and the span length to draw object.
MGDrawParam(
	const MGObject& obj,///< Target object to get box.
	double span_length_wire=20///< span length.
);

///Constructor.
MGDrawParam(
	const MGContext& contx,///< Context to get from.
	double span_length_wire=20///< span length.
);

///Conversion constructor.
MGDrawParam(
	const mgTLInputParam& tlpara
);

//////////// Operator overload.////////

MG_DLL_DECLR friend std::ostream& operator<< (std::ostream& out, const MGDrawParam& para);

// Serialization fucntion.
MG_DLL_DECLR friend MGOfstream& operator<< (MGOfstream& buf, const MGDrawParam& para);
MG_DLL_DECLR friend MGIfstream& operator>> (MGIfstream& buf, MGDrawParam& para);

double span_length_wire()const{return m_span_length_wire;};
int line_desity_wire_face()const{return m_line_desity_wire_face;};
double curve_tolerance_tess()const{return m_curve_tolerance_tess;};
double surface_tolerance_tess()const{return m_surface_tolerance_tess;};
double maximum_edge_length_tess()const{return m_maximum_edge_length_tess;};

void set_line_density(int line_density=1){m_line_desity_wire_face=line_density;};
void set_span_length(double span_length){m_span_length_wire=span_length;};
void set_maximum_edge_length_tess(double maximum_edge_length){
	m_maximum_edge_length_tess=maximum_edge_length;};

private:

//Build m_curve_tolerance_tess and m_surface_tolerance_tess.
void build_crv_srf_tolerance(
	double curve_tolerance,//Maximum deviation allowed to approximate curves.
	double surface_tolerance//Maximum deviation allowed to approximate surfaces.
);

////////////Member Data//////////
double m_span_length_wire;//The maximum length of the one line segment to approximate
		//curves by a polyline.
int m_line_desity_wire_face;//Indicates how many wires be drawn for the inner lines of
		///the wire representation of the face.
double m_curve_tolerance_tess;//The maximum deviation from a curve allowed for the tessellaiton
	/// of a surface. The length of the one line segment of this approximation
	/// is limitted by m_maximum_edge_length.
double m_surface_tolerance_tess;//The maximum deviation allowed of a traiangle from surfaces
	/// forn the tessellation.
double m_maximum_edge_length_tess;//The maximum length of a triangle of the tessellation.

};

///@cond
double compute_curve_tolerance(
	const MGObject& obj,	///<Target object 
	double span_length		///<Span length of MGOpenglView.
);
///@endcond

/** @} */ // end of BASE group

#endif
