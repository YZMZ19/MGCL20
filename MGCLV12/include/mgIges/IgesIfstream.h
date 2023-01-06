/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno             */
/* All rights reserved.                                             */
/********************************************************************/

#pragma once

/// @file
///	@brief  Declaration for class MGIgesIfstream.
///	@author System fugen

#include <fstream>
#include "mg/Straight.h"
#include "mg/LBRep.h"
#include "mg/RLBRep.h"
#include "mg/SBRep.h"
#include "mg/RSBRep.h"
#include "mg/Plane.h"
#include "mg/Sphere.h"
#include "mgIges/IgesFstream.h"
#include "mgIges/IgesFstream.h"
#include "mgIges/IgesVertexListMap.h"
#include "mgIges/Iges504EdgeListMap.h"

// forward declarations
class MGObject;
class MGCurve;
class MGSurface;
class MGGroup;
class MGGel;
class MGColor;
class MGBVertex;
class MGEdge;
class MGLoop;
class MGFace;
class MGShell;

/** @addtogroup FileInputOutput
 *  @{
 */

///MGIgesIfstream read in *.iges; *.igs file, transforming IGES objects to MGCL objects.
class MG_DLL_DECLR MGIgesIfstream: public MGIgesFstream{
	
	friend class MGIges504Edge;
	friend class MGIges504EdgeListMap;
	friend class MGIgesPD142;
	friend class MGIges508Edge;
	friend class MGIgesPD510;

	///VERTEX list of the directory entry of type 502 vertex.
	mutable MGIgesVertexListMap m_vertexListMap;

	///EDGE list of the directory entry of type 504 edge.
	mutable MGIges504EdgeListMap m_edgeListMap;

public:

/// Constructors.

	/// Creates an object of class MGIgesIfstream with a filename.
	MGIgesIfstream(const TCHAR* filename);

	///Open the file. When this is opened already, this will be closed, then will
	///be opened.
	void open(const TCHAR* filename);
	void close();

	///Test if input stream is open.
	bool is_open(){ return m_ifstream.is_open(); };

	///Test if input stream is good().
	///True will be returned when good.
	bool good()const{return m_ifstream.good();};

/// Operator overload.

	/// Reads all objects of the IGES stored in the file
	/// as MGCL objects.
	MGIgesIfstream& operator>>(MGGroup&);

private:
/// Member functions.
	
	///Convert i-th MGIgesDirectoryEntry object(m_DirectoryEntries[i])
	///to MGObject that is a newed object.
	///When de was not an independent object, null will be returned.
	///The 2nd form convert from MGIgesDirectoryEntry.
	MGGel* convert_to_gel(
		int i
	)const;
	MGGel* convert_to_gel(
		const MGIgesDirectoryEntry& de
	)const;


	///From the current stream position, get one line data.
	void get_one_line(
		char* lineData,	///<line data without ID letter and sequence number
			///<(that is, from column 1 to nchar) will be output.
			///<buffer length must be >=(nchar+1).
		char& sectionID_letter,	///<section identification letter of the line.
		int& sequence,	///<ascending sequence number of the line.
		int nchar=72	///<number of characters of one line
			///<(When Parameter Data section nchar=64. Otherwise, nchar=72)
	);

	///Convert de(type=314: Color definition entity) to MGColor.
	///Returned is a newed MGColor.
	MGColor* convert_color(
		const MGIgesDirectoryEntry& de
	)const;

	///Convert de to MGEllipse(a newed object). de must be of type 100(Circular Arc).
	MGEllipse* convert_arc(
		const MGIgesDirectoryEntry& de
	)const;

	///Convert de to MGCompositeCurve(a newed object). de must be of type 102
	///(Composite curve).
	MGCompositeCurve* convert_composite(
		const MGIgesDirectoryEntry& de
	)const;

	///Convert de to MGObject(a newed object). de must be of type 104(conic arc).
	///Output MGObject is a MGEllipse, or MGRLBRep(when the arc is parabola or hyperbola).
	MGObject* convert_conic_arc(
		const MGIgesDirectoryEntry& de
	)const;

	///Convert de(type=108:Plane) to MGObject(MGPlane or MGFace).
	///Returned is a newed object.
	MGObject* convert_plane(
		const MGIgesDirectoryEntry& de
	)const;

	///Convert de to MGObject(a newed object). de must be of type 190(Plane Surface).
	///Output MGObject is a MGPlane(unbounded infinite plane).
	MGPlane* convert_planeSurface(
		const MGIgesDirectoryEntry& de
	)const;

	///Convert de to MGObject(a newed object). de must be of type 192
	///(Right circular cylindrical surface).
	///Output MGObject is a MGCylinder(unbounded infinite cylinder).
	MGCylinder* convert_cylinder(
		const MGIgesDirectoryEntry& de
	)const;

	///Convert de to MGObject(a newed object). de must be of type 196
	///(Sphere surface).
	///Output MGObject is a MGSphere.
	MGSphere* convert_sphere(
		const MGIgesDirectoryEntry& de
	)const;

	///Convert de to MGObject(a newed object). de must be of type 158(SPHERE).
	MGSphere* convert_sphere158(
		const MGIgesDirectoryEntry& de
	)const;

	///Convert de(type=110: LINE) to MGStraigh.
	///Returned is a newed object.
	MGStraight* convert_line(
		const MGIgesDirectoryEntry& de
	)const;

	///Convert de(type=112: spline curve) to MGLBRep.
	///Returned is a newed object.
	MGLBRep* convert_spline(
		const MGIgesDirectoryEntry& d
	)const;

	///Convert de(type=116: point) to MGPoint.
	///Returned is a newed object.
	MGPoint* convert_point(
		const MGIgesDirectoryEntry& de
	)const;

	///Convert de(type=118: Ruled surface) to MGSBRep.
	///Returned is a newed object.
	MGSurface* convert_ruled_surface(
		const MGIgesDirectoryEntry& de
	)const;

	///Convert de(type=120: Revolution surface) to MGSurface.
	///Returned is a newed object.
	MGSurface* convert_revolution_surface(
		const MGIgesDirectoryEntry& de
	)const;

	///Convert de(type=122: Tabulated cylinder) to MGSurface.
	///Returned is a newed object.
	MGSurface* convert_tab_cyl(
		const MGIgesDirectoryEntry& de
	)const;

	///Convert de(type=126: Rational b-spline) to MGCurve(MGLBRep or MBRLBRep).
	///Returned is a newed object.
	MGCurve* convert_nurbs(
		const MGIgesDirectoryEntry& de
	)const;

	///Convert de(type=128: NURBS surface) to MGSurface(MGSBRep or MGRSBRep).
	///Returned is a newed object.
	MGSurface* convert_nurbs_surface(
		const MGIgesDirectoryEntry& de
	)const;

	///Convert de(type=143: bounded surface) to MGFace.
	///Returned is a newed object.
	MGFace* convert_bounded_surface(
		const MGIgesDirectoryEntry& de
	)const;

	///Convert de(type=144: Trimmed surface) to MGFace.
	///Returned is a newed object.
	MGFace* convert_trimmed_surface(
		const MGIgesDirectoryEntry& de
	)const;

	///Convert de(type=186: Manifold Solid B-Rep Object) to MGGell.
	///Returned is a newed MGGel.
	///When MSBO has one or more void shells, MSBO is represented as a MGGroup.
	///When MSBO does not have a void, MSBO is represented as a MGShell.
	MGGel* convert_MSBO(
		const MGIgesDirectoryEntry& de
	)const;

	///Convert de(type=402: Associatibity entries) to MGGroup.
	///Accepted type numbers are only group associativity: form number 1,7,14, and 15
	///Returned is a newed MGGroup.
	MGGroup* convert_group(
		const MGIgesDirectoryEntry& de
	)const;

	///Convert de(type=502: VERTEX List entries) to MGGroup.
	MGBVertex* convert_vertex(
		int vertex_list,///<the DE index of the VERTEX List Entry.
		int vertex  ///<List index of the vertices in ertex List Entry DE.
	)const;

	///Convert an edge in type=504 edge list entries to MGEdge.
	MGEdge* convert_edge(
		int edge_list,///<the DE index of the EDGE List Entry.
		int edge  ///<List index of the start vertex in the EDGE List DE.
	)const;

	///Convert de(type=508: LOOP entry) to MGLoop.
	///Returned is a newed MGLoop object.
	MGLoop* convert_loop(
		const MGIgesDirectoryEntry& de,	///<directory entry of type 508.
		const MGSurface& srf ///<Base surface whose boundary loo this loop will make.
	)const;

	///Convert de(type=510: FACE) to MGFace.
	///Returned is a newed MGFace object.
	MGFace* convert_face(
		const MGIgesDirectoryEntry& de///<pointer to the DE of the FACE Entry.
	)const;

	///Convert de(type=514: SHELL) to MGShell.
	///Returned is a newed MGShell object.
	MGShell* convert_shell(
		const MGIgesDirectoryEntry& de///<pointer to the DE of the SHELL Entry.
	)const;

	///initialization
	void initialize();

	///Transform obj if de has the transformation matrix.
	void transform(
		const MGIgesDirectoryEntry& de,	///<de of the object obj.
		MGGel& obj					///<Object to transform.
	)const;
	
private:
	 std::ifstream m_ifstream;///<Input stream.

};

/** @} */ // end of FileInputOutput group
