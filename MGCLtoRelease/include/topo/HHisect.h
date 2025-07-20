/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#ifndef _MGHHisect_HH_
#define _MGHHisect_HH_
/** @file */
/** @addtogroup IsectContainer
 *  @{
 */

#include <deque>
#include "mg/isect.h"
#include "mg/CompositeCurve.h"
#include "mg/FPline.h"

class MGSSisect;
class MGFSurface;
class MGHHisects;

///MGHHisect is to represent one continuous intersection line of shells.

///Intersection lines a shell with a shell, a face, or a surface.
///(MGCompositeCurve* iline, deque<MGFPline> uvl1, deque<MGFPline> uvl2)
///where iline is a world coordinate rep of the line, uvl1 is a deque of
///1st shell's face parameter rep, uvl2 is a deque of the 2nd shell's or
///face's parameter rep of the intersection line.
///uvl1[i] corresponds to uvl2[i] one by one for all i.
///The parameter ranges of all the uvl1[i] are continuous and the total of them is
///equal to the parameter range of iline. For uvl2, the same.
///Let uvl1[i]'s start parameter be t1, and end parameter t2, then
///uvl2[i]'s parameter range is also from t1 to t2.
///Let sf1 be MGSurfCurve(f1's surface, uvl1[i]), then sf1 is the same curve as
///the iline's part of the parameter range t1 to t2. And sf1 is also equal to
///MGSurfCurve(f2's surface, uvl2[i]).
///MGHHisect uses MGFPline to represent the intersection lines.
///MGHHisect is also used to represent a pojection curve. In this case the size of
///uvline2 is zero.
///**** Projection line rep and intersection line rep cannot be mixed. ****

class MG_DLL_DECLR MGHHisect:public MGisect{
public:
	typedef std::deque<MGFPline> container_type;
	typedef container_type::iterator iterator;
	typedef container_type::const_iterator const_iterator;

private:
	mutable std::unique_ptr<MGCompositeCurve> m_iline;
			///<World coordinates (x,y,z) rep of the intersection line.
	container_type m_uvlines1;
			///<1st (u,v) parameter rep line of the intersection line.
	container_type m_uvlines2;
			///<2nd (u,v) parameter rep line of the intersection line.
			///<This vector lenght can be 0 even m_uvlines1's vector length is more than 0,
			///<which means this isect does not have 2nd face data.

public:

////////Special member functions/////////
MGHHisect()=default;///void constructor.
~MGHHisect()=default;
MGHHisect(const MGHHisect& rhs)=delete;
MGHHisect& operator=(const MGHHisect& rhs)=delete;
MGHHisect(MGHHisect&& rhs);
MGHHisect& operator=(MGHHisect&& rhs);
	
///Construct from MGSSisect.
MGHHisect(
	const MGFSurface* face1,///<face1. This must not be null.
	const MGFSurface* face2,///<face2. This may be null
							///<(e.g. for face2 that is actually a surface).
	std::unique_ptr<MGSSisect>&& ssi	///<intersection line of face1 and face2 expressed as MGSSisect.
);

///Construct two faces intersection lines.
///uvline1 and 2 makes uvlines whose vector length is 1.
///MGHHisect takes the ownership of iline, uvline1, and uvline2
///(These must be newed objects).
///When face2 is null, and uvlines2!=null, it indicates face2 is actually a surface.
///When both of face2 and uvlines2 are null, it indivates (face2, uvline2) are not
///used. This case occurs when MGHHisect is used to represent projection lines.
MGHHisect(
	MGCurve* iline,	///<Intersection line of world coordinates.
			///<iline must be a newed object pointer and MGHHisect takes the ownership.
	const MGFSurface* face1,///<face1. This must not be null.
	MGCurve* uvline1,	///<Intersection line of face1's (u,v) coordinates.
			///<uvline1 must be a newed object pointer and MGHHisect takes the ownership.
	const MGFSurface* face2=0,///<When face2 is null, and uvlines2!=null, it indicates
						///<face2 is actually a surface.
	MGCurve* uvline2=0  ///<Intersection line of face2's (u,v) coordinates.
						///<Takes the ownership.
);

///Ordering functions.
bool operator== (const MGisect& is)const;
bool operator< (const MGisect& is)const;
bool operator< (const MGCCisect& is)const{return false;};
bool operator< (const MGCSisect& csi)const{return false;};
//bool operator< (const MGCFisect& is)const{return false;};
bool operator< (const MGSSisect& is)const{return false;};
//bool operator< (const MGFFisect& is)const{return false;};

///Comparison operator.
bool operator== (const MGHHisect& hhi2)const;
bool operator< (const MGHHisect& hhi2)const;
bool operator> (const MGHHisect& hhi2)const{return hhi2<(*this);};
bool operator<= (const MGHHisect& hhi2)const{return !(hhi2<(*this));};
bool operator>= (const MGHHisect& hhi2)const{return !((*this)<hhi2);};
bool operator!= (const MGHHisect& hhi2)const{return !operator==(hhi2);};


///Extract a connected line from hhivec one by one and build one continuous
///line. Function's return value is:
///true: if line(s) are extracted and released from hhivec.
///false: no line are extracted.
bool build_one(MGHHisect& hhi);
bool build_one(MGHHisects& hhivec);

///Change parameter range, be able to change the direction by providing
///t0 greater than t1.
void change_range(
	double t0,	///<Parameter value for the start of original. 
	double t1	///<Parameter value for the end of original. 
);

///Connect a line to this HHisect.
///When both of face2 and uvlines2 are null, it indivates (face2, uvline2) are not
///used. This case occurs when MGHHisect is used to represent projection lines.
///****Projection lines rep. and intersection lines rep cannot be mixed.
///iline, uvline1, and uvline2(or hhi2) must have the same direction.
///iline's direction must be equal to this HHisect's.
///MGHHisect takes the ownership of iline, uvline1, and uvline2
///(These must be newed objects).
void connect_line_to_end(
	MGCurve* iline,	///<Intersection line of world coordinates.
	const MGFSurface* face1,///<face1. This must not be null.
	MGCurve* uvline1,///<Intersection line of face1's (u,v) coordinates.
	const MGFSurface* face2=0,///<face2. This may be null
					///<(e.g. for face2 that is actually a surface).
	MGCurve* uvline2=0///<Intersection line of face2's (u,v) coordinates.
					///<takes the ownership of all the curves of ssi.
);
void connect_line_to_end(
	MGHHisect&& hhi2		///<After connected, this hhi2's member data's ownership
						///<will be transfered to this MGHHisect,
						///<just like std::auto_ptr's assignment.
);

///Connect a line to this HHisect.
///When both of face2 and uvlines2 are null, it indivates (face2, uvline2) are not
///used. This case occurs when MGHHisect is used to represent projection lines.
///****Projection lines rep. and intersection lines rep cannot be mixed.
///iline, uvline1, and uvline2(or hhi2) must have the same direction.
///iline's direction must be opposite to this HHisect's.
///MGHHisect takes the ownership of iline, uvline1, and uvline2
///(These must be newed objects).
void connect_line_to_start(
	MGCurve* iline,		///<Intersection line of world coordinates.
	const MGFSurface* face1,///<face1. This must not be null.
	MGCurve* uvline1,	///<Intersection line of face1's (u,v) coordinates.
	const MGFSurface* face2=0,///<face2. This may be null
						///<(e.g. for face2 that is actually a surface).
	MGCurve* uvline2=0///<Intersection line of face2's (u,v) coordinates.
						///<takes the ownership of all the curves of ssi.
);
void connect_line_to_start(
	MGHHisect&& hhi2		///<After connected, this hhi2's member data's ownership
						///<will be transfered to this MGHHisect,
						///<just like std::auto_ptr's assignment.
);

///Test if this has face2's information.
bool has_face2_data()const{return !m_uvlines2.empty();};

///Return the world coordinate isect data.
const MGCurve& iline()const{return *m_iline;};
MGCurve& iline(){return *m_iline;};

///Test if this is a null HHisect or not.
bool is_null()const{return m_iline.get()==nullptr;};

///Get the start/end parameter value of the isect line.
double param_e()const{return m_iline->param_e();};
double param_s()const{return m_iline->param_s();};

///Reverse the direction of this intersection line.
void reverse_direction();

///return number of uvlines.
int num_of_uvline() const{return int(m_uvlines1.size());};

///Release the pointer of the last curve.
///Returned will be the released MGCurve pointer.
void release_back(
	MGCurve*& ilineLast,
	MGFPline& uvline1Last,
	MGFPline& uvline2Last
);

///Release the pointer of the 1st curve.
///Returned be the released MGCurve pointer.
void release_front(
	MGCurve*& iline1st,
	MGFPline& uvline11st,
	MGFPline& uvline21st
);

///Release the pointer of the iline curve.
///Returned will be the released MGCurve pointer.
MGCompositeCurve* release_line(){return m_iline.release();};

///Return uvline1.
const std::deque<MGFPline>& uvlines1() const{return m_uvlines1;};

///Return uvline2.
const std::deque<MGFPline>& uvlines2() const{return m_uvlines2;};

///Return i-th uvline.
const MGFPline& uvline1(int i)const;
const MGFPline& uvline2(int i)const;

///exchange12 1st and 2nd lines.
///This can be used only for intersection line rep.
void exchange12();

///Return the manifold dimension of the intersection, i.e.
///0: when the intersection is a point,
///1: when                  is a curve,
///2: when                  is a surface.
int manifold_dimension()const{return 1;};

/// Output virtual function.
std::ostream& toString(std::ostream& ostrm)const;

///Return the object of the intersection(world coordinates representation).
const MGObject& isect()const{return *m_iline;};

};

/** @} */ // end of IsectContainer group
#endif
