/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#ifndef _MGFPline_HH_
#define _MGFPline_HH_

class MGCurve;
class MGFSurface;

/** @file */
/** @addtogroup GEORelated
 *  @{
 */

///Face's (u,v) parameter value line.

///MGFPline is to represent an parameter (u,v) line of a face.
///(MGFSurface* f, MGCurve uvline) where f is a face pointer, and uvline is
///a parameter (u,v) line of the face f.
///
///MGFPline is used to express Shell's intersection lines.
///The behavior of MGFPline is like an auto_ptr. Copy or assignment
///of MGFPline means transfer of the ownership of all the included curves
///to copied or assigned MGFPline and original MGFPline does not have the
///ownership of the curves any more. Users should be aware of it.
class MG_DLL_DECLR MGFPline{

public:
	
///String stream Function
MG_DLL_DECLR friend std::ostream& operator<< (std::ostream& ostrm, const MGFPline& fpl);

////////Special member functions/////////
MGFPline():m_face(nullptr){;};	///void constructor.
~MGFPline()=default;			///Destructor.
MGFPline(const MGFPline& fpl)=delete;///Copy constructor.
MGFPline& operator= (const MGFPline& fpl)=delete;///Copy assignment.
MGFPline(MGFPline&& fpl);		///Move constructor.
MGFPline& operator= (MGFPline&& fpl);///Move assignment.

///Construct from all the necessary data.
MGFPline(
	const MGFSurface* face,	///<face1.
	MGCurve* uvline///<(u,v) line of the face, takes the ownership of the curve.
					///<That is uvline must be a newed object pointer.
):m_face(face), m_uvline(uvline){;};

///Comparison operator.
bool operator< (const MGFPline& fpl2)const;
bool operator> (const MGFPline& fpl2)const{return fpl2<(*this);};
bool operator<= (const MGFPline& fpl2)const{return !(fpl2<(*this));};
bool operator>= (const MGFPline& fpl2)const{return !((*this)<fpl2);};
bool operator== (const MGFPline& fpl2)const;
bool operator!= (const MGFPline& fpl2)const{return !operator==(fpl2);};

////////Member function////////

///Change parameter range, be able to change the direction by providing
///t1 greater than t2.
void change_range(
	double t0,	///<Parameter value for the start of original. 
	double t1	///<Parameter value for the end of original. 
);

///Get face's pointer.
const MGFSurface* face()const{return m_face;};

///Reverse the direction of this line.
void reverse_direction();

///Release the uvline curve pointer from this.
MGCurve* release_line();

///Return face's (u,v) parameter representation line.
const MGCurve& uvline() const{return *m_uvline;}
MGCurve& uvline() {return *m_uvline;}

private:
	const MGFSurface* m_face;	///<Face pointer.
	mutable std::unique_ptr<MGCurve> m_uvline;///<2D line whose coordinates are
	///(u, v) of the m_face.
};

/** @} */ // end of GEORelated group

#endif
