/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#ifndef _mgCPointsVec_HH_
#define _mgCPointsVec_HH_

#include <memory>
#include <iosfwd>
#include <vector>
#include <deque>
#include <algorithm>
#include "topo/LEPoint.h"
#include "topo/Loop.h"
#include "topo/Face.h"

///@cond

class MGTrimLoop;
class mgCPoints;
class mgCPointsVec;
class mgTL2PlBridge;
using UniqueTrimLoop = std::unique_ptr<MGTrimLoop>;

class MG_DLL_DECLR mgCONNECT_POINT{
public:
	enum OUTIN{coming_in, going_out};
private:
	MGTrimLoop* m_trim_loop;//MGTrimLoop that touch the inner loop.
	OUTIN m_outin;
		//indicates if the contact point is coming_in of m_trim_loop(that is,
		//the end point of m_trim_loop), or going_out(that is, the start point
		//of m_trim_loop).

public:
	mgCONNECT_POINT():m_trim_loop(0){;};
	mgCONNECT_POINT(MGTrimLoop* loop, OUTIN outin)
		:m_trim_loop(loop), m_outin(outin){;}

	bool operator< (const mgCONNECT_POINT& linf2)const;
	bool operator== (const mgCONNECT_POINT& linf2)const;
	bool operator!= (const mgCONNECT_POINT& linf2)const{return !operator==(linf2);};
	MGVector eval_deriv()const;
	bool is_null()const{return m_trim_loop->is_null();};
	bool is_going_out(){return m_outin==going_out;};
	bool is_coming_in(){return m_outin==coming_in;};
	MGLEPoint lep()const{
		if(m_outin==coming_in)
			return m_trim_loop->end_lep();
		else
			return m_trim_loop->start_lep();
	};
	std::unique_ptr<MGLoop> loop_clone(){
		return std::unique_ptr<MGLoop>(new MGLoop(*(m_trim_loop->loop())));
	};
	OUTIN outin()const{return m_outin;};
	MGTrimLoop* trim_loop(){return m_trim_loop;};
	const MGTrimLoop* trim_loop()const{return m_trim_loop;};

	//obtain the next point of this connect point.
	//The next point may be in the same loop, or in a different loop.
	//function's return value is the mgCONNECT_POINT of the next point.
	mgCONNECT_POINT next_point(
		mgCPointsVec& cpointsVec//container of MGTrimLoop.
	);

	//obtain the previous point of this connect point.
	//The next point may be in the same loop, or in a different loop.
	//function's return value is the mgCONNECT_POINT of the previous point.
	mgCONNECT_POINT prev_point(
		mgCPointsVec& cpointsVec//container of MGTrimLoop.
	);

	//Debug Function
	MG_DLL_DECLR friend std::ostream& operator<< (std::ostream& out, const mgCONNECT_POINT& cpoint);

};

class MG_DLL_DECLR mgCPoints{
public:
	std::vector<mgCONNECT_POINT> m_cpoints;

	mgCONNECT_POINT& operator[](int i){return m_cpoints[i];};

	int size()const{return (int)m_cpoints.size();}
	void sort(){std::sort(m_cpoints.begin(), m_cpoints.end());}
	void push_back(const mgCONNECT_POINT& cpoint){m_cpoints.push_back(cpoint);};

	int next_id(int j)const{return (j+1)%m_cpoints.size();};
	int prev_id(int j)const{
		int n=(int)m_cpoints.size();
		return (j+n-1)%n;
	};
	mgCONNECT_POINT next_point(int j){
		return m_cpoints[next_id(j)];
	}
	mgCONNECT_POINT prev_point(int j){
		return m_cpoints[prev_id(j)];
	}

	//get the id of m_cpoints that includes tloop
	int find_trim_loop(const MGTrimLoop* tloop, mgCONNECT_POINT::OUTIN out_in)const;

	//Debug Function
	MG_DLL_DECLR friend std::ostream& operator<< (std::ostream& out, const mgCPoints& cpoints);
};

class MG_DLL_DECLR mgCPointsVec{

const MGFace& m_face;//Target face to trim
const MGPosition& m_uv;//trimming positioning data.
std::vector<UniqueTrimLoop> m_trim_loops;//Trimloop's
std::vector<mgCPoints> m_cpoints_vec;//vector of mgCPoints.

public:
	
////////Special member functions/////////
mgCPointsVec()=delete;
~mgCPointsVec()=default;//Destructor.
mgCPointsVec(const mgCPointsVec& rhs)=delete;//Copy constructor.
mgCPointsVec& operator=(const mgCPointsVec& gel2)=delete;//Copy assignment.
mgCPointsVec(mgCPointsVec&& rhs)=default;//Move constructor.
mgCPointsVec& operator=(mgCPointsVec&& rhs)=default;//Move assignment.

explicit mgCPointsVec(
	const MGFace& face, //target face.
	const MGPosition& uv=MGPosition()
):m_face(face), m_uv(uv), m_cpoints_vec(face.get_number_of_boundaries()){;};

//Extract MGTrimLoop & mgCPoints info, and build this mgCPointsVec data.
void build_CPointsVec(
	mgTL2PlBridge& bridge	//original loops to trim.
);

//Extract MGTrimLoop & mgCPoints info, and build this mgCPointsVec data.
//While processing, a closed outer loop or closed inner loop is found, store it in
//the oloops or iloops.
//When m_uv is not null, oloops contains the only closed outer loop that contain m_uv,
//and return function's return code as 1. If the closed loop was not found,
//generate m_trim_loops(the vector of MGTrimLoop).
//When m_uv is null, put all the outer loops into oloops, and inner loops into iloops.
//Function's return value is
//  1: when m_uv is not null, a closed loop that includes m_uv was extracted into oloops.
//  0: when m_uv is not null, no closed outer loop that includes m_uv was found
//     and trim loop info were output into m_trim_loops and
//     mgCONNECT_POINT's into m_cpoints_vec. Closed inner loops may be output to iloops.
//When m_uv is null, function's return value is always 0.
int extract_trim_loop_info(
	const MGLoop& network,	//original loops to trim.
	std::vector<UniqueLoop>& oloops,//When m_uv is not nul, the closed outerboundary loop
		//that includes m_uv is output and the return code ==1.
		//When m_uv is null, all the detected closed outerboundary loops will be output,
		//and the function's return value is 0.
	std::vector<UniqueLoop>& iloops//closed inner boundary loop that
							//includes uv is output when return code =0 or 1.
);

//Extract a boundary out of networks that includes m_uv.
//m_uv is a parameter of m_face.
//Function's return value is
//  0: no loops were extracted,
//  2: Loops to trim face were output into used_tloops.
//     If trimming is performed by used_tloops[i] for i=0,...,n-1 one by one,
//     face will be trimmed by the part of the face that includes uv.
int extract_uv_boundary(
	std::unique_ptr<MGLoop>& loop,	//the loop that includes m_uv will be output
		//when return code=2.
	std::deque<MGTrimLoop*>& used_tloops
);

//extract the loop whose start point is first_point.
//Function's return value is:
//0: no loop was extracted since first_point was null.
//1: a loop was extracted into loop.
int extract_loop(
	mgCONNECT_POINT& first_point,
	std::unique_ptr<MGLoop>& loop, //the loop will be output.
	std::deque<MGTrimLoop*>& tloops//used MGTrimLoop to extract loop will be output.
);

mgCPoints& operator[](int i){return m_cpoints_vec[i];};
void push_at(int i, const mgCONNECT_POINT& cpoint){m_cpoints_vec[i].push_back(cpoint);};
int size()const{return (int)m_cpoints_vec.size();};
bool no_uv()const{return m_uv.is_null();};

MG_DLL_DECLR friend std::ostream& operator<< (std::ostream& out, const mgCPointsVec& cpvec);

};

///@endcond
#endif
