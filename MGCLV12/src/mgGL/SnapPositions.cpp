/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/

#include "StdAfx.h"
#include <set>
#include "mg/Tolerance.h"
#include "mg/DNameControl.h"
#include "mg/Plane.h"
#include "mg/Curve.h"
#include "mg/Point.h"
#include "mg/Group.h"
#include "mg/PickObjects.h"
#include "topo/BVertex.h"
#include "topo/Loop.h"
#include "topo/Face.h"
#include "topo/Shell.h"
#include "mgGL/GLSLProgram.h"
#include "mgGL/SnapPositions.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//
//Implements MGSnapPositions Class.
//MGSnapPositions is a class to store array(vector) of MGPosition's,
//used for MGLocateState specifically.
//

MGSnapPositions::MGSnapPositions(snap_kind kind)
:m_snap_kind(kind){
}

///Copy constructor.
MGSnapPositions::MGSnapPositions(const MGSnapPositions& sp)
:mgVBO(sp),m_obj_nums(sp.m_obj_nums),m_snap_kind(sp.m_snap_kind)
,m_positions(sp.m_positions){
}

///Assignment
MGSnapPositions& MGSnapPositions::operator=(
	const MGSnapPositions& sp
){
	m_obj_nums=sp.m_obj_nums;
	m_snap_kind=sp.m_snap_kind;
	m_positions=sp.m_positions;
	return *this;
}

///Virtual Destructor
MGSnapPositions::~MGSnapPositions(){
}

//append this positions in m_positions into points.
void MGSnapPositions::append_position_data(std::vector<MGPosition>& points)const{
	points.insert(points.end(), m_positions.begin(), m_positions.end());
}

void MGSnapPositions::clear(){
	mgVBO::clearElements();
	m_positions.clear();
	m_obj_nums.clear();
}

//Extract position data.
void MGSnapPositions::extract(
	const MGPoint& point	//the curve to extract.
){
	switch(m_snap_kind){
	case endpos:
		m_obj_nums.push_back(obj_num(&point,1));
		m_positions.push_back(point.position());
		return;
	default:
		;
	}
}

//Extract position data.
void MGSnapPositions::extract(
	const MGCurve& crv	//the curve to extract.
){
	switch(m_snap_kind){
	case endpos:
		m_obj_nums.push_back(obj_num(&crv,2));
		m_positions.push_back(crv.start_point());
		m_positions.push_back(crv.end_point());
		return;
	case knotpos:
		{
		const MGKnotVector& t=crv.knot_vector();
		int k=t.order(), n=t.bdim();
		int numpoint=0;
		for(int i=k; i<n; i++){
			if(t[i]==t[i-1]) continue;
			m_positions.push_back(crv.eval(t[i]));
			numpoint++;
		}
		m_obj_nums.push_back(obj_num(&crv,numpoint));
		return;
		}
	case centerpos:
		m_obj_nums.push_back(obj_num(&crv,1));
		m_positions.push_back(crv.center());
		return;
	default:
		;
	}
}

///Get curve crv's parameter value from its snap_kind and point.
void get_param(
	const MGCurve& crv,	//the curve objective.
	MGSnapPositions::snap_kind snap_kind,
	const MGPosition& point,
	MGPosition& param
){
	param.resize(1);
	double& t=param(0);
	double error=MGTolerance::mach_zero();

	t=crv.param_s();
	switch(snap_kind){
	case MGSnapPositions::endpos:
		{
		MGVector Ps=crv.start_point(), Pe=crv.end_point();
		MGVector dif1=Ps-point, dif2=Pe-point;
		if(dif1%dif1<= dif2%dif2)
			break;
		else{
			t=crv.param_e();
			break;
		}
		}
	case MGSnapPositions::centerpos:
		param=crv.center_param();
		break;
	case MGSnapPositions::knotpos:
		{
		const MGKnotVector& tv=crv.knot_vector();
		int k=tv.order(), n=tv.bdim();
		for(int i=k; i<n; i++){
			t=tv[i];
			MGVector dif=crv.eval(t)-point;
			if(dif%dif<=error)
				break;
		}
		break;
		}
	default:
		;
	}
}

//Extract position data.
void MGSnapPositions::extract(
	const MGSurface& srf//the surface to extract.
){
	double u0=srf.param_s_u(), u1=srf.param_e_u();
	double v0=srf.param_s_v(), v1=srf.param_e_v();
	int numObj=0;
	switch(m_snap_kind){
	case vertexpos:
		{
		const MGPlane* pl=dynamic_cast<const MGPlane*>(&srf);
		if(pl) return;
		m_positions.push_back(srf.eval(u0,v0));
		m_positions.push_back(srf.eval(u1,v0));
		m_positions.push_back(srf.eval(u0,v1));
		m_positions.push_back(srf.eval(u1,v1));
		numObj=4;
		break;
		}

	case centerpos:
		m_positions.push_back(srf.center());
		numObj=1;
		break;

	default:
		;
	}
	if(numObj)
		m_obj_nums.push_back(obj_num(&srf,numObj));
}

///Get surface srf's parameter value from its snap_kind and point.
void get_param(
	const MGSurface& srf,	//the MGSurface objective.
	MGSnapPositions::snap_kind snap_kind,
	const MGPosition& point,
	MGPosition& param
){
	double u0=srf.param_s_u(), u1=srf.param_e_u();
	double v0=srf.param_s_v(), v1=srf.param_e_v();
	param.resize(2);
	double error=MGTolerance::mach_zero();

	MGVector dif;
	MGPosition Q;
	switch(snap_kind){
	case MGSnapPositions::vertexpos:
		{
		const MGPlane* pl=dynamic_cast<const MGPlane*>(&srf);
		if(pl)
			return;//This does not take place, since MGPlane does not have vertices.

		param(0)=u0; param(1)=v0;
		Q=srf.eval(param);
		dif=Q-point;
		if(dif%dif<=error){
			break;
		}
		param(0)=u1;
		Q=srf.eval(param);
		dif=Q-point;
		if(dif%dif<=error){
			break;
		}
		param(1)=v1;
		Q=srf.eval(param);
		dif=Q-point;
		if(dif%dif<=error){
			break;
		}
		param(0)=u0;
		Q=srf.eval(param);
		dif=Q-point;
		if(dif%dif<=error){
			break;
		}
		}

	case MGSnapPositions::centerpos:
		param=srf.center_param();
		break;
	default:
		;
	}
}

//Extract vertex of a topology.
void MGSnapPositions::extract(
	const MGFace& face	//the face to extract.
){
	int numObj=0;
	switch(m_snap_kind){
	case vertexpos:{
		int nlp=face.number_of_boundaries();
		for(int i=0; i<nlp; i++){
			const MGLoop& lp=*(face.loop(i));
			std::vector<MGBCell*> bvrtices;
			lp.get_all_boundary_binders(bvrtices);
			for(auto& j:bvrtices){
				const MGBVertex* v=dynamic_cast<const MGBVertex*>(j);
				MGPosition uv=v->position();
				m_positions.push_back(face.eval(uv));
				numObj++;
			}
		}
	}
		break;

	case centerpos:
		m_positions.push_back(face.center());
		numObj=1;
		break;

	default:
		;
	}
	if(numObj)
		m_obj_nums.push_back(obj_num(&face,numObj));
}

///Get face f's parameter value from its snap_kind and point.
void get_param(
	const MGFace& f,	//the MGSurface objective.
	MGSnapPositions::snap_kind snap_kind,
	const MGPosition& point,
	MGPosition& param
){
	param.resize(2);
	double error=MGTolerance::mach_zero();
	int nlp=f.number_of_boundaries();

	switch(snap_kind){
	case MGSnapPositions::vertexpos:
		for(int i=0; i<nlp; i++){
			const MGLoop& lp=*(f.loop(i));
			std::vector<MGBCell*> bvrtices;
			lp.get_all_boundary_binders(bvrtices);
			for(auto& j : bvrtices){
				const MGBVertex* v=dynamic_cast<const MGBVertex*>(j);
				param=v->position();
				MGVector dif=f.eval(param)-point;
				if(dif%dif<=error){
					return;
				}
			}
		}
		param.set_null();//Hopefully this does not happen.
		break;
		
	case MGSnapPositions::centerpos:
		param=f.center_param();
		break;
	default:
		param.set_null();//Hopefully this does not happen.
		;
	}
}

void MGSnapPositions::extract(
	const MGShell& shell	//the surface to extract.
){
	switch(m_snap_kind){
	case vertexpos:
		{
		MGShell::const_iterator i=shell.pcell_begin(), ie=shell.pcell_end();
		for(; i!=ie; i++){
			const MGFace* fi=shell.face(i);
			extract(*fi);
		}
		return;
		}
	default:
		;
	}
}

void MGSnapPositions::extract(
	const MGGel& gel	//the gel to extract.
){
	const MGGel* gelP = &gel;
	const MGCurve* crv=dynamic_cast<const MGCurve*>(gelP);
	if(crv){ extract(*crv); return;}

	const MGFace* f = dynamic_cast<const MGFace*>(gelP);
	if(f){ extract(*f); return;}

	const MGSurface* srf = dynamic_cast<const MGSurface*>(gelP);
	if(srf){ extract(*srf); return;}

	const MGPoint* P = dynamic_cast<const MGPoint*>(gelP);
	if(P){	extract(*P); return;}

	const MGShell* shl = dynamic_cast<const MGShell*>(gelP);
	if(shl){  extract(*shl); return;}

	const MGGroup* group = dynamic_cast<const MGGroup*>(gelP);
	if(group){
		MGGroup::const_iterator i=group->begin(), ie=group->end();	
		for(; i!=ie; i++)
			extract(**i);
		return;
	}
}

void MGSnapPositions::extract(
	const std::list<const MGGel*>& gel_list	//the group to extract.
){
	std::list<const MGGel*>::const_iterator i=gel_list.begin(), ie=gel_list.end();	
	for(; i!=ie; i++)
		extract(**i);
}

void MGSnapPositions::extract(
	const MGPickObjects& pobjs	//array of pick objects to extract.
){
	MGPickObjects::const_iterator i=pobjs.begin(), ie=pobjs.end();
	for(; i!=ie; i++)
		extract(*((**i).leaf_object()));
}

///Get the object of the position posID of m_positions.
const MGObject* MGSnapPositions::object(int posID)const{
	std::vector<obj_num>::const_iterator i=m_obj_nums.begin(), ie=m_obj_nums.end();
	int numTotal=0;
	for(; i!=ie; i++){
		int n=i->second;//number of points in this object.
		if(!n) continue;

		numTotal+=n;
		if(posID<numTotal)
			return i->first;
	}
	return 0;
}

void MGSnapPositions::make_display_list(
	MGCL::VIEWMODE vmode
){
	clearElements(MGCL::WIRE_AND_SHADING);

	setStaticAttribPointSize(5.f);
	std::vector<obj_num>::const_iterator i=m_obj_nums.begin(), ie=m_obj_nums.end();
	int j=0;//Point id counter.
	for(; i!=ie; i++){
		int n=i->second;//number of points in this object.
		if(!n) continue;
		for(int k=0; k<n; k++){
			Begin(GL_POINTS);Vertex3dv(m_positions[j++].data()); End();
		}
	}
	setDirty(false);
}

void MGSnapPositions::selectionDraw(MGCL::VIEWMODE viewMode){
	if(!is_made())
		make_display_list();

	mgGLSLProgram* glsl=mgGLSLProgram::getCurrentGLSLProgram();
	glsl->setFuncType(mgGLSL::Select);

	mgCoordinateTypeSwitcher coordType(m_coordinateType, getAnchorPosition());//save the coordinate type.
	size_t n=m_elements.size();
	for(unsigned i=0; i<n; i++){
		mgGLSL::setColorAsSelectionName(i+1);
		UniqueVBOElement& elmi=m_elements[i];
		elmi->selectionDraw(MGCL::WIREVIEW);
	}
}

void MGSnapPositions::get_pick_data(
	const std::set<unsigned>& selected,///Selected data of pick_to_select_buf.
	MGPosition& point,	//point data will be output.
	const MGObject*& obj,//When function's return value is nearpos, end, or knot,
				//the point's parameter value of the object will be returned.
	MGPosition& param	///<When obj is an MGCurve and
		///<this snap kind is nearpos, endpos, knotpos, or centerpos,
		///<the point's parameter value of the curve be returned: param.sdim()=1.
		///<When obj is an MGFSurface and this snap kind is centerpos, or vertexpos
		///<the point's parameter value(u,v) be returned:param.sdim()=2.
)const{

	const std::set<unsigned>::const_iterator
		i=selected.begin(), iend=selected.end();
	assert(i!=iend);

	param.set_null(), point.set_null();
	unsigned int pos=*i-1;//Other name is the id of the m_points.
				///-1 because selectionDraw put name as position+1.
	if(pos>=m_positions.size())
		return;

	double error=MGTolerance::mach_zero();
	point=m_positions[pos];
	obj=object(pos);
	const MGCurve* crv=dynamic_cast<const MGCurve*>(obj);
	if(crv){
		get_param(*crv,m_snap_kind,point,param);
		return;
	}

	const MGSurface* srf = dynamic_cast<const MGSurface*>(obj);
	if(srf){
		get_param(*srf,m_snap_kind,point,param);
		return;
	}
		
	const MGFace* f = dynamic_cast<const MGFace*>(obj);
	if(f){
		get_param(*f,m_snap_kind,point,param);
	}
}
