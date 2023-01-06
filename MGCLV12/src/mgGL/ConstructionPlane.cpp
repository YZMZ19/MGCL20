/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#include <bitset>
#include "StdAfx.h"
#include "mg/Box.h"
#include "mg/Position.h"
#include "mg/CSisect.h"
#include "mg/Straight.h"
#include "mg/Ofstream.h"
#include "mg/Ifstream.h"
#include "mg/LBRep.h"
#include "mgGL/Color.h"
#include "mgGL/glViewAttrib.h"
#include "mgGL/ConstructionPlane.h"
#include "mgGL/glslprogram.h"
#include "mgGL/Context.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

MGConstructionPlane::MGConstructionPlane(
):m_disabled(false),m_bind_to_grid(false){
	//m_lineColor.set_color(m_lColorV);
}

MGConstructionPlane::MGConstructionPlane(
	const MGPlane& plane,//construction plane.
	double uspan,		//span length along u axis.
	double vspan,		//span length along v axis.
	int uline_num,		//number of lines along u axis.
	int vline_num,		//number of lines along v axis.
	double nspan		//span length along normal axis.
):m_disabled(false),m_bind_to_grid(false)
,m_plane(MGUnit_vector(plane.u_deriv()),MGUnit_vector(plane.v_deriv()),MGPosition(plane.root_point()))
,m_uspan(uspan), m_vspan(vspan), m_nspan(nspan)
,m_unum(uline_num), m_vnum(vline_num)
{
	assert(uline_num>0 && vline_num>0);
	assert(uspan>0. && vspan>0. && nspan>0.);
	//m_lineColor.set_color(m_lColorV);
}

MGConstructionPlane::MGConstructionPlane(
	double origin[3],	//Origin's coordinate value.
	double uaxis[3],	//A vector value of the horizontal direction.
	double vaxis[3],	//A vector value of the vertical direction.
	double uspan,		//span length along u axis.
	double vspan,		//span length along v axis.
	int uline_num,		//number of lines along u axis.
	int vline_num,		//number of lines along v axis.
	double nspan		//span length along normal axis.
):m_disabled(false),m_bind_to_grid(false)
{
	assert(uline_num>0 && vline_num>0);
	assert(uspan>0. && vspan>0. && nspan>0.);
	m_plane=MGPlane(
		MGUnit_vector(MGVector(3,uaxis)),MGUnit_vector(MGVector(3,vaxis)),MGPosition(3,origin));
}

///Obtain the grid data of this.
void MGConstructionPlane::get_grid_data(
	double& uspan,		///<span length along u axis.
	double& vspan,		///<span length along v axis.
	int& uline_num,		///<number of lines along u axis.
	int& vline_num,		///<number of lines along v axis.
	double& nspan		///<span length along normal axis.
)const{
	uspan=m_uspan;
	vspan=m_vspan;
	nspan=m_nspan;
	uline_num=m_unum;
	vline_num=m_vnum;
}

void MGConstructionPlane::set_plane(const MGPlane& plane){
	m_plane=MGPlane(
		MGUnit_vector(plane.u_deriv()),MGUnit_vector(plane.v_deriv()),MGPosition(plane.root_point()));
	setDirty(true);
}


//Compute grid data and the plane from the plane and the grid span data.
void MGConstructionPlane::set_grid_data(
	const MGPlane& plane,//construction plane.
	double uspan,		//span length along u axis.
	double vspan,		//span length along v axis.
	int uline_num,		//number of lines along u axis.
	int vline_num,		//number of lines along v axis.
	double nspan		//span length along normal axis.
){
	assert(uline_num>0 && vline_num>0);
	assert(uspan>0. && vspan>0. && nspan>0.);

	set_plane(plane);
	m_uspan=uspan;
	m_vspan=vspan;
	m_nspan=nspan;
	m_unum=uline_num;
	m_vnum=vline_num;
	mgVBO::setDirty(true);
}

///Construct grid and the plane data from box and view.
///When view_num=2(x,y):(u,v)=(X-axis, Y-axis) where (u,v) is the direction of this plane,
///     view_num=3(y,z):(u,v)=(Y-axis, Z-axis),
///     view_num=4(x,z):(u,v)=(X-axis, Z-axis).
///The colors of the axes are set according to the kind of the axis.
///X's color=gridColors[1], Y's color=gridColors[2], Z's color=gridColors[3].
///The other line color is gridColor[0].
void MGConstructionPlane::setGridDataByBox(
	const MGBox& box,
	int view_num,	///<Standard view number:
	///< =1: 3D perspective view, 2:(x,y) 2D view, 3:(y,z) 2D view, 4:(z,x) 2D view
	///< =0: non standard view.
	const MGColor* gridColors,
	double spanIn///< span length of the grid. If spanIn<=0., default value from the box is used.
){
	double span;
	int lnum;
	MGPosition mid;

	int sdid;//maximum area coordinate pair will be output.
				//0:(x,y), 1:(y,z), 2:(z,x)
	MGcplane_parameter(box*INITIAL_SCALE,span,lnum,sdid,mid);
	bool view_is_fixed=(2<=view_num && view_num<=4);
	if(view_is_fixed)
		sdid=int(view_num-2);
	if(spanIn>0.)
		span=spanIn;
	set_span(span);
	set_num(lnum);
	MGVector uderi(0.,0.,0.),vderi(0.,0.,0.);
	int sdidp1=(sdid+1)%3;
	uderi(sdid)=1.; vderi(sdidp1)=1.;
	m_plane=MGPlane(uderi,vderi,mid);
	if(gridColors){
		set_colorsByViewID(view_num,gridColors);
	}
	mgVBO::setDirty(true);
}

//Bind the input point uv(this construction plane's parameter value)
//to the nearest grid point of this construction plane.
void MGConstructionPlane::bind_to_grid(
	const MGPosition& uv, MGPosition& uvout
)const{
	uvout.resize(2);
	double x=double(int(fabs(uv[0])/m_uspan+.5))*m_uspan;
	if(uv[0]<0.) uvout(0)=-x; else uvout(0)=x;
	double y=double(int(fabs(uv[1])/m_vspan+.5))*m_vspan;
	if(uv[1]<0.) uvout(1)=-y; else uvout(1)=y;
}

///Change origin point.
void MGConstructionPlane::change_origin(
	const MGPosition& new_origin
){
	m_plane.change_root_point(new_origin);
	mgVBO::setDirty(true);
}

//Convert cplane coordinates to the normal world coordinates.
//Function's return is the world coordinates.
MGPosition MGConstructionPlane::convert_to_world(const MGPosition& cplane_coord)const{
	MGPosition world(m_plane.root_point());
	world+=m_uspan*cplane_coord[0]*m_plane.u_deriv();
	world+=m_vspan*cplane_coord[1]*m_plane.v_deriv();
	world+=m_nspan*cplane_coord[2]*m_plane.normal();
	return world;
}

//Convert to cplane coordinates from the normal world coordinates.
//Function's return is the cplane coordinates.
MGPosition MGConstructionPlane::convert_from_world(const MGPosition& world_coord)const{
	MGPosition cplaneP(world_coord-m_plane.root_point());
	cplaneP(0)/=m_uspan;
	cplaneP(1)/=m_vspan;
	cplaneP(2)/=m_nspan;
	return cplaneP;
}

//Get line and axis colors
void MGConstructionPlane::get_colors(
	MGColor& lineColor,		//Grid line color
	MGColor& uaxisColor,	//u axis
	MGColor& vaxisColor		//v axis
)const{
	lineColor=m_lineColor;
	uaxisColor=m_uaxisColor;
	vaxisColor=m_vaxisColor;
}

//locate a point on this plane, given straight line.
//the located point will be the intersection(or nearest bound) point
//of sl and the plane.
//Function's return value is the world coordinate located.
MGPosition MGConstructionPlane::locate(
	const MGStraight& sl,//input the ray straight line of the cursor.
	MGPosition& uv		//the plane's parameter value (u,v) will be output.
)const{
	MGCSisect is;
	sl.relation(plane(), is);
	if(is_bind_to_grid()){
		bind_to_grid(is.param_surface(),uv);
	}else
		uv=is.param_surface();
	return plane().eval(uv);
}

//Draw this plane using OpenGL.
void MGConstructionPlane::make_display_list(MGCL::VIEWMODE vmode){
	if(disabled()||!valid())
		return;

	initializeVBO();
	const MGPosition& origin=m_plane.root_point();
	const double* origind=origin.data();
	MGVector udire=m_plane.u_deriv()*m_uspan;
	MGVector vdire=m_plane.v_deriv()*m_vspan;
	MGPosition u0(origin-udire*m_unum), u1(origin+udire*m_unum);
	MGPosition v0(origin-vdire*m_vnum), v1(origin+vdire*m_vnum);

	const int UNITNUM=10;
	//int i,j;
	MGPosition vS1(v0), vE1(v1), vS2(v0), vE2(v1);
	MGPosition uS1(u0), uE1(u1), uS2(u0), uE2(u1);
	LineWidth(1.f);//glLineWidth(1.f);
	setStaticAttribColor(m_lineColor);//m_lineColor.exec();
	Begin(GL_LINES);
		int i=std::min(UNITNUM,m_unum);
		int im1=i-1;
		do{
			for(int i1=1; i1<=im1; i1++){
				vS1+=udire; vE1+=udire;
				Vertex3dv(vS1.data());Vertex3dv(vE1.data());
				vS2-=udire; vE2-=udire;
				Vertex3dv(vS2.data());Vertex3dv(vE2.data());
			}
			vS1+=udire; vE1+=udire;
			vS2-=udire; vE2-=udire;
			i+=UNITNUM;
		}while(i<=m_unum);

		int j= std::min(UNITNUM,m_vnum);
		int jm1=j-1;
		do{
			for(int j1=1; j1<=jm1; j1++){
				uS1+=vdire; uE1+=vdire;
				Vertex3dv(uS1.data());Vertex3dv(uE1.data());
				uS2-=vdire; uE2-=vdire;
				Vertex3dv(uS2.data());Vertex3dv(uE2.data());
			}
			uS1+=vdire; uE1+=vdire;
			uS2-=vdire; uE2-=vdire;
			j+=UNITNUM;
		}while(j<=m_vnum);
	End();

	LineWidth(2.f);//glLineWidth(2.f);
	////u minus axis and v minus axis.
	Begin(GL_LINE_STRIP);
		Vertex3dv(u0.data());Vertex3dv(origind);Vertex3dv(v0.data());
	End();

	vS1=vS2=v0; vE1=vE2=v1;
	uS1=uS2=u0; uE1=uE2=u1;
	Begin(GL_LINES);
		MGVector udire10=udire*UNITNUM;
		for(i=UNITNUM; i<=m_unum; i+=UNITNUM){
			vS1+=udire10; vE1+=udire10;
			Vertex3dv(vS1.data());Vertex3dv(vE1.data());
			vS2-=udire10; vE2-=udire10;
			Vertex3dv(vS2.data());Vertex3dv(vE2.data());
		}
		MGVector vdire10=vdire*UNITNUM;
		for(j=UNITNUM; j<=m_vnum; j+=UNITNUM){
			uS1+=vdire10; uE1+=vdire10;
			Vertex3dv(uS1.data());Vertex3dv(uE1.data());
			uS2-=vdire10; uE2-=vdire10;
			Vertex3dv(uS2.data());Vertex3dv(uE2.data());
		}
	End();

	LineWidth(2.f);//glLineWidth(2.f);
	setStaticAttribColor(m_uaxisColor);//m_uaxisColor.exec();////u plus axis.
	Begin(GL_LINES);
		Vertex3dv(origind);Vertex3dv(u1.data());
	End();
	setStaticAttribColor(m_vaxisColor);//m_vaxisColor.exec();////v plus axis.
	Begin(GL_LINES);
		Vertex3dv(origind);Vertex3dv(v1.data());
	End();

	setDirty(false);
}

//set the colors to default ones.
void MGConstructionPlane::set_colors(
	const MGColor colors[3]///< [0]=line, [1]=u-axis, [2]=vaxis)
){
	m_lineColor=colors[0];
	m_uaxisColor=colors[1];
	m_vaxisColor=colors[2];
	setDirty(true);
}

///set colors by VID. colors[0] is the one of grid lines.
///colors[i] are the ones for axis lines for 1<=i<=3.
///As a standard one, colors[1]=x-axis, [2]=y-axis, [3]=z-axis.
///When 2<=vid<=4, (u,v) colors are:
///vid=2:(u,v)=(x,y), vid=3:(u,v)=(y,z), vid=4:(u,v)=(z,x).
///Other vid is treated as vid=2.
void MGConstructionPlane::set_colorsByViewID(
	int vid,
	const MGColor colors[4]
){
	if(vid<=1 || 4<vid)
		vid=2;
	int sdid=vid-1;
	int sdidp1=sdid%3+1;
	m_lineColor=colors[0];
	m_uaxisColor=colors[sdid];
	m_vaxisColor=colors[sdidp1];
	setDirty(true);
}

//set grid line color.
void MGConstructionPlane::set_line_color(const MGColor& color){
	m_lineColor=color;
	setDirty(true);
}

//set uaxis color.
void MGConstructionPlane::set_uaxis_color(const MGColor& color){
	m_uaxisColor=color;
	setDirty(true);
}

void MGConstructionPlane::set_span(double span){
	m_uspan=m_vspan=m_nspan=span;
	setDirty(true);
}
void MGConstructionPlane::set_uspan(double span){
	m_uspan=span;
	setDirty(true);
}
void MGConstructionPlane::set_vspan(double span){
	m_vspan=span;
	setDirty(true);
}
void MGConstructionPlane::set_num(int line_num){
	m_vnum=m_unum=line_num;
	setDirty(true);
}
void MGConstructionPlane::set_unum(int unum){
	m_unum=unum;
	setDirty(true);
}
void MGConstructionPlane::set_vnum(int vnum){
	m_vnum=vnum;
	setDirty(true);
}

//set vaxis color.
void MGConstructionPlane::set_vaxis_color(const MGColor& color){
	m_vaxisColor=color;
	setDirty(true);
}

///Import grid data from a MGContext, which are
///  1)(u,v) grid numbers, 2)(u,v) grid spans. Colors are not imported.
void MGConstructionPlane::importGridAttrib(const MGContext& ctx){
	// グリッド密度の設定し直し
	const int* gridNum = ctx.gridNum();
	set_unum(gridNum[0]);
	set_vnum(gridNum[1]);

	// 細グリッド線間隔を設定
	const double* gridSpan = ctx.gridSpan();
	set_uspan(gridSpan[0]);
	set_vspan(gridSpan[1]);
}

#define MESH_NUM 500.
//Compute cplane parameter from 3D box.
void MGcplane_parameter(
	const MGBox& box,
	double& span,	//span length will be output.
	int& lnum,	//number of lines along vertical and horizontal will be output.
	int& sdid,	//maxmum area coordinate pair will be output.
					//0:(x,y), 1:(y,z), 2:(z,x)
	MGPosition& mid)//rounded mid point will be output.
{
	double len[3];
	mid.resize(3);
	int i;
	for(i=0; i<3; i++){
		const MGInterval& intrvl=box[i];
		len[i]=intrvl.length();
		mid(i)=intrvl.mid_point();
	}
	double maxlen=len[0];
	if(maxlen<len[1]){
		if(len[1]<len[2]) maxlen=len[2];
		else maxlen=len[1];
	}else{
		if(maxlen<len[2]) maxlen=len[2];
	}
	span=maxlen*1.1/MESH_NUM;
	if(span>10000.) span=100000.;
	else if(span>1000.) span=10000.;
	else if(span>100.) span=1000.;
	else if(span>10.) span=100.;
	else if(span>1.) span=10.;
	else if(span>.1) span=1.;
	else if(span>.01) span=.1;
	else span=.01;
	lnum=int(MESH_NUM*.5);

	double area0=len[0]*len[1], area1=len[1]*len[2], area2=len[2]*len[0];
	sdid=0;
	if(area0<area1){
		if(area1<area2) sdid=2;
		else sdid=1;
	}else{
		if(area0<area2) sdid=2;
	}
//	double span10=span*10.;
//	double hspan10=span10*.5;
	double hspan=span*.5;
	for(i=0; i<3; i++)
		mid(i)=int((mid(i)+hspan)/span)*span;
}

// Serialization.
MGOfstream& operator<< (MGOfstream& buf, const MGConstructionPlane& cpl){
	std::bitset<32> boolData;
	boolData.set(0,cpl.m_disabled);
	boolData.set(1,cpl.m_bind_to_grid);

	buf << boolData.to_ulong();
	buf.WritePointer(&(cpl.m_plane));
	buf << cpl.m_vspan;
	buf << cpl.m_uspan;
	buf << cpl.m_vnum;
	buf << cpl.m_unum;

	cpl.m_lineColor.WriteMembers(buf);
	cpl.m_uaxisColor.WriteMembers(buf);
	cpl.m_vaxisColor.WriteMembers(buf);
	return buf;
}
MGIfstream& operator>> (MGIfstream& buf, MGConstructionPlane& cpl){
	long boollong;
	buf >> boollong;
	std::bitset<32> boolData(boollong);
	cpl.m_disabled=boolData[0];
	cpl.m_bind_to_grid=boolData[1];

	MGGel* gel=	buf.ReadPointer();
	MGPlane* pl=dynamic_cast<MGPlane*>(gel);
	cpl.m_plane=*pl; delete pl;
	buf >> cpl.m_vspan;
	buf >> cpl.m_uspan;
	buf >> cpl.m_vnum;
	buf >> cpl.m_unum;

	cpl.m_lineColor.ReadMembers(buf);
	cpl.m_uaxisColor.ReadMembers(buf);
	cpl.m_vaxisColor.ReadMembers(buf);
	return buf;
}
	
//Debug Function.
std::ostream& operator<< (std::ostream& out, const MGConstructionPlane& pln){
	out<<"MGConstructionPlane="<<&pln<<",m_disabled="<<pln.m_disabled;
	out<<",m_bind_to_grid="<<pln.m_bind_to_grid<<std::endl;
	out<<",m_plane="<<pln.m_plane<<",span=("<<pln.m_uspan<<","<<pln.m_vspan<<")";
	out<<",num=("<<pln.m_unum<<","<<pln.m_vnum<<")";
	return out;
}
