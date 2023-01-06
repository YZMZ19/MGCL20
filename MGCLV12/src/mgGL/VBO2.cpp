/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/

#include "StdAfx.h"
#include "mg/Tolerance.h"
#include "mg/AttribedGel.h"
#include "mg/DNameControl.h"
#include "mg/BPointSeq.h"
#include "mg/Curve.h"
#include "mg/LBRep.h"
#include "mg/SPointSeq.h"
#include "mg/Surface.h"
#include "mg/Plane.h"
#include "mg/MGStl.h"
#include "topo/Complex.h"
#include "topo/Cell.h"
#include "topo/Edge.h"
#include "topo/Loop.h"
#include "topo/Face.h"
#include "tl2/TL2Triangles.h"
#include "mgGL/glslprogram.h"
#include "mgGL/VBO.h"
#ifndef _CONSOLE
#include "mgGL/MGStringWriter.h"
#endif

/////////////////////////////////////////////////////////////////////////////
/// Vertex Buffer Object Class.

///mgVBOに対して描画データ作成後の情報を保持するためのクラス。
///mgVBOはmgVBOLeafとmgVBO自身をelementとして保持する。
///mgVBOをelementとするのはMGGroupまたはMGShellのメンバーデータの表示情報を
///保持するため

//Draw an arrow symbol with implementation of OpenGL.
//data[0] is the origin of the arrow, data[1] is the top of the arrow,
//data[2], [3] are two bottoms of arrowhead.
void mgVBO::drawArrow(const MGPosition pos[4]){
	Begin(GL_LINES);
		Vertex3dv(pos[0].data());
		Vertex3dv(pos[1].data());
		Vertex3dv(pos[1].data());
		Vertex3dv(pos[2].data());
		Vertex3dv(pos[1].data());
		Vertex3dv(pos[3].data());
	End();
}

/// Draw an object of class MGBox, by wireframe.
void mgVBO::drawBox(const MGBox& box){
	int sd = box.sdim();
	const MGInterval& x = box[0];
	const MGInterval& y = box[1];
	double x0=x.low_point(), x1=x.high_point();
	double y0=y.low_point(), y1=y.high_point();

	if(sd == 2){
		Begin(GL_LINE_LOOP);
			Vertex3d(x0, y0, 0.); Vertex3d(x1, y0, 0.);
			Vertex3d(x1, y1, 0.); Vertex3d(x0, y1, 0.);
		End();
		return;
	}else if(sd == 3){
		const MGInterval& z = box[2];
		double z0=z.low_point(), z1=z.high_point();
		Begin(GL_LINES);
			Vertex3d(x0, y0, z0);
			Vertex3d(x1, y0, z0);

			Vertex3d(x1, y0, z0);
			Vertex3d(x1, y1, z0);
			
			Vertex3d(x1, y1, z0);
			Vertex3d(x0, y1, z0);
			
			Vertex3d(x0, y1, z1);
			Vertex3d(x0, y0, z1);
	//
			Vertex3d(x0, y0, z1);
			Vertex3d(x1, y0, z1);

			Vertex3d(x1, y0, z1);
			Vertex3d(x1, y1, z1);
			
			Vertex3d(x1, y1, z1);
			Vertex3d(x0, y1, z1);
			
			Vertex3d(x0, y1, z1);
			Vertex3d(x0, y0, z1);
	//
			Vertex3d(x0, y0, z0);
			Vertex3d(x0, y0, z1);

			Vertex3d(x1, y0, z0);
			Vertex3d(x1, y0, z1);

			Vertex3d(x1, y1, z0);
			Vertex3d(x1, y1, z1);

			Vertex3d(x1, y0, z0);
			Vertex3d(x1, y0, z1);
		End();
	}
}

/// Draw a control points, dotted lines shall be drawn
/// between point[i-1] and point[i], for i = 1, .., length()-1.
void mgVBO::drawPointSeq(
	const MGBPointSeq& bp,
	bool draw_points		//True if points be drawn.
){
	drawPolyline(bp);
	if(draw_points){
		int n=bp.length();
		for(int i = 0; i < n; i++){
			drawPoint(bp(i));
		}
	}
}
void mgVBO::drawPointSeq(
	const MGSPointSeq& sp,
	bool draw_points		//True if points be drawn.
){
	int i, nu = sp.length_u(), nv = sp.length_v();
	for(i = 0; i < nu; i++){
		Begin(GL_LINE_STRIP);
		for(int j = 0; j < nv; j++){
			Vertex3d(sp(i,j,0),sp(i,j,1),sp(i,j,2));
		}
		End();
	}

	for(i = 0; i < nv; i++){
		Begin(GL_LINE_STRIP);
		for(int j = 0; j < nu; j++){
			Vertex3d(sp(j,i,0),sp(j,i,1),sp(j,i,2));
		}
		End();
	}

	if(draw_points){
		for(int i = 0; i < nu; i++){
			for(int j = 0; j < nv; j++){
				drawPoint(sp(i, j));
			}
		}
	}
}

//Draw 3D curve in the topology's star cell world coordinates.
///obj is a boundary of the star cell, and the curves are extracted from the
///boundary of the star cell and drawn.
void mgVBO::drawWire_in_star(const MGLoop& loop){
	loop.drawWire_in_star(*this);
}

///Draw curvature variation graph so-called Hige.
void mgVBO::drawCurvaGraph(
	const MGCurve& curve,///<The target curve.
	int density,///<Dinsity of the graph.
	bool use_radius,///<Indicates if curvature is used(=false) or curvature radius(=true). 
	double scaleRelative///<Scale of the graph. =1. is default length.
){
	assert(density > 0);

	double lengthBase = curve.curvatureLengthDisplay(use_radius);
	if(!use_radius)
		lengthBase *=-1.;
	double len = lengthBase*scaleRelative;

	MGVector v(curve.start_point());  // v : vertex of the hige's polyline
	auto f = [&](double tin){
		Vertex3dv(v.data());
		MGVector pos = curve.eval(tin);     // pos : on the curve
		MGVector T, N, B;
		double curva, torsion;
		curve.Frenet_frame(tin, T, N, B, curva, torsion);
		// N is from pos to the center of the osculating circle at pos.
		if(use_radius)
			curva = 1./curva;

		v = pos+len*curva*N;
		Vertex3dv(v.data());
		Vertex3dv(pos.data());
		Vertex3dv(v.data());
	};

	double oneOverDensity = 1. / density;
	const MGKnotVector& kv = curve.knot_vector();
Begin(GL_LINES);
	for(int i = kv.order()-1, n=kv.bdim(); i<n; i++){
		double t = kv[i];
		double dt = (kv[i+1]-t)*oneOverDensity;
		for(int j = 0; j < density; j++, t += dt)
			f(t);
	}
	f(curve.param_e());
End();
}

///Draw a point using openGL functions.
void mgVBO::drawPoint(double x,double y,double z,double size){
	if(size<=0.)
		size=mgVBOElement::getDefaultPointSize();
	float innerSize=(float)size-2.f;
	if(innerSize>=1.f)
		drawPointWithColor(x,y,z,size,innerSize);
}
void mgVBO::drawPoint(const MGPosition& pos,double size){
	drawPoint(pos[0], pos[1], pos[2],size);
}
void mgVBO::drawPointInverseColor(double x,double y,double z,double size){
	const MGColor* inner=&MGColor::get_instance(MGColor::White);
	float outer[4];
	glGetFloatv(GL_CURRENT_COLOR,outer);
	MGColor outclr(outer[0],outer[1],outer[2],outer[3]);
	if(size<=0.)
		size=mgVBOElement::getDefaultPointSize();
	float innerSize=(float)size-2.f;
	if(innerSize>=1.f)
		drawPointWithColor(x,y,z,size,innerSize,&outclr,inner);
}
void mgVBO::drawPointInverseColor(const MGPosition& pos,double size){
	drawPointInverseColor(pos[0],pos[1],pos[2],size);
}
void mgVBO::drawPointWithColor(double x,double y,double z,
	double outerSize, double innerSize,
	const MGColor* colorInner, const MGColor* colorOuter
){
	Begin(GL_POINTS);
	if(colorOuter)
		setStaticAttribColor(*colorOuter);//Outer color.
	setStaticAttribPointSize((float)outerSize);
	Vertex3d(x,y,z);
	End();

	if(innerSize>=1.f){
		const MGColor* inner=colorInner;
		if(!inner)
			inner=&MGColor::get_instance(MGColor::White);
		Begin(GL_POINTS);
		setStaticAttribColor(*inner);//Inner color.
		setStaticAttribPointSize((float)innerSize);
		Vertex3d(x,y,z);
		End();
	}
}
void mgVBO::drawPointWithColor(const MGPosition& pos,
	double outerSize, double innerSize,
	const MGColor* colorInner, const MGColor* colorOuter
){
	drawPointWithColor(pos[0], pos[1], pos[2],
		outerSize,innerSize,
		colorInner,colorOuter);
}

///draw points sequence ipos with 2 colors, inner and outer.
void mgVBO::drawPoints(
	const MGColor& boundary_color,
	const MGColor& inner_color,
	const std::vector<MGPosition>& ipos,
	double size
){
	if(ipos.size()<=0)
		return;

	if(size<=0.)
		size=mgVBOElement::getDefaultPointSize();
	std::vector<MGPosition>::const_iterator i,ibegin=ipos.begin(), iend=ipos.end();
	Begin(GL_POINTS);
		setStaticAttribColor(boundary_color);
		setStaticAttribPointSize((float)size);//
		for(i=ibegin; i!=iend; i++)
			Vertex(*i);
	End();

	Begin(GL_POINTS);
		setStaticAttribColor(inner_color);
		setStaticAttribPointSize((float)size-2.f);//
		for(i=ibegin; i!=iend; i++)
			Vertex(*i);
	End();
}

///Draw a polyline using openGL functions.
///The last argument must end with nullptr.
void mgVBO::drawOpenPolyline(const MGPosition* P0, ...){
	const MGPosition* Pi=P0;

	va_list args;
    va_start(args, P0);
	Begin(GL_LINE_STRIP);
	while(Pi){
		const MGPosition& P=*Pi;
		Vertex3d(P[0],P[1],P[2]);
		Pi = va_arg(args, const MGPosition*);
	}
	End();
    va_end(args);
}

///Draw a polyline using openGL functions.
///The last argument must end with nullptr.
void mgVBO::drawClosedPolyline(const MGPosition* P0, ...){
	const MGPosition* Pi=P0;

	va_list args;
    va_start(args, P0);
	Begin(GL_LINE_LOOP);
	while(Pi){
		const MGPosition& P=*Pi;
		Vertex3d(P[0],P[1],P[2]);
		Pi = va_arg(args, const MGPosition*);
	}
	End();
    va_end(args);
}

///Draw a polyline using openGL functions.
///When clodes=true, 1st and last points will be connected.
void mgVBO::drawPolyline(const MGBPointSeq& line, bool closed){
	int n = line.length();
	if(n<=0)
		return;

	int kind=GL_LINE_STRIP;
	if(closed) kind=GL_LINE_LOOP;
	Begin(kind);
	for(int i = 0; i < n; i++){
		Vertex3d(line(i,0),line(i,1),line(i,2));
	}
	End();
}

///Draw a polyline using openGL functions.
///When cloded=true, 1st and last points will be connected.
void mgVBO::drawPolyline(const std::vector<MGPosition>& line, bool closed){
	int kind=GL_LINE_STRIP;
	if(closed) kind=GL_LINE_LOOP;
	size_t n(line.size());
	if(n<=0)
		return;

	Begin(kind);
	for(size_t i = 0; i < n; i++){
		const MGPosition& Pi=line[i];
		Vertex3d(Pi[0],Pi[1],Pi[2]);
	}
	End();
}

///Draw a polyline using openGL functions.
///When cloded=true, 1st and last points will be connected.
void mgVBO::drawPolyline(size_t nPoints, const MGPosition line[], bool closed){
	int kind=GL_LINE_STRIP;
	if(closed) kind=GL_LINE_LOOP;
	if(nPoints<=0)
		return;

	Begin(kind);
	for(size_t i = 0; i < nPoints; i++){
		const MGPosition& Pi=line[i];
		Vertex3d(Pi[0],Pi[1],Pi[2]);
	}
	End();
}

///Draw a line from start to end.
void mgVBO::drawStraight(const MGPosition& end, const MGPosition& start){
	Begin(GL_LINE_STRIP);
	Vertex3d(start[0],start[1],start[2]);
	Vertex3d(end[0],end[1],end[2]);
	End();
}

//Draw an object in its parameter space(MGDraw_in_parameter_space).
//This is valid only for Surface, Face, Loop, Edge.
void mgVBO::drawObjInParameterSpace(const MGObject& obj){
	const MGSurface* sf=dynamic_cast<const MGSurface*>(&obj);
	if(sf){
		const MGPlane* pl=dynamic_cast<const MGPlane*>(sf);
		if(pl) return;//Plane will not be drawn.
		MGBox uv=sf->param_range();
		double u0=uv[0].low_point(), u1=uv[0].high_point();
		double v0=uv[1].low_point(), v1=uv[1].high_point();
		Begin(GL_LINE_STRIP);
		Vertex3d(u0, v0, 0.); Vertex3d(u1, v0, 0.);
		Vertex3d(u1, v1, 0.); Vertex3d(u0, v1, 0.);
		Vertex3d(u0, v0, 0.);
		End();
		drawPoint(u0, v0, 0.);
		drawPoint(u1, v0, 0.);
		drawPoint(u1, v1, 0.);
		drawPoint(u0, v1, 0.);
		return;
	}
	const MGFace* f=dynamic_cast<const MGFace*>(&obj);
	if(f){
		if(f->hasOuterBoundaryLoop()){
			const MGLoop& ol=*(f->loop(int(0)));
			ol.drawWire(*this);
		}else{
			std::vector<UniqueCurve> crvs=f->outer_boundary_param();
			size_t n=crvs.size();
			for(size_t i=0; i<n; i++)
				crvs[i]->drawWire(*this);
		}
		int i,if0;
		int nib=f->number_of_inner_boundaries(if0);
		for(i=0; i<nib; i++,if0++)
			(f->loop(if0))->drawWire(*this);

		int nlp=f->number_of_boundaries();
		for(; if0<nlp; if0++)
			(f->loop(if0))->drawWire(*this);

		for(i=0; i<nlp; i++){
			const MGLoop& lp=*(f->loop(i));
			lp.drawVertex(*this);
		}
		return;
	}
	const MGLoop* lp=dynamic_cast<const MGLoop*>(&obj);
	if(lp){
		lp->drawWire(*this);
		lp->drawVertex(*this);
		return;
	}
	const MGEdge* edg=dynamic_cast<const MGEdge*>(&obj);
	if(edg){
		edg->drawWire(*this);
		edg->drawVertex(*this);
		return;
	}
}

///Draw a rectangle.
///(param_rectangle)
void mgVBO::drawRectangle(
	const MGBox& box	//Box to draw.
){
	double u0=box[0].low_point(), u1=box[0].high_point();
	double v0=box[1].low_point(), v1=box[1].high_point();
	Begin(GL_LINE_LOOP);
		Vertex3d(u0, v0, 0.); Vertex3d(u1, v0, 0.);
		Vertex3d(u1, v1, 0.); Vertex3d(u0, v1, 0.);
	End();
}

///
void getCurvature(
	const MGSurface& surf,
	const MGPosition& uv,
	MGCL::SURFACE_CURVATURE_KIND kind,
	MGPosition& curvaData///<	(0)=curvature, (1-3)=normal, (4-6)=point data.

){
	double curvature[4];
	MGUnit_vector N;
	surf.curvatures(uv,curvature,N);
	curvaData.resize(7);
	curvaData(0)=curvature[kind];
	curvaData.store_at(1,N,0,3);
	curvaData.store_at(4,surf.eval(uv),0,3);
}

/// Renders curvatures mapping that comes into colorful image.
/// A point whose curvature is within [lower, upper], the color varies.
void mgVBO::drawSurfaceCurvature(
	const mgTL2Triangles& tld,
	MGCL::SURFACE_CURVATURE_KIND kind,
	double lower, double upper //minimum and maximum value of the curvatures of the kind.
){
	assert(tld.get_kind()==MGCL::UV);
	const MGSurface* surfP=tld.surface();
	assert(surfP);
	const MGSurface& surf=*surfP;

	double mzero=MGTolerance::mach_zero();
	if(upper-lower<=2.*mzero){
		upper+=mzero;
		lower-=mzero;
	}

	const double R[3] = {1.0, 0.0, 0.0};
	const double Y[3] = {1.0, 1.0, 0.0};
	const double G[3] = {0.0, 1.0, 0.0};
	const double C[3] = {0.0, 1.0, 1.0};
	const double B[3] = {0.0, 0.0, 1.0};

	MGBPointSeq bp(5, 3);
	bp.store_at(4, R);
	bp.store_at(3, Y);
	bp.store_at(2, G);
	bp.store_at(1, C);
	bp.store_at(0, B);
	int err = 0;
	MGLBRep color; color.buildByInterpolation(bp, 2);
	color.change_range(lower, upper);

	//const mgTLTriangles& tris=tld.triangles();
	mgTL2Triangles::const_iterator i=tld.begin(), ie=tld.end();
	for(; i!=ie; ++i){
		const mgTL2Triangle& tri=**i;
		if(tri.size()<3) continue;

		mgTL2Triangle::const_iterator j=tri.begin(), je=tri.end();
		GLenum triType=(tri.getGeometryType()== mgTESTRIANG::mgTESTRIANG_FAN)?
			GL_TRIANGLE_FAN : GL_TRIANGLE_STRIP;

		Begin(triType, MGCL::WIRE);
		for(;j!=je; ++j){
			const MGPosition& uv=*j;
			MGPosition curvaDataj;
			getCurvature(surf,uv,kind,curvaDataj);

			//Nomal.
			Normal3d(curvaDataj[1],curvaDataj[2],curvaDataj[3]);

			//Curvature color.
			double curvature = curvaDataj[0];
			if(curvature > upper)
				curvature = upper;
			else if(curvature < lower)
				curvature = lower;
			Color3dv(color.eval(curvature).data());

			//Position data.
			Vertex3d(curvaDataj[4],curvaDataj[5],curvaDataj[6]);
		}
		End(GL_FILL);
	}
}

/// Renders curvatures mapping that comes into colorful image.
/// A point whose curvature is within [lower, upper], the color varies.
void mgVBO::drawSurfaceCurvature(
	const std::vector<mgTL2Triangles>& tldvec,///<target triangulated data to draw.
	MGCL::SURFACE_CURVATURE_KIND kind,
	double lower, double upper //minimum and maximum value of the curvatures of the kind.
){
	std::vector<mgTL2Triangles>::const_iterator datai=tldvec.begin(), dataie=tldvec.end();
	for(; datai!=dataie; datai++){
		drawSurfaceCurvature(*datai,kind,lower,upper);
	}
}

/// MGStlオブジェクトを描画する
void mgVBO::drawSTL(
	const MGStl& stl, // 描画するMGStlオブジェクト
	MGCL::DRAW_TARGET target,///<When target=WIRE, built elements are
			///stored as wire mode display, else as shading mode display
	GLenum polygonMode//Polygon mode to draw, GLPOINT, GL_LINE, or GLFILL.
){
	const std::vector<MGPosition>& vertices = stl.positions();
	const std::vector<MGUnit_vector>& normals = stl.normals();

	Begin(GL_TRIANGLES,target);// 描画を開始
	size_t nTriang(normals.size()); // 三角形の個数だけループ
	int indices[3];
	for(int j = 0; size_t(j) < nTriang; j++){
		stl.GetVertIndices(j, indices);

		if(polygonMode==GL_FILL){// シェーディングを行う場合
			const double* normalP=normals[j].data();
			Normal3dv(normalP);
			Normal3dv(normalP);
			Normal3dv(normalP);
		}

		// 頂点のインデックスを元に頂点の座標を取得する
		const MGPosition& position1 = vertices[indices[0]];
		const MGPosition& position2 = vertices[indices[1]];
		const MGPosition& position3 = vertices[indices[2]];
		Vertex3dv(position1.data());
		Vertex3dv(position2.data());
		Vertex3dv(position3.data());
	}
	End(polygonMode);
}

///OpenGL shading display of a tesselated data tris.
void mgVBO::drawShade(
	const mgTL2Triangles& tris,
	MGCL::DRAW_TARGET target,///<When target=WIRE, built elements are
			///stored as wire mode display, else as shading mode display
	GLenum polygonMode//Polygon mode to draw, GLPOINT, GL_LINE, or GLFILL.
){
	mgTL2Triangles::const_iterator i=tris.begin(), ie=tris.end();
	for(; i!=ie; ++i){
		const mgTL2Triangle& tri=**i;
		if(tri.size()<3) continue;

		GLenum type=GL_TRIANGLE_STRIP;
		if(tri.getGeometryType()== mgTESTRIANG::mgTESTRIANG_FAN){
			type=GL_TRIANGLE_FAN;
		}

		Begin(type,target);
		mgTL2Triangle::const_iterator j=tri.begin(), je=tri.end();
		for(int k=0;j!=je; ++j,k++){
			const MGPosition& P = *j;
			const double* Pdata=P.data();
			int nsd=P.sdim();
			if(nsd>3){
				Normal3dv(Pdata+3);
			}
			if(nsd==2){
				Vertex2dv(Pdata);
			}else{
				Vertex3dv(Pdata);
			}
		}
		End(polygonMode);
	}
}

///OpenGL shading display of a tesselated data tris.
void mgVBO::drawShade(
	const std::vector<mgTL2Triangles>& trisVector,///<target triangulated data to draw.
	MGCL::DRAW_TARGET target,///<When target=WIRE, built elements are
			///stored as wire mode display, else as shading mode display
	GLenum polygonMode//Polygon mode to draw, GLPOINT, GL_LINE, or GLFILL.
){
	std::vector<mgTL2Triangles>::const_iterator i=trisVector.begin(), iend=trisVector.end();
	for(; i!=iend; i++)
		drawShade(*i,target,polygonMode);
}

#ifdef _CONSOLE
//World coordinates
void mgVBO::drawString(const CString& str, const MGPosition& P, const MGColor* colr) {
	std::cout << "mgVBO::drawString:str=" << str << ", P=" << P << ", colr=" << colr << std::endl;
}
void mgVBO::drawString(const std::string& str, const MGPosition& P, const MGColor* colr) {
	std::cout << "mgVBO::drawString:str=" << str << ", P=" << P << ", colr=" << colr << std::endl;
}
void mgVBO::drawString(const char* str, const MGPosition& P, const MGColor* colr) {
	std::cout << "mgVBO::drawString:str=" << *str << ", P=" << P << ", colr=" << colr << std::endl;
}
void mgVBO::drawString(const wchar_t* str, const MGPosition& P, const MGColor* colr) {
	std::cout << "mgVBO::drawString:str=" << *str << ", P=" << P << ", colr=" << colr << std::endl;
}

//Screen coordinates
void mgVBO::drawStringByScreen(const CString& str, const MGPosition& P, const MGColor* colr) {
	std::cout << "mgVBO::drawStringByScreen:str=" << str << ", P=" << P << ", colr=" << colr << std::endl;
}
void mgVBO::drawStringByScreen(const std::string& str, const MGPosition& P, const MGColor* colr) {
	std::cout << "mgVBO::drawStringByScreen:str=" << str << ", P=" << P << ", colr=" << colr << std::endl;
}
void mgVBO::drawStringByScreen(const char* str, const MGPosition& P, const MGColor* colr) {
	std::cout << "mgVBO::drawStringByScreen:str=" << str << ", P=" << P << ", colr=" << colr << std::endl;
}
void mgVBO::drawStringByScreen(const wchar_t* str, const MGPosition& P, const MGColor* colr) {
	std::cout << "mgVBO::drawStringByScreen:str=" << str << ", P=" << P << ", colr=" << colr << std::endl;
}

#else //_CONSOLE

//World coordinates
void mgVBO::drawString(const CString & str, const MGPosition & P, const MGColor * colr){
	push_back_element(MGStringWriter::Draw(str, P,colr));
}
void mgVBO::drawString(const std::string& str, const MGPosition & P, const MGColor * colr){
	push_back_element(MGStringWriter::Draw(str.c_str(), P, colr));
}
void mgVBO::drawString(const char * str, const MGPosition & P, const MGColor * colr){
	push_back_element(MGStringWriter::Draw(str, P, colr));
}
void mgVBO::drawString(const wchar_t * str, const MGPosition & P, const MGColor * colr){
	push_back_element(MGStringWriter::Draw(str, P, colr));
}

//Screen coordinates
void mgVBO::drawStringByScreen(const CString & str, const MGPosition & P, const MGColor * colr){
	push_back_element(MGStringWriter::DrawByScreen(str, P, colr));
}
void mgVBO::drawStringByScreen(const std::string& str, const MGPosition & P, const MGColor * colr){
	push_back_element(MGStringWriter::DrawByScreen(str.c_str(), P, colr));
}
void mgVBO::drawStringByScreen(const char * str, const MGPosition & P, const MGColor * colr){
	push_back_element(MGStringWriter::DrawByScreen(str, P, colr));
}
void mgVBO::drawStringByScreen(const wchar_t * str, const MGPosition & P, const MGColor * colr){
	push_back_element(MGStringWriter::DrawByScreen(str, P, colr));
}
#endif //_CONSOLE