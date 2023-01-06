/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#include "StdAfx.h"
#include "cskernel/bler.h"
#include "cskernel/bpval2.h"
#include "cskernel/bk2fli.h"
#include "cskernel/Blcbpn.h"
#include "mg/Point.h"
#include "mg/Ellipse.h"
#include "mg/Straight.h"
#include "mg/Straight.h"
#include "mg/LBRep.h"
#include "mg/RLBRep.h"
#include "mg/TrimmedCurve.h"
#include "mg/CompositeCurve.h"
#include "mg/SurfCurve.h"
#include "mg/Plane.h"
#include "mg/MGStl.h"
#include "mgGL/VBO.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//Implements the drawWire functions of all the classes.

//Draw 3D point(vertex) in world coordinates.
//The object is converted to point(s) and is drawn.
//This is valid only for topology objects or MGPoint.
void MGPoint::drawVertex(
	mgVBO& vbo
) const{vbo.drawPoint(ref(0), ref(1), ref(2));}

void MGPoint::drawWire(
	mgVBO& vbo,
	int line_density	//line density to draw a surface in wire mode.
)const{
	vbo.drawPoint(ref(0), ref(1), ref(2));
}

#define MINIMUMANGLE 0.392
void MGEllipse::drawSE(
	mgVBO& vbo,
	double t0,			//Start parameter value of the curve.
	double t1			//End parameter value of the curve.
						//Draw will be performed from t0 to t1.
)const{
	t0=gp_to_radian(t0); t1=gp_to_radian(t1);
	if(t1<t0){double save=t1; t1=t0; t0=save;}

	const MGDrawParam& para=mgVBOElement::getDrawParam();
	double span_length=para.span_length_wire();

	double r=major_len();
	double whole_span=t1-t0;
	int m=int(r*whole_span/span_length);m++;
	int mmin=int(10.*whole_span/mgPAI)+4;
	if(m<mmin)
		m=mmin;
	double t=t0, dt=whole_span/double(m);

	vbo.Begin(GL_LINE_STRIP);

	MGVector V;
	for(int i=0; i<m; i++){
		V=eval_in_radian2(t);
		vbo.Vertex3d(V[0], V[1], V[2]);
		t+=dt;
	}
	V=eval_in_radian2(t1);
	vbo.Vertex3d(V[0], V[1], V[2]);

	vbo.End();
}

//////////////////////////////////////////////

void MGStraight::drawSE(
	mgVBO& vbo,
	double t0,			//Start parameter value of the curve.
	double t1			//End parameter value of the curve.
						//Draw will be performed from t0 to t1.
)const{
	MGPosition Vs=eval(t0);
	MGPosition Ve=eval(t1);
	vbo.drawStraight(Ve,Vs);
}

#define INFINITE_LINE_LENGTH 50.
void MGStraight::drawWire(
	mgVBO& vbo,
	int line_density	//line density to draw a surface in wire mode.
) const{
	double t0=param_s(),t1=param_e();
	if(infinite_above()) t1=INFINITE_LINE_LENGTH;
	if(infinite_below()) t0=-INFINITE_LINE_LENGTH;

	double len=direction_len()*INFINITE_LINE_LENGTH/20.;
	MGUnit_vector X(direction()); MGVector Y, Z;
	X.orthonormal(X,Y,Z);
	MGVector Y2(Y*len);
	MGVector P0=eval(t0), P1=eval(t1), V;

	vbo.Begin(GL_LINE_STRIP);
	if(infinite_below()){
		V=P0+Y2;
		vbo.Vertex3d(V[0], V[1], V[2]);
	}
	vbo.Vertex3d(P0[0], P0[1], P0[2]);
	vbo.Vertex3d(P1[0], P1[1], P1[2]);
	if(infinite_above()){
		V=P1-Y2-X*len;
		vbo.Vertex3d(V[0], V[1], V[2]);
	}
	vbo.End();
}

//////////////////////////////////////////////

void MGLBRep::drawSE(
	mgVBO& vbo,
	double t0,			//Start parameter value of the curve.
	double t1			//End parameter value of the curve.
						//Draw will be performed from t0 to t1.
)const{
	int k=order(), n=bdim();
	if(!n)
		return;
	if(t1<t0){
		double save=t1; t1=t0; t0=save;
	}

	vbo.Begin(GL_LINE_STRIP);
	if(k==2){
		const MGKnotVector& t=knot_vector();
		if(t0<t[1]) t0=t[1];
		if(t1>t[n]) t1=t[n];
		MGVector P=eval(t0);
		vbo.Vertex3d(P[0],P[1],P[2]);
		int i0=t.locate(t0), i1=t.locate(t1);
		const MGBPointSeq& bp=line_bcoef();
		for(int i=i0; i<i1; i++)
			vbo.Vertex3d(bp(i,0),bp(i,1),bp(i,2));
		P=eval(t1);
		vbo.Vertex3d(P[0],P[1],P[2]);
	}else{
		drawgl(vbo,t0, t1);
	}
	vbo.End();
}

//////////////////////////////////////////////
void MGRLBRep::drawSE(
	mgVBO& vbo,
	double t0,			//Start parameter value of the curve.
	double t1			//End parameter value of the curve.
						//Draw will be performed from t0 to t1.
)const{
	vbo.Begin(GL_LINE_STRIP);
	drawgl(vbo, t0, t1);
	vbo.End();
}

//////////////////////////////////////////////
void MGTrimmedCurve::drawSE(
	mgVBO& vbo,
	double t0,			//Start parameter value of the curve.
	double t1			//End parameter value of the curve.
						//Draw will be performed from t0 to t1.
)const{
	double ts=param_s(), te=param_e();
	if(t0<ts) t0=ts;
	if(t1>te) t1=te;
	m_curve->drawSE(vbo,t0,t1);
}

//////////////////////////////////////////////
void MGCompositeCurve::drawSE(
	mgVBO& vbo,
	double t0,			//Start parameter value of the curve.
	double t1			//End parameter value of the curve.
						//Draw will be performed from t0 to t1.
)const{
	if(m_composite.size()==0) return;

	double s0,s1;
	if(t0>t1){ s0=t1; s1=t0;}
	else{ s0=t0; s1=t1;}
	int i0=find(s0), i1=find(s1);
	if(i0==i1){
		m_composite[i0]->drawSE(vbo,s0,s1);
	}else{
		m_composite[i0]->drawSE(vbo,s0,m_composite[i0]->param_e());
		for(int i=i0+1; i<i1; i++)
			m_composite[i]->drawSE(vbo,
			m_composite[i]->param_s(),m_composite[i]->param_e());
		m_composite[i1]->drawSE(vbo,
			m_composite[i1]->param_s(), s1);
	}
}

#define NDIVIDE 3
#define NMAX 150
void MGCurve::drawSE(
	mgVBO& vbo,
	double t0,			//Start parameter value of the curve.
	double t1			//End parameter value of the curve.
						//Draw will be performed from t0 to t1.
)const{
	int i;
	double s0,s1,s[NDIVIDE+1];
	if(t0>t1){ s0=t1; s1=t0;}
	else{ s0=t0; s1=t1;}
	s[0]=range(s0); s[NDIVIDE]=range(s1);
	double ds=(s[NDIVIDE]-s[0])/NDIVIDE;
	for(i=1; i<NDIVIDE; i++) s[i]=s[i-1]+ds;
	vbo.Begin(GL_LINE_STRIP);
	MGVector P=eval(s0);
	vbo.Vertex3d(P[0], P[1], P[2]);

	MGVector dfdt=eval(s[0]); double vlen1=dfdt.len();
	for(int m=0; m<NDIVIDE; m++){
	
	s0=s[m]; s1=s[m+1];
	double s2=(s0+s1)*.5;
	double vlen0=vlen1;
	dfdt=eval(s2,1); double vlen2=dfdt.len();
	dfdt=eval(s1,1); vlen1=dfdt.len();
	double vlen;
	if(vlen0<=vlen1){
		if(vlen1<=vlen2) vlen=vlen0+vlen2;
		else if(vlen0>vlen2) vlen=vlen1+vlen2;
		else vlen=vlen1+vlen0;
	}else{
		if(vlen0<=vlen2) vlen=vlen2+vlen1;
		else if(vlen1>vlen2) vlen=vlen0+vlen2;
		else vlen=vlen0+vlen1;
	}
	
	const MGDrawParam& para=mgVBOElement::getDrawParam();
	double span_length=para.span_length_wire();

	vlen*=.5;
	double span=(s1-s0);
	double df=vlen*span;
	int n=int(df/span_length); n+=2;
	if(n>NMAX) n=NMAX;
	double dt=span/double(n);
	double t=s0;
	for(int i=1; i<n; i++){
		t+=dt;
		P=eval(t);
		vbo.Vertex3d(P[0], P[1], P[2]);
	}
	P=eval(s1);
	vbo.Vertex3d(P[0], P[1], P[2]);

	}
	vbo.End();
}

//Draw this by converting straight line segments.
void MGLBRep::drawgl(
	mgVBO& vbo,
	double tstart, double tend	//start and end parameter value of this.
)const{
	const int c0=0, c1=1;
	const int k=order(), n=bdim(), irc=line_bcoef().capacity();
	const double* t=knot_data();
	const double* rcoef=coef_data();

    double x[3]={0.,0.,0.}, dx[3];

// ***** INITIAL SET. 
	int sd=sdim(); if(sd>3) sd=3;
	double t0=t[k-1], t1=t[n];
	if(tstart>=tend){
		tstart=t0; tend=t1;
	}else{
		if(tstart<t0) tstart=t0;
		if(tend>t1) tend=t1;
	}
	if(t0>=t1)
		return;

// ******* CONVERT TO PP AND DRAW EACH LINE *****
	int np1=n+1, is1,is2;
	is1=bk2fli_(np1, t, tstart);
	is2=bk2fli_(np1, t, tend);	while(t[is2-1]==tend) is2--;

// ***** DRAW LINE FROM RW(I) TO RW(I+1) 
//   *** MOVE TO THE FIRST POSITION 
	int i;
	for(i=0;i<sd;i++) x[i]=bler_(k,n,t,rcoef+i*irc,tstart,c0);
	vbo.Vertex3d(x[0],x[1],x[2]);

//   *** DRAW LINE ONE KNOT SPAN BY ONE
	const MGDrawParam& para=mgVBOElement::getDrawParam();
	double span_length=para.span_length_wire();
	double *wk1=new double[3*k*k+3*k];
	double* wk1p3k=wk1+3*k;
	for(int j=is1; j <= is2; ++j){
		double ts, te;
		double tnow=ts=t[j-1]; te=t[j]; if(ts>=te) continue;
		if(j==is1) ts=tnow=tstart;
		if(j==is2) te=tend;
//     ...CONVERT THE ONE SPAN TO PP-REP
		int jmk=j-k; int lpp;
		double brk[2];
		blcbpn_(k,k,t+jmk,rcoef+jmk,irc,sd,c1,wk1p3k,brk,wk1,&lpp);
		double tm=(ts+te)*0.5;
		double vlen=0.;
		for(i=0;i<sd;i++){
			dx[i]=bpval2_(brk,wk1+i*k,k,tm,c1,c1);
			vlen+=dx[i]*dx[i];
		}
	    vlen = sqrt(vlen);//length of 1st deriv at the middle point of the span.
		double span=te-ts;
	    int mj = (int)(vlen*span/span_length) + 1;
		if(mj<4) mj=4;
		double dt=span/(double)mj;

		//Draw middle points between the knots BY INCREMENTING DT.
		for(int l=1; l<mj; l++){
			tnow+=dt;
			for(i=0;i<sd;i++) x[i]=bpval2_(brk,wk1+i*k,k,tnow,c1,c0);
			vbo.Vertex3d(x[0],x[1],x[2]);
		}
		
		//Draw to the end point of the knot.
		for(i=0;i<sd;i++) x[i]=bpval2_(brk,wk1+i*k,k,te,c1,c0);
		vbo.Vertex3d(x[0],x[1],x[2]);

	}
	delete[] wk1;
}

//Draw this by converting straight line segments.
void MGRLBRep::drawgl(
	mgVBO& vbo,
	double tstart, double tend	//start and end parameter value of this.
)const{
	const int c0=0, c1=1;
	const int k=order(), n=bdim(), irc=line_bcoef().capacity();
	const double* t=knot_data();
	const double* rcoef=coef_data();
    double x[3]={0.,0.,0.}, dx[3];

// ***** INITIAL SET.
	double t0=t[k-1], t1=t[n];
	if(tstart>=tend){
		tstart=t0; tend=t1;
	}else{
		if(tstart<t0) tstart=t0;
		if(tend>t1) tend=t1;
	}

// ******* CONVERT TO PP AND DRAW EACH LINE *****
	int is1, is2;
	int np1=n+1;
	is1=bk2fli_(np1, t, tstart);
	is2=bk2fli_(np1, t, tend);
	while(t[is2-1]==tend) is2--;

    int i;
	int ncd=sdim();
	int sd=ncd; if(sd>3) sd=3;

	// ***** DRAW LINE FROM RW(I) TO RW(I+1)
//   *** MOVE TO THE FIRST POSITION
	double w=bler_(k,n,t,rcoef+ncd*irc,tstart,c0);
	for(i=0;i<sd;i++) x[i]=bler_(k,n,t,rcoef+i*irc,tstart,c0)/w;
	vbo.Vertex3d(x[0],x[1],x[2]);

//   *** DRAW LINE ONE KNOT SPAN BY ONE
	const MGDrawParam& para=mgVBOElement::getDrawParam();
	double span_length=para.span_length_wire();
	int ncdp1=ncd+1;
	double*  pcoef=new double[ncdp1*(k*k+k)];
	double* wpcoef=pcoef+k*ncd;//pcoef for weight of rational b-coef.
	double* scratch=wpcoef+k;//work area.
	for(int j=is1; j<=is2; ++j){
		double ts,te;
		double tnow=ts=t[j-1]; te=t[j]; if(ts>=te) continue;
		if(j==is1) ts=tnow=tstart;
		if(j==is2) te=tend;
//     ...CONVERT THE ONE SPAN TO PP-REP
		int jmk=j-k; int lpp;
		double brk[2];
	    blcbpn_(k,k,&t[jmk],rcoef+jmk,irc,ncdp1,c1,scratch,brk,pcoef,&lpp);
		double tm=(ts+te)*0.5;
		double vlen=0.;
		for(i=0;i<sd;i++){
			x[i]=bpval2_(brk,pcoef+i*k,k,tm,c1,c0);
			dx[i]=bpval2_(brk,pcoef+i*k,k,tm,c1,c1);
		}
		w=bpval2_(brk,wpcoef,k,tm,c1,c0);
		double dw=bpval2_(brk,wpcoef,k,tm,c1,c1);
		vlen=0.;
		for(i=0;i<sd;i++){
			double vln=(dx[i]-dw*x[i]/w)/w; vlen+=vln*vln;
		}
	    vlen = sqrt(vlen);//length of 1st deriv at the middle point of the span.
		double span=te-ts;
	    int mj = int(vlen*span/span_length) + 1;
		if(mj<4) mj=4;
		double dt=span/double(mj);

		//Draw middle points between the knots BY INCREMENTING DT
		for(int l=1; l<mj; l++){
			tnow+=dt;
			w=bpval2_(brk,wpcoef,k,tnow,c1,c0);
			for(i=0;i<sd;i++) x[i]=bpval2_(brk,pcoef+i*k,k,tnow,c1,c0)/w;
			vbo.Vertex3d(x[0],x[1],x[2]);
		}

		//Draw to the end point of the knot.
		w=bpval2_(brk,wpcoef,k,te,c1,c0);
		for(i=0;i<sd;i++) x[i]=bpval2_(brk,pcoef+i*k,k,te,c1,c0)/w;
		vbo.Vertex3d(x[0],x[1],x[2]);
	}
    delete[] pcoef;
}

void MGPlane::get_uv_display_vector(
	MGVector& u,
	MGVector& v
)const{
	const double L = 10.;
	u = u_deriv(); u *= L / u.len();
	v = v_deriv(); v *= L / v.len();
}

void MGPlane::drawWirePlane(
	mgVBO& vbo,
	int line_density,	//line density to draw a surface in wire mode.
	MGCL::DRAW_TARGET target
)const{

struct DrawCross{
	DrawCross(mgVBO* vbo,const MGVector& u, const MGVector& v,
		MGCL::DRAW_TARGET target) :
		m_vbo(vbo),m_u(u), m_v(v),m_target(target){}

	void operator()(const MGVector& start, const MGVector& end,	bool is_u) const{
		MGVector v1, v2;
		if(is_u){
			v1 = start - m_u + m_v;
			v2 = start - m_u - m_v;
		}else{
			v1 = start - m_u - m_v;
			v2 = start + m_u - m_v;
		}

		m_vbo->Begin(GL_LINE_STRIP,m_target);
			m_vbo->Vertex3d(v1[0], v1[1], v1[2]);
			m_vbo->Vertex3d(start[0], start[1], start[2]);
			m_vbo->Vertex3d(v2[0], v2[1], v2[2]);
		m_vbo->End();

		m_vbo->Begin(GL_LINE_STRIP,m_target);
			m_vbo->Vertex3d(start[0], start[1], start[2]);
			m_vbo->Vertex3d(end[0], end[1], end[2]);
		m_vbo->End();

		MGVector tmp(end - start);
		v1 += tmp;
		v2 += tmp;
		m_vbo->Begin(GL_LINE_STRIP,m_target);
			m_vbo->Vertex3d(v1[0], v1[1], v1[2]);
			m_vbo->Vertex3d(end[0], end[1], end[2]);
			m_vbo->Vertex3d(v2[0], v2[1], v2[2]);
		m_vbo->End();
	};

	mgVBO* m_vbo;
	const MGVector& m_u;
	const MGVector& m_v;
	MGCL::DRAW_TARGET m_target;
};

	MGVector u,v;
	get_uv_display_vector(u,v);

	const MGPosition& cen = center();
	vbo.Begin(GL_LINE_STRIP,target);
		MGVector V(cen + u + v);  // (L, L).
		vbo.Vertex3d(V[0], V[1], V[2]);
		V -= 2 * u;
		vbo.Vertex3d(V[0], V[1], V[2]);  // (-L, L)
		V -= 2 * v;
		vbo.Vertex3d(V[0], V[1], V[2]);  // (-L, -L)
		V += 2 * u;
		vbo.Vertex3d(V[0], V[1], V[2]);  // (L, -L)
		V += 2 * v;
		vbo.Vertex3d(V[0], V[1], V[2]);  // (L, L)
	vbo.End();

	const double ratio = 4.11;
	u /= ratio;
	v /= ratio;
	DrawCross da(&vbo,u, v,target);
	da(cen-3*u, cen+3*u, true);
	da(cen-3*v, cen+3*v, false);
}
