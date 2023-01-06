/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#include "StdAfx.h"
#include "mg/Box.h"
#include "mg/Matrix.h"
#include "mg/SPointSeq.h"
#include "mg/Straight.h"
#include "mg/Ellipse.h"
#include "mg/RLBRep.h"
#include "mg/RSBRep.h"
#include "mg/Tolerance.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

// Implements MGRSBRep Class
//
// Defines Rational Surface B-Representation.
// This NURBS is of homogeneous form, i.e., B-Coefficients have
// weight included values. 
// When usual NURBS form is (xij, yij, zij, wij) ,
// MGRSBRep form is (xij*wij, yij*wij, zij*wij, wij)
//				 for i=0,..., m-1, and j=0,..., n-1.

//Approximate an original B-Rep old by a new knot configuration.
//The new knot vectors are input from the member data.
//Use setKnotVector() before buildByNewKnotVectorWithKTV() to set the knotvector.
//The new knot config must be inside the range of the original B-Rep
//parameter. However new knots may be coarse or fine.
//Function's return value is Error flag.
//Error is detected only when ut (=2) or vt(=12) is illegal.
//When error!=0, the original old_brep is copied to this.
int MGRSBRep::buildByNewKnotVectorWithKTV(
	const MGRSBRep& old//Original B-Rep.
){
	invalidateBox();
	copy_appearance(old);
	return  m_surface.buildByNewKnotVectorWithKTV(old.m_surface);
}

//Gets new B-Rep by adding knots to an original B-Rep.
void MGRSBRep::addKnots(
	const MGKnotArray& uknots,	//Knots to add for u-direction
	const MGKnotArray& vknots)	//Knots to add for v-direction.
{
	m_surface.addKnots(uknots, vknots);
	invalidateBox();
}

// Gets new NURBS Surface by computing a part of the original.
//New one is exactly the same as the original except that it is partial.
//If multiple==true(!=0), knot_u(i)=t1 and knot_u(n+i)=t2 for i=0,..., k-1
//will be guaranteed. Here, n=bdim_u(), k=order_u(),
//t1=uvrange(0).low_point(), and t2=uvrange(0).high_point().
//About knot_v(j), the same.
//Both u-range and v-range must be inside the range of old.
void MGRSBRep::shrinkToParameters(
	const MGBox& uvrange,		//u and v parameter range.
	MGRSBRep& newBrep,///<Shrinked surface is output, which can be this.
	int multiple //Indicates if start and end knot multiplicities
					//are necessary. =0:unnecessary, !=0:necessary.
)const{
	m_surface.shrinkToParameters(uvrange, newBrep.m_surface, multiple);
	newBrep.copy_appearance(*this);
}

//Convert the bspline coefficients to homogeneous form.
void MGRSBRep::convertToHomogeneous(){
	int sd=sdim();
	int bdu=m_surface.bdim_u();
	int bdv=m_surface.bdim_v();

	//Multiply each weight to all control polygons.
	for(int i=0;i<bdu; i++){
		for(int j=0;j<bdv; j++){
			double weight=m_surface.coef(i,j,sd);
			for(int r=0; r<sd; r++)
				m_surface.coef(i,j,r)*=weight;
		}
	}
}

///Construct this MGRSBRep, given weights and bcoef's.
///The knot vectors must be set before use this by setKnotVeector().
///MGSPointSeq bcoef is non homogeneous form(do not include weights).
void MGRSBRep::buildRSBRepFromMemberDataWithKTV(
	const MGSPointSeq& bcoef,	
	//Control Vertex of rational surface B-Rep that does not includes weights.
	const MGSPointSeq& weights	//weights, weights(i,j,0) is for bcoef(i,j,.)
){
	assert(bcoef.length_u()==knot_vector_u().bdim()
		   && bcoef.length_v()==knot_vector_v().bdim());

	double weight;
	int dim=bcoef.sdim(), bdimu=knot_vector_u().bdim(), bdimv=knot_vector_v().bdim();
	m_surface.m_surface_bcoef.resize(bdimu,bdimv,dim+1);

	//Multiply each weight to all control polygons.
	for(int i=0;i<bdimu; i++){
		for(int j=0;j<bdimv; j++){
			weight=weights(i,j,0);
			for(int r=0; r<dim; r++)
				m_surface.coef(i,j,r)=bcoef(i,j,r)*weight;
			m_surface.coef(i,j,dim)=weight;
		}
	}
}

//Construct surface of revolution, given planar MGRLBRep.
//Parameterization of the surface is:
//	u=const parameter line generates given rlb(when u=0.).
//  v=const parameter line generates circle.
void MGRSBRep::buildRevolutionSurface(
		const MGRLBRep& rlb,	//Planar MGRLBRep to rotate.
		const MGStraight& sl,	//Rotation axis. This is treated as infinite
								//one, even if it is not.
		double angle			//Rotation angle in radian.
	//-2*pai<=angle<=2*pai.
	//If angle is positive, circle is anti-clockwise around direction Vector N
	//of sl. If negative, circle is clockwise around N.
){
	invalidateBox(); copy_appearance(rlb);
	MGUnit_vector N(sl.direction());
	MGVector V=sl.nearest_to_origin();
		//V is to translate sl to pass through origin.
	MGMatrix mat_toZ;
	mat_toZ.to_axis(N, 2);	//mat_toZ is to rotate sl to be z axis.

	MGBPointSeq cp=rlb.line_bcoef().non_homogeneous();
	cp-=V; cp*=mat_toZ;

	int n=cp.length(); int mid=n/2;
	MGVector P_sample(2,cp(0)+cp(mid)+cp(n-1));
	double ulen=P_sample.len()/3.;

	MGUnit_vector U_sample(P_sample);
	MGRLBRep circle(MGEllipse(MGPosition(0.,0.), MGPosition(U_sample), angle));
	MGBPointSeq circle_cp(circle.non_homogeneous_bcoef());
	int m=circle.bdim();
	MGSPointSeq sfcpR(m,n,4);
	MGSPointSeq sfcp(m,n,3);

	//Compute coefficient in normalized space where
	//rotation axis is z-axis.
	int i,j;
	int wid=rlb.sdim();			//Weight id of rlb.
	for(j=0; j<n; j++){
	//Get non homogeneous coef in sfcp(.,.,0-2) and weight in sfcpR(.,.,3).
		double weightv=rlb.coef(j,wid);
		double len=U_sample%cp(j);
		for(i=0; i<m; i++){
			double weightu=circle.coef(i,2);
			sfcpR(i,j,3)=weightu*weightv;
			sfcp(i,j,0)=circle_cp(i,0)*len;
			sfcp(i,j,1)=circle_cp(i,1)*len;
			sfcp(i,j,2)=cp(j,2);
		}
	}
	MGMatrix mat_toN; mat_toN.from_axis(N, 2);
	//mat_toN is to rotate z axis to be parallel to sl.

	//Get non homogeneous coef in original space.
	sfcp*=mat_toN; sfcp+=V;
	for(j=0; j<n; j++)
		for(i=0; i<m; i++) sfcpR.store_at(i,j,sfcp(i,j));
	m_surface.m_surface_bcoef=std::move(sfcpR);
	m_surface.m_uknot=circle.knot_vector()*ulen;
	m_surface.m_vknot=rlb.knot_vector();
	convertToHomogeneous();
}

//**** Conversion Constructor.****

// Convert Non ratoinal to Rational form.
MGRSBRep::MGRSBRep(
	const MGSBRep& brep,	//Original SBRep. This can be ordinary SBRep, or 
		//homogeneous form of MGRSBRep. When homogeneous form,
		//the last space dimension elements are weights.
	int homogeneous)		//true(non zero): homogeneous form,
							//false(zero):ordinary SBRep.
:MGSurface(brep), m_surface(brep){
	if(!homogeneous){
		int m=brep.bdim_u(), n=brep.bdim_v(), dimm1=brep.sdim();
		MGSPointSeq cp(m,n,dimm1+1);
		int i,j,r;
		for(i=0; i<m; i++){
			for(j=0; j<n; j++){
				for(r=0; r<dimm1; r++) cp(i,j,r)=brep.coef(i,j,r);
				cp(i,j,dimm1)=1.;
			}
		}
		m_surface.m_surface_bcoef=std::move(cp);
		m_surface.m_uknot=brep.knot_vector_u();
		m_surface.m_vknot=brep.knot_vector_v();
	}
}

// Construct a Line NURBS by changing space dimension and ordering of
//coordinates.
MGRSBRep::MGRSBRep(
	int dim,			// New space dimension.
	const MGRSBRep& rsb,// Original Surface B-rep.
	int start1, 		// Destination order of new line.
	int start2) 		// Source order of original line.
:MGSurface(rsb){
	invalidateBox();
	int dim0=rsb.sdim();
	MGSPointSeq cp1(dim0,rsb.surface_bcoef());//Exclude weights.
	MGSPointSeq cp2(dim,cp1,start1,start2);  //Change order of coordinates.
	MGSPointSeq cp3(dim+1,cp2);			     //Get area for weights.
	for(int i=0; i<cp3.length_u(); i++)
	for(int j=0; j<cp3.length_v(); j++)
		cp3(i,j,dim)=rsb.surface_bcoef()(i,j,dim0);//Set weights.
	m_surface.m_surface_bcoef=std::move(cp3);
	m_surface.m_uknot=rsb.knot_vector_u();
	m_surface.m_vknot=rsb.knot_vector_v();
}

//Member Function

//Return minimum box that includes the whole line.
void MGRSBRep::compute_box(MGBox& bx) const{
	m_surface.surface_bcoef().non_homogeneous().compute_box(bx);
}
	
	
// ���͂̃p�����[�^�͈͂̋Ȗʕ������͂ރ{�b�N�X��Ԃ��B
//Return minimum box that includes the partial line.
MGBox MGRSBRep::box_limitted(
	const MGBox& uvbox
)const{
	double u1=uvbox(0).low_point(), u2=uvbox(0).high_point();
	double us=param_s_u(), ue=param_e_u();
	if(u1<us) u1=us; if(u2>ue) u2=ue;

	double v1=uvbox(1).low_point(), v2=uvbox(1).high_point();
	double vs=param_s_v(), ve=param_e_v();
	if(v1<vs) v1=vs; if(v2>ve) v2=ve;

	if(MGREqual_base(u1,u2,knot_vector_u().param_span())){
		MGRLBRep line=parameter_line(1,(u1+u2)*.5);
		return line.box_limitted(uvbox(1));
	}else if(MGREqual_base(v1,v2,knot_vector_v().param_span())){
		MGRLBRep line=parameter_line(0,(v1+v2)*.5);
		return line.box_limitted(uvbox(0));
	}
	MGSBRep temp;
	m_surface.shrinkToParameters(uvbox,temp);
	return temp.surface_bcoef().non_homogeneous().box();
}

//Changing this object's space dimension.
void MGRSBRep::change_dimension(
	int dim,		// new space dimension
	int start1, 		// Destination order of new object.
	int start2) 		// Source order of this object.
{
	int dim0=sdim();
	MGSPointSeq cp1(dim0,surface_bcoef());	//Exclude weights.
	MGSPointSeq cp2(dim,cp1,start1,start2); //Change order of coordinates.
	MGSPointSeq cp3(dim+1,cp2);			    //Get area for weights.
	const MGSPointSeq& sp=surface_bcoef();
	for(int i=0; i<cp3.length_u(); i++)
	for(int j=0; j<cp3.length_v(); j++)
		cp3(i,j,dim)=sp(i,j,dim0);//Set weights.
	m_surface.m_surface_bcoef=std::move(cp3);
	invalidateBox();
}

//Change parameter range, be able to change the direction by providing
//t1 greater than t2.
MGRSBRep& MGRSBRep::change_range(
	int is_u,				//if true, (t1,t2) are u-value. if not, v.
	double t1,				//Parameter value for the start of original. 
	double t2)				//Parameter value for the end of original. 
{
	m_surface.change_range(is_u,t1,t2);
	invalidateBox();
	return *this;
}

//Construct new surface object by copying to newed area.
//User must delete this copied object by "delete".
MGRSBRep* MGRSBRep::clone() const{return new MGRSBRep(*this);}

//Construct new surface object by changing
//the original object's space dimension.
//User must delete this copied object by "delete".
MGRSBRep* MGRSBRep::copy_change_dimension(
	int sdim,		// new space dimension
	int start1, 		// Destination order of new line.
	int start2) 		// Source order of this line.
	const{
	return new MGRSBRep(sdim,*this,start1,start2);
}
//Evaluate right continuous ndu'th and ndv'th derivative data.
//Function's return value is (d(ndu+ndv)f(u,v))/(du**ndu*dv**ndv).
// ndu=0 and ndv=0 means positional data evaluation.
MGVector MGRSBRep::eval(
	double u, double v,	// Parameter value of the surface.
	int ndu,			// Order of Derivative along u.
	int ndv			// Order of Derivative along v.
)const{
	int m;
	int dim=sdim();
	MGVector result(dim);
	if(ndu==0 && ndv==0){
		MGVector data=m_surface.eval(u,v);
		double weight=data.ref(dim);
		for(m=0; m<dim; m++) result(m)=data.ref(m)/weight;
	}else{
		double* deriv; double deriva[27];
		int len=dim*(ndu+1)*(ndv+1);
		if(len<=27) deriv=deriva; else deriv=new double[len];
		eval_all(u,v,ndu,ndv,deriv);
		result=MGVector(dim,deriv+len-dim);
		if(len>27) delete[] deriv;
	}
	return result;
}

//Evaluate surface data.
MGVector MGRSBRep::eval(
	const MGPosition& uv	// Parameter value of the surface.
	, int ndu			// Order of derivative along u.
	, int ndv			// Order of derivative along v.
) const{	return eval(uv.ref(0),uv.ref(1),ndu,ndv);}

//Compute position, 1st and 2nd derivatives.
// �p�����[�^�l��^���Ĉʒu�A�ꎟ�����l�A�񎟔����l�����Ƃ߂�B
//Evaluate right continuous surface data.
//Evaluate all positional data and 1st and 2nd derivatives.
void MGRSBRep::eval_all(
	double u, double v,		// Parameter value of the surface.
	MGPosition& f,			// Positional data.
	MGVector&   fu,			// df(u,v)/du
	MGVector&   fv,			// df/dv
	MGVector&   fuv,		// d**2f/(du*dv)
	MGVector&   fuu,		// d**2f/(du**2)
	MGVector&   fvv			// d**2f/(dv**2)
	) const
{
	int dim=sdim(); int iid=3*dim;
	double data_area[27];
	double* data;
	if(dim<=3) data=data_area; else data=new double[dim*3*3];
	eval_all(u,v,2,2,data);
	f=MGPosition(dim,data);
	fu=MGVector(dim,data+iid);
	fv=MGVector(dim,data+dim);
	fuv=MGVector(dim,data+dim+iid);
	fuu=MGVector(dim,data+2*iid);
	fvv=MGVector(dim,data+2*dim);
	if(dim>3) delete[] data;
}

//Evaluate all of i'th derivative data for 0<=i<=nderiv.
//Output will be put on deriv[j+i*sdim()]
//for 0<=i<=nderiv and 0<=j<sdim(), i.e. 
//deriv[j+i*sdim()] is i-th derivative data for 0<=j<sdim(). 
void MGRSBRep::eval_all(
	double u, double v,		//Parameter value to evaluate.
	int ndu, int ndv,	//Order of Derivative along u and v direction.
	double* deriv	//Output. (d(i+j)f(u,v))/(du**i*dv**j) in
					//deriv[r+j*dim+i*(ndv+1)*dim] for 0<= r <dim=sdim().
					//for 0<=i<=ndu and 0<=j<=ndv.
					//deriv is an array of deriv[ndu+1][ndv+1][r].
	) const
{
	int dim=sdim(); int dimp1=dim+1;
	int mdu=ndu+1, mdv=ndv+1;
	int md=mdu; if(md<mdv) md=mdv;
	int idderi=mdv*dim; //deriv[r+j*dim+i*idderi] for deriv[i][j][r].

	//Prepare necessary binominal coefficient data.
	double *bc; double bca[9];
	if(md<=3) bc=bca; else bc=new double[md*md];
	MGCL::Binominal(md-1, bc);
	//Actually bc is an array of bc[md][md], and contains
	//binominal(i,j).

	double data_area[36];// 36=(2+1)*(2+1)*(3+1), ie ndu=2, ndv=2, dim=3.
	double *data;
	int data_len=dimp1*mdu*mdv;
	if(data_len<=36) data=data_area;
	else data=new double[data_len];
	//Actually data is 3 dimensional array of data[mdu][mdv][dimp1],
	//i.e., data[r+dimp1*j+mdv*dimp1*i]=data[i][j][r];
	int iddata=mdv*dimp1; //data[r+dimp1*j+iddata*i]=data[i][j][r];

	m_surface.eval_all(u,v,ndu,ndv,data);
	// data[r+dimp1*j+iddata*i]=data[i][j][r] is the data of
	//i-th along u and j-th along v derivative
	// of r-th space dimension element. For r=dim, wieght.

	double val, weight=data[dim];
	int r, i, j, m,n;
	for(r=0; r<dim; r++){
		for(m=0; m<=ndu; m++){
		int mmd=m*md;
		for(n=0; n<=ndv; n++){

			int nmd=n*md;
			val=data[r+dimp1*n+iddata*m];
			int id1=r+m*idderi;
			for(j=1; j<=n; j++)
				val-= bc[nmd+j]*data[dim+dimp1*j]*deriv[id1+(n-j)*dim];

			int id2=r+n*dim;
			for(i=1; i<=m; i++){
				val-=bc[mmd+i]*data[dim+iddata*i]*deriv[id2+(m-i)*idderi];
				double val2=0.;
				int id3=r+(m-i)*idderi, id4=dim+iddata*i;
				for(j=1; j<=n; j++)
					val2+=bc[nmd+j]*data[id4+dimp1*j]*deriv[id3+(n-j)*dim];
				val-=bc[mmd+i]*val2;
			}
			deriv[r+n*dim+m*idderi]=val/weight;

		}
		}
	}

	if(data_len>36) delete[] data;
	if(md>3) delete[] bc;
}

//Modify the original Surface by extrapolating the specified perimeter.
//The extrapolation is C2 continuous if the order >=4.
//The extrapolation is done so that extrapolating length is "length"
//at the position of the parameter value "param" of the perimeter.
MGRSBRep& MGRSBRep::extend(
	int perimeter,	//perimeter number of the Surface.
					// =0:v=min, =1:u=max, =2:v=max, =3:u=min.
	double param,	// parameter value of above perimeter.
	double length,	//chord length to extend at the parameter param of the perimeter.
	double dk){  //Coefficient of how curvature should vary at
//    extrapolation start point. When dk=0, curvature keeps same, i.e.
//    dK/dS=0. When dk=1, curvature becomes zero at length extrapolated point,
//    i.e. dK/dS=-K/length at extrapolation start point.
//    (S=parameter of arc length, K=Curvature at start point)
//    That is, when dk reaches to 1 from 0, curve changes to flat.

	assert(sdim()<=3);
	assert(perimeter>=0 && perimeter<4);
	
	invalidateBox();
	const int ncd=surface_bcoef().sdim();
	int at_start=1;//starting perimeter
	int nu, nv;
	int order; int n,m; MGKnotVector* t;
	if(perimeter==1 || perimeter==3){	// Extrapolate to u-direction
		order=order_u();
		n=bdim_u();
		t=&(knot_vector_u());
		if(perimeter==1)
			at_start=0;//ending perimeter
		m=nv=bdim_v();
	}else{
		// Extrapolate to v-direction
		order=order_v();
		n=bdim_v();
		t=&(knot_vector_v());
		if(perimeter==2)
			at_start=0;//ending perimeter
		m=nu=bdim_u();
	}
	//(nu,nv) are new surface B-Rep dimensions of u and v.
	//(order,n,t) is line B-rep to extrapolate.
	//m is the number of line B-reps to extrapolate.

	MGSPointSeq surf;
	MGRLBRep lbtemp;
	MGKnotVector& t1=lbtemp.knot_vector();
	MGBPointSeq& coeftemp=lbtemp.line_bcoef();
	coeftemp.resize(n,ncd);
	double tse;
	if(at_start)
		tse=t->param_s();
	else
		tse=t->param_e();

	MGPosition uv=perimeter_uv(perimeter,param);//Surface parameter value of param.
	int ndu=0,ndv=0;
	if(perimeter==0 || perimeter==2) ndv=1;
	else                             ndu=1;
	double slen=length/(eval(uv,ndu,ndv)).len();

	int nnew; double firstd_len,dlen;
	for(int i=0; i<m; i++){
		if(perimeter==0 || perimeter==2){
			for(int j=0; j<n; j++)
				for(int k=0; k<ncd; k++) coeftemp(j,k)=coef(i,j,k);
		}else{
			for(int j=0; j<n; j++)
				for(int k=0; k<ncd; k++) coeftemp(j,k)=coef(j,i,k);
		}
		coeftemp.set_length(n);

	//Compute first derivative length at the end of the extrapolating line.
		t1=*t; 
		firstd_len=lbtemp.eval(tse,1).len();
		dlen=firstd_len*slen;
		lbtemp.extend(at_start,dlen,dk);
		nnew=lbtemp.bdim();
		if(perimeter==0 || perimeter==2){
			if(i==0){
				nv=nnew;
				surf.resize(nu,nv,ncd);
			}
			for(int j=0; j<nnew; j++)
				for(int k=0; k<ncd; k++) surf(i,j,k)=coeftemp(j,k);
		}else{
			if(i==0){
				nu=nnew;
				surf.resize(nu,nv,ncd);
			}
			for(int j=0; j<nnew; j++)
				for(int k=0; k<ncd; k++) surf(j,i,k)=coeftemp(j,k);
		}
	}

	*t=std::move(t1);
	surf.set_length(nu,nv);
	surface_bcoef()=std::move(surf);
	return *this;
}

///Compute binominal coefficients.	
///Let i=degree, then bc(i,j) contains j-th coefficient of the degree i.
///bc(i,j) for 0<=i<=m and 0<=j<=i in bc[(m+1)*i+j].
///bc is an arrary of length (m+1)*(m+1).
void MGCL::Binominal(int m, double* bc){
	assert(m>=1);

	int mp1=m+1;
	bc[0]=bc[mp1]=bc[mp1+1]=1.;
	int i,j,n,n_1=mp1;
	for(i=2; i<=m; i++){
		n=mp1+n_1;
		bc[n]=1.;
		for(j=1;j<i; j++) bc[n+j]=bc[n_1+j-1]+bc[n_1+j];
		bc[n+i]=1.;
		n_1=n;
	}
}

//Return non_homogeneous B-Coefficients with weights of
//the rational B-Spline. This MGSPointSeq includes weights.
MGSPointSeq MGRSBRep::non_homogeneous_bcoef() const{
	int i,j,r, sd=sdim(), m=bdim_u(), n=bdim_v();
	MGSPointSeq cp(m,n,sd+1);
	for(i=0; i<m; i++){
	for(j=0; j<n; j++){
		double weight=coef(i,j,sd);
		for(r=0; r<sd; r++) cp(i,j,r)=coef(i,j,r)/weight;
		cp(i,j,sd)=weight;
	}
	}
	return cp;
}

//Test if this is actually non_rational, i.e. , all of the weights are
//same values.
int MGRSBRep::non_rational()const{
	int sd=sdim();
	double weight=coef(0,0,sd);
	for(int i=0; i<bdim_u(); i++)
	for(int j=0; j<bdim_v(); j++)
		if(!MGREqual2(weight,coef(i,j,sd))) return 0;
	return 1;
}

MGBox MGRSBRep::param_range() const{return m_surface.param_range();}

// Compute parameter curve.
//Returned is newed area pointer, and must be freed by delete.
MGCurve* MGRSBRep::parameter_curve(
	int is_u				//Indicates x is u-value if is_u is true.
	, double x				//Parameter value.
							//The value is u or v according to is_u.
)const{
	MGRLBRep* line=new MGRLBRep;
	line->homogeneous()=m_surface.parameter_line(is_u,x);
	return line;
}

// Compute parameter line.
MGRLBRep MGRSBRep::parameter_line(
		int is_u			//Indicates x is u-value if is_u is true.
		, double x			//Parameter value.
							//The value is u or v according to is_u.
)const{
	MGRLBRep line;
	line.homogeneous()=m_surface.parameter_line(is_u,x);
	return line;
}

//Compute part of the surface limitted by the parameter range bx.
//bx(0) is the parameter (us,vs) and bx(1) is (ue,ve).
//That is u range is from us to ue , and so on.
MGRSBRep* MGRSBRep::part(const MGBox& bx,int multiple)const{
	MGRSBRep* rsb=new MGRSBRep;
	shrinkToParameters(bx, *rsb, multiple);
	return rsb;
}

// Compute perimeter line B-Rep.
MGRLBRep MGRSBRep::perimeter(int i) const
// i is perimeter number:
// =0: v=min line, =1: u=max line, =2: v=max line, =3: u=min line
{
	assert(i<4);

	int is_u; double x;
	switch(i){
		case 0:  is_u=0; x=param_s_v(); break;
		case 1:  is_u=1; x=param_e_u(); break;
		case 2:  is_u=0; x=param_e_v(); break;
		default: is_u=1; x=param_s_u(); break;
	}
	return parameter_line(is_u, x);
}

//Operator overload

//Assignment.
//When the leaf object of this and srf2 are not equal, this assignment
//does nothing.
MGRSBRep& MGRSBRep::operator=(const MGGel& gel2){
	const MGRSBRep* gel2_is_this=dynamic_cast<const MGRSBRep*>(&gel2);
	if(gel2_is_this)
		operator=(*gel2_is_this);
	return *this;
}
MGRSBRep& MGRSBRep::operator=(MGGel&& gel2){
	MGRSBRep* gel2_is_this=dynamic_cast<MGRSBRep*>(&gel2);
	if(gel2_is_this)
		operator=(std::move(*gel2_is_this));
	return *this;
}

// �Ȑ��̕��s�ړ����s���I�u�W�F�N�g�𐶐�����B
//Translation of the curve.
MGRSBRep MGRSBRep::operator+ (const MGVector& v) const{
	MGSPointSeq bc(m_surface.surface_bcoef());
	bc.homogeneous_transform(v);
	MGRSBRep rsb;
	rsb.buildRSBRepFromMemberData(std::move(bc), knot_vector_u(), knot_vector_v());
	return rsb;
}

// �^�x�N�g�������Ȑ��𕽍s�ړ����Ď��g�Ƃ���B
//Translation of the curve.
MGRSBRep& MGRSBRep::operator+= (const MGVector& v){
	m_surface.surface_bcoef().homogeneous_transform(v);
	m_box+=v;
	return *this;
}

// �Ȑ��̋t�����ɕ��s�ړ����s���I�u�W�F�N�g�𐶐�����B
//Translation of the curve.
MGRSBRep MGRSBRep::operator- (const MGVector& v) const{
	MGSPointSeq bc(m_surface.surface_bcoef());
	bc.homogeneous_transform(-v);
	MGRSBRep rsb;
	rsb.m_surface.buildSBRepFromMemberData(
		std::move(bc), knot_vector_u(), knot_vector_v());
	return rsb;
}

// �^�x�N�g�������Ȑ����}�C�i�X�����ɕ��s�ړ����Ď��g�Ƃ���B
//Translation of the curve.
MGRSBRep& MGRSBRep::operator-= (const MGVector& v){
	m_surface.surface_bcoef().homogeneous_transform(-v);
	m_box-=v;
	return *this;
}

// �^����ꂽ�X�P�[���������I�u�W�F�N�g�𐶐�����B
//generate line by scaling.
MGRSBRep MGRSBRep::operator* (double s) const{
	MGSPointSeq bc(m_surface.surface_bcoef());
	bc.homogeneous_transform(s);
	MGRSBRep rsb;
	rsb.m_surface.buildSBRepFromMemberData(
		std::move(bc), knot_vector_u(), knot_vector_v());
	return rsb;
}

// �^����ꂽ�X�P�[���������I�u�W�F�N�g�𐶐�����B
//generate line by scaling.
MGRSBRep operator* (double scale, const MGRSBRep& rsb){
	return rsb*scale;
}

// ���g�̋Ȑ��ɗ^����ꂽ�X�P�[����������B
//Scale the curve.
MGRSBRep& MGRSBRep::operator*= (double scale){
	m_surface.surface_bcoef().homogeneous_transform(scale);
	invalidateBox();
	return *this;
}

// �^����ꂽ�ϊ��ŋȐ��̕ϊ����s���I�u�W�F�N�g�𐶐�����B
//Matrix transformation of the curve.
MGRSBRep MGRSBRep::operator* (const MGMatrix& mat) const{
	MGSPointSeq bc(m_surface.surface_bcoef());
	bc.homogeneous_transform(mat);
	MGRSBRep rsb;
	rsb.m_surface.buildSBRepFromMemberData(
		std::move(bc), knot_vector_u(), knot_vector_v());
	return rsb;
}

// �^����ꂽ�ϊ��ŋȐ��̕ϊ����s�����g�̋Ȑ��Ƃ���B
//Matrix transformation of the curve.
MGRSBRep& MGRSBRep::operator*=(const MGMatrix& mat){
	m_surface.surface_bcoef().homogeneous_transform(mat);
	invalidateBox();
	return *this;
}

// �^����ꂽ�ϊ��ŋȐ��̃g�����X�t�H�[�����s���I�u�W�F�N�g�𐶐�����B
//General transformation of the curve.
MGRSBRep MGRSBRep::operator* (const MGTransf& tr) const{
	MGSPointSeq bc(m_surface.surface_bcoef());
	bc.homogeneous_transform(tr);
	MGRSBRep rsb;
	rsb.m_surface.buildSBRepFromMemberData(
		std::move(bc), knot_vector_u(), knot_vector_v());
	return rsb;
}

// �^����ꂽ�ϊ��ŋȐ��̃g�����X�t�H�[�����s�����g�Ƃ���B
//General transformation of the curve.
MGRSBRep& MGRSBRep::operator*= (const MGTransf& tr){
	m_surface.surface_bcoef().homogeneous_transform(tr);
	invalidateBox();
	return *this;
}

// �_�����Z�q�̑��d��`
// �^�Ȑ��Ǝ��g�����������̔�r������s���B
//RSB and SB comparison.
bool MGRSBRep::operator==(const MGSBRep& sb)const{
	if(sdim()!=sb.sdim())
		return 0;//Check of space dimension.
	if(order_u()!=sb.order_u())
		return 0;	//Check of order.
	if(order_v()!=sb.order_v())
		return 0;	//Check of order.

	int bdu=bdim_u(), bdv=bdim_v();
	if(bdu!=sb.bdim_u())
		return 0;		//Check of B-Rep dimension.
	if(bdv!=sb.bdim_v())
		return 0;		//Check of B-Rep dimension.
	if(bdu<=0 && bdv<=0)
		return 1;
	if(!non_rational())
		return 0;		//Check of rationality.
	if(knot_vector_u() != sb.knot_vector_u())
		return 0;//Check of knot vector.
	if(knot_vector_v() != sb.knot_vector_v())
		return 0;//Check of knot vector.

	//Finally, check of control polygon.
	return m_surface.surface_bcoef().non_homogeneous()==sb.surface_bcoef();
}

bool MGRSBRep::operator==(const MGRSBRep& rsb2)const{
	int sd1=sdim(), sd2=rsb2.sdim();
	if(sd1!=sd2)
		return 0;
	if(order_u()!=rsb2.order_u())
		return 0;	//Check of order.
	if(order_v()!=rsb2.order_v())
		return 0;	//Check of order.
	int bdu1=bdim_u(), bdu2=rsb2.bdim_u();
	if(bdu1!=bdu2)
		return 0;
	int bdv1=bdim_v(), bdv2=rsb2.bdim_v();
	if(bdv1!=bdv2)
		return 0;
	if(bdu1<=0 && bdv1<=0)
		return 1;
	if(knot_vector_u() != rsb2.knot_vector_u())
		return 0;
	if(knot_vector_v() != rsb2.knot_vector_v())
		return 0;

	double ratio=rsb2.coef(0,0,sd2)/coef(0,0,sd1);
	//Check if weights are equal.
	for(int i=0; i<bdu1; i++)
	for(int j=0; j<bdv1; j++)
		if(!MGREqual2(ratio,rsb2.coef(i,j,sd2)/coef(i,j,sd1)))
			return 0;

	return m_surface.surface_bcoef().non_homogeneous()
			==rsb2.m_surface.surface_bcoef().non_homogeneous() ;
}
bool MGRSBRep::operator<(const MGRSBRep& gel2)const{
	return m_surface<gel2.m_surface;
}
bool MGRSBRep::operator==(const MGGel& gel2)const{
	const MGRSBRep* gel2_is_this=dynamic_cast<const MGRSBRep*>(&gel2);
	if(gel2_is_this)
		return operator==(*gel2_is_this);
	else{
		const MGSBRep* gel2_is_sb=dynamic_cast<const MGSBRep*>(&gel2);
		if(gel2_is_sb)
			return operator==(*gel2_is_sb);
	}
	return false;
}
bool MGRSBRep::operator<(const MGGel& gel2)const{
	const MGRSBRep* gel2_is_this=dynamic_cast<const MGRSBRep*>(&gel2);
	if(gel2_is_this)
		return operator<(*gel2_is_this);
	return identify_type() < gel2.identify_type();
}
