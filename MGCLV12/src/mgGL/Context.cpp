/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#include "StdAfx.h"
#include "mg/Tolerance.h"
#include "mg/Ofstream.h"
#include "mg/Ifstream.h"
#include "mgGL/Context.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

////////////////////////////////////////////////////////////////////
// MGContext defines the attributes of a document.

MGContext::MGContext():m_line_density(1),m_smooth(float(.01)),
m_pick_aperture(float(5.)),m_appearance(0)
{
	//m_view.cplane().get_colors(
	//	m_gridColor[0],m_gridColor[1],m_gridColor[2]);
	////setDefaultTolerance();
}

MGContext::MGContext(const MGBox& bx)
:m_line_density(1),m_smooth(float(.01)),
m_pick_aperture(float(5.)),m_appearance(0),
m_view(bx){

	// グリッド定義の作成
	double span;
	int lnum;
	MGPosition mid;

	int sdid;//maximum area coordinate pair will be output.
				//0:(x,y), 1:(y,z), 2:(z,x)
	MGcplane_parameter(bx*INITIAL_SCALE,span,lnum,sdid,mid);
	// Contextにセット
	m_gridNum[0]=m_gridNum[1]=lnum;
	m_gridSpan[0]=m_gridSpan[1]=span;

	//m_view.cplane().get_colors(
	//	m_gridColor[0],m_gridColor[1],m_gridColor[2]);
	//setDefaultTolerance();
}

MGContext::MGContext(
	const MGSnapAttrib& snap_attrib,//Snap data
	
	int line_density,		//line density for a surface to draw in wire mode.
	const MGColor& Bcolor,	//Background color.
	const MGColor& Gcolor,	//Object lines color.
	const MGColor& Hcolor,	//Object highlight color.
	float smooth,	//Smoothness of the curves to draw.
		// 1/smooth is the division number of a curve whose length is the window width.
		// When smooth becomes small, smoothness increases.	
	float pick_aperture,//Pick aperture. Number of pixels to allow picking.
	const MGglViewAttrib& view, ///<MGglViewAttrib of the mainView of the document.
	const MGColor gridColor[4],///< Construction plane's grid colors.
	const int gridNum[2], ///< Construction plane's grid number.
	const double gridSpan[2], ///< Construction plane's grid  span.
	const double torelance[6],//MGCL Tolerance.
		//[0]=wc_zero;
		//[1]=rc_zero;
		//[2]=mach_zero;
		//[3]=line_zero;
		//[4]=angle_zero;
		//[5]=max_knot_ratio;
	const mgTLInputParam& tessellate_param,//tessellation parameter.
	MGAppearance* appearance//must be a newed object of MGAppearance.
		//the ownership will be transfered to this MGContext.
):m_snap_attrib(snap_attrib),m_line_density(line_density),m_smooth(smooth),
m_pick_aperture(pick_aperture),m_tessellate_param(tessellate_param),
m_appearance(appearance),
m_Bcolor(Bcolor),	//Background color.
m_Gcolor(Gcolor),	//Object lines color.
m_Hcolor(Hcolor),	//Object highlight color.
m_view(view)
{
	for(int i=0; i<4; i++){
		m_gridColor[i]=gridColor[i];
	}

	for(int i =0; i<2; ++i){
		m_gridNum[i]=gridNum[i];
	}

	for(int i =0; i<2; ++i){
		m_gridSpan[i]=gridSpan[i];
	}

	for(int i=0; i<6; i++){
		m_tolerance[i]=m_tolerance[i];
	}
}

//copy constructor
MGContext::MGContext(
	const MGContext& context
):m_snap_attrib(context.m_snap_attrib),m_line_density(context.m_line_density),
m_smooth(context.m_smooth),m_pick_aperture(context.m_pick_aperture),
m_tessellate_param(context.m_tessellate_param)
,m_appearance(0),
m_Bcolor(context.m_Bcolor),	//Background color.
m_Gcolor(context.m_Gcolor),	//Object lines color.
m_Hcolor(context.m_Hcolor),	//Object highlight color.
m_view(context.m_view){
	for(int i=0; i<4; i++){
		m_gridColor[i]=context.m_gridColor[i];
	}

	for(int i =0; i<2; ++i){
		m_gridNum[i]=context.m_gridNum[i];
	}

	for(int i =0; i<2; ++i){
		m_gridSpan[i]=context.m_gridSpan[i];
	}

	for(int i=0; i<6; i++){
		m_tolerance[i]=context.m_tolerance[i];
	}

	if(context.m_appearance){
		m_appearance=context.m_appearance->clone();
	}

}

//Assignment operator
MGContext& MGContext::operator=(const MGContext& context){
	if(this==&context)
		return *this;

	m_snap_attrib=context.m_snap_attrib;
	m_line_density=context.m_line_density;
	m_smooth=context.m_smooth;
	m_pick_aperture=context.m_pick_aperture;
	m_tessellate_param=context.m_tessellate_param;

	delete m_appearance; m_appearance=0;
	if(context.m_appearance)
		m_appearance=context.m_appearance->clone();

	m_Bcolor=context.m_Bcolor;	//Background color.
	m_Gcolor=context.m_Gcolor;	//Object lines color.
	m_Hcolor=context.m_Hcolor;	//Object highlight color.

	m_view=context.m_view;
	for(int i=0; i<4; i++)
		m_gridColor[i]=context.m_gridColor[i];

	for(int i=0; i<6; i++)
		m_tolerance[i]=context.m_tolerance[i];
	return *this;
}
MGContext& MGContext::operator=(const MGGel& gel2){
	const MGContext* gel2_is_this=dynamic_cast<const MGContext*>(&gel2);
	if(gel2_is_this)
		operator=(*gel2_is_this);
	return *this;
}

bool MGContext::operator<(const MGContext& gel2)const{
	return this<&gel2;
}
bool MGContext::operator<(const MGGel& gel2)const{
	const MGContext* gel2_is_this=dynamic_cast<const MGContext*>(&gel2);
	if(gel2_is_this)
		return operator<(*gel2_is_this);
	return false;
}

//////////// Destructor.////////

MGContext::~MGContext(){
	delete m_appearance;
}

void MGContext::set_Bcolor(const MGColor& Bcolor){
	m_Bcolor=Bcolor;	//Background color.
}
void MGContext::set_Gcolor(const MGColor& Gcolor){
	m_Gcolor=Gcolor;	//Background color.
}
void MGContext::set_Hcolor(const MGColor& Hcolor){
	m_Hcolor=Hcolor;	//Background color.
}
void MGContext::set_gridColors(const MGColor gridColor[4]){
	for(int j=0; j<4; j++){
		m_gridColor[j]=gridColor[j];
	}
}

void MGContext::set_gridNum(const int gridNum[2]){
	for(int j=0; j<2; ++j){
		m_gridNum[j]=gridNum[j];
	}
}

void MGContext::set_gridSpan(const double gridSpan[2]){
	for(int j=0; j<2; ++j){
		m_gridSpan[j]=gridSpan[j];
	}
}

//Set the view data of the view view_num
void MGContext::set_view(
	const MGglViewAttrib& view
){
	m_view=view;
};

//Execution of MGCL tolerance.
void MGContext::exec_tolerance()const{
	MGTolerance::set_wc_zero(m_tolerance[0]);//
	MGTolerance::set_rc_zero(m_tolerance[1]);//
	MGTolerance::set_mach_zero(m_tolerance[2]);//
	MGTolerance::set_line_zero(m_tolerance[3]);//
	MGTolerance::set_angle_zero(m_tolerance[4]);//
	MGTolerance::set_max_knot_ratio(m_tolerance[5]);//
}

//Tolerance
void MGContext::set_tolerance(
	double wc_zero,
	double rc_zero,
	double mach_zero,
	double line_zero,
	double angle_zero,
	double max_knot_ratio
){
	m_tolerance[0]=wc_zero;
	m_tolerance[1]=rc_zero;
	m_tolerance[2]=mach_zero;
	m_tolerance[3]=line_zero;
	m_tolerance[4]=angle_zero;
	m_tolerance[5]=max_knot_ratio;
}

//Appearance data.
//set appearance.
void MGContext::set_appearance(
	MGAppearance* appearance	//appearance must be a newed object. The ownership will
								//be transfered to this MGContext.
){
	delete m_appearance;
	m_appearance=appearance;
}

std::unique_ptr<MGAppearance> MGContext::remove_appearance(){
	std::unique_ptr<MGAppearance> app(m_appearance);
	m_appearance=nullptr;
	return app;
}

//Execution of drawing OpenGL attributes.
//Valid OpenGL rendering context must be made current.
void MGContext::exec_draw_attributes(mgVBO& vbo)const{
	//glDisable(GL_LIGHTING);
	if(!m_appearance) return;
	m_appearance->drawAttrib(vbo);
}

//Execution of rendering OpenGL attributes.
//Valid OpenGL rendering context must be made current.
void MGContext::exec_render_attributes(mgVBO& vbo)const{
	if(!m_appearance) return;
	m_appearance->render(vbo);
}


// Output virtual function.
std::ostream& MGContext::toString(std::ostream& ostrm) const{
	const float* Bcolor=m_Bcolor.color();
	const float* Gcolor=m_Gcolor.color();
	const float* Hcolor=m_Hcolor.color();
	ostrm<<"MGContext="<<this<<","<<m_snap_attrib<<std::endl;
	ostrm<<",Bcolor["<<Bcolor[0]<<","<<Bcolor[1]<<","<<Bcolor[2]<<"]";
	ostrm<<",Gcolor["<<Gcolor[0]<<","<<Gcolor[1]<<","<<Gcolor[2]<<"]";
	ostrm<<",Hcolor["<<Hcolor[0]<<","<<Hcolor[1]<<","<<Hcolor[2]<<"]";
	ostrm<<",smooth="<<m_smooth;
	ostrm<<",pick_aperture="<<m_pick_aperture;

	ostrm<<std::endl<<",view=";
	ostrm<<m_view;

	const float* gridColor= m_gridColor[0].color();
	const float* axisXcolor = m_gridColor[1].color();
	const float* axisYcolor = m_gridColor[2].color();
	const float* axisZcolor = m_gridColor[3].color();

	ostrm<<",gridColor["<<gridColor[0]<<","<<gridColor[1]<<","<<gridColor[2]<<"]";
	ostrm<<",axisXcolor["<<axisXcolor[0]<<","<<axisXcolor[1]<<","<<axisXcolor[2]<<"]";
	ostrm<<",axisYcolor["<<axisYcolor[0]<<","<<axisYcolor[1]<<","<<axisYcolor[2]<<"]";
	ostrm<<",axisZcolor["<<axisZcolor[0]<<","<<axisZcolor[1]<<","<<axisZcolor[2]<<"]";

	ostrm<<std::endl<<",gridNum=("<< m_gridNum[0] << "," << m_gridNum[1]<<")";
	ostrm<<",gridSpan=("<< m_gridSpan[0] << "," << m_gridSpan[1]<<")";

	ostrm<<",wc_zero="<<m_tolerance[0];
	ostrm<<",rc_zero="<<m_tolerance[1];
	ostrm<<",mach_zero="<<m_tolerance[2];
	ostrm<<",line_zero="<<m_tolerance[3];
	ostrm<<",angle_zero="<<m_tolerance[4];
	ostrm<<",max_knot_ratio="<<m_tolerance[5];

	ostrm<<std::endl<<",tessellate_param="<<m_tessellate_param;

	ostrm<<",appearance=";
	if(m_appearance){
		ostrm<<(*m_appearance);
	}else{
		ostrm<<"NULL";
	}
	return ostrm;
}

//Write all member data.
// ヘッダでの宣言順にRead/Write
void MGContext::WriteMembers(MGOfstream& buf)const{
	buf<<m_snap_attrib;
	buf<<m_line_density;
	m_Bcolor.WriteMembers(buf);
	m_Gcolor.WriteMembers(buf);
	m_Hcolor.WriteMembers(buf);
	buf<<m_smooth;
	buf<<(double)m_pick_aperture;
	buf<<m_view;
	for(int i=0; i<4; i++){
		m_gridColor[i].WriteMembers(buf);
	}

	for(int i=0;i<2;++i){
		buf << m_gridNum[i];
	}

	for(int i=0;i<2;++i){
		buf << m_gridSpan[i];
	}

	for(int i=0; i<6; i++){
		buf<<m_tolerance[i];
	}

	buf<<m_tessellate_param;
	buf.WritePointer(m_appearance);
}

//Write all member data
// ヘッダでの宣言順にRead/Write
void MGContext::ReadMembers(MGIfstream& buf){
	buf>>m_snap_attrib;
	buf>>m_line_density;
	m_Bcolor.ReadMembers(buf);
	m_Gcolor.ReadMembers(buf);
	m_Hcolor.ReadMembers(buf);

	buf>>m_smooth;
	double aperture;
	buf>>aperture;m_pick_aperture=float(aperture);
	buf>>m_view;
	for(int i=0; i<4; i++){
		m_gridColor[i].ReadMembers(buf);
	}

	for(int i=0;i<2;++i){
		buf >> m_gridNum[i];
	}

	for(int i=0;i<2;++i){
		buf >> m_gridSpan[i];
	}

	for(int i=0; i<6; i++){
		buf >> m_tolerance[i];
	}

	buf>>m_tessellate_param;

	m_appearance=static_cast<MGAppearance*>(buf.ReadPointer());
}
