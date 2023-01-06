/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
//! @file
//!	@brief  Declaration for class MGIgesPD142.
//!	@author System fugen

#include "StdAfx.h"
#include "mg/Tolerance.h"
#include "mg/CompositeCurve.h"
#include "topo/Loop.h"
#include "topo/Face.h"
#include "mgiges/IgesGsec.h"
#include "mgiges/IgesPD142.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

//!	@brief MGIgesPD142 is the class for Iges parameter data type 142(Curve on parameteric space).
using namespace MGIges;

// Constructors.

//! Constructs an object of class MGIgesPD142.
MGIgesPD142::MGIgesPD142(MGIgesDirectoryEntry* DEpointer)
:MGIgesPD(CURVE_ON_PARAMETRIC_SURFACE,DEpointer),m_created_way(0),m_prefered(0),
m_surface_DE(0),m_param_curve_DE(0),m_model_curve_DE(0){
}

//! Constructs an object of class MGIgesPD142.
MGIgesPD142::MGIgesPD142(
	const MGLoop& loop,	//loop to make PD142. This is a loop of the face.
	int surface_DE,	//the base surface. The surface must be output to IGES file first.
	MGIgesOfstream& igesfile//Iges file to output.
):MGIgesPD(CURVE_ON_PARAMETRIC_SURFACE),m_created_way(0),m_prefered(1),
m_surface_DE(surface_DE),
m_param_curve_DE(0),m_model_curve_DE(0){
	MGIgesDEStatusNumber::SESwitch pld=MGIgesDEStatusNumber::PLDependent;

	std::vector<UniqueCurve> uvcurves=loop.curves();
	int i, nuvcurves=(int)uvcurves.size();
	if(!nuvcurves)
		return;
	MGCurve* c0=uvcurves[0].release();
	MGCompositeCurve uvcomposites(c0);
	for(i=1; i<nuvcurves; i++){
		MGCurve* ci=uvcurves[i].release();
		uvcomposites.connect_to_end(ci);
	}
	m_param_curve_DE=uvcomposites.out_to_IGES(igesfile,pld);
	MGIgesFstream::UniqueDE& pcDE=igesfile.directoryEntry(m_param_curve_DE);
	pcDE->set_SubordinateEntitySwitch(MGIgesDEStatusNumber::PLDependent);
								//Phisically and logically dependent.

	std::vector<UniqueCurve> curves=loop.curves_world();
	int ncurves=(int)curves.size();
	if(!ncurves)
		return;
	MGCurve* cv0=curves[0].release();
	MGCompositeCurve composites(cv0);
	for(i=1; i<ncurves; i++){
		MGCurve* cvi=curves[i].release();
		composites.connect_to_end(cvi);
	}
	m_model_curve_DE=composites.out_to_IGES(igesfile,pld);
}

//Obtain both the parametric space curve of the surface and the model space curve.
void MGIgesPD142::trim_face(
	const MGIgesIfstream& igesifstrm,
	std::unique_ptr<MGFace>& face,//Face to be trimmed by this boundary MGIgesPD142.
	bool outer	//True if this be the outer boundary.
)const{
	const MGSurface& srf=*(face->surface());
	
	std::unique_ptr<MGCurve> param_curve;
	MGCurve* param_curve2=0;
	std::unique_ptr<MGGel> uvcrvobj;
	if(m_param_curve_DE){
		mgTolSetWCZero wczeroSet(srf.param_error());//Set&save the error.
		uvcrvobj=std::unique_ptr<MGGel>(igesifstrm.convert_to_gel(m_param_curve_DE));
		param_curve2=dynamic_cast<MGCurve*>(uvcrvobj.get());
		if(param_curve2){
			param_curve=std::unique_ptr<MGCurve>(dynamic_cast<MGCurve*>(uvcrvobj.release()));
			param_curve->change_dimension(2);//Since the original is 3D curve.
		}
	}

	std::unique_ptr<MGCurve> model_curve;
	MGCurve* model_curve2=0;
	std::unique_ptr<MGGel> modelcrvobj;
	if(m_model_curve_DE){
		modelcrvobj=std::unique_ptr<MGGel>(igesifstrm.convert_to_gel(m_model_curve_DE));
		model_curve2=dynamic_cast<MGCurve*>(modelcrvobj.get());
		if(model_curve2)
			model_curve=std::unique_ptr<MGCurve>(
				dynamic_cast<MGCurve*>(modelcrvobj.release()));
	}

	std::unique_ptr<MGLoop> loop;
	if(param_curve2){
		MGCompositeCurve* ccrv=dynamic_cast<MGCompositeCurve*>(param_curve2);
		if(ccrv){
			//If ccrv is composite, model curve is not used.
			std::unique_ptr<MGCurve> model_curveDummy;
			loop=std::unique_ptr<MGLoop>(new MGLoop(param_curve,model_curveDummy));
		}else{
			if(model_curve2){
				loop=std::unique_ptr<MGLoop>(new MGLoop(param_curve,model_curve));
			}
		}
	}else{
		if(!model_curve.get())
			return;
		loop=face->build_loop(*model_curve);
	}

	if(outer){
		if(loop->area()<0.)
			loop->negate();
		loop->make_close();
		face->prepend_boundary(loop.release());
	}else{
		if(loop->area()>0.)
			loop->negate();
		loop->make_close();
		face->append_boundary(loop.release());
	}

}

//Read in parameter data from string stream data.
void MGIgesPD142::read_in(
	char pDelimeter,
	std::istringstream& pdstream
){
	MGCL::get_integer(pDelimeter,pdstream,m_created_way);
	get_DEpointer(pDelimeter,pdstream,m_surface_DE);
	get_DEpointer(pDelimeter,pdstream,m_param_curve_DE);
	get_DEpointer(pDelimeter,pdstream,m_model_curve_DE);
	MGCL::get_integer(pDelimeter,pdstream,m_prefered);
}

//Write out this PD as MGIgesParamLine's(into plines).
//Except for string data, one integer or double data is output
//into one MGIgesParamLine, not striding over more than one line.
//Only when string data is output(to Holleris string), the data
//may stride over more than one lines.
//plines[i] for 0<=i<plines.size() are valid.
void MGIgesPD142::write_out_into_string(
	const MGIgesGSec& gsec,	
	std::vector<std::string>& plines ///<output plines.
)const{
	put_integer(m_created_way,gsec,plines);
	put_DEpointer(m_surface_DE,gsec,plines);
	put_DEpointer(m_param_curve_DE,gsec,plines);
	put_DEpointer(m_model_curve_DE,gsec,plines);
	put_integer(m_prefered,gsec,plines);
}
