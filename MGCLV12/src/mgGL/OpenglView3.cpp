/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
// MGOpenGLView.cpp : インプリメンテーション ファイル
//
#include "StdAfx.h"
#include "mg/Tolerance.h"
#include "mg/Position.h"
#include "mg/AttribedGel.h"
#include "mg/Group.h"
#include "mg/GelPositions.h"
#include "mg/CParam_list.h"
#include "mg/DnameControl.h"
#include "topo/Edge.h"
#include "topo/Loop.h"
#include "topo/Face.h"
#include "topo/Shell.h"
#include "mgGL/OpenGLView.h"
#include "mgGL/GLAttrib.h"
#include "mgGL/SysGLList.h"
#include "mgGL/glViewAttrib.h"
#include "mgGL/VBO.h"
#include "mgGL/glslprogram.h"

using namespace std;

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

namespace {
	GLuint glErr;
	const glm::vec3 ORIGIN(0., 0., 0.);
};

///Set up drawing environment.
void MGOpenGLView::setupDrawEnv(
	const MGColor& backColor,//When pick mode, backColor is not used.
	const float* centrApertr //centrApertu = nullptr means standard draw,
	//and centrApertu != null means selection mode.
){
///// initialize OpenGL

	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LEQUAL);
	const float* Bcolr = backColor.color();
	glClearColor(Bcolr[0], Bcolr[1], Bcolr[2], Bcolr[3]);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	if (centrApertr) {//When selection
		glDisable(GL_BLEND);
	}else{//When non selection
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	}
	glClearDepth(1.0);
	glEnable(GL_POLYGON_OFFSET_FILL);// ポリゴンオフセットフィルを設定
	glPolygonOffset(1., 1.);

	mgGLSLProgram* glsl = mgGLSLProgram::getCurrentGLSLProgram();

///Shader program's uniform set.

//(1) Projection matrix
	glm::mat4 projMat;///Projection matrix to set into uniform.
	get_projection_matrix(m_viewPort, projMat);
	if(centrApertr){//When selection
		glm::vec2 c2(centrApertr[0], centrApertr[1]);
		glm::vec2 delta2(centrApertr[2], centrApertr[3]);
		glm::ivec4 vpOld(m_viewPort[0], m_viewPort[1], m_viewPort[2], m_viewPort[3]);
		projMat = glm::pickMatrix(c2, delta2, vpOld)*projMat;
	}
	glsl->setUniform(mgGLSLProgram::projMatrix, projMat);

//(2) Model_View matrix.
	glm::mat4 modelViewMat;
	get_model_matrix(modelViewMat);
	glsl->setUniform(mgGLSLProgram::modelViewMatrix, modelViewMat);

//(3) Normal Matrix ModelViewMatrixとNormalMatrixはセットです。
	glm::mat3 normalMatrix = glm::mat3(modelViewMat);
	glsl->setUniform(mgGLSLProgram::normalMatrix, normalMatrix);

//(4) ProjModelView matrix (to save matrix computation in shader program)
	glm::mat4 ProjModelView = projMat * modelViewMat;
	glsl->setUniform(mgGLSLProgram::modelViewProjMatrix, ProjModelView);

//(5) Anchor Point変換用の係数設定
	float widthHalf = float(m_viewPort[2])*0.5f;
	float heightHalf = float(m_viewPort[3])*0.5f;
	glm::vec3 scaleFactor(1.0f / widthHalf, 1.0f / heightHalf, 1.0f);
	glm::mat4 ndcMtx = glm::scale(
		glm::translate(glm::mat4(), glm::vec3(-widthHalf, -heightHalf, 0.0f)),
		scaleFactor);
	glsl->setUniform(mgGLSLProgram::ndcMarix, ndcMtx);

	glm::mat3 ndcScaleMtx = glm::mat3(ndcMtx);
	glsl->setUniform(mgGLSLProgram::ndcScaleMatrix, ndcScaleMtx);

// (6) dpiFactor 
	float func = m_dpi / 72.0f;
	glsl->setUniform(mgGLSLProgram::dpiFactor, func);
}

void MGOpenGLView::execDefaultStaticAttrib() {
	glDisable(GL_CULL_FACE);// 両面を描画
	mgGLSL::execStaticColorAttrib(Gcolor());
	mgGLSL::execStaticLineWidth(1.f);
	mgGLSL::execStaticLineStipple(0, 0);
	mgGLSL::execLightMode(0);//Light Off
}

///Constructglm::lookAtMatrix from eye_position() and view_up_vector().
void MGOpenGLView::setLookAtMat() {
	const MGPosition& eyeP = eye_position(); glm::vec3 eye(eyeP[0], eyeP[1], eyeP[2]);
	const MGVector& upP = view_up_vector(); glm::vec3 up(upP[0], upP[1], upP[2]);
	m_lookAtMat = glm::lookAt(eye, ORIGIN, up);
}

//get ModelView matrix of OpenGL.
void MGOpenGLView::get_model_matrix(
	glm::mat4& modelMat	//double modelMat[16] ///<OpenGL's model matrix
)const {
	glm::mat4 viewMat = m_lookAtMat * m_viewAttrib.m_modelViewMat;
	const MGPosition& cntr = center(); glm::vec3 cntr2(-cntr[0], -cntr[1], -cntr[2]);
	modelMat = glm::translate(viewMat, cntr2)*m_viewAttrib.m_PreCenterMat;//glTranslated;
}

///get projection matrix, given the viewport data 
void MGOpenGLView::get_projection_matrix(
	const int vp[4],///<viewport data ={left, bottom, widht, height}
	glm::mat4& projMat	///<OpenGL's projection matrix
)const {
	double height2 = view_volume_height()*.5;
	double wide2 = height2 * double(vp[2]) / double(vp[3]);
	float left = float(m_viewAttrib.m_cx - wide2), right = float(m_viewAttrib.m_cx + wide2),
		bottom = float(m_viewAttrib.m_cy - height2), top = float(m_viewAttrib.m_cy + height2);

	float znear = (float)view_volume_near(), zfar = (float)view_volume_far();
	projMat = is_perspective() ?
				glm::frustum(left, right, bottom, top, znear, zfar)
			  : glm::ortho(left, right, bottom, top, znear, zfar);
}

void MGOpenGLView::getModelViewProjectionMatrices(
	glm::mat4& modelM,//double modelM[16],
	glm::mat4& projM,// double projM[16],
	int vp[4]		//={left, bottom, width, height}.
)const{
	get_model_matrix(modelM);
	get_viewport(vp);
	get_projection_matrix(&vp[0], projM);
}

//Update the center and the scale of the view.
///The pespectiveness and the cplane are unchanged.
void MGOpenGLView::updateCenterScalle(
	const MGPosition& center,
	double diameter///<diameter of the view. This is set to m_diameter.
		///<diameter of the sphere that sorround the model.
		///<If diameter<=0. the current diameter is not updated.
) {
	for (int i = 0; i < 3; i++)
		m_center_current[i] = (float)center[i];
	if (diameter <= 0.)
		diameter = m_viewAttrib.diameter();
	m_viewAttrib.compute_viewing_environment(center, diameter);
	setLookAtMat();
}

static const MGColor NullColor(0., 0., 0., 0.);
//Function's return value is the number of hit objects.
int MGOpenGLView::pick_to_select_buf(
	const float centrApertr[4],
	///<Screen coordinates. (left, bottom) is (0,0) and (aperturex, aperturey).
	mgVBO* display_list,	//display list that includes pick objects.
	std::set<unsigned>& selected///Selected data will be returned. This data consist of
			///the data set by selectName.
){
	float width=centrApertr[2], height=centrApertr[3];
	int viewport[4] ={int(centrApertr[0]-width/2.+.5),int(centrApertr[1]-height/2.+.5),
					  int(width), int(height)};
	int& x = viewport[0]; int& y = viewport[1];//(left, bottom)
	int& w = viewport[2]; int& h = viewport[3];//(width, height)

	int numHit = 0;
	if(h > 0 && w > 0){
		int xOld = m_viewPort[0]; int yOld = m_viewPort[1];
		int wOld = m_viewPort[2]; int hOld = m_viewPort[3];

		if (x < xOld)
			x = xOld;
		if ((x + w) > (xOld + wOld))
			w = xOld + wOld - x;
		if (y < yOld)
			y = yOld;
		if ((y + h) > (yOld + hOld))
			h = yOld + hOld - y;
		if(h&&w){
			//Set the target selection viewport.
			glViewport(x, y, w, h);
			glDrawBuffer(GL_BACK);
			float cAper[4]={float(x+w*.5), float(y+h*.5), float(w), float(h)};
			setupDrawEnv(NullColor, cAper);

			//Target objects drawing
			display_list->selectionDraw(viewMode());
			glReadBuffer(GL_BACK);
			extractSelected(viewport, selected);
			numHit = (int)selected.size();
		}
	}
	return numHit;
}

//Pick objects in the display list generated by make_display_list.
//Function's return value is MGPickObject vector.
//All the objects which were inside the pick aperture will be output.
//This data can be accessed using current_object(), or current_PickObject().
///pick_glv needs a current view.
MGPickObjects MGOpenGLView::pick_glv(
	const float centrApertr[4],///<specifies pick center and aperture.
	const MGAbstractGels& objtypes
) {
	MGPickObjects pobjs;
	std::set<unsigned> selected;
	mgVBO* vboGroup = display_list();
	int objnum = vboGroup ? pick_to_select_buf(centrApertr, vboGroup, selected):0;
	if (!objnum)
		return pobjs;

	MGDNameControl& dnc = getDNameControlInstance();
	for (auto& i: selected) {
		mgVBO* vboi = dnc.VBO_from_dlistName(i);
		if (!vboi)
			continue;
		MGPickObject pobj;
		vboi->buildGelPosition(vboGroup, pobj);

		MGObject* obj = pobj.leaf_object();
		if (pobj.is_shell_face()) {//When is_shell_face.
			MGShell* shelli = pobj.get_shell_of_shell_face();
			if (obj->type_is(objtypes) || shelli->type_is(objtypes))
				pobjs.push_back(pobj);
		}else if(pobj.leaf_isObject() && obj->type_is(objtypes))
			pobjs.push_back(pobj);
	}

	return pobjs;
}

class mgPerimeterSelection : public mgVBO {
public:
	const MGSurface& m_surf;
	const MGDrawParam& m_dparam;
	mgPerimeterSelection(const MGSurface& surf, const MGDrawParam& para)
		:m_surf(surf), m_dparam(para) {};

	///m_gelの描画データ作成のみをおこなう。
	///すでに作成済みであっても強制的に再作成を行う。
	///m_gel=0のときはなにもしない。
	void make_display_list(MGCL::VIEWMODE viewMode = MGCL::DONTCARE);

	///描画関数selectionDraw()は、Object選択のための表示処理をする。
	///通常のdrawとの相違：///Colorとしてm_bufferIDを用い、size処理以外の
	///attributesの処理（normal, texture, color)をしない。
	void selectionDraw(MGCL::VIEWMODE viewMode = MGCL::DONTCARE);
};

void mgPerimeterSelection::make_display_list(MGCL::VIEWMODE viewMode) {
	clearElements(MGCL::WIRE_AND_SHADING);

	int nperi = m_surf.perimeter_num();
	int ldensity = m_dparam.line_desity_wire_face();
	for (int i = 0; i < nperi; i++) {
		std::unique_ptr<MGCurve> peri(m_surf.perimeter_curve(i));
		peri->drawWire(*this, ldensity);
	}
	setDirty(false);
}
void mgPerimeterSelection::selectionDraw(MGCL::VIEWMODE viewMode) {
	if (!is_made())
		make_display_list();

	mgGLSLProgram* glsl = mgGLSLProgram::getCurrentGLSLProgram();
	glsl->setFuncType(mgGLSL::Select);

	mgCoordinateTypeSwitcher coordType(m_coordinateType, getAnchorPosition());//save the coordinate type.
	size_t n = m_elements.size();
	for (unsigned i = 0; i < n; i++) {
		mgGLSL::setColorAsSelectionName(i + 1);
		UniqueVBOElement& elmi = m_elements[i];
		elmi->selectionDraw(MGCL::WIREVIEW);
	}
}

//Pick a perimeter of the surface surf. That is, obtain the perimeter number
//that passes input (sx,sy) when drawn in the current view matrix.
//Function's return value is perimeter number picked.
//When no perimeters are picked, -1 will be returned.
int MGOpenGLView::pick_perimeter_glv(
	const MGSurface& surf,
	int sx, int sy,	///<Screen coordinates. (left, bottom) is (0,0).
	MGPosition* uv,	//surface parameter (u,v) nearest to (sx,sy) will be returned.
	float aperturex,//specifies pick aperture of x and y.
	float aperturey//When <=0. value is specified, default value(the value
			//obtained by pick_aperture() will be used.
) {
	mgPerimeterSelection periSel(surf, draw_param());

	int perimeter = -1;
	if (aperturex <= 0.) aperturex = pick_aperture();
	if (aperturey <= 0.) aperturey = pick_aperture();
	std::set<unsigned> selected;
	float centrApertr[4] = { (float)sx,(float)sy, aperturex,aperturey };
	pick_to_select_buf(centrApertr, &periSel, selected);
	int objnum = (int)selected.size();
	if (objnum > 0) {
		perimeter = *(selected.begin()) - 1;
		if (uv) { // if parameter is needed
			double t;
			std::unique_ptr<MGCurve> peri(surf.perimeter_curve(perimeter));
			get_near_position(peri.get(), centrApertr, t);
			*uv = surf.perimeter_uv(perimeter, t);
		}
	}
	return perimeter;
}

class mgEdgeSelection : public mgVBO {
public:
	const MGLoop& m_loop;
	const MGDrawParam& m_dparam;
	mgEdgeSelection(const MGLoop& loop, const MGDrawParam& para)
		:m_loop(loop), m_dparam(para) {};

	///m_gelの描画データ作成のみをおこなう。
	///すでに作成済みであっても強制的に再作成を行う。
	///m_gel=0のときはなにもしない。
	void make_display_list(MGCL::VIEWMODE viewMode = MGCL::DONTCARE);

	///描画関数selectionDraw()は、Object選択のための表示処理をする。
	///通常のdrawとの相違：///Colorとしてm_bufferIDを用い、size処理以外の
	///attributesの処理（normal, texture, color)をしない。
	void selectionDraw(MGCL::VIEWMODE viewMode = MGCL::DONTCARE);
};
void mgEdgeSelection::make_display_list(MGCL::VIEWMODE viewMode) {
	int ldensity = m_dparam.line_desity_wire_face();
	int nedge = m_loop.number_of_edges();
	for (int j = 0; j < nedge; j++) {
		const MGEdge* edge2 = m_loop.edge(j);
		MGEdge& be = *(edge2->make_binder_with_curve());
		MGTrimmedCurve cij = be.trimmed_curve();
		cij.drawWire(*this, ldensity);
	}
	setDirty(false);
}
void mgEdgeSelection::selectionDraw(MGCL::VIEWMODE viewMode) {
	if (!is_made())
		make_display_list();

	mgGLSLProgram* glsl = mgGLSLProgram::getCurrentGLSLProgram();
	glsl->setFuncType(mgGLSL::Select);

	mgCoordinateTypeSwitcher coordType(m_coordinateType, getAnchorPosition());//save the coordinate type.
	size_t n = m_elements.size();
	for (unsigned i = 0; i < n; i++) {
		mgGLSL::setColorAsSelectionName(i + 1);
		UniqueVBOElement& elmi = m_elements[i];
		elmi->selectionDraw(MGCL::WIREVIEW);
	}
}

//Pick an edge of the face f. That is, obtain the edge number
//that passes input (sx,sy) when drawn in the current view matrix.
//Function's return value is the edge pointer picked.
//When no edges are picked, null will be returned.
const MGEdge* MGOpenGLView::pick_edge_glv(
	const MGFace& f,
	int sx, int sy,	///<Screen coordinates. (left, bottom) is (0,0).
	MGPosition* uv,	//surface parameter (u,v) nearest to (sx,sy) will be returned.
	float aperturex,//specifies pick aperture of x and y.
	float aperturey//When <=0. value is specified, default value(the value
			//obtained by pick_aperture() will be used.
) {
	const MGEdge* edge = 0;
	if (aperturex <= 0.) aperturex = pick_aperture();
	if (aperturey <= 0.) aperturey = pick_aperture();
	float centrApertr[4] = { float(sx),float(sy),aperturex, aperturey};

	int nloop = f.number_of_loops();
	for (int i = 0; i < nloop; i++) {
		const MGLoop& li = *(f.loop(i));
		mgEdgeSelection edgeSel(li, draw_param());
		std::set<unsigned> selected;
		pick_to_select_buf(centrApertr, &edgeSel, selected);
		size_t objnum = selected.size();
		if (objnum) {
			edge = li.edge(*selected.begin() - 1);
			if (uv) {
				MGEdge& be = *(edge->make_binder_with_curve());
				MGTrimmedCurve cij = be.trimmed_curve();
				double t;
				get_near_position(&cij, centrApertr, t);
				*uv = edge->eval(be.param_pcell(t));
			}
			break;
		}
	}
	return edge;
}

//Determine if screen coordinate (sx,sy) is closer to the start point or to the end
//of the curve curve.
//Functin's return value is 0: if start point, 1: if end point.
int MGOpenGLView::pick_start_end_glv(
	const MGCurve& curve,
	int sx, int sy	//Screen coordinates. (left, bottom) is (0,0).
) {
	MGStraight sl;
	unproject_to_sl_glv(sx, sy, sl);
	MGPosition P0 = curve.start_point(), P1 = curve.end_point();
	if (sl.distance(P0) <= sl.distance(P1))
		return 0;
	return 1;
}

///Extract selectionName data from the frame buffer drawn by selectionDraw();
void MGOpenGLView::extractSelected(
	const int viewport[4],///Viewport of the selection target window.
	std::set<unsigned>& selected///Selected name data will be returned.
		/// This data consist of the data set by selectionDraw.
){
	GLint x = viewport[0], y = viewport[1];
	GLsizei width = viewport[2], height = viewport[3];
	int numPixels = width * height;
	unsigned* pixels = new unsigned[numPixels];
	glReadPixels(x, y, width, height, GL_RGBA, GL_UNSIGNED_BYTE, pixels);
	//if((glErr=glGetError()) != GL_NO_ERROR){
	//	CString msg(gluErrorString(glErr));
	//	COUT<<"MGOpenGLView::extractSelected::glReadPixels::"<<(TCAST)msg<<std::endl;
	//}
	for (int i = 0; i < numPixels; i++) {
		mgGLSLProgram::SELECT_NAME nub;
		nub.uiName = pixels[i];
		unsigned pixeli = nub.uiName;//unsigned pixeli=pixels[i];
		if (pixeli != 0) {
			selected.insert(pixeli);
		}
	}
	delete[] pixels;
}

//Project world coordinates to OpenGL's screen coordinates.
//If modelMat, projMat, or vp is not input, project will ask OpenGL to get them.
//Generally, users of project are recommended to get modelMat, projlMat, or
//vp, and input them to project if continuous multiple use of project will take place.
void MGOpenGLView::project(
	const MGPosition& world,
	MGPosition& screen,
	const glm::mat4* modelMat,	//OpenGL's model matrix
	const glm::mat4* projlMat	//OpenGL's projection matrix
) const {
	const glm::mat4* model = modelMat;
	const glm::mat4* proj = projlMat;
	glm::mat4 modelM, projM;

	if (!proj) {
		proj = &projM;
		get_projection_matrix(m_viewPort, projM);
	}
	if (!model) {
		model = &modelM;
		get_model_matrix(modelM);
	}

	glm::vec3 obj(world[0], world[1], world[2]);
	glm::ivec4 vp2(m_viewPort[0],m_viewPort[1],m_viewPort[2],m_viewPort[3]);
	glm::vec3 v = glm::project(obj, *model, *proj, vp2);
	screen = MGPosition(v[0], v[1], v[2]);
}
