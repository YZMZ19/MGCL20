/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/

#include "StdAfx.h"
#include "mg/AttribedGel.h"
#include "mg/DNameControl.h"
#include "mg/Curve.h"
#include "mg/GelPosition.h"
#include "topo/Shell.h"
#include "mgGL/StaticGLAttrib.h"
#include "mgGL/OpenGLView.h"
#include "mgGL/glslprogram.h"
#include "mgGL/VBOLeafBuilder.h"
#include "mgGL/VBOPointer.h"
#include "mgGL/VBO.h"
#include "mgGL/VBOLeaf.h"

static const float white[4]={1.,1.,1.,1.};//Highlight back color.

/////////////////////////////////////////////////////////////////////////////
/// Vertex Buffer Object Class.

//static const float edgeColor[4]={1.,.5,.5,0.};//Edge color to hilight.
//static const float white[4]={1.,1.,1.,0.};//Highlight back color.
//static const float endPointColor[4]={.5,1.,.5,0.};//Start/End point color to hilight.

/// mgVBO��OpenGL 4�@�p�`��̂��߂̃N���X
/// ��{�I�ȕ`��̂��߂̃c�[����OpenGL4���g�p�����`��@�\��񋟂���B
///
///******���p��̎菇******
///(1) MGAttribedGel�p��mgVBO�́AmgVBO(const MGAttribedGel&)��constructor�𗘗p���č\�z���A
///   �K�v�ɉ�����drawXXX�܂���begin-end�ɂ��MGAttribedGel�ɒǉ���mgVBOLeaf���쐬��
///   make_display_list()(make_display_list()��m_gel�ɑ΂�make_display_list�𔭐M��m_gel�̕`��f�[�^���쐬����j����B
///�@ ��MGAttribedGel��mgVBO��drawXXX, begin()-end()�𗘗p���ĕ`��f�[�^�쐬����(make_display_list())���s��
///(2) draw()�ɂ��`�揈��������B
///   MGAttribedGel�p��mgVBO�ɑ΂��Ă�make_display_list�����������Ȃ�draw()���Ăяo���Ă��悢�B
///   ���̏ꍇmgVBO��m_gel�ɑ΂���make_display_list()�𔭐M���Adraw���s���B
///
///   draw()/make_display_list���ĂԑO�ɂ͎�ނ̈قȂ镡���̐}�`�̍쐬���s���Ă�
///   �悢(drawXXX, begin()-end())���Adraw()/make_display_list�܂łɍ쐬�����}�`��
///   draw()/make_display_list�ɂ��ЂƂ܂Ƃ߂̈����ƂȂ�B
///   �`��֐�draw()�̓����o�[�f�[�^�̊emgVBOElement�ɑ΂���draw()�𔭐M����B
///
/// begin()-end()�ň͂܂ꂽvertex(����т��̕t��attribute�j�A�����
/// ���ׂĂ�drawXXX�́A��܂��͈�ȏ��mgVBOLeaf��draw()�̂��߂�
/// �쐬����imake_display_list�j�����Ŏ��ۂ̕`��͍s��Ȃ��B
///
///mgVBO��mgVBOLeaf��mgVBO���g��element�Ƃ��ĕێ�����B
///mgVBO��element�Ƃ���̂�MGGroup�܂���MGShell�̃����o�[�f�[�^�̕\������
///�ێ����邽�߁B����ȊO��MGObject��mgVBOLeaf(��ʂɂ͕����j��element�Ƃ��ĕێ�����
namespace{GLenum glErr;};

// �I�y���[�V����

///��MGAttribedGel�p��constructor.
mgVBO::mgVBO(mgGLSL::CoordinateType coordType):m_gel(0),m_lineWidthStatic(-1.f), m_pointSizeStatic(-1.f)
,m_stippleFactor(-1),m_LineStipplePattern(0),m_lightMode(-1)
,m_coordinateType(coordType)
,m_elementsDirty(true),m_elementsShadeDirty(true)
{
	MGDNameControl& dnc=getDNameControlInstance();
	m_dname=dnc.getNewNamebyVBO(this);
	m_target_elements=&m_elements;
}

///MGAttribedGel�p��constructor.
mgVBO::mgVBO(const MGAttribedGel& gel)
:m_gel(const_cast<MGAttribedGel*>(&gel))
,m_lineWidthStatic(-1.f), m_pointSizeStatic(-1.f)
,m_stippleFactor(-1),m_LineStipplePattern(0),m_lightMode(-1)
,m_coordinateType(mgGLSL::World)
,m_elementsDirty(true),m_elementsShadeDirty(true)
{
	MGDNameControl& dnc=getDNameControlInstance();
	m_dname=dnc.getNewNamebyVBO(this);
	m_target_elements=&m_elements;
}

///Copy constructor.
mgVBO::mgVBO(const mgVBO& vbo):m_gel(0)
,m_colorStatic(vbo.m_colorStatic)
,m_lineWidthStatic(vbo.m_lineWidthStatic), m_pointSizeStatic(vbo.m_pointSizeStatic)
,m_stippleFactor(vbo.m_stippleFactor),m_LineStipplePattern(vbo.m_LineStipplePattern)
,m_lightMode(vbo.m_lightMode)
,m_coordinateType(mgGLSL::World)
,m_elementsDirty(true),m_elementsShadeDirty(true)
{
	assert(!vbo.m_gel);///Copy constructor not allowed for gel VBO.
	MGDNameControl& dnc=getDNameControlInstance();
	m_dname=dnc.getNewNamebyVBO(this);
	if(vbo.m_target_elements==&(vbo.m_elements))
		m_target_elements=&m_elements;
	else
		m_target_elements=&m_elementsShade;
}

///Assignment.
mgVBO& mgVBO::operator=(const mgVBO& vbo){
	m_no_display=vbo.m_no_display;

	clearElements();
	m_colorStatic=vbo.m_colorStatic;
	m_lineWidthStatic=vbo.m_lineWidthStatic;
	m_pointSizeStatic=vbo.m_pointSizeStatic;
	m_stippleFactor=vbo.m_stippleFactor;
	m_LineStipplePattern=vbo.m_LineStipplePattern;
	m_lightMode=vbo.m_lightMode;
	m_gel=0;
	m_coordinateType = vbo.m_coordinateType;
	setDirty(true);
	if(vbo.m_target_elements==&(vbo.m_elements))
		m_target_elements=&m_elements;
	else
		m_target_elements=&m_elementsShade;

	return *this;
}

mgVBO::~mgVBO(){
	if(m_gel){
		mgVBO* vbo=m_gel->m_VBO.release();
		assert(vbo==0 || vbo==this);
	}
	if(m_dname){
		MGDNameControl& dnc=getDNameControlInstance();
		dnc.deleteDlistMap(m_dname);
	}
}

///VBO data������������B
///(1) display mode��\���ɂ���
///(2) viewMode�ɏ]���AclearElements()������
///(3) clearStaticAttributes()���M
///(4) gel�p��vbo�ł���΁AdrawAttrib()���M
void mgVBO::initializeVBO(MGCL::VIEWMODE viewMode){
	set_display();
	if(viewMode!=MGCL::SHADING){
		clearElements(MGCL::WIRE);
	}
	if(viewMode!=MGCL::WIRE){
		clearElements(MGCL::SHADING);
	}
	clearStaticAttributes();
	if(m_gel)
		m_gel->drawAttrib(*this);
	m_target_elements=&m_elements;
}

///����mgVBOElement��null(���܂�draw/make_display_list()��������Ă��Ȃ��j���𔻒�
///mgVBO��m_gel�ɑ΂���make_display_list()�������Ȃ���Ă��Ȃ��Ƃ�false���Ԃ����
bool mgVBO::is_made(MGCL::VIEWMODE viewMode){
	if(viewMode==MGCL::HIGHLIGHT)
		viewMode=MGCL::WIREVIEW;
	else if(m_gel){
		if(m_gel->manifold_dimension()<=1)
			viewMode=MGCL::WIREVIEW;
	}

	switch(viewMode){
		case MGCL::WIRE: return !m_elementsDirty;
		case MGCL::SHADING: return !m_elementsShadeDirty;
		default:
			return !m_elementsDirty && !m_elementsShadeDirty;
	}
}

/// <summary>
/// Given top Group(File MGGroup) mgVBO, build MGGelPosition of this,
/// This must be built by MGOpenGLView.
/// </summary>
void mgVBO::buildGelPosition(mgVBO* vboGroup, MGGelPosition& gelp){
	MGObject* obj = dynamic_cast<MGObject*>(gel());//Lowest object pointer.
	if (!obj)
		return;//This must not happen(Lowest name mst be MGObject).

	std::vector<mgVBO*> vbos;
	if (!buildVBOHierarchy(*vboGroup, vbos))
		return;//This must not happen.

	size_t n = vbos.size(); assert(n >= 2);
	MGAttribedGel* grpA = vboGroup->gel(); assert(grpA);
	MGGroup* grp = static_cast<MGGroup*>(grpA);//Top MGGroup.
	gelp=MGGelPosition(grp, obj);
	for (size_t k = 1; k <= n - 2; k++) {
		mgVBO* vbok = vbos[k];
		MGAttribedGel* gl = vbok->gel();//Lower gel is MGGroup or MGShell.
		gelp.append_lower_gel(gl);
	}
}

///Clear all the data.
void mgVBO::clearElements(MGCL::DRAW_TARGET target){
	if(target!=MGCL::SHADING)
		m_elements.clear();
	if(target!=MGCL::WIRE)
		m_elementsShade.clear();
	m_builder.reset(0);
}

///type��mgVBOLeaf���ЂƂ�(Begin-End�łЂƂj�쐬���Ēǉ�����B
///target��Begin()-End()�ō쐬�����mgVBOLeaf���ǂ�element�ɓ���邩�������B
///target=WIRE:Wire(highlight�p��element(m_elements)
///target=SHADING:Shading�p��element(m_elementsShade).
///type�͎��̂��̂��Ƃ���(GL_xxxx��xxxx�������j
///POINTS,LINES,LINE_STRIP,LINE_LOOP,TRIANGLE_FAN,TRIANGLE_STRIP,QUAD_STRIP.
///Begin()--End()�łЂƂ�mgVBOLeaf�Ƃ���邱�Ƃ��ۏ؂����B
void mgVBO::Begin(GLenum type, MGCL::DRAW_TARGET target){
	setElementTarget(target);
	if(m_builder.get())
		assert(m_builder->is_null());

	float size=m_lineWidthStatic;
	if(type==GL_POINTS)
		size=m_pointSizeStatic;
	m_builder.reset(new mgVBOLeafBuilder(type,m_colorStatic,size));
	if(type==GL_LINES || type==GL_LINE_STRIP || type==GL_LINE_LOOP)
		m_builder->setLineStipple(m_stippleFactor,m_LineStipplePattern);
	m_builder->setLightMode(m_lightMode);
}

///Begin()�Ŏn�߂�VBO���ЂƂ�mgVBOLeaf�Ƃ��č쐬�A�ǉ�����B
///polygonMode��Begin()��type��GL_QUAD_STRIP, GL_TRIANGLES, GL_TRIANGLE_STRIP, GL_TRIANGLE_FAN
///�̍ۂɂ̂ݗL���ŁA����PolygonMode�iGL_POINT, GL_LINE, GL_FILL)���w�肷��B
///�쐬���ꂽmgVBOElement*���Ԃ���邪����mgVBOElement��mgVBO�����L���Ă���B
///null���ԋp���ꂽ�Ƃ��A�쐬�Ɏ��s�B
mgVBOLeaf* mgVBO::End(GLenum polygonMode){
	mgVBOLeaf* leaf=0;
	if(m_builder.get()){
		if(!m_builder->is_null()){
			leaf=new mgVBOLeaf(*m_builder);
			m_target_elements->emplace_back(leaf);
			leaf->setPolygonMode(polygonMode);
		}
	}
	m_builder.reset(0);
	return leaf;
}

///����Begin()����Ă��邩�ǂ����𒲂ׂ�BBegin()����AEnd()�̑O�ł����
///true��Ԃ�
bool mgVBO::is_InBegin(){
	if(m_builder.get())
		return true;
	return false;
}

///���_���Ƃ̐F���w�肷��
void mgVBO::Color(const MGColor& colr){
	if(!is_InBegin())
		return;
	m_builder->push_backColor(colr);
}
void mgVBO::Color3fv(const float colr[3]){
	if(!is_InBegin())
		return;
	float c[4];
	for(int i=0; i<3; i++) c[i]=colr[i];
	c[3]=1.f;
	vboColor vc(c);
	m_builder->push_backColor(vc);
}
void mgVBO::Color3dv(const double colr[3]){
	if(!is_InBegin())
		return;
	float c[4];
	for(int i=0; i<3; i++) c[i]=(float)colr[i];
	c[3]=1.f;
	vboColor vc(c);
	m_builder->push_backColor(vc);
}
void mgVBO::Color4fv(const float colr[4]){
	if(!is_InBegin())
		return;
	vboColor vc(colr);
	m_builder->push_backColor(vc);
}
void mgVBO::Color4ubv(const unsigned char rgba[4]){
	float r=float(rgba[0])/255.f;
	float g=float(rgba[1])/255.f;
	float b=float(rgba[2])/255.f;
	float a=float(rgba[3])/255.f;
	Color(MGColor(r,g,b,a));
}

///���_���Ƃ�normal���w�肷��
void mgVBO::Normal(const MGVector& norml){
	if(!is_InBegin())
		return;
	m_builder->push_backNormal(norml);
}
void mgVBO::Normal(float x, float y, float z){
	if(!is_InBegin())
		return;
	vboFPoint fp(x,y,z);
	m_builder->push_backNormal(fp);
}
void mgVBO::Normal3d(double x, double y, double z){
	if(!is_InBegin())
		return;
	vboFPoint fp(x,y,z);
	m_builder->push_backNormal(fp);
}
void mgVBO::Normal3fv(const float norml[3]){
	if(!is_InBegin())
		return;
	vboFPoint fp(norml[0],norml[1],norml[2]);
	m_builder->push_backNormal(fp);
}
void mgVBO::Normal3dv(const double norml[3]){
	if(!is_InBegin())
		return;
	vboFPoint fp(norml[0],norml[1],norml[2]);
	m_builder->push_backNormal(fp);
}

///���_�̍��W�l���w�肷��
void mgVBO::Vertex(const MGPosition& v){
	if(!is_InBegin())
		return;
	m_builder->push_backVertex(v);
}
void mgVBO::Vertex(float x, float y, float z){
	if(!is_InBegin())
		return;
	vboFPoint fp(x,y,z);
	m_builder->push_backVertex(fp);
}
void mgVBO::Vertex3d(double x, double y, double z){
	if(!is_InBegin())
		return;
	vboFPoint fp(x,y,z);
	m_builder->push_backVertex(fp);
}
void mgVBO::Vertex2fv(const float v[2]){
	if(!is_InBegin())
		return;
	vboFPoint fp(v[0],v[1],0.f);
	m_builder->push_backVertex(fp);
}
void mgVBO::Vertex3fv(const float v[3]){
	if(!is_InBegin())
		return;
	vboFPoint fp(v[0],v[1],v[2]);
	m_builder->push_backVertex(fp);
}
void mgVBO::Vertex2dv(const double v[2]){
	if(!is_InBegin())
		return;
	vboFPoint fp(v[0],v[1],0.);
	m_builder->push_backVertex(fp);
}
void mgVBO::Vertex3dv(const double v[3]){
	if(!is_InBegin())
		return;
	vboFPoint fp(v[0],v[1],v[2]);
	m_builder->push_backVertex(fp);
}

///Texture�̍��W�l���w�肷��
void mgVBO::TexCoord(const MGPosition& v){
	if(!is_InBegin())
		return;
	m_builder->push_backTexture(v);
}
void mgVBO::TexCoord(float s, float t){
	if(!is_InBegin())
		return;
	vboFP2D fp(s,t);
	m_builder->push_backTexture(fp);
}
void mgVBO::TexCoord2d(double s, double t){
	if(!is_InBegin())
		return;
	vboFP2D fp(s,t);
	m_builder->push_backTexture(fp);
}
void mgVBO::TexCoord2fv(const float v[2]){
	if(!is_InBegin())
		return;
	vboFP2D fp(v[0],v[1]);
	m_builder->push_backTexture(fp);
}
void mgVBO::TexCoord2dv(const double v[2]){
	if(!is_InBegin())
		return;
	vboFP2D fp(v[0],v[1]);
	m_builder->push_backTexture(fp);
}

///Static attribute��ݒ肷��B
///begin() - end()�̊Ԃł���΂���mgVBOLeaf�ɑ΂��Ă����K�p����
///begin() - end()�̊O�ł���΂��̌��mgVBOLeaf���ׂĂɓK�p�����B
void mgVBO::setStaticAttribColor(const MGColor& color){
	if(is_InBegin()){
		//Begin()�̊Ԃł�setStaticAttribColor.
		m_builder->setStaticAttribColor(color);
	}else{
		//Begin()�̊O�ł�setStaticAttribColor.
		m_colorStatic=color;
	}
}
void mgVBO::setStaticAttribColor(MGColor::ColorID id){
	const MGColor& color=MGColor::get_instance(id);
	setStaticAttribColor(color);
}
void mgVBO::setStaticAttribColor(const float color[4]){
	MGColor clr(color[0],color[1],color[2],color[3]);
	if(m_builder.get()){
		//Begin()�̊Ԃł�setStaticAttribColor.
		m_builder->setStaticAttribColor(clr);
	}else{
		//Begin()�̊O�ł�setStaticAttribColor.
		m_colorStatic=clr;
	}
}
void mgVBO::setStaticAttribColor(float r, float g, float b){
	MGColor clr(r,g,b);
	if(is_InBegin()){
		//Begin()�̊Ԃł�setStaticAttribColor.
		m_builder->setStaticAttribColor(clr);
	}else{
		//Begin()�̊O�ł�setStaticAttribColor.
		m_colorStatic=clr;
	}
}
void mgVBO::setStaticAttribLineWidth(GLfloat size){
	if(is_InBegin()){
		//Begin()�̊Ԃł�setStaticAttribLineWidth.
		if(m_builder->typeBegin()!=GL_POINTS)
			m_builder->setStaticAttribSize(size);
	}else{
		//Begin()�̊O�ł�setStaticAttribSize.
		m_lineWidthStatic=size;
	}
}
void mgVBO::setStaticAttribPointSize(GLfloat size){
	if(is_InBegin()){
		//Begin()�̊Ԃł�setStaticAttribPointSize.
		if(m_builder->typeBegin()==GL_POINTS)
			m_builder->setStaticAttribSize(size);
	}else{
		//Begin()�̊O�ł�setStaticAttribSize.
		m_pointSizeStatic=size;
	}
}

///Static Attributes �����ׂ�default�ɂ��ǂ�(display/noDisplay�͑ΏۊO�j
void mgVBO::clearStaticAttributes(){
	m_colorStatic.set_undefined();
	m_lineWidthStatic=-1.f;
	m_pointSizeStatic=-1.f;
	m_stippleFactor=-1;
	m_lightMode=-1;
}


///Line stipple�������Z�b�g����B
///When factor=0 is input, line pattern is disabled. This means the line is solid.
///When factor<0, the stipple attribute is undefined. This means the attribute
///is defined by the environment.
///When factor<=0, pattern is unnecessary.
void mgVBO::setLineStipple(short int factor, GLushort pattern){
	if(is_InBegin()){
		//Begin()�̊Ԃł�setLineStipple.
		unsigned type=m_builder->typeBegin();
		if(type==GL_LINES || type==GL_LINE_STRIP || type==GL_LINE_LOOP)
			m_builder->setLineStipple(factor,pattern);
	}else{
		m_stippleFactor=factor;
		m_LineStipplePattern=pattern;
	}
}
void mgVBO::disableLinePattern(){
	if(is_InBegin()){
		//Begin()�̊Ԃł�setLineStipple.
		unsigned type=m_builder->typeBegin();
		if(type==GL_LINES || type==GL_LINE_STRIP || type==GL_LINE_LOOP)
			m_builder->setLineStipple(0,0);
	}else{
		m_stippleFactor=0;
		m_LineStipplePattern=0;
	}
}

///m_gel�̕`��f�[�^�쐬�݂̂������Ȃ��B
///���łɍ쐬�ς݂ł����Ă������I�ɍč쐬���s���B
void mgVBO::make_display_list(
	MGCL::VIEWMODE vmode
){
	if(m_gel){
		if(vmode!=MGCL::SHADING){
			m_gel->make_display_list(MGCL::WIREVIEW);
		}
		if(vmode!=MGCL::WIRE){
			if(!m_gel->displayList_is_made(MGCL::WIREVIEW))
				m_gel->make_display_list(MGCL::WIREVIEW);
			m_gel->make_display_list(MGCL::SHADINGVIEW);
		}
	}else
		initializeVBO(vmode);
}

void mgVBO::execStaticGLAttrib(){
	if(m_colorStatic.defined()){
		mgGLSL::execStaticColorAttrib(m_colorStatic);
	}
	if(m_lineWidthStatic>0.)
		mgGLSL::execStaticLineWidth(m_lineWidthStatic);

	if(m_stippleFactor>0)
		mgGLSL::execStaticLineStipple((GLint)m_stippleFactor,m_LineStipplePattern);

	if(m_pointSizeStatic>0.){
		mgGLSLProgram* glslP=mgGLSLProgram::getCurrentGLSLProgram();
		glEnable(GL_PROGRAM_POINT_SIZE);
		glslP->setUniform(mgGLSLProgram::pointSize,m_pointSizeStatic);
	}
	if(m_lightMode>=0){
		mgGLSL::execLightMode(m_lightMode);
	}
}

///�`��֐�draw()�́A!is_made()�ł���΁A�쐬���A�\������������B
///is_made()(�`��f�[�^�쐬�ς݁j�ł���΁A���łɍ쐬���ꂽmgVBOElement�̕`����s���B
void mgVBO::draw(MGCL::VIEWMODE viewMode){
	if(is_no_display())
		return;

	if(viewMode!=MGCL::SHADING)
		if(!is_made(MGCL::WIREVIEW))
			make_display_list(MGCL::WIREVIEW);
	if(viewMode==MGCL::SHADING || viewMode==MGCL::WIRE_AND_SHADING)
		if(!is_made(MGCL::SHADINGVIEW))
			make_display_list(MGCL::SHADINGVIEW);

	//save the coordinate type and update to m_coordinateType.
	mgCoordinateTypeSwitcher coordType(m_coordinateType, getAnchorPosition());

mgGLSL::pushStaticGLAttrib();
	if(viewMode!=MGCL::HIGHLIGHT)
		execStaticGLAttrib();

	///Process for mgVBOLeaf in m_elementsShade.
	if(viewMode==MGCL::SHADING || viewMode==MGCL::WIRE_AND_SHADING){
		mgGLSL::execLightMode(m_lightMode);
		size_t nShade=m_elementsShade.size();
		for(size_t i=0; i<nShade; i++){
			UniqueVBOElement& elmi=m_elementsShade[i];
			elmi->draw(viewMode);
		}
	}
	
	///Process for mgVBOPointer and mgVBOLeaf in m_elements.
	size_t n=m_elements.size();
	for(size_t i=0; i<n; i++){
		UniqueVBOElement& elmi=m_elements[i];
		mgVBO* vboi=elmi->vboPointer();
		if(vboi)//if mgVBOPointer.
			vboi->draw(viewMode);
		else{//If mgVBOLeaf.
			if(!m_gel || m_gel->manifold_dimension()<2 || viewMode!=MGCL::SHADING)
				elmi->draw(viewMode);
			else{
				if(dynamic_cast<MGGroup*>(m_gel) || dynamic_cast<MGShell*>(m_gel))
					elmi->draw(viewMode);
				else
					continue;
			}
		}
	}
mgGLSL::popStaticGLAttrib();

}

///draw()��mgVBOLeaf���쐬�ς݁inot null)�ł���΍쐬�������s��Ȃ����A
///redraw()�͋����I�ɍč쐬���s���`�揈���������Ȃ��B
void mgVBO::redraw(MGCL::VIEWMODE viewMode){
	make_display_list(viewMode);
	draw(viewMode);
}

///gel��mgVBOPointer���쐬��mgVBOElement�Ƃ��Ēǉ�����B
///gel must be valid when draw event happens since gel is referenced at that time.
void mgVBO::drawGel(const MGAttribedGel& gel){
	mgVBO* vbo=gel.dlist_name();
	if(vbo)
		m_elements.emplace_back(new mgVBOPointer(*vbo));
}

///gel��mgVBOLeafPointer���쐬��mgVBOElement�Ƃ��Ēǉ�����B
void mgVBO::drawVBOLeaf(
	const mgVBOLeaf& leaf,
	MGCL::DRAW_TARGET target///<When target=WIRE, built elements are
			///stored as wire mode display, else as shading mode display
){
	if(target==MGCL::SHADING || target==MGCL::WIRE_AND_SHADING)
		m_elementsShade.emplace_back(new mgVBOLeafPointer(leaf));
	if(target==MGCL::WIRE || target==MGCL::WIRE_AND_SHADING)
		m_elements.emplace_back(new mgVBOLeafPointer(leaf));
}

void mgVBO::setGel(const MGAttribedGel* gel){
	MGAttribedGel* gel2=const_cast<MGAttribedGel*>(gel);
	m_gel=gel2;
}

///gel��mgVBOPointer�������o�[����O��
void mgVBO::deleteGel(const MGAttribedGel& gel){
	const mgVBO* vboGel=gel.m_VBO.get();
	if(!vboGel)
		return;

	int n=int(m_elements.size());
	for(int i=n-1; i>=0; i--){
		UniqueVBOElement& elmi=m_elements[i];
		mgVBOPointer* vboptr=dynamic_cast<mgVBOPointer*>(elmi.get());
		if(!vboptr)
			continue;
		if(vboptr->vboPointer()==vboGel){
			m_elements.erase(m_elements.begin()+i);
			break;
		}
	}
}

///hilight�����ŕ\������
void mgVBO::highlight(){
	draw(MGCL::HIGHLIGHT);
}

///Selection�ɐݒ肷�閼�O�����߂�B=0�̂Ƃ��A���O�̐ݒ菈�������Ȃ��B
///0(null) �͂���mgVBOElement��mgVBOLeaf�ł���A���O�������Ȃ����Ƃ������B
///���O��mgVBO����������
GLuint mgVBO::getSelectionName()const{
	return m_dname;
}

///�`��֐�selectionDraw()�́AObject�I���̂��߂̕\������������B
///�ʏ��draw�Ƃ̑���F///Color�Ƃ���m_bufferID��p���Asize�����ȊO��
///attributes�̏����inormal, texture, color)�����Ȃ��B
void mgVBO::selectionDraw(MGCL::VIEWMODE viewMode){
	if(is_no_display())
		return;

	if(!is_made(viewMode))
		return;

	unsigned nameSelect=getSelectionName();
	mgGLSL::setColorAsSelectionName(nameSelect);
	mgGLSLProgram* glsl=mgGLSLProgram::getCurrentGLSLProgram();
	glsl->setFuncType(mgGLSL::Select);

	//save the coordinate type.
	mgCoordinateTypeSwitcher coordType(m_coordinateType, getAnchorPosition());
	///Process for mgVBOLeaf.
	if(viewMode==MGCL::SHADING || viewMode==MGCL::WIRE_AND_SHADING){
		size_t nShade=m_elementsShade.size();
		for(size_t i=0; i<nShade; i++){
			UniqueVBOElement& elmi=m_elementsShade[i];
			elmi->selectionDraw(viewMode);
		}
	}
	
	///Process for mgVBOPointer and mgVBOLeaf.
	size_t n=m_elements.size();
	for(size_t i=0; i<n; i++){
		UniqueVBOElement& elmi=m_elements[i];
		mgVBO* vboi=elmi->vboPointer();
		if(vboi)//if mgVBOPointer.
			vboi->selectionDraw(viewMode);
		else{//If mgVBOLeaf.
			if(!m_gel || m_gel->manifold_dimension()<2 || viewMode!=MGCL::SHADING)
				elmi->selectionDraw(viewMode);
			else{

				if(dynamic_cast<MGGroup*>(m_gel) || dynamic_cast<MGShell*>(m_gel))
					elmi->selectionDraw(viewMode);
				else
					continue;
			}
		}
	}
}

///Build VBO hierarchy.
///Let n=vbos.size(), then 
///vbos[i] includes vbos[i+1] as mgVBOPointer for i=0,...,n-2.
///vbos[0]=&parent, and vbos[n-1] = this mgVBO.
bool mgVBO::buildVBOHierarchy(
	mgVBO& parent,//Parent mgVBO that may include this vbo as its child
		/// or the descendant.
	std::vector<mgVBO*>& vbos
){
	vbos.push_back(&parent);
	size_t n=parent.m_elements.size();
	for(size_t i=0; i<n; i++){
		UniqueVBOElement& elmi=parent.m_elements[i];
		mgVBO* vboi=elmi->vboPointer();
		if(!vboi)
			continue;
		if(vboi==this){
			vbos.push_back(vboi);
			return true;
		}else{
			if(buildVBOHierarchy(*vboi,vbos))
				return true;
		}
	}
	vbos.pop_back();
	return false;
}

GLfloat mgVBO::getPointSize()const{
	return m_pointSizeStatic;
}

GLfloat mgVBO::getLineWidth()const{
	return m_lineWidthStatic;
}

///Set dirty flag(s) of m_elements or m_elementsShade.
void mgVBO::setDirty(bool is_dirty, MGCL::DRAW_TARGET target)const{
	if(target!=MGCL::WIRE)
		m_elementsShadeDirty=is_dirty;
	if(target!=MGCL::SHADING)
		m_elementsDirty=is_dirty;
}

///Function ID���Z�b�g����.
///setFunctionID()��Begin()-End()�̊Ԃł̂ݗL���B
//void mgVBO::setFunctionID(int functionID){
//	assert(is_InBegin());
//	if(!is_InBegin())
//		return;
//	m_builder->setFunctionID(functionID);
//}
void mgVBO::setDrawType(mgGLSL::DrawType drawType){
	assert(is_InBegin());
	if(!is_InBegin())
		return;
	m_builder->setDrawType(drawType);
}

///Texture��set/get����B
///texture�͎Q�Ƃ���̂�
///setTexture()��Begin()-End()�̊Ԃł̂ݗL���B
void mgVBO::setTexture(mgTexture* texture){
	assert(is_InBegin());
	if(!is_InBegin())
		return;
	m_builder->setTexture(texture);
}

///Set Elements target.
void mgVBO::setElementTarget(MGCL::DRAW_TARGET target){
	if(target==MGCL::SHADING)
		m_target_elements=&m_elementsShade;
	else
		m_target_elements=&m_elements;
}
