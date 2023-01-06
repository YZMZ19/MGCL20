#pragma once

#ifndef _CONSOLE

#include <string>
#include <glm/glm.hpp>
#include <iostream>
#include "mg/MGCL.h"
#include "mg/Position.h"

class mgStaticGLAttrib;
class MGColor;
class mgGLSLProgram;

/** @file */

/** @addtogroup DisplayHandling
 *  @{
 */

///mgGLSL is a namespace for OpenGL handling, maily used by mgGLSLProgram.
namespace mgGLSL{

/// function ID
/// �`��^�C�v(FuncType��standard�̂Ƃ��̂ݗL��)
/// Primitive:�ʏ�̃v���~�e�B�u
/// Texture:�e�N�X�`���`��
/// gridTexture:�_�̃X�v���C�g�`��
typedef enum {
	Primitive = 0,
	Texture = 11,
	gridTexture=12,
} DrawType;

/// �`��^�C�v
/// standard:�ʏ�`��
/// Select:Select�����p�`��
/// Analysis:��͕\���p�`��
typedef enum {
	standard=0,
	Select = 1,
	Analysis = 5,
} FuncType;

/// ���W�ϊ��^�C�v(Vertex�̍��W�n)
/// World:���[���h���W�n
/// NdcScreen:NDC��ԍ��W�n
/// AnchorPoint:�A���J�[�|�C���g����̑��΋���(�A���J�[�|�C���g�̓��[���h���W�n�Ŏw��)
/// AnchorPointScreen:�A���J�[�|�C���g����̑��΋���(�A���J�[�|�C���g��NDC���W�n�Ŏw��)
typedef enum {
	World = 0,
	NdcScreen=2,
	AnchorPoint = 3,
	AnchorPointScreen = 4,
}CoordinateType;

/// �ƌ�����
/// NoShading:�ƌ��Ȃ�
/// Shading:�ƌ�����
typedef enum {
	NoShading=0,
	Shading=1,
} ShadeMode;

/// �[�u�������^�C�v
/// ZebraVert:����������
/// ZebraHorizon:����������
typedef enum {
	ZebraVert = 0,
	ZebraHorizon = 1,
}ZebraType;

///change error code to string and print it.
MG_DLL_DECLR void printOpenGLError(int errorCode);

/// Returns 1 if an OpenGL error occurred, 0 otherwise.
MG_DLL_DECLR int checkForOpenGLError(const char*, int);

/// Dump OpenGL version info.
MG_DLL_DECLR void dumpGLInfo(bool dumpExtensions = false);

/// get LightEnable (Shade On/Off) data.
MG_DLL_DECLR bool LightEnabled();

MG_DLL_DECLR void CALLBACK debugCallback(GLenum source,
                                           GLenum type,
                                           GLuint id,
                                           GLenum severity,
                                           GLsizei length,
                                           const GLchar* message,
                                           void* userParam);

///Initialize StaticGLAttribStack.
MG_DLL_DECLR void initializeStaticGLAttribStack();

///Set static color as selection name.
void setColorAsSelectionName(unsigned name);

std::stack<mgStaticGLAttrib>&  getStaticGLAttribStack();
mgStaticGLAttrib& getCurrentStaticGLAttrib();

void execStaticGLAttrib(const mgStaticGLAttrib& attrib);
void execStaticColorAttrib(const MGColor& color);
void execStaticColorAttrib(const float color[4]);
void execStaticLineWidth(float lineWidth);
void execStaticLineStipple(short int factor, GLuint pattern);
void execLightMode(// set Light(Shading On/Off)
	int mode/// <0: undefined, =0:Light is disabled, >0:Light is enabled.
);
void pushStaticGLAttrib();
void popStaticGLAttrib();

};

///mgGLSLProgram��OpenGL Shader Program��compile , link���Ă���uniform�ϐ��̊Ǘ����s���܂�.

///���p���@�F
///�P�AmgGLSLProgram�I�u�W�F�N�g�̍쐬
///�Q�AcompileShaderFromFile�A�܂���compileShaderFromString�ŕK�v��shader��
///    �R���p�C��
///�R�A�K�v�ł����bindAttribLocation(), bindFragDataLocation()���Ă񂾌�Alink()
///�@use()�ɂ��shader�v���O������I�����AsetUniform()�ɂ��uniform�ϐ��ɒl���Z�b�g
///
///***** MGCL�ł�shader�ɑ΂��Ď���location�ƕϐ�����\�񂵂Ă���F
///location�ɑ΂��Ă����̒l���͎g�p���Ă͂Ȃ�Ȃ��B
/// �܂��A�����̕ϐ��͂��ׂĂ�vertex shader�ɉ��L�̂悤�ɐ錾����Ȃ���΂Ȃ�Ȃ��B
///uniform int functionID;//Function ID of the shaders:
///   0: standard(no shading functions);
///  11: Texture that uses texture1 sampler2D.(MGPlaneImage Texture uses this one.)
///	
/// uniform mat4 modelViewProjMatrix;//=projMatrix*modelViewMatrix.
/// uniform mat3 normalMatrix;
/// uniform mat4 modelViewMatrix;
/// uniform mat4 projMatrix;
/// uniform sampler2D texture1;
///
/// layout(location = 0) in vec3 vPosition;
/// layout(location = 1) in vec4 vColor;
/// layout(location = 2) in vec3 vNormal;
/// layout(location = 3) in vec2 vTextureCoord;
class MG_DLL_DECLR mgGLSLProgram{

public:
	//m_uniform_locations(UniformName) array size.
	static const size_t MAX_ATTRIB_LOC = 20;

	// Light�܂��
	static const GLint LIGHT_NUM = 10;
	static const GLint PROP_NUM = 11;

	/// CreateShader�̃^�C�v�񋓎q
     typedef enum{
        VERTEX, FRAGMENT, GEOMETRY,
        TESS_CONTROL, TESS_EVALUATION
    } GLSLShaderType;

	///Vertex attrib�ϐ���`�p
	typedef enum {
		vPosition=0,
		vColor,
		vNormal,
		vTextureCoord
	} VertexAttribId;

	union SELECT_NAME{
		unsigned uiName;
		unsigned short usName[2];
		unsigned char ucName[4];//
	};

	///Uniform name(id).
	//Last name of UNIFORM_ATTRIB_ID must be less than MAX_ATTRIB_LOC.
	typedef enum{
		modelViewProjMatrix = 0,	// mat4
		modelViewMatrix,		// mat4
		projMatrix,				// mat4
		normalMatrix,			// mat3 

		ndcMarix,				// mat4
		ndcScaleMatrix,			// mat3
		dpiFactor,				// float

		FuncType,               // draw function(Draw/Select/Zebra Analysis) FuncType��Enum��`
		DrawType,				// func type��Draw�̂Ƃ��̂ݗL�� DrawType��Enum��`
		ShaderMode,				// int 0:NoShading 1:Shading

		CoordinateType,			// CoordinateType (World/NDC/AnchorPoint/AnchorPointScreen) CoordinateType��Enum��`

		texture2D,				// texture.
		pointSize,				// Vertex�̕\���T�C�Y�B�J�[�\���Ȃǂ�POINT_SPLITE�\���Ɏg��

		anchorPoint,			// AnchorPoint�BBillboard�\��(�g��B�k�����Ă��T�C�Y���ς��Ȃ��`��)�̊�_

		LightTwoSides,			// boolean True:TwoSide, False:�Ж�
		ForceLight,				// false�̂Ƃ��́A�L����Light�������ꍇ�́ANO_SHADING���[�h�ŕ`��B

		ZebraAxis,				// �[�u������ ZebraType�Ŏw��
		ZebraSize,				// �[�u���̕��B �[�u���̃X�e�b�v���݂� step(0.5f, fract(1.0f / ZebraSize * atan(r.y, r.x)));
	}UniformName;

	///Light property.
	typedef enum {
		isEnabled = 0,	// bool
		ambientColor,	// vec4 ambient Color
		diffuseColor,	// vec4 diffulse color
		specularColor,	// vec4 specular color
		position,		// vec4 position
		spotDirection,	// vec3
		spotExponent,	// float
		spotCutoff,		// float
		constantAttenuation,	// float
		linearAttenuation,		// float
		quadraticAttenuation	// float
	}LightProps;

private://Member data.
	static mgGLSLProgram* m_CurrrentGLSL;//The current glslProgram to use.

	int  m_handle;///<Shader program handler.
	GLint m_VersionMajor, m_VersionMinor;//OpenGL version.
	bool m_linked;//If the program is linked or not.
	std::string m_logString;

	//Tehese values are set when mgGLSLProgram::link() is invoked.
	GLint m_uniform_locations[MAX_ATTRIB_LOC];//uniform loc of UniformName.
	GLint m_LightPropsIds[LIGHT_NUM][PROP_NUM];//uniform loc of lights.

public:
	mgGLSLProgram();
	~mgGLSLProgram();

	///compile�����܂��s��ꂽ�Ƃ���true���A���s�����Ƃ���false��Ԃ��܂�
    bool compileShaderFromFile(const char* fileName, GLSLShaderType type);

    ///compile�����܂��s��ꂽ�Ƃ���true���A���s�����Ƃ���false��Ԃ��܂�
	bool compileShaderFromString(const std::string& source, GLSLShaderType type);

	///Delete this program and initialize this.
	void freeProgram();

    ///link�����܂��s��ꂽ�Ƃ���true���A���s�����Ƃ���false��Ԃ��܂�
	bool link();

	///����mgGLSLProgram�̗��p���J�n���܂��B
    void use();

	///�G���[���N���������̃��O���e�𕶎��ŋ��߂܂�
    std::string& log();

	///Program�n���h����Ԃ��܂�
    int getHandle();

	///Link�ς��ۂ���q�˂�
    bool isLinked();

	///name�Ɋ��蓖�Ă�location���w�肷��Blocation�͂��̌��link()�ŗL���ƂȂ�B
    void bindAttribLocation( GLuint location, const char * name);

	///Current mgGLSLProgram�����߂�B
	static void setCurrentGLSLProgram(mgGLSLProgram* glsl);

	///Current mgGLSLProgram�����߂�B
	static mgGLSLProgram* getCurrentGLSLProgram();

	///name�Ɋ��蓖�Ă�location���w�肷��Blocation�͂��̌��link()�ŗL���ƂȂ�B
	void bindFragDataLocation( GLuint location, const char * name );

	///�K���vetex attrib ��location�����߂�B
	int getvPositionLocation()const;///vPosition
	int getvColorLocation()const;///vColor
	int getvNormalLocation()const;///vNormal
	int getvTextureCoordLocation()const;///vTexture

	/// Get Uniform location.
	GLint getUniformLocation(const char* name)const;
	GLint getUniformLocation(const std::string& name)const{
		return getUniformLocation(name.c_str());
	};

	/// Uniform�ɒl���Z�b�g(invoke glUniform)����֐��Q�B

	/// loc is a uniform location obtained by getUniformLocation,
	/// which is stored in m_uniform_locations
    void setUniform( GLint loc, float x, float y, float z);
    void setUniform( GLint loc, const glm::vec3& v);
    void setUniform( GLint loc, const glm::vec4& v);
    void setUniform( GLint loc, const glm::mat4& m);
    void setUniform( GLint loc, const glm::mat3& m);
    void setUniform( GLint loc, float val);
    void setUniform( GLint loc, int val);
    void setUniform( GLint loc, bool val);

	///Functions to store data in uniform variables(invoke glUniform) through UniformName.
	void setUniform(UniformName name, float x, float y, float z);
	void setUniform(UniformName name, const glm::vec3& v);
	void setUniform(UniformName name, const glm::vec4& v);
	void setUniform(UniformName name, const glm::mat4& m);
	void setUniform(UniformName name, const glm::mat3& m);
	void setUniform(UniformName name, float val);
	void setUniform(UniformName name, int val);
	void setUniform(UniformName name, bool val);

	///Enable lights when bEnabled=true, else disable.
	void EnableLights(bool bEnabled = true);

	///Test if light is enabled.
	bool LightEnabled();

	///Functions to set light data in uniform variables(invoke glUniform).
    void setUniformLights( GLint lightNo, LightProps name, float x, float y, float z);
    void setUniformLights( GLint lightNo, LightProps name, const glm::vec3& v);
    void setUniformLights( GLint lightNo, LightProps name, const glm::vec4& v);
    void setUniformLights( GLint lightNo, LightProps name, const glm::mat4& m);
    void setUniformLights( GLint lightNo, LightProps name, const glm::mat3& m);
    void setUniformLights( GLint lightNo, LightProps name, float val);
    void setUniformLights( GLint lightNo, LightProps name, int val);
    void setUniformLights( GLint lightNo, LightProps name, bool val);

	///Get version by glGetInteger(), and set the info in this program.
	void setOpenGLVersion();

	///Get OpendGL version obtained by setOpenGLVersion.
	void getOpenGLVerion(GLint& major, GLint& minor)const;

	///Print each info.
	void printActiveUniforms();
	void printActiveAttribs();
	void printPrjMatrix();
	void printLightProps();
	
	///Set & Get function type.
	void setFuncType(mgGLSL::FuncType type);
	int getFuncType();

	///Set & Get coordinate type.
	void setCoordinateType(mgGLSL::CoordinateType type);
	mgGLSL::CoordinateType getCoordinateType();

	///Get Anchoe position.
	void getAnchor(MGPosition& anchor);

private:
	//m_uniform_locations[]�ɒl���Z�b�g����(link()���̏����j
	void build_attribLocations();

	//m_uniform_locationsLights[]�ɒl���Z�b�g����(ling()���̏����j
	void build_attribLocationsLight();

    int  getAttribLocation(const char* name)const;
	bool fileExists(const std::string& fileName)const;
};

///Utility class to invoke glsl's setFuncType.

///mgFuncTypeSwitcher saves the current function type and invoke input function type.
///When mgFuncTypeSwitcher is destructed, the saved original type is restored.
class MG_DLL_DECLR mgFuncTypeSwitcher{
public:
	mgFuncTypeSwitcher(mgGLSL::FuncType type);

	~mgFuncTypeSwitcher();

private:
	mgGLSLProgram* m_pGLSL;//The current glsl program is saved.
	int m_orgFuncType;//The original function type is saved.
};

///mgCoordinateTypeSwitcher saves the current coordinate type and invoke input coordinate type.

///When mgCoordinateTypeSwitcher is destructed, the saved original type is restored.
class MG_DLL_DECLR mgCoordinateTypeSwitcher{
public:
	mgCoordinateTypeSwitcher(mgGLSL::CoordinateType coordinateType, const MGPosition* anchorP = nullptr);

	~mgCoordinateTypeSwitcher();

private:
	mgGLSLProgram* m_pGLSL;//The current glsl program is saved.
	mgGLSL::CoordinateType m_coordinateType;
	MGPosition m_anchorPoint;//Valid only when m_coordinateType=AnchorPoint or AnchorPointScreen.
};

#else //_CONSOLE

class MGColor;
class mgStaticGLAttrib;
namespace mgGLSL {
	typedef enum {NoShading = 0,Shading = 1,} ShadeMode;
	void execStaticColorAttrib(const MGColor& color);
	void execStaticGLAttrib(const mgStaticGLAttrib& attrib);
	void execStaticColorAttrib(const float color[4]);
	void execStaticLineWidth(float lineWidth);
	void execStaticLineStipple(short int factor, GLuint pattern);
	void execLightMode(int mode);
};

#endif //_CONSOLE

/** @} */ // end of DisplayHandling group
