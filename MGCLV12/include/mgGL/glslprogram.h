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
/// 描画タイプ(FuncTypeがstandardのときのみ有効)
/// Primitive:通常のプリミティブ
/// Texture:テクスチャ描画
/// gridTexture:点のスプライト描画
typedef enum {
	Primitive = 0,
	Texture = 11,
	gridTexture=12,
} DrawType;

/// 描画タイプ
/// standard:通常描画
/// Select:Select処理用描画
/// Analysis:解析表示用描画
typedef enum {
	standard=0,
	Select = 1,
	Analysis = 5,
} FuncType;

/// 座標変換タイプ(Vertexの座標系)
/// World:ワールド座標系
/// NdcScreen:NDC空間座標系
/// AnchorPoint:アンカーポイントからの相対距離(アンカーポイントはワールド座標系で指定)
/// AnchorPointScreen:アンカーポイントからの相対距離(アンカーポイントはNDC座標系で指定)
typedef enum {
	World = 0,
	NdcScreen=2,
	AnchorPoint = 3,
	AnchorPointScreen = 4,
}CoordinateType;

/// 照光処理
/// NoShading:照光なし
/// Shading:照光あり
typedef enum {
	NoShading=0,
	Shading=1,
} ShadeMode;

/// ゼブラ処理タイプ
/// ZebraVert:垂直方向縞
/// ZebraHorizon:水平方向縞
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

///mgGLSLProgramはOpenGL Shader Programをcompile , linkしてそのuniform変数の管理を行います.

///利用方法：
///１、mgGLSLProgramオブジェクトの作成
///２、compileShaderFromFile、またはcompileShaderFromStringで必要なshaderを
///    コンパイル
///３、必要であればbindAttribLocation(), bindFragDataLocation()を呼んだ後、link()
///　use()によりshaderプログラムを選択し、setUniform()によりuniform変数に値をセット
///
///***** MGCLではshaderに対して次のlocationと変数名を予約している：
///locationに対してこれらの値をは使用してはならない。
/// また、これらの変数はすべてのvertex shaderに下記のように宣言されなければならない。
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

	// Lightまわり
	static const GLint LIGHT_NUM = 10;
	static const GLint PROP_NUM = 11;

	/// CreateShaderのタイプ列挙子
     typedef enum{
        VERTEX, FRAGMENT, GEOMETRY,
        TESS_CONTROL, TESS_EVALUATION
    } GLSLShaderType;

	///Vertex attrib変数定義用
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

		FuncType,               // draw function(Draw/Select/Zebra Analysis) FuncTypeでEnum定義
		DrawType,				// func typeがDrawのときのみ有効 DrawTypeでEnum定義
		ShaderMode,				// int 0:NoShading 1:Shading

		CoordinateType,			// CoordinateType (World/NDC/AnchorPoint/AnchorPointScreen) CoordinateTypeでEnum定義

		texture2D,				// texture.
		pointSize,				// Vertexの表示サイズ。カーソルなどのPOINT_SPLITE表示に使う

		anchorPoint,			// AnchorPoint。Billboard表示(拡大。縮小してもサイズが変わらない描画)の基準点

		LightTwoSides,			// boolean True:TwoSide, False:片面
		ForceLight,				// falseのときは、有効なLightが無い場合は、NO_SHADINGモードで描画。

		ZebraAxis,				// ゼブラ方向 ZebraTypeで指定
		ZebraSize,				// ゼブラの幅。 ゼブラのステップ刻みは step(0.5f, fract(1.0f / ZebraSize * atan(r.y, r.x)));
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

	///compileがうまく行われたときはtrueを、失敗したときはfalseを返します
    bool compileShaderFromFile(const char* fileName, GLSLShaderType type);

    ///compileがうまく行われたときはtrueを、失敗したときはfalseを返します
	bool compileShaderFromString(const std::string& source, GLSLShaderType type);

	///Delete this program and initialize this.
	void freeProgram();

    ///linkがうまく行われたときはtrueを、失敗したときはfalseを返します
	bool link();

	///このmgGLSLProgramの利用を開始します。
    void use();

	///エラーが起こった時のログ内容を文字で求めます
    std::string& log();

	///Programハンドルを返します
    int getHandle();

	///Link済か否かを尋ねる
    bool isLinked();

	///nameに割り当てるlocationを指定する。locationはこの後のlink()で有効となる。
    void bindAttribLocation( GLuint location, const char * name);

	///Current mgGLSLProgramを求める。
	static void setCurrentGLSLProgram(mgGLSLProgram* glsl);

	///Current mgGLSLProgramを求める。
	static mgGLSLProgram* getCurrentGLSLProgram();

	///nameに割り当てるlocationを指定する。locationはこの後のlink()で有効となる。
	void bindFragDataLocation( GLuint location, const char * name );

	///規定のvetex attrib のlocationを求める。
	int getvPositionLocation()const;///vPosition
	int getvColorLocation()const;///vColor
	int getvNormalLocation()const;///vNormal
	int getvTextureCoordLocation()const;///vTexture

	/// Get Uniform location.
	GLint getUniformLocation(const char* name)const;
	GLint getUniformLocation(const std::string& name)const{
		return getUniformLocation(name.c_str());
	};

	/// Uniformに値をセット(invoke glUniform)する関数群。

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
	//m_uniform_locations[]に値をセットする(link()内の処理）
	void build_attribLocations();

	//m_uniform_locationsLights[]に値をセットする(ling()内の処理）
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
