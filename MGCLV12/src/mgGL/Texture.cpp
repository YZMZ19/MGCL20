/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno             */
/* All rights reserved.                                             */
/********************************************************************/

/**
 * @file Texture.h
 * @brief クラス MGGridCursor の宣言
 */

#include "StdAfx.h"
#include "mgGL/Texture.h"
#include "mgGL/GLSLProgram.h"
#include "mgGL/Image.h"
namespace{int glErr;};

///Texture definition class.
///Before use of mgTexture, set_image must be invoked.
///set_sampler is optionally invoked. If set_sampler is not invoked,
///uniform of spampler2D  "texture1" is assumed.
///When set_textureUnit is not invoked, unit number 0 is assumed.

const std::string& mgTexture::getDefaultSamplerVariable()const{
	static const std::string m_DefaultSamplerVariable="texture1";
	return m_DefaultSamplerVariable;
}

/// デフォルトコンストラクター
///Before use of mgTexture, set_image and set_sampler must be invoked.
mgTexture::mgTexture(
	GLenum target///Currently only GL_TEXTURE_2D is allowed.
):m_target(target),m_width(0), m_height(0),m_numberOfTextureLevels(1),
m_textureUnit(0),m_textureID(0),m_samplerLocation(-1){
}

mgTexture::~mgTexture(){
	if(m_textureID)
		glDeleteTextures(1,&m_textureID);
}

///Set the glsl program and sampler variable name of 
void mgTexture::set_sampler(mgGLSLProgram* glsl, const std::string& samplerName){
	//glsl->use();
	m_samplerLocation=glsl->getUniformLocation(samplerName.c_str());
}

///image データをtextureにセットする
/// @note pixelsピクセルデータのフォーマットについて説明すると、
/// GLuint 型の配列を利用し、各 4 バイトを 1 ピクセルの色に対応させている。
/// 整数の上位バイトから A, B, G, R の順に格納されていることを前提とする。
/// 例えば値 0xFF000000, 0xFFFF0000, 0xFF00FF00, 0xFF0000FF はそれぞれ
/// 黒、青、緑、赤を意味する。
///Set the image data. set_image() invokes:
///1. glGenTextures(if m_textureID was 0),
///2. glBindTexture
///3. glTexImage2D or glTexStorege2D
///4. glTexSubImage2D
///*******After set_image is invoked, image data is unnecessary since the data 
///is transfered to GPU.
void mgTexture::set_image(
	GLsizei width, ///< ビットマップデータの横幅。ピクセル幅と考えて差し支えない。
	GLsizei height,///< ビットマップデータの縦幅。
	const GLuint* pixels,///< RGBA バイト列。各整数の上位バイトが A 値。
	bool mutableTexture
	,bool isPointSprite
	,GLint wrap///GL_TEXTURE_WRAP_S&GL_TEXTURE_WRAP_Tに指定するparameterを指定
		///GL_REPEAT, GL_CLAMP_TO_EDGE, GL_CLAMP_TO_BORDER, GL_MIRRORED_REPEAT,
		///or GL_MIRROR_CLAMP_TO_EDGEGL_REPEAT.
	,GLint magminFilter///GL_TEXTURE_MAG_FILTER&GL_TEXTURE_MIN_FILTERに
		///指定するparameterを指定する:GL_NEAREST,GL_LINEAR,GL_NEAREST_MIPMAP_NEAREST,
		///GL_LINEAR_MIPMAP_NEAREST,GL_NEAREST_MIPMAP_LINEAR,GL_LINEAR_MIPMAP_LINEAR
		///いずれか。GL_TEXTURE_MAG_FILTER&GL_TEXTURE_MIN_FILTER双方に同じ値がセットされる。
){
	if(!m_textureID)
		glGenTextures(1,&m_textureID);
	mgGLSLProgram* glsl=mgGLSLProgram::getCurrentGLSLProgram();
	GLint major, minor;
	glsl->getOpenGLVerion(major,minor);
	if(minor<=1)
		mutableTexture=true;

	glBindTexture(m_target, m_textureID);	assert(glIsTexture(m_textureID));

	if(isPointSprite){
		glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
		glTexEnvi(GL_POINT_SPRITE, GL_COORD_REPLACE, GL_TRUE);
	}
 
	GLint textureLevel=0;
	if(mutableTexture)
		glTexImage2D(m_target,textureLevel,	GL_RGBA8,//internal format
			width,height,
			0,GL_RGBA,GL_UNSIGNED_BYTE,0);
	else
		glTexStorage2D(m_target,m_numberOfTextureLevels,GL_RGBA8,width,height);//not usable in 4.1

	glTexParameteri(m_target,GL_TEXTURE_WRAP_S,wrap);
	glTexParameteri(m_target,GL_TEXTURE_WRAP_T,wrap);
	glTexParameteri(m_target,GL_TEXTURE_MAG_FILTER,magminFilter);
	glTexParameteri(m_target,GL_TEXTURE_MIN_FILTER,magminFilter);
	glPixelStorei(GL_UNPACK_ALIGNMENT, 4);
	glPixelStorei(GL_UNPACK_SWAP_BYTES, GL_TRUE);
	glTexSubImage2D(m_target,textureLevel,
		0,0,width,height,//xoffset, yoffset,  width,hwight
		GL_RGBA,GL_UNSIGNED_BYTE,pixels);
	glBindTexture(m_target, 0);
	m_width=width; m_height=height;
}

///Set the image data. set_image() invokes:
///1. glGenTextures(if m_textureIDwas 0),
///2. glBindTexture
///3. glTexImage2D or  glTexStorege2D
///4. glTexSubImage2D
///*******After set_image is invoked, image data is unnecessary since teh data 
///is transfered to GPU.
void mgTexture::set_image(
	const MGImage& image,
	bool mutableTexture
	,bool isPointSprite
	,GLint wrap///GL_TEXTURE_WRAP_S&GL_TEXTURE_WRAP_Tに指定するparameterを指定
		///GL_REPEAT, GL_CLAMP_TO_EDGE, GL_CLAMP_TO_BORDER, GL_MIRRORED_REPEAT,
		///or GL_MIRROR_CLAMP_TO_EDGEGL_REPEAT.
	,GLint magminFilter///GL_TEXTURE_MAG_FILTER&GL_TEXTURE_MIN_FILTERに
		///指定するparameterを指定する:GL_NEAREST,GL_LINEAR,GL_NEAREST_MIPMAP_NEAREST,
		///GL_LINEAR_MIPMAP_NEAREST,GL_NEAREST_MIPMAP_LINEAR,GL_LINEAR_MIPMAP_LINEAR
		///いずれか。GL_TEXTURE_MAG_FILTER&GL_TEXTURE_MIN_FILTER双方に同じ値がセットされる。
){
	set_image(image.width(),image.height(),(GLuint*)image.image()
		,mutableTexture, isPointSprite,wrap,magminFilter);
}

///Use this texture after set_image and set_sampler are invoked.
///Use() will invoke OpenGL:
///1. glActiveTexture()
///2. glBindTexture()
///4. glSetUniform1i(m_samplerLocation).
///*****glUseProgram() is not invoked, m_samplerLocation must be 
///the current program's.
void mgTexture::use()const{
	int TUnitEnum=GL_TEXTURE0+m_textureUnit;
	glActiveTexture(TUnitEnum);

	assert(m_textureID);
	glBindTexture(m_target,m_textureID);

	mgGLSLProgram* glsl=mgGLSLProgram::getCurrentGLSLProgram();
	glsl->setUniform(mgGLSLProgram::texture2D, m_textureUnit);
}
