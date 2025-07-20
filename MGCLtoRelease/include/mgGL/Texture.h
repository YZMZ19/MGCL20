/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/

/**
 * @file Texture.h
 * @brief �N���X MGGridCursor �̐錾
 */
#ifndef _mgTexture_HH
#define _mgTexture_HH

#include "mg/MGCL.h"

class MGImage;
class mgGLSLProgram;
/** @addtogroup GLAttrib
 *  @{
 */

/// @class mgTexture Texture.h "mgGL/Texture.h"

///Texture definition class.

///Before use of mgTexture, set_image must be invoked, and
///set_sampler is optionally invoked. If set_sampler is not invoked,
///uniform of spampler2D  "texture1" is assumed.
///Since mgTexture does not invoke glUseProgram(mgGLSLProgram's use() function)
///the current program's uniform variable texture1(or the variable set by set_sampler)
///must be spamler2D.
///When set_textureUnit is not invoked, unit number 0 is assumed.
class MG_DLL_DECLR mgTexture{

public:

	/// �f�t�H���g�R���X�g���N�^�[
	///Before use of mgTexture, set_image and set_sampler must be invoked.
	mgTexture(
		GLenum target=GL_TEXTURE_2D///Currently only GL_TEXTURE_2D is allowed.
	);

	~mgTexture();

	///Set the glsl program and sampler variable name of 
	///*****glUseProgram() is not invoked, glsl must be the current program.
	void set_sampler(mgGLSLProgram* glsl, const std::string& samplerName);

	const std::string& getDefaultSamplerVariable()const;


	///image �f�[�^��texture�ɃZ�b�g����
	/// @note pixels�s�N�Z���f�[�^�̃t�H�[�}�b�g�ɂ��Đ�������ƁA
	/// GLuint �^�̔z��𗘗p���A�e 4 �o�C�g�� 1 �s�N�Z���̐F�ɑΉ������Ă���B
	/// �����̏�ʃo�C�g���� A, B, G, R �̏��Ɋi�[����Ă��邱�Ƃ�O��Ƃ���B
	/// �Ⴆ�Βl 0xFF000000, 0xFFFF0000, 0xFF00FF00, 0xFF0000FF �͂��ꂼ��
	/// ���A�A�΁A�Ԃ��Ӗ�����B
	///Set the image data. set_image() invokes:
	///1. glGenTextures(if m_textureID was 0)
	///2. glBindTexture
	///3. glTexImage2D or glTexStorege2D
	///4. glTexSubImage2D
	///*******After set_image is invoked, image data is unnecessary since the data 
	///is transfered to GPU.
	void set_image(
		GLsizei width ///< �r�b�g�}�b�v�f�[�^�̉����B�s�N�Z����
		,GLsizei height///< �r�b�g�}�b�v�f�[�^�̏c���B
		,const GLuint* pixels///< RGBA �o�C�g��B�e�����̏�ʃo�C�g�� A �l�B
		,bool mutableTexture=true///<true if mutable, false, if immutable.
		,bool isPointSprite = false///<true if point sprite.
		,GLint wrap=GL_REPEAT///<GL_TEXTURE_WRAP_S&GL_TEXTURE_WRAP_T�Ɏw�肷��parameter���w��
			///<GL_REPEAT, GL_CLAMP_TO_EDGE, GL_CLAMP_TO_BORDER, GL_MIRRORED_REPEAT,
            ///<or GL_MIRROR_CLAMP_TO_EDGEGL_REPEAT.
		,GLint magminFilter=GL_LINEAR///<GL_TEXTURE_MAG_FILTER&GL_TEXTURE_MIN_FILTER��
			///<�w�肷��parameter���w�肷��:GL_NEAREST,GL_LINEAR,GL_NEAREST_MIPMAP_NEAREST,
			///<GL_LINEAR_MIPMAP_NEAREST,GL_NEAREST_MIPMAP_LINEAR,GL_LINEAR_MIPMAP_LINEAR
			///<�����ꂩ�BGL_TEXTURE_MAG_FILTER&GL_TEXTURE_MIN_FILTER�o���ɓ����l���Z�b�g�����B
	);

	///Set the image data from other MGImage. This type of set_image uses image.width(), .height(),
	///and .image(), and invoke above set_image.
	///set_image() invokes:
	///1. glGenTextures(if m_textureIDwas 0),
	///2. glBindTexture
	///3. glTexImage2D or  glTexStorege2D
	///4. glTexSubImage2D
	///*******After set_image is invoked, image data is unnecessary since teh data 
	///is transfered to GPU.
	void set_image(
		const MGImage& image///<Original image.
		,bool mutableTexture=true///<true if mutable, false, if immutable.
		,bool isPointSprite=false///<true if point sprite.
		,GLint wrap=GL_REPEAT///GL_TEXTURE_WRAP_S&GL_TEXTURE_WRAP_T�Ɏw�肷��parameter���w��
			///GL_REPEAT, GL_CLAMP_TO_EDGE, GL_CLAMP_TO_BORDER, GL_MIRRORED_REPEAT,
            ///or GL_MIRROR_CLAMP_TO_EDGEGL_REPEAT.
		,GLint magminFilter=GL_LINEAR///GL_TEXTURE_MAG_FILTER&GL_TEXTURE_MIN_FILTER��
			///�w�肷��parameter���w�肷��:GL_NEAREST,GL_LINEAR,GL_NEAREST_MIPMAP_NEAREST,
			///GL_LINEAR_MIPMAP_NEAREST,GL_NEAREST_MIPMAP_LINEAR,GL_LINEAR_MIPMAP_LINEAR
			///�����ꂩ�BGL_TEXTURE_MAG_FILTER&GL_TEXTURE_MIN_FILTER�o���ɓ����l���Z�b�g�����B
	);

	///When set_textureUnit is not invoked, m_textureUnit=0 is assumed.
	void set_textureUnit(int textureUnit){m_textureUnit=textureUnit;};

	///get texture ID
	GLuint getTextureID()const{return m_textureID;};

	///Use this texture after set_image and set_sampler are invoked.
	///Use() will invoke OpenGL:
	///1. glActiveTexture()
	///2. glBindTexture()
	///4. glSetUniform1i(m_samplerLocation).
	///*****glUseProgram() is not invoked, m_samplerLocation must be 
	///the current program's.
	void use()const;

	GLsizei width()const{return m_width;};
	GLsizei height()const{return m_height;};

private:
	GLsizei m_numberOfTextureLevels;//Currently this is set to 1.
	GLsizei m_width, m_height;///Width and height of the texture.
	GLenum m_target;///Specifies texture target.
		///Must be GL_TEXTURE_XX(GLTEXTURE_1D, GL_TEXTURE_2D, GL_TEXTURE_3D, or etc.)

	int m_textureUnit;///Texture unit numger(Input of setUnifrom).
		///Input of glActiveTexture is GL_TEXTURE0+m_textureUnit(note that this is int).
	
	GLuint m_textureID;///Texture object id generated by glGenTextures.

	mutable int m_samplerLocation;///m_glsl's uniform sampler(2D?) location..
		///mgTexture will invoke glSetUniform1i(m_samplerLocation);
};///_mgTexture_HH

/** @} */ // end of GLAttrib group
#endif