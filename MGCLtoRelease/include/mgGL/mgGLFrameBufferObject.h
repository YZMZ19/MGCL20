/********************************************************************/
#ifndef _MGGLFBO_HH_
#define _MGGLFBO_HH_


#include <GL/glew.h>
#include <GL/wglew.h>
#include <mg/MGCL.h>

/** @addtogroup GLAttrib
 *  @{
 */

/// �I�t�X�N���[�������_�����O�����Ŏg�p����t���[���o�b�t�@�I�u�W�F�N�g�N���X.

/// OpenGL��FBO�ɑ������܂��B
///  ���p�菇�͎��̂Ƃ���F
/// * (1) CreateBuffer�ɂ�菉��������
/// * (2) beginRenderTexture()���M
/// * (3) �`�揈��(MGOpenGLView::drawScene()�Ȃǁj
/// * (4) CopyRenderTexture()
/// * (5) endRenderTexture()���M
class MG_DLL_DECLR mgGLFramebufferObject {

//////////////////////////////////////////////////////////////////////////
// �p�u���b�N ���\�b�h
public:
	mgGLFramebufferObject();
	~mgGLFramebufferObject(){;};
	
	/// ����Ԃ��܂��B
	int getWidth(){	return m_Width;	};

	/// ������Ԃ��܂��B
	int getHeight(){return m_Height;};

	/// �e�N�X�`�����ʎq��Ԃ��܂��B
	GLuint getTexName(){return m_textureId;};


	/**
	 * �w�肳�ꂽ���ƍ����Ńt���[���o�b�t�@�I�u�W�F�N�g (FBO) ���\�����܂��B<p>
	 * ���Ƀt���[���o�b�t�@�I�u�W�F�N�g (FBO) ���\������Ă���ꍇ�́A
	 * ���݂̃t���[���o�b�t�@�I�u�W�F�N�g (FBO) ���폜���ĐV�����t���[���o�b�t�@�I�u�W�F�N�g (FBO) ���\�����܂��B
	 * 
	 * @param width ��
	 * @param height ����
	 * @param format �����t�H�[�}�b�g
	 * @throws RuntimeException �t���[���o�b�t�@�̍\���Ɏ��s�����ꍇ�B
	 */
	///Buffer �� creat. format��CopyRenderTexture�Ŏ��o���Ƃ���format���w��
	///����fromat��GL_RGBA8�Ƃ����
	bool CreateBuffer(int width, int height, GLenum format = GL_RGB);

	/**
	 * �N���[���A�b�v���s���܂��B
	 */
	void ReleaseBuffer();

	/**
	 * ���̃t���[���o�b�t�@�I�u�W�F�N�g���o�C���h���ėL���ɂ��܂��B
	 */
	void beginRenderTexture();

	/**
	 * ���̃t���[���o�b�t�@�I�u�W�F�N�g���A���o�C���h���Ė����ɂ��܂��B
	 */
	void endRenderTexture();

	///FBO.beginRenderTexture(), draw(), ..., ....
	///�ō쐬����screen image�f�[�^��CreateBuffer�Ŏw�肵���`����pixels�ɓǂ݂���
	void CopyRenderTexture(
		GLvoid* pixels	///format�ɏ]���ė̈�m�ی�̗̈�pointer����͂���
	);

	///FBO.beginRenderTexture(), draw(), ..., ....
	///�ō쐬����screen image�f�[�^��Gdiplus::Bitmap�`���œǂ݂���
	///function�̖߂�l��new�Ŋm�ۂ��ꂽ�̈�ł���A���p�҂�delete����K�v������
	Gdiplus::Bitmap* readViewAsBitmap();
	
private:
	
	int m_Width;/// �����_�����O�C���[�W�̕�
	int m_Height;/// �����_�����O�C���[�W�̍���
	GLuint m_pixFormat;///CopyRenderTexture�Ŏ��o��image�̌`�����w�肷��
	GLuint m_framebufferId;///�t���[���o�b�t�@���ʎq��ێ����܂��B	
	GLuint m_renderbufferId;///�����_�[�o�b�t�@���ʎq��ێ����܂��B	
	GLuint m_textureId;///�e�N�X�`�����ʎq��ێ����܂��B

	///���̃t���[���o�b�t�@�̃o�b�N�A�b�v
	///beginRenderTexture�ŋL������AendRenderTexture�Ŗ߂����
	GLint m_oldBufferId;

	/// ���̃h���[�o�b�t�@
	///beginRenderTexture�ŋL������AendRenderTexture�Ŗ߂����
	GLint m_oldDrawBuf;

};

/** @} */ // end of GLAttrib group
#endif //_MGGLFBO_HH_