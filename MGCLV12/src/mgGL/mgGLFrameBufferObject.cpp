#include "StdAfx.h"
#include "mgGL/mgGLFrameBufferObject.h"

mgGLFramebufferObject::mgGLFramebufferObject():
m_Width(0),m_Height(0),
	m_framebufferId(0),m_renderbufferId(0),m_textureId(0)
	,m_oldBufferId(0)
	,m_oldDrawBuf(GL_BACK)
	,m_pixFormat(GL_RGB)
{;}

bool mgGLFramebufferObject::CreateBuffer(int width, int height, GLenum format)
{
	ASSERT(width>0 && height>0);

	m_Width = width;
	m_Height = height;
	m_pixFormat = format;

	// �t���[���o�b�t�@�̐����ƃo�C���h
	// ���݂̃t���[���o�b�t�@��ޔ�
	glGetIntegerv(GL_FRAMEBUFFER_BINDING, &m_oldBufferId);
	glGenFramebuffers(1, &m_framebufferId);
	glBindFramebuffer(GL_FRAMEBUFFER, m_framebufferId);

	// �e�N�X�`���̐�����ColorAttachment0�ւ̊֘A�t��
	glGenTextures(1, &m_textureId);
	glBindTexture(GL_TEXTURE_2D, m_textureId);
	glTexImage2D(GL_TEXTURE_2D,0,GL_RGBA8,width,height,0,m_pixFormat,GL_UNSIGNED_BYTE,NULL);
	glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, m_textureId, 0);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

	// �����_�[�o�b�t�@�̐�����DepthAttachment�ւ̊֘A�t��
	glGenRenderbuffers(1, &m_renderbufferId);
	glBindRenderbuffer(GL_RENDERBUFFER, m_renderbufferId);
	glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, m_renderbufferId);
	glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT24, width, height);

	// �t���[���o�b�t�@�̏�Ԃ��`�F�b�N
	if(glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE) {
		ReleaseBuffer();
		return false;
	}

	// �t���[���o�b�t�@�����ɖ߂�
	glBindFramebuffer(GL_FRAMEBUFFER, m_oldBufferId);
	m_oldBufferId = 0;

	return true;
}

/**
* �N���[���A�b�v���s���܂��B
*/
void mgGLFramebufferObject::ReleaseBuffer() {
	// �Â��t���[���o�b�t�@������΃o�C���h����
	if(m_oldBufferId) {
		glBindFramebuffer(GL_FRAMEBUFFER, m_oldBufferId);
		m_oldBufferId = 0;
	}

	// �����_�[�o�b�t�@�̊J��
	if(m_renderbufferId) {
		glDeleteRenderbuffers(1, &m_renderbufferId);
		m_renderbufferId = 0;
	}

	// �e�N�X�`���̊J��
	if(m_textureId) {
		glDeleteTextures(1, &m_textureId);
		m_textureId = 0;
	}

	// �t���[���o�b�t�@�̊J��  
	if(m_framebufferId) {
		glDeleteFramebuffers(1, &m_framebufferId);
		m_framebufferId = 0;
	}
}


// �����_�[�e�N�X�`���J�n�֐�
void mgGLFramebufferObject::beginRenderTexture() {

	// ���݂̃t���[���o�b�t�@��ޔ�
	glGetIntegerv(GL_FRAMEBUFFER_BINDING, &m_oldBufferId);

	// �����_�[�e�N�X�`���p�̃t���[���o�b�t�@���o�C���h
	glBindFramebuffer(GL_FRAMEBUFFER, m_framebufferId);
	glGetIntegerv(GL_DRAW_BUFFER, &m_oldDrawBuf);
	GLenum fboBuffers[]={GL_COLOR_ATTACHMENT0};
	glDrawBuffers(1, fboBuffers);
}

// �����_�[�e�N�X�`���I���֐�
void mgGLFramebufferObject::endRenderTexture() {
	glDrawBuffer(m_oldDrawBuf);
	// �t���[���o�b�t�@�����ɖ߂�
	glBindFramebuffer(GL_FRAMEBUFFER, m_oldBufferId);
	m_oldBufferId = 0;
}

// �����_�[�e�N�X�`���R�s�[�֐�
void mgGLFramebufferObject::CopyRenderTexture(GLvoid* pixels) {
	glFlush();
	GLint oldReadBuf=GL_BACK;
	glGetIntegerv(GL_READ_BUFFER, &oldReadBuf);
	glReadBuffer(GL_COLOR_ATTACHMENT0);
	glReadPixels(0,0, m_Width, m_Height,m_pixFormat, GL_UNSIGNED_BYTE, pixels);
	glReadBuffer(oldReadBuf);
}

Gdiplus::Bitmap* mgGLFramebufferObject::readViewAsBitmap(){
	//// Change the pixel data format
	int width = getWidth();
	int height = getHeight();
	Gdiplus::Bitmap* bitmap = 0;//��������bitmap

	try{
		bitmap = new Gdiplus::Bitmap(width, height, PixelFormat32bppARGB);
	}catch(std::bad_alloc){
		delete bitmap;
		return 0;
	}

	Gdiplus::BitmapData bitmapData;
	Gdiplus::Rect rect(0, 0, width, height);
	Gdiplus::Status sts=bitmap->LockBits(&rect,
		Gdiplus::ImageLockModeRead|Gdiplus::ImageLockModeWrite,PixelFormat32bppARGB,&bitmapData);

	if(sts==Gdiplus::Ok){
		GLvoid* bgra = (GLvoid*)bitmapData.Scan0;
		glFlush();
		
		GLint oldReadBuf=GL_BACK;
		glGetIntegerv(GL_READ_BUFFER, &oldReadBuf);
		glReadBuffer(GL_COLOR_ATTACHMENT0);
		glReadPixels(0,0, width, height, GL_BGRA_EXT, GL_UNSIGNED_BYTE, bgra);
		glReadBuffer(oldReadBuf);
		glFinish();
		bitmap->UnlockBits(&bitmapData);
		bitmap->RotateFlip( Gdiplus::RotateNoneFlipY );

		return bitmap;
	}else{
		delete bitmap;
		return 0;
	}
}
