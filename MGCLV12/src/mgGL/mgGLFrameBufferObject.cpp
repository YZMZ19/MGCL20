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

	// フレームバッファの生成とバインド
	// 現在のフレームバッファを退避
	glGetIntegerv(GL_FRAMEBUFFER_BINDING, &m_oldBufferId);
	glGenFramebuffers(1, &m_framebufferId);
	glBindFramebuffer(GL_FRAMEBUFFER, m_framebufferId);

	// テクスチャの生成とColorAttachment0への関連付け
	glGenTextures(1, &m_textureId);
	glBindTexture(GL_TEXTURE_2D, m_textureId);
	glTexImage2D(GL_TEXTURE_2D,0,GL_RGBA8,width,height,0,m_pixFormat,GL_UNSIGNED_BYTE,NULL);
	glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, m_textureId, 0);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

	// レンダーバッファの生成とDepthAttachmentへの関連付け
	glGenRenderbuffers(1, &m_renderbufferId);
	glBindRenderbuffer(GL_RENDERBUFFER, m_renderbufferId);
	glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, m_renderbufferId);
	glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT24, width, height);

	// フレームバッファの状態をチェック
	if(glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE) {
		ReleaseBuffer();
		return false;
	}

	// フレームバッファを元に戻す
	glBindFramebuffer(GL_FRAMEBUFFER, m_oldBufferId);
	m_oldBufferId = 0;

	return true;
}

/**
* クリーンアップを行います。
*/
void mgGLFramebufferObject::ReleaseBuffer() {
	// 古いフレームバッファがあればバインドする
	if(m_oldBufferId) {
		glBindFramebuffer(GL_FRAMEBUFFER, m_oldBufferId);
		m_oldBufferId = 0;
	}

	// レンダーバッファの開放
	if(m_renderbufferId) {
		glDeleteRenderbuffers(1, &m_renderbufferId);
		m_renderbufferId = 0;
	}

	// テクスチャの開放
	if(m_textureId) {
		glDeleteTextures(1, &m_textureId);
		m_textureId = 0;
	}

	// フレームバッファの開放  
	if(m_framebufferId) {
		glDeleteFramebuffers(1, &m_framebufferId);
		m_framebufferId = 0;
	}
}


// レンダーテクスチャ開始関数
void mgGLFramebufferObject::beginRenderTexture() {

	// 現在のフレームバッファを退避
	glGetIntegerv(GL_FRAMEBUFFER_BINDING, &m_oldBufferId);

	// レンダーテクスチャ用のフレームバッファをバインド
	glBindFramebuffer(GL_FRAMEBUFFER, m_framebufferId);
	glGetIntegerv(GL_DRAW_BUFFER, &m_oldDrawBuf);
	GLenum fboBuffers[]={GL_COLOR_ATTACHMENT0};
	glDrawBuffers(1, fboBuffers);
}

// レンダーテクスチャ終了関数
void mgGLFramebufferObject::endRenderTexture() {
	glDrawBuffer(m_oldDrawBuf);
	// フレームバッファを元に戻す
	glBindFramebuffer(GL_FRAMEBUFFER, m_oldBufferId);
	m_oldBufferId = 0;
}

// レンダーテクスチャコピー関数
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
	Gdiplus::Bitmap* bitmap = 0;//生成するbitmap

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
