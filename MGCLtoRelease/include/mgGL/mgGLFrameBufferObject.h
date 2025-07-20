/********************************************************************/
#ifndef _MGGLFBO_HH_
#define _MGGLFBO_HH_


#include <GL/glew.h>
#include <GL/wglew.h>
#include <mg/MGCL.h>

/** @addtogroup GLAttrib
 *  @{
 */

/// オフスクリーンレンダリング処理で使用するフレームバッファオブジェクトクラス.

/// OpenGLのFBOに相当します。
///  利用手順は次のとおり：
/// * (1) CreateBufferにより初期化処理
/// * (2) beginRenderTexture()発信
/// * (3) 描画処理(MGOpenGLView::drawScene()など）
/// * (4) CopyRenderTexture()
/// * (5) endRenderTexture()発信
class MG_DLL_DECLR mgGLFramebufferObject {

//////////////////////////////////////////////////////////////////////////
// パブリック メソッド
public:
	mgGLFramebufferObject();
	~mgGLFramebufferObject(){;};
	
	/// 幅を返します。
	int getWidth(){	return m_Width;	};

	/// 高さを返します。
	int getHeight(){return m_Height;};

	/// テクスチャ識別子を返します。
	GLuint getTexName(){return m_textureId;};


	/**
	 * 指定された幅と高さでフレームバッファオブジェクト (FBO) を構成します。<p>
	 * 既にフレームバッファオブジェクト (FBO) が構成されている場合は、
	 * 現在のフレームバッファオブジェクト (FBO) を削除して新しいフレームバッファオブジェクト (FBO) を構成します。
	 * 
	 * @param width 幅
	 * @param height 高さ
	 * @param format 内部フォーマット
	 * @throws RuntimeException フレームバッファの構成に失敗した場合。
	 */
	///Buffer を creat. formatはCopyRenderTextureで取り出すときのformatを指定
	///内部fromatはGL_RGBA8とされる
	bool CreateBuffer(int width, int height, GLenum format = GL_RGB);

	/**
	 * クリーンアップを行います。
	 */
	void ReleaseBuffer();

	/**
	 * このフレームバッファオブジェクトをバインドして有効にします。
	 */
	void beginRenderTexture();

	/**
	 * このフレームバッファオブジェクトをアンバインドして無効にします。
	 */
	void endRenderTexture();

	///FBO.beginRenderTexture(), draw(), ..., ....
	///で作成したscreen imageデータをCreateBufferで指定した形式でpixelsに読みだす
	void CopyRenderTexture(
		GLvoid* pixels	///formatに従って領域確保後の領域pointerを入力する
	);

	///FBO.beginRenderTexture(), draw(), ..., ....
	///で作成したscreen imageデータをGdiplus::Bitmap形式で読みだす
	///functionの戻り値はnewで確保された領域であり、利用者はdeleteする必要がある
	Gdiplus::Bitmap* readViewAsBitmap();
	
private:
	
	int m_Width;/// レンダリングイメージの幅
	int m_Height;/// レンダリングイメージの高さ
	GLuint m_pixFormat;///CopyRenderTextureで取り出すimageの形式を指定する
	GLuint m_framebufferId;///フレームバッファ識別子を保持します。	
	GLuint m_renderbufferId;///レンダーバッファ識別子を保持します。	
	GLuint m_textureId;///テクスチャ識別子を保持します。

	///元のフレームバッファのバックアップ
	///beginRenderTextureで記憶され、endRenderTextureで戻される
	GLint m_oldBufferId;

	/// 元のドローバッファ
	///beginRenderTextureで記憶され、endRenderTextureで戻される
	GLint m_oldDrawBuf;

};

/** @} */ // end of GLAttrib group
#endif //_MGGLFBO_HH_