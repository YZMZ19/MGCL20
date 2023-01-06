#include "StdAfx.h"
#include "FTGL/ftgl.h"
#include <tchar.h>
#include <atlstr.h>
#include "mgGL/Appearance.h"
#include "mgGL/VBObyAnchorPt.h"
#include "mgGL/MGStringWriter.h"


MGStringWriter::MGStringWriter(void)
{
	// TODO フォントファイルは同梱すること。
	const char* fontPath = "C:\\Windows\\Fonts\\arial.ttf";
	unsigned int faceSize = 12;
	unsigned int resolution =96; 

	m_pFont = new FTPolygonFont(fontPath);
    m_pFont->FaceSize(faceSize,resolution);
    m_pFont->CharMap(ft_encoding_unicode);
	m_pFont->UseDisplayList(false);
}

///Set font data.
///The default font is "C:\\Windows\\Fonts\\arial.ttf".
///If this is not the case, setFont must be invoked.
void MGStringWriter::setFont(
	const char* fontPath,	///< font file path.
	unsigned int faceSize,	///< the face size in points(1/72 inch).
	unsigned int resolution	///<　the resolution of the target device.
){
	MGStringWriter& writer=*(getInstance());
	delete writer.m_pFont;

	writer.m_pFont = new FTPolygonFont(fontPath);
    writer.m_pFont->FaceSize(faceSize,resolution);
    writer.m_pFont->CharMap(ft_encoding_unicode);
	writer.m_pFont->UseDisplayList(false);
}

MGStringWriter::~MGStringWriter(void)
{
	delete m_pFont;
}
void MGStringWriter::Init(){
	MGStringWriter::getInstance();
}

MGStringWriter* MGStringWriter::getInstance(){
	static MGStringWriter instance;
	return &instance;
}

VBObyAnchorPt* MGStringWriter::Draw(
	const char *str,
	const MGPosition& pos,
	const MGColor* color)
{
	return getInstance()->DrawString(mgGLSL::AnchorPoint, str, pos, color);
}

VBObyAnchorPt*  MGStringWriter::Draw(
	const wchar_t *str,
	const MGPosition& pos,
	const MGColor* color)
{
	return getInstance()->DrawString(mgGLSL::AnchorPoint,str,pos,color);
}

VBObyAnchorPt* MGStringWriter::DrawByScreen(
	const char *str,
	const MGPosition& pos,
	const MGColor* color)
{
	return getInstance()->DrawString(mgGLSL::AnchorPointScreen, str, pos, color);
}

VBObyAnchorPt* MGStringWriter::DrawByScreen(
	const wchar_t *str,
	const MGPosition& pos,
	const MGColor* color)
{
	return getInstance()->DrawString(mgGLSL::AnchorPointScreen, str, pos, color);
}
