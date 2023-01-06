/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/

#pragma once

#include "mgGL/GLAttrib.h"
#include "mgGL/Pixel.h"

class MGOfstream;
class MGIfstream;
class MGColors;
class mgVBO;

//
//Define MGColor Class.
/** @file */

/** @addtogroup GLAttrib
 *  @{
 */

///MGColor defines the OpenGL color (R,G,B,A).
class MG_DLL_DECLR MGColor:public MGGLAttrib{

public:

/// Color enumuration.
/// This is completely compatible to GDI+ color enumuration.

///The color id from 1 to 8(black to white) is the id of IGES.
enum ColorID{
	dummy                = 0,
    Black                = 1,
    Red                  = 2,
    Green                = 3,
    Blue                 = 4,
    Yellow               = 5,
    Magenta              = 6,
    Cyan                 = 7,
    White                = 8,
    AliceBlue            ,
    AntiqueWhite         ,
    Aqua                 ,
    Aquamarine           ,
    Azure                ,
    Beige                ,
    Bisque               ,
    BlanchedAlmond       ,
    BlueViolet           ,
    Brown                ,
    BurlyWood            ,
    CadetBlue            ,
    Chartreuse           ,
    Chocolate            ,
    Coral                ,
    CornflowerBlue       ,
    Cornsilk             ,
    Crimson              ,
    DarkBlue             ,
    DarkCyan             ,
    DarkGoldenrod        ,
    DarkGray             ,
    DarkGreen            ,
    DarkKhaki            ,
    DarkMagenta          ,
    DarkOliveGreen       ,
    DarkOrange           ,
    DarkOrchid           ,
    DarkRed              ,
    DarkSalmon           ,
    DarkSeaGreen         ,
    DarkSlateBlue        ,
    DarkSlateGray        ,
    DarkTurquoise        ,
    DarkViolet           ,
    DeepPink             ,
    DeepSkyBlue          ,
    DimGray              ,
    DodgerBlue           ,
    Firebrick            ,
    FloralWhite          ,
    ForestGreen          ,
    Fuchsia              ,
    Gainsboro            ,
    GhostWhite           ,
    Gold                 ,
    Goldenrod            ,
    Gray                 ,
    GreenYellow          ,
    Honeydew             ,
    HotPink              ,
    IndianRed            ,
    Indigo               ,
    Ivory                ,
    Khaki                ,
    Lavender             ,
    LavenderBlush        ,
    LawnGreen            ,
    LemonChiffon         ,
    LightBlue            ,
    LightCoral           ,
    LightCyan            ,
    LightGoldenrodYellow ,
    LightGray            ,
    LightGreen           ,
    LightPink            ,
    LightSalmon          ,
    LightSeaGreen        ,
    LightSkyBlue         ,
    LightSlateGray       ,
    LightSteelBlue       ,
    LightYellow          ,
    Lime                 ,
    LimeGreen            ,
    Linen                ,
    Maroon               ,
    MediumAquamarine     ,
    MediumBlue           ,
    MediumOrchid         ,
    MediumPurple         ,
    MediumSeaGreen       ,
    MediumSlateBlue      ,
    MediumSpringGreen    ,
    MediumTurquoise      ,
    MediumVioletRed      ,
    MidnightBlue         ,
    MintCream            ,
    MistyRose            ,
    Moccasin             ,
    NavajoWhite          ,
    Navy                 ,
    OldLace              ,
    Olive                ,
    OliveDrab            ,
    Orange               ,
    OrangeRed            ,
    Orchid               ,
    PaleGoldenrod        ,
    PaleGreen            ,
    PaleTurquoise        ,
    PaleVioletRed        ,
    PapayaWhip           ,
    PeachPuff            ,
    Peru                 ,
    Pink                 ,
    Plum                 ,
    PowderBlue           ,
    Purple               ,
    RosyBrown            ,
    RoyalBlue            ,
    SaddleBrown          ,
    Salmon               ,
    SandyBrown           ,
    SeaGreen             ,
    SeaShell             ,
    Sienna               ,
    Silver               ,
    SkyBlue              ,
    SlateBlue            ,
    SlateGray            ,
    Snow                 ,
    SpringGreen          ,
    SteelBlue            ,
    Tan                  ,
    Teal                 ,
    Thistle              ,
    Tomato               ,
    Transparent          ,
    Turquoise            ,
    Violet               ,
    Wheat                ,
    WhiteSmoke           ,
    YellowGreen          ,
	endID
};

MGColor():MGGLAttrib(static_cast<int>(mgGLMode::UNDEFINED)){;};

/// Construct a color from float values.
MGColor(float red, float green, float blue, float alpha=1.);
MGColor(const MGPixel& pixel);

/// Construct a color from ARGB value. In argb, each 8bits of (a,r,g,b) is the value of
///the range 0 to 255 .
explicit MGColor(unsigned int argb);

///Assignment
MGColor& operator=(const MGGel& gel2);
MGColor& operator=(const MGColor& gel2);

///Scaling the color values by the factor scale.
///All of the elements except trancparency element will be multiplied
///by the scale and be clamped between 0. and 1.
MGColor& operator*=(float scale);

///Add a color values to RGB data..
///All of the elements except trancparency element will be added by value
///and be clamped between 0. and 1.
MGColor& operator+=(float value);

///comparison
bool operator<(const MGColor& gel2)const;
bool operator<(const MGGel& gel2)const;

////////Member Function////////

/// helper function
void argb_to_float(unsigned int argb, float out[4]) const;

///Generate a newed clone object.
MGColor* clone()const;

///Invoke appropriate OpenGL fucntion to the drawing environment.
void exec()const;

///vboに対して色属性をセットする
void exec(mgVBO& vbo)const;

float* color(){return m_color;};
const float* color()const{return m_color;};

///get the color id of GDI+ color enumeration.
///If not found, 0 will be returned.
///maxID is the max id of ColorID enumeration to search.
///Searching is done from id 0 to maxID in the table.
///NOTE if not founed, -1 is returned, not 0.
int get_ColorID(int maxID)const;

///get the color id of GDI+ color enumeration as the return value.
///If not found, -1 be returned, and the the color value of GDI+ is set in ColorValue.
///NOTE if not founed, -1 is returned, not 0.
int getColorIdOrValue(
	int maxID,//Maximum enumeration id of GDI+ to search in the table.
			//Searching is done from id 0 to maxID in the table.
	unsigned long& ColorValue
)const;

///Get the color instance reference by the color id.
static const MGColor& get_instance(ColorID id);

///Test if this is highlight attrib or not.
bool is_highlight_attrib()const{return true;};

///Get unsigned integer valur of a color id.
///Function's return value is (A, R, G, B) data of GDI+.
static int get_ARGBinstance(MGColor::ColorID id);

///get unsigned integer value of this color
///return 0xAARRGGBB
unsigned int get_colorAsUInt()const;

void set_color(const float color[4]);
void set_color(float red, float green, float blue, float alpha=1.);
void get_color(float color[4])const;
void get_color(float& red, float& green, float& blue, float& alpha)const;

///draw GLAttribute process.
void drawAttrib(
	mgVBO& vbo,///<Target graphic object.
	bool no_color=false	///<if true, color attribute will be neglected.
)const;

///render GLAttribute process.
void render(mgVBO& vbo)const{exec(vbo);};

///Turn on the appropriate mask bit for this attribute. See glPushAttrib().
void set_draw_attrib_mask(unsigned int& mask)const{set_Amask(mask,CURRENT_BIT);};

///Turn off the appropriate mask bit for this attribute. See glPushAttrib().
void reset_draw_attrib_mask(unsigned int& mask)const{reset_Amask(mask,CURRENT_BIT);};

///Set this color be enabled.
///set_disabled() is provided by MGGLAttrib.
///This can be used display/undisplay.
void set_enabled(){m_flag=static_cast<int>(mgGLMode::ENABLED);};

///Turn on the appropriate mask bit for this attribute. See glPushAttrib().
void set_render_attrib_mask(unsigned int& mask)const{set_Amask(mask,CURRENT_BIT);};

///Turn off the appropriate mask bit for this attribute. See glPushAttrib().
void reset_render_attrib_mask(unsigned int& mask)const{reset_Amask(mask,CURRENT_BIT);};

/// Return This object's typeID
long identify_type() const{return MGCOLOR_TID;};

///Get the name of the class.
std::string whoami()const{return "Color";};

///Read all member data.
void ReadMembers(MGIfstream& buf);
///Write all member data
void WriteMembers(MGOfstream& buf)const;

/// Output function.
std::ostream& toString(std::ostream&) const;

///Output to IGES stream file(Color=PD314).
///Function's return value is the directory entry id created.
int out_to_IGES(
	MGIgesOfstream& igesfile,
	int SubordinateEntitySwitch=0
)const;

private:

    float m_color[4] = { 0.,0.,0.,1. };///<color data.
	static MGColors m_colors;///<color instance;

	friend class MGColors;

};

///@cond

///MGColors defines standard MGColor array that can be accessed through ColorID.
class MG_DLL_DECLR MGColors{
public:
	const MGColor& color(int i){return *(m_colors[i]);};
	~MGColors();
private:
	MGColor** m_colors;
	MGColors();
	friend class MGColor;
};

///@endcond

#ifndef _CONSOLE

namespace{

///Convert Windows GDI color(0x00bbggrr) to MGColor.
MGColor FromCOLORREF(COLORREF color){
	return MGColor(
		GetRValue(color)/255.f,
		GetGValue(color)/255.f,
		GetBValue(color)/255.f);
}

}
	
///Convert dwColor(0xaarrggbb) to COLORREF(0x00bbggrr).
inline COLORREF ARGBtoCOLORREF(DWORD dwColor){
	BYTE r = (BYTE)((dwColor & 0x00FF0000) >> 16);
	BYTE g = (BYTE)((dwColor & 0x0000FF00) >> 8);
	BYTE b = (BYTE)(dwColor & 0x000000FF);
	return RGB(r, g, b);
}

///Convert COLORREF(0x00bbggrr) to dwColor(0x00rrggbb).
inline DWORD COLORREFtoARGB(COLORREF cl){
	BYTE r = GetRValue(cl);
	BYTE g = GetGValue(cl);
	BYTE b = GetBValue(cl);
	return (r << 16) | (g << 8) | b;
}

///Convert MGColor to Windows GDI color(0x00bbggrr).
inline COLORREF MGCOLORtoCOLORREF(const MGColor& c){
	const float* color=c.color();
	return RGB(color[0]*255, color[1]*255, color[2]*255);
}

#endif //_CONSOLE

/** @} */ // end of GLAttrib group
