/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno             */
/* All rights reserved.                                             */
/********************************************************************/
#ifndef _MGPixel_HH_
#define _MGPixel_HH_

#include "mg/MGCL.h"

/** @addtogroup DisplayHandling
 *  @{
 */

///Define MGPixel Class of (R,G,B,A) pixel data.

///MGPixel defines one pixel data of MGImage.
///That is, one (R,G,B,A) data
class MG_DLL_DECLR MGPixel{

public:

MGPixel(unsigned data=0){m_pixel.uint=data;};

///Conversion constructor from MGColor to MGPixel.
MGPixel(const MGColor& color);

unsigned char getRed()const{return m_pixel.uchar[0];};
unsigned char getGreen()const{return m_pixel.uchar[1];};
unsigned char getBlue()const{return m_pixel.uchar[2];};
unsigned char getAlpha()const{return m_pixel.uchar[3];};

unsigned int getPixel()const{return m_pixel.uint;};
const unsigned char* getPixelAsChar()const{return m_pixel.uchar;};

void setRed(unsigned char R){m_pixel.uchar[0]=R;};
void setGreen(unsigned char G){m_pixel.uchar[1]=G;};
void setBlue(unsigned char B){m_pixel.uchar[2]=B;};
void setAlpha(unsigned char A){m_pixel.uchar[3]=A;};

///Set only RGB of pixel2 without updating Alpha data.
void setRGB(const MGPixel& pixel2);

void setPixel(unsigned int pixel){m_pixel.uint=pixel;};

unsigned char& operator[](unsigned i){return m_pixel.uchar[i];};
const unsigned char& operator[](unsigned i)const{return m_pixel.uchar[i];};
MGPixel& operator&=(unsigned int data){m_pixel.uint &=data;return *this;};
MGPixel& operator|=(unsigned int data){m_pixel.uint |=data;return *this;};

///Test if this pixel is null pixel(that is, RGBA=(0,0,0,0)).
bool isNull()const{return m_pixel.uint==0;};

private:

///Union to treat one pixel as a unit and RGBA elements.
///One pixel consists of 4 bytes, which are (red, green , blue, alpha).
///If uint is used, the image's pixel is treated as a image,
///and if uchar[i] is used each element of the pixel can be handled.
union pixelData{
	unsigned int uint;
	unsigned char uchar[4];//[0-2]=(R,G,B), [3]=A.
};

pixelData m_pixel;//One pixel data.

};

/** @} */ // end of DisplayHandling group
#endif
