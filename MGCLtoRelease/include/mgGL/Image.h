/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#pragma once

#include "mg/MGCL.h"
#include "mgGL/Pixel.h"
#include <gdiplus.h>

//class MGOfstream;
//class MGIfstream;

/** @addtogroup DisplayHandling
 *  @{
 */

///MGImage defines bit map image data.

///MGImage defines the attributes of Image.
///MGImage's behavior is just like std::auto_ptr<>.
///That is , newed area of m_image is transfered if copied, of assined.
class MG_DLL_DECLR MGImage{

public:

////////Special member functions/////////
virtual ~MGImage();
MGImage(const MGImage& image2);///Copy constructor.
MGImage(MGImage&& image2)noexcept;//Move constructor.
MGImage& operator=(const MGImage& image2);//Assignment.
MGImage& operator=(MGImage&&)noexcept;

MGImage():m_width(0), m_height(0),m_image(0){;};
MGImage(int wdth, int hght);///Garbage image data constructor.

///Extract a part of image2.
MGImage(
	const MGImage& image2,///<The target image.
	int x,///<left bottom address of image2, x.
	int y,	///< y.
	int width,///< Width.
	int height///< Height.
);

///Conversion constructor from Gdiplus::Bitmap.
MGImage(Gdiplus::Bitmap& bitmap);

///Extract a part of bitmap.
MGImage(
	Gdiplus::Bitmap& bitmap,///<The target bitmap image.
	int x,///<left bottom address of image2, x.
	int y,	///< y.
	int width,///< Width.
	int height///< Height.
);

///Extract a part of bitmap with transparent
MGImage(
	Gdiplus::Bitmap& bitmap,///<The target bitmap image.
	int x,///<left bottom address of image2, x.
	int y,	///< y.
	int width,///< Width.
	int height,///< Height.
	double alpha///< 0. <= alpha <= 1.
);

///imageのcontrast変換 
MGImage(
	Gdiplus::Bitmap& bitmap,
	const Gdiplus::BrightnessContrast& bc,
	double alpha
);

///Generate a cloned MGImage.
///Returned is newed one, must be deleted.
///clone() does not affect this image data.
MGImage* clone()const;

int width()const{return m_width;};
int height()const{return m_height;};
MGPixel* image(){return m_image;};
const MGPixel* image()const{return m_image;};
MGPixel& operator()(int i, int j);
const MGPixel& operator()(int i, int j)const;

///Test if pixel at (i,j) has zero_alpha value.
///Returns true if alpha value is zero.
bool is_zero_alpha(int i, int j)const;

///Fill all the pixel of this with the input color pdata.
void fill_color(
	const MGPixel& pdata
);

///Fill all the pixels of the range (j, i1,i2) with the input color pdata.
//The color of Pixel(i,j) for i=i1,...i2 is set to pdata.
void fill_color(
	const MGPixel& pdata,
	int j,
	int i1,
	int i2
);

///Fill all the pixels of ranges with the input color pdata.
//In ranges, range(j,i1,i2) are stored.
//Let m=ranges.size(), then m=3n(always a mutiple of 3) where n is the number of ranges.
//Let (j,i1,i2)=(ranges[3*k],ranges[3*k+1], ranges[3*k+2]) for k=0,...,n-1,
//then PixelData(i,j) for i=i1,...i2 is one range for the height j of this mesh.
void fill_color(
	const MGPixel& pdata,
	const std::vector<int>& ranges// Ranges(j,i1,i2) are input.
);

///全てのpixelについて、alphaはそのままでinputされた色に変更する。
///inputされた色のalphaは無視される。
void fill_color_NoChangeAlpha(
	const MGPixel& pdata
);
///(j, i1,i2) で指定された全てのピクセルについて、
///alphaはそのままでinputされた色に変更する。
void fill_color_NoChangeAlpha(
	const MGPixel& pdata,
	int j,
	int i1,
	int i2
);

///Copy all the pixels of image2 into this.
///This image's width and height are not changed. And part of this and image2
///are copied into this.
void copy_color(
	const MGImage& image2,//Source image data.
	int x, int y,	///<left bottom address of image2.
	int wdth, int hght
);

///Test if any one of the four pixels of (i,j) to (i+1,j+1) is non zero or not.
///i must be < width()-1, and j must be < hieght()-1.
///If any one of them is nonzero, return true.
bool includeNonZeroAlpha(int i, int j)const;

///Copy all the pixels of ranges in image2 into this.
///In ranges, range(j,i1,i2) are stored.
///Let m=ranges.size(), then m=3n(always a mutiple of 3) where n is the number of ranges.
///Let (j,i1,i2)=(ranges[3*k],ranges[3*k+1], ranges[3*k+2]) for k=0,...,n-1,
///then PixelData(i,j) for i=i1,...i2 is one range for the height j of this mesh.
void copy_color(
	const MGImage& image2,///<Source image data.
	const std::vector<int>& ranges///< Ranges(j,i1,i2) are input.
		///<ranges indicate the places of both this and image2.
);

///resize the image size to (width,height).
///Scaling the whole image to the size (width, height).
void resize(
	int width,
	int height
);

///resize the image size to (width,height) filling the color to
///the extra part for the size(width,height).
///resize_with_fill_color() does not perform scaling to the image.
///width and height can be less than the original length. In this case,
///image trimming will be done.
void resize_with_fill_color(
	int width,
	int height,
	const MGPixel& pdata
);

///Add border color(0,0,0,0)=(alfa=0) of pixel size2 for each perimeter.
void resize_and_add_zero_border(int nwidth2, int nheight2);

private:
	int m_width, m_height;		///<number of pixels about the width and the height.
	MGPixel* m_image;///<array of pixels of size m_width*m_height;

///Extract a part of bitmap into this, from(x,y) to (x+width, y+height).
void extract(
	Gdiplus::Bitmap& bitmap,
	int x,  int y,	///<left bottom address of bitmap.
	int width, int height,
	double alpha=-1.0
);

};

///Compute 2's power of width and height.
void MGImageCompute_2spower(
	int width,///< Width
	int height,///< Height
	int& width2,///<The smallest 2's power of width will be output.
	int& height2///<The smallest 2's power of height will be output.
);

/** @} */ // end of DisplayHandling group
