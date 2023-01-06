/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/

#include "StdAfx.h"

#include <Gdiplusgraphics.h >
#include <GdiPlusEffects.h>
#include <Gdiplusheaders.h>

#include "mgGL/Color.h"
#include "mgGL/Image.h"

using namespace Gdiplus;

//Set only RGB of pixel2 without updating Alpha data.
void MGPixel::setRGB(const MGPixel& pixel2){
	unsigned char A=getAlpha();
	*this=pixel2;
	setAlpha(A);
}

///Conversion conrtructor from MGColor to MGPixel.
MGPixel::MGPixel(const MGColor& color){
	float rgba[4];
	color.get_color(rgba[0],rgba[1],rgba[2],rgba[3]);
	for(int i=0; i<4; i++){
		if(rgba[i]>1.f)
			rgba[i]=1.f;
		else if(rgba[i]<0.f)
			rgba[i]=0.f;
	}
	setRed(unsigned char(rgba[0]*255.));
	setGreen(unsigned char(rgba[1]*255.));
	setBlue(unsigned char(rgba[2]*255.));
	setAlpha(unsigned char(rgba[3]*255.));
}


//Garbage image data constructor.
MGImage::MGImage(int wdth, int hght)
:m_width(wdth),m_height(hght),m_image(nullptr){
	if(wdth>0 && hght>0)
		m_image=new MGPixel[wdth*hght];
}

MGImage::MGImage(MGImage&& image2)noexcept
	:m_width(image2.width()),m_height(image2.height()),m_image(image2.m_image){
	image2.m_image=nullptr;
}

///Extract a part of image2.
MGImage::MGImage(
	const MGImage& image2,
	int x,  int y,	///<left bottom address of image2.
	int width, int hight
):m_width(width),m_height(hight),m_image(new MGPixel[width*hight]){
	assert((x + width) <= image2.width() && (y + hight) <= image2.height());
	copy_color(image2, x, y, width, hight);
}
MGImage::MGImage(const MGImage& image2)
	:MGImage(image2,0,0,image2.width(), image2.height()){
}

///Conversion constructor from Gdiplus::Bitmap.
MGImage::MGImage(
	Gdiplus::Bitmap& bitmap
):m_width(bitmap.GetWidth()),m_height(bitmap.GetHeight()),
m_image(new MGPixel[m_width*m_height]){
	extract(bitmap,0,0,m_width,m_height);
}

//imageのcontrast変換 
MGImage::MGImage(
	Gdiplus::Bitmap& bitmap,
	const Gdiplus::BrightnessContrast& bc,
	double alpha
):m_width(bitmap.GetWidth()),m_height(bitmap.GetHeight()),
m_image(new MGPixel[m_width*m_height]){

	// 四角形の左上隅と右下隅の座標を定義
	RECT rect = {0, 0, m_width, m_height};

    // Increase the brightness in a portion of the image.
	Gdiplus::Effect& effect = (Gdiplus::Effect&)bc;
	bitmap.ApplyEffect(&effect, &rect);

	extract(bitmap,0,0,m_width,m_height,alpha);
}

//Extract a part of bitmap.
MGImage::MGImage(
	Gdiplus::Bitmap& bitmap,
	int x,  int y,	//left bottom address of bitmap.
	int width, int height
):m_width(width),m_height(height),
m_image(new MGPixel[width*height]){
	extract(bitmap,x,y,width,height);
}

//Extract a part of bitmap with transparent
MGImage::MGImage(
	Gdiplus::Bitmap& bitmap,
	int x,  int y,	//left bottom address of bitmap.
	int width, int height,
	double alpha
):m_width(width),m_height(height),
m_image(new MGPixel[width*height]){
	extract(bitmap,x,y,width,height,alpha);
}

MGImage::~MGImage(){
	delete[] m_image;
}

//Assignment.
MGImage& MGImage::operator=(MGImage&& image2)noexcept {
	delete[] m_image;
	m_width=image2.width();
	m_height=image2.height();
	m_image=image2.m_image;
	image2.m_image=nullptr;
	return *this;
}
MGImage& MGImage::operator=(const MGImage& image2){
	delete[] m_image;
	m_width = image2.width();
	m_height = image2.height();
	m_image = new MGPixel[m_width * m_height];
	copy_color(image2, 0, 0, m_width, m_height);
	return *this;
}

MGPixel& MGImage::operator()(int i, int j){
	MGPixel* pixelP=image();
	return pixelP[j*m_width+i];
}
const MGPixel& MGImage::operator()(int i, int j)const{
	const MGPixel* pixelP=image();
	return pixelP[j*m_width+i];
}

///Generate a cloned MGImage.
///Returned is newed one, must be deleted.
///clone() does not affect this image data.
MGImage* MGImage::clone()const{
	MGImage* image2=new MGImage;
	image2->m_width=m_width;
	image2->m_height=m_height;
	size_t tlen=m_width*m_height;
	image2->m_image=new MGPixel[tlen];
	for(size_t i=0; i<tlen; i++){
		image2->m_image[i]=m_image[i];
	}
	return image2;
}

//Extract a part of bitmap into this, from(x,y) to (x+width, y+height).
void MGImage::extract(
	Gdiplus::Bitmap& bitmap,
	int x,  int y,	//left bottom address of bitmap.
	int width, int height,
	double alpha
){
	int bmW=bitmap.GetWidth(), bmH=bitmap.GetHeight();
	assert(x+width<=bmW);
	assert(y+height<=bmH);


	// Change the MGPixel data format
	MGPixel* pixelP=m_image;
	if(!pixelP)
		return;

	Rect rect(0,0,bmW,bmH);
	BitmapData bmData;
	bitmap.LockBits(&rect,ImageLockModeRead,PixelFormat32bppARGB,&bmData);
	int totalhm1=bmH-1;//total height -1.
	UINT* gdiPixels=(UINT*)bmData.Scan0;;
	for(int j=0; j<height; j++){
		int jpy=j+y;
		int jrow_s=j*m_width;
		int ybitmap=totalhm1-(j+y);
			//This is because MFC's row address is reverse order to MGImage's.
		for(int i=0; i<width; i++){
			MGPixel& pixelij=pixelP[jrow_s+i];
			Gdiplus::Color gdipixel=gdiPixels[ybitmap*bmData.Stride/4 + (i+x)];
			pixelij[0]=gdipixel.GetR();
			pixelij[1]=gdipixel.GetG();
			pixelij[2]=gdipixel.GetB();
			if (0 <= alpha && alpha <= 1.0){
				pixelij[3]=(unsigned char)(gdipixel.GetA() * alpha);
			} else {
				pixelij[3]=gdipixel.GetA();
			}
		}
	}
	bitmap.UnlockBits(&bmData);
}

//resize the image size to (width,height).
void MGImage::resize(
	int width,
	int height
){
	MGPixel* glimage2ptr(new MGPixel[width*height]);
	if(!glimage2ptr)
		return;

	int error=gluScaleImage(GL_RGBA,
		m_width,m_height,GL_UNSIGNED_BYTE,m_image,
		width,height,GL_UNSIGNED_BYTE,glimage2ptr);
	if(error){
		delete[] glimage2ptr;
	}else{
		delete[] m_image;
		m_image=glimage2ptr;
		m_width=width;
		m_height=height;
	}
}

///Fill all the pixel of this with the input color.
void MGImage::fill_color(
	const MGPixel& pdata
){
	for(int j=0; j<m_height; j++){
		int jw1=j*m_width;
		for(int i=0; i<m_width; i++)
			m_image[jw1+i]=pdata;
	}
}

///Fill all the pixels of the range (j, i1,i2) with the input color pdata.
//The color of Pixel(i,j) for i=i1,...i2 is set to pdata.
void MGImage::fill_color(
	const MGPixel& pdata,
	int j,
	int i1,
	int i2
){
	assert(0<=j && j<m_height);
	assert(0<=i1 && i1<m_width);
	assert(0<=i2 && i2<m_width);

	int jw1=j*m_width;
	for(int i=i1; i<=i2 ; i++)
		m_image[jw1+i]=pdata;
}

///Fill all the pixels of ranges with the input color pdata.
//In ranges, range(j,i1,i2) are stored.
//Let m=ranges.size(), then m=3n(always a mutiple of 3) where n is the number of ranges.
//Let (j,i1,i2)=(ranges[3*k],ranges[3*k+1], ranges[3*k+2]) for k=0,...,n-1,
//then PixelData(i,j) for i=i1,...i2 is one range for the height j of this mesh.
void MGImage::fill_color(
	const MGPixel& pdata,
	const std::vector<int>& ranges// Ranges(j,i1,i2) are input.
){
	int n=(int)ranges.size();
	for(int k3=0; k3<n;){
		int j=ranges[k3++];
		int i1=ranges[k3++];
		int i2=ranges[k3++];
		int jw1=j*m_width;
		for(int i=i1; i<=i2; i++)
			m_image[jw1+i]=pdata;
	}
}

void MGImage::fill_color_NoChangeAlpha(const MGPixel& pdata){
	for(int j=0; j<m_height; j++){
		int jw1=j*m_width;
		for(int i=0; i<m_width; i++){
			m_image[jw1+i].setRGB(pdata);
		}
	}
}
void MGImage::fill_color_NoChangeAlpha(const MGPixel& pdata, int j, int i1, int i2){
	assert(0<=j && j<m_height);
	assert(0<=i1 && i1<m_width);
	assert(0<=i2 && i2<m_width);
	int jw1=j*m_width;
	for(int i=i1; i<=i2 ; i++)
		m_image[jw1+i].setRGB(pdata);
}

//Copy all the pixels of image2 into this.
//This image's width and height are not changed. And part of this and image2
//are copied into this.
void MGImage::copy_color(
	const MGImage& image2,//Source image data.
	int x, int y,	///<left bottom address of image2.
	int width, int hight
){
	assert((x + width) <= image2.width() && (y + hight) <= image2.height());

	for (int j = 0; j < hight; j++) {
		int jrow_s = j * m_width;
		int jrow_s2 = (j + y) * image2.m_width;
		for (int i = 0; i < width; i++) {
			MGPixel& pixelij = m_image[jrow_s + i];
			MGPixel& pixel2ij = image2.m_image[jrow_s2 + x + i];
			pixelij = pixel2ij;
		}
	}

}

//Copy all the pixels of ranges in image2 into this.
//In ranges, range(j,i1,i2) are stored.
//Let m=ranges.size(), then m=3n(always a mutiple of 3) where n is the number of ranges.
//Let (j,i1,i2)=(ranges[3*k],ranges[3*k+1], ranges[3*k+2]) for k=0,...,n-1,
//then PixelData(i,j) for i=i1,...i2 is one range for the height j of this mesh.
void MGImage::copy_color(
	const MGImage& image2,//Source image data.
	const std::vector<int>& ranges// Ranges(j,i1,i2) are input.
		///ranges indicate the places of both this and image2.
){
	int n=(int)ranges.size();
	for(int k3=0; k3<n;){
		int j=ranges[k3++];
		int i1=ranges[k3++];
		int i2=ranges[k3++];
		int jw1=j*m_width;
		int jw2=j*image2.m_width;
		for(int i=i1; i<=i2; i++)
			m_image[jw1+i]=image2.m_image[jw2+i];
	}
}

///Test if pixel at (i,j) has zero_alpha value.
///Returns true if alpha value is zero.
bool MGImage::is_zero_alpha(int i, int j)const{
	const MGPixel& pelij=operator()(i,j);
	return pelij.getAlpha()==0;
}

//Test if any one of the four pixels of (i,j) to (i+1,j+1) is non zero or not.
//i must be < width()-1, and j must be < hieght()-1.
//If any one of them is nonzero, return true.
bool MGImage::includeNonZeroAlpha(int i, int j)const{
	const MGImage& image=*this;
	const MGPixel& pelij=image(i,j);
	if(pelij.getAlpha())
		return true;
	const MGPixel& pelip1j=image(i+1,j);
	if(pelip1j.getAlpha())
		return true;
	const MGPixel& pelip1jp1=image(i+1,j+1);
	if(pelip1jp1.getAlpha())
		return true;
	const MGPixel& pelijp1=image(i,j+1);
	if(pelijp1.getAlpha())
		return true;

	return false;
}

//Resize the image size to (width,height) filling the color to
//the extra part for the size(width,height).
//resize_with_fill_color() does not perform scaling to the image.
//width and height can be less than the original length. In this case,
//image trimming will be done.
void MGImage::resize_with_fill_color(
	int width,
	int height,
	const MGPixel& pdata
){
	MGPixel* glimage2ptr(new MGPixel[width*height]);
	if(!glimage2ptr)
		return;

	int w2=m_width;
	if(width<w2) w2=width;
	int h2=m_height;
	if(height<h2) w2=height;

	int i,j;
	MGPixel* pixelP1=m_image;
	for(j=0; j<h2; j++){
		int jw1=j*m_width, jw2=j*width;
		for(i=0; i<w2; i++)
			glimage2ptr[jw2+i]=pixelP1[jw1+i];
		for(;i<width; i++)
			glimage2ptr[jw2+i]=pdata;

	}
	for(; j<height; j++){
		int jw2=j*width;
		for(i=0; i<width; i++)
			glimage2ptr[jw2+i]=pdata;

	}
	delete[] m_image;
	m_image=glimage2ptr;
	m_width=width;
	m_height=height;
}

//Add border color(0,0,0,0)=(alfa=0) of pixel size2 for each perimeter.
void MGImage::resize_and_add_zero_border(int nwidth2,int nheight2){
	resize(nwidth2-4,nheight2-4);
	int i, j;
	MGPixel* image3ptr(new MGPixel[nwidth2*nheight2]);
	if(!image3ptr)
		return;

	int width2by2=nwidth2*2;
	int width2Mborder=nwidth2-4;

	for(i=0; i<width2by2; i++) image3ptr[i]=0;
	for(j=0; j<nheight2-4; j++){
		int jrow_s3=(j+2)*nwidth2;
		int jrow_s2=j*width2Mborder;
		image3ptr[jrow_s3]=0;image3ptr[jrow_s3+1]=0;
		int jrow_s3p2=jrow_s3+2;
		for(i=0; i<width2Mborder; i++){
			image3ptr[jrow_s3p2+i]=m_image[jrow_s2+i];
		}
		image3ptr[jrow_s3p2+width2Mborder]=0;
		image3ptr[jrow_s3p2+width2Mborder+1]=0;
	}
	int erow_s3=(nheight2-2)*nwidth2;
	for(i=0; i<width2by2; i++)
		image3ptr[erow_s3+i]=0;

	delete[] m_image;
	m_image=image3ptr;
	m_width=nwidth2;
	m_height=nheight2;
}

//Compute 2's power of width and height.
void MGImageCompute_2spower(
	int width, int height,
	int& width2, int& height2//The smallest 2's power of width and height will be output.
){
	int nshift=5;
	width2=64;
	while(width2<width && nshift<32){
		width2=width2<<1;nshift++;
	}

	nshift=5;
	height2=64;
	while(height2<height && nshift<32){
		height2=height2<<1 ;nshift++;
	}
}
