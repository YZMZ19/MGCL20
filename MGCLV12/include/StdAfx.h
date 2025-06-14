/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
// StdAfx.h : Include File to define standard indespencible include files
#pragma once

#include "targetver.h"

// C/C++ standard library
#include <assert.h>
#include <malloc.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <algorithm>
#include <bitset>
#include <deque>
#include <functional>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <list>
#include <memory>
#include <sstream>
#include <stack>
#include <string>
#include <utility>
#include <vector>

// Windows specific library
#define VC_EXTRALEAN // Windows �w�b�_�[����w�ǎg�p����Ȃ��X�^�b�t�����O���܂��B
#include <afxwin.h> // MFC �̃R�A����ѕW���R���|�[�l���g
#include <afxext.h> // MFC �̊g������
#include <afxdisp.h> // MFC �̃I�[�g���[�V���� �N���X
#include <afxdtctl.h> // MFC �� Internet Explorer 4 �R���� �R���g���[�� �T�|�[�g

#ifndef _AFX_NO_AFXCMN_SUPPORT
#include <afxcmn.h> // MFC �� Windows �R���� �R���g���[�� �T�|�[�g
#endif // _AFX_NO_AFXCMN_SUPPORT

#ifdef _CONSOLE

using GLenum = unsigned int;
using GLuint = unsigned int;
using GLint = int;
using GLsizei = int;

#define GL_FILL 0x1B02
#define GL_LINES 0x0001
#define GL_LINE_STRIP 0x0003
#define GL_FRONT 0x0404
#define GL_BACK 0x0405
#define GL_FRONT_AND_BACK 0x0408
#define GLsizei 0x0404
#define GL_REPEAT 0x2901
#define GL_LINEAR 0x2601
#define GL_TEXTURE_2D 0x0DE1

#else

#include <gdiplus.h>
#include <Gdiplusimaging.h>

// OpenGL
#include <gl/glew.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#endif //_CONSOLE
