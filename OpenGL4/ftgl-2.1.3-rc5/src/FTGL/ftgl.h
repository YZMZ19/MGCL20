/*
 * FTGL - OpenGL font library
 *
 * Copyright (c) 2001-2004 Henry Maddocks <ftgl@opengl.geek.nz>
 * Copyright (c) 2008 Sam Hocevar <sam@zoy.org>
 * Copyright (c) 2008 Sean Morrison <learner@brlcad.org>
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#ifndef __ftgl__
#define __ftgl__

/* We need the Freetype headers */
#include <ft2build.h>
#include FT_FREETYPE_H
#include FT_GLYPH_H
#include FT_OUTLINE_H

/* Floating point types used by the library */
typedef double   FTGL_DOUBLE;
typedef float    FTGL_FLOAT;

/* Macros used to declare C-linkage types and symbols */
#   define FTGL_BEGIN_C_DECLS extern "C" { namespace FTGL {
#   define FTGL_END_C_DECLS } }

namespace FTGL
{
    typedef enum
    {
        RENDER_FRONT = 0x0001,
        RENDER_BACK  = 0x0002,
        RENDER_SIDE  = 0x0004,
        RENDER_ALL   = 0xffff
    } RenderMode;

    typedef enum
    {
        ALIGN_LEFT    = 0,
        ALIGN_CENTER  = 1,
        ALIGN_RIGHT   = 2,
        ALIGN_JUSTIFY = 3
    } TextAlignment;
}

// Compiler-specific conditional compilation
#ifdef _MSC_VER // MS Visual C++

    // Disable various warning.
    // 4786: template name too long
    #pragma warning(disable : 4251)
    #pragma warning(disable : 4275)
    #pragma warning(disable : 4786)


	// Define FTGL_LIB_PRAGMAS to 1 to include library
	// pragmas or to 0 to exclude library pragmas.
	// The default behavior depends on the compiler/platform.
	//
	#ifndef FTGL_LIB_PRAGMAS
		#if ( defined(_MSC_VER) || defined(__WATCOMC__) ) && !defined(_WIN32_WCE)
			#define FTGL_LIB_PRAGMAS 1
		#else
			#define FTGL_LIB_PRAGMAS 0
		#endif
	#endif


    // The following definitions control how symbols are exported.
    // If the target is a static library ensure that FTGL_LIBRARY_STATIC
    // is defined. If building a dynamic library (ie DLL) ensure the
    // FTGL_LIBRARY macro is defined, as it will mark symbols for
    // export. If compiling a project to _use_ the _dynamic_ library
    // version of the library, no definition is required.
    #ifdef FTGL_LIBRARY_STATIC      // static lib - no special export required
		#define FTGL_EXPORT
        /* Link with static lib */
		//#if FTGL_LIB_PRAGMAS
		//	#pragma comment (lib, "ftgl_static.lib")
		//#endif
		
    #elif FTGL_LIBRARY              // dynamic lib - must export/import symbols appropriately.
		#define FTGL_EXPORT   __declspec(dllexport)
    #elif FTGL_DLL
		#define FTGL_EXPORT   __declspec(dllimport)
        /* Link with ftgl link lib */
		/* #if FTGL_LIB_PRAGMAS
			#pragma comment (lib, "ftgl.lib")
		#endif */
	#endif

/* Drag in other Windows libraries as required by FreeGLUT */
#   if FTGL_LIB_PRAGMAS
#       pragma comment (lib, "glu32.lib")    /* link OpenGL Utility lib     */
#       pragma comment (lib, "opengl32.lib") /* link Microsoft OpenGL lib   */
#   endif

#else
    // Compiler that is not MS Visual C++.
    // Ensure that the export symbol is defined (and blank)
    #define FTGL_EXPORT
#endif

#include <FTGL/FTPoint.h>
#include <FTGL/FTBBox.h>
#include <FTGL/FTBuffer.h>

#include <FTGL/FTGlyph.h>
#include <FTGL/FTPolyGlyph.h>

#include <FTGL/FTFont.h>
#include <FTGL/FTGLPolygonFont.h>

#endif  //  __ftgl__
