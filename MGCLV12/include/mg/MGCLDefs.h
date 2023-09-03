#pragma once

// 標準ライブラリのストリームの文字コード対応
#include <tchar.h>
#include <sstream> 
#include <vector>
#include <memory>

//
// The following macros are used to enable DLL export/import.
// MG_DLL_DECLR for class, global functions, and global valiable values in declaration.
#ifdef NO_MGCLDLL
	// Not DLL.
#	define MG_DLL_DECLR

#else

#	ifdef MGCLDLL_EXPORTS
#		define MG_DLL_DECLR		__declspec(dllexport)
#	else                     // Import DLL.
#		define MG_DLL_DECLR		__declspec(dllimport)
#	endif// MGCL_IMPORTS

#endif	//NO_MGCLDLL

class MGGel;
class MGObject;
class MGGeometry;
class MGStl;
class MGGLAttrib;
class MGCurve;
class MGLBRep;
class MGRLBRep;
class MGFSurface;
class MGSurface;
class MGSBRep;
class MGRSBRep;
class MGCell;
class MGBCell;
class MGPVertex;
class MGBVertex;
class MGFace;
class MGLoop;
class MGEdge;
class MGLight;
class mgSysGL;
class MGPickObject;

//To treat UNICODE
using tstring=std::basic_string<TCHAR>;
using tostringstream=std::basic_ostringstream<TCHAR>;
using tistringstream=std::basic_istringstream<TCHAR>;
using tostream=std::basic_ostream<TCHAR>;
using tistream=std::basic_istream<TCHAR>;
using tifstream=std::basic_ifstream<TCHAR>;
using tofstream=std::basic_ofstream<TCHAR>;

using UniqueGel = std::unique_ptr<MGGel>;
using UniqueObject = std::unique_ptr<MGObject>;
using UniqueStl = std::unique_ptr<MGStl>;
using UniqueGLAttrib = std::unique_ptr<MGGLAttrib>;
using UniqueGLAttribVec = std::vector<UniqueGLAttrib>;
using UniqueGeometry = std::unique_ptr<MGGeometry>;
using UniqueCurve = std::unique_ptr<MGCurve>;
using UniqueLBRep = std::unique_ptr<MGLBRep>;
using UniqueRLBRep = std::unique_ptr<MGRLBRep>;
using UniqueFSurface = std::unique_ptr<MGFSurface>;
using UniqueSurface = std::unique_ptr<MGSurface>;
using UniqueSBRep = std::unique_ptr<MGSBRep>;
using UniqueRSBRep = std::unique_ptr<MGRSBRep>;
using UniqueCell = std::unique_ptr<MGCell>;
using UniquePVertex = std::unique_ptr<MGPVertex>;
using SharedBCell = std::shared_ptr<MGBCell>;
using SharedBVertex = std::shared_ptr<MGBVertex>;
using UniqueFace = std::unique_ptr<MGFace>;
using UniqueEdge = std::unique_ptr<MGEdge>;
using SharedEdge = std::shared_ptr<MGEdge>;
using UniqueLoop = std::unique_ptr<MGLoop>;
using UniqueLight = std::unique_ptr<MGLight>;
using UniqueSysGL = std::unique_ptr<mgSysGL>;
using UniquePickObject = std::unique_ptr<MGPickObject>;

//Extract const T* from originalVec's first to last into o.
//OutputIterator o must have appropriate size to store the pointer.
template<class InputIterator, class OutputIterator>
void extractConstPointerVec(
	InputIterator first,
	InputIterator last,
	OutputIterator o
){
	for (;first!=last; first++, o++)
		*o = &(**first);
};

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ は前行の直前に追加の宣言を挿入します。

#ifdef _UNICODE
  #define TCAST const wchar_t*
  #define COUT std::wcout
  #define CERR std::wcerr
#else
  #define TCAST const char*
  #define COUT std::cout
  #define CERR std::cerr
#endif

