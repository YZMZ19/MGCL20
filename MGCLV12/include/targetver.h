#pragma once

// ビルドターゲットのWindowsバージョン設定
#include <winsdkver.h>

#ifndef _WIN32_WINNT
#define _WIN32_WINNT _WIN32_WINNT_MAXVER
#endif

#include <SDKDDKVer.h>

// 利用するGDI+のバージョン設定
// ここではGDI+1.1を利用

#ifdef GDIPVER
  #if (GDIPVER < 0x0110)
    #error This application needs GDI+1.1! 
  #endif
#else
  #if (_WIN32_WINNT >= 0x0600)
    #define GDIPVER 0x0110
  #else
    #error GDI+1.1 needs windows vista or later.
  #endif
#endif

