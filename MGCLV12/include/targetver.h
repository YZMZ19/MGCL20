#pragma once

// �r���h�^�[�Q�b�g��Windows�o�[�W�����ݒ�
#include <winsdkver.h>

#ifndef _WIN32_WINNT
#define _WIN32_WINNT _WIN32_WINNT_MAXVER
#endif

#include <SDKDDKVer.h>

// ���p����GDI+�̃o�[�W�����ݒ�
// �����ł�GDI+1.1�𗘗p

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

