MGCL20ソリューション

本ソリューションは、結果としてMGCLtoReleaseフォルダー配下の include, lib, dllの
各フォルダー配下のファイルを更新する。ビルドではlib, binフォルダーのファイルを直接更新し、
ビルド後のイベントとして、MGCLの開発で変更されている可能性のあるファイル (includeファイルすべて)を
MGCLV12/includeからMGCLtoReleaseへコピーする：

MGCLのDLLはMGCLtoReaseフォルダー全体をコピーして利用する。
MGCLtoReleasは下記の構成となる：

MGCLtoRelease
 |
 +- bin(MGCLDLL、OpenGL関連のDLL、OpenGL Shaderプログラム)
 |
 +- lib(MGCLとOpenGL関連のDLL用のlib群)
 |
 +- include(MGCLとOpenGLのincludeファイル群)

１）includeファイル(MGCLtoRelease/include)
includeフォルダー配下に次のサブフォルダーがある：
 GL(glew）, glm, mg(MGCL geometry and other basic classes), 
topo(MGCL topology classes), mgGL(MGCL Graphic classes), Tl2(MGCL tesellation classes), 
mgiges(MGCL iges classes)

2) libファイル(MGCLtoRelease/lib)
DLLとして提供される、glewとMGCLのリンク用ファイル
glew32.lib　とMGCLVxx.lib(xxはMGCLのバージョン番号)

3) dllファイルとOpenGL Shaderプログラム(MGCLtoRelease/bin)
　1) MGCLVxx.dll, 2) glew32.dll, 3) mgclShade.frag, 4) mgclShader.vert
　
