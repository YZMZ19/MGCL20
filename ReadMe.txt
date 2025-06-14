
１）MGCL_DIR
下記のフォルダー構成を前提として、
MGCL_DIRマクロのみを定義すれば残りのマクロはすべて定義されるfolder構造となっている：
提供されるsolution.propsではfugenのソリューションディレクトリ
（fugen3.slnが配置されるところ）は以下にMGCLtoReleaseフォルダーがあり、
このフォルダーがMGCL_DIRとして定義されている。

MGCLV11toRelease :最上位のMGCL Solutionのfolder.
 |
 +- bin-debug-win32(MGCLDLLとOpenGL関連のDLL)
 |
 +- lib-debug-win32(MGCLとOpenGL関連のDLL用のlib群)
 |
 +- MGCLinclude(MGCLのincludeファイル群)
 |
 +- OpenGLinclude(OpenGL：Freetype, Ftgl, Glew) のincludeファイル群

２）MGCL_INC_DIR　MGCLのincludeファイルディレクトリー
３）OPENGL_INC_DIR　OpenGLのincludeファイルディレクトリー
４）MGCL_LIB_DIR　MGCLのdll用libファイルの作成先
５）MGCL_BIN_DIR　MGCLのdllファイルの作成先
