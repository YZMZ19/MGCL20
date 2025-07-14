MGCLを利用した開発にはincludeファイル(.hファイル）、リンクのための.libファイル、
および実行ファイル(DLLとOpenGL Shaderプログラム)が必要とされる。

１）includeファイル
includeフォルダー配下に次のサブフォルダーがある：
 GL(glew）, glm, mg(MGCL geometry and other basic classes), 
topo(MGCL topology classes), mgGL(MGCL Graphic classes), Tl2(MGCL tesellation classes), 
mgiges(MGCL iges classes)

2) libファイル
DLLとして提供される、glewとMGCLのリンク用ファイル
glew32.lib　とMGCLVxx.lib(xxはMGCLのバージョン番号)

3) dllファイルとOpenGL Shaderプログラム
　1) MGCLVxx.dll, 2) glew32.dll, 3) mgclShade.frag, 4) mgclShader.vert

MGCL20のプロジェクトではビルド後のイベントとして、MGCLの開発で変更されている可能性のあるファイルを
fugenソリューションのMGCLtoReleaseへコピーするためのバッチコマンドを走らせている：

（１）作成されたdllとその利用のためのlibファイルのコピー
 $(FUGEN_MGCL_DIR)bin\　と　$(FUGEN_MGCL_DIR)lib\　へ
（２）MGCLをdllで利用するためのincludeファイル
mg, mgGL, tl2, topo, mgigesの各フォルダーを$(FUGEN_MGCL_DIR)includeへ

ここで、FUGEN_MGCL_DIRは MGCLSolution.propsで定義されているマクロで、
下記のフォルダー構成を前提としている：

MGCLtoRelease :fugenプロジェクトで利用するMGCLのfolder.
 |
 +- bin(MGCLDLLとOpenGL関連のDLL)
 |
 +- lib(MGCLとOpenGL関連のDLL用のlib群)
 |
 +- include(MGCLとOpenGLのincludeファイル群)

