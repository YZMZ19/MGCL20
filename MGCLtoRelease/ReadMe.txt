MGCL　Version12 Releaseファイル

1. ビルドについて
・OpenGLはV4.2以上が稼働する必要があります。
・利用しているC++のバージョンはC++20です。

２．MGCL V12の主な特徴
(1) C++ Move Semantics対応
MGVLV12ではムーブセマンティクスを利用しています。IGES関連プログラムではまだムーブセマンティクス対応の
見直しが必要ですが、他では、特にB-Spline関連クラス（MGNDDArray, MGKnotVector, MGBPointSeq, MGSPointSeqなど）は、
この機能を利用するよう改良されています。

(2) unique_ptrの利用
MGCLV12では積極的にunique_ptrを利用しています。このためにMGCLDefs.hで、多くのUniqueXxxxxを定義しています。
これらを利用してunique_ptrのvector, listなどを利用しています。
従来使用していたauto_ptrの利用はすべてunique_ptrの利用に変更しています。

(3) Tessellationプログラム改良
Tessellationプログラムを大幅に改良しています。今までうまくいかなかった多くの例も可となっています。

(4) B-Splineクラス MGLBRep, MGSBRepのctor改良
今までコンストラクター(ctor)として提供していたB-Splineの生成（Curve, Surfaceともに）を
buildxxxxxとしてメンバー関数化してわかりやすくした。

３．使用している外部ソフト
OpenGL V4.2のため、glm, glew, freetype, ftgl を利用しています。
MGCLの利用はMITライセンスで、利用に際して制約がありませんが、商用で利用する場合上記のライセンス状況をよく調べて利用してください。

freetype, glew, glmは全く手を加えることなくオリジナルなものをそのまま利用していますが、
ftglにはMGCLのVBOクラスを利用できるよう手を加えています。
freetype, glewのdll, .libはMGCLのdownloadで入手できますが、ソースプログラムは各サイトから入手してください。
MGCLtoReleaseにはDLLとそのためのlibファイル、およびコンパイルに必要なincludeファイルのみを梱包してあります。
ftglはMGCLV12.dllとしてMGCLの一部に組み込まれています。

(1) glm
glmはコンパイル時にのみ使用します。includeファイルだけを利用し、cppファイル（実行ファイル）はありません。
Linkのためのlibファイルおよび実行のためのdllはありません。

(2) glew
"The OpenGL Extension Wrangler Library (GLEW) is a cross-platform open-source C/C++ extension loading library"
とあるようにglewは種々のOpenGLモジュールのアドレス解決をしてくれます。
OpenGLのモジュール利用の橋渡しをしてくれ、Linkのためのlibファイルと実行のためのdllを必要とします。
.libファイルとdllが梱包されています。

(3) freetype
freetypeはftglが利用しているため必要で、MGCLは直接利用していません。
libファイルが提供されます。dllはありません。

(4) ftgl
ftglは文字描画のためのプログラムです。
ftglは各種フォントを利用できますが、MGCLは、ftglのフォントすべてをサポートしているわけではなく、
一部しかサポートしていません。今後の課題です。

４．OpenGLのshaderプログラム
MGCLは、一部Windowsに依存しています。
画像処理（.pngなどの画像ファイルの読み込み）、OpenGLのWindows用contextがWindows対応となっています。
OpenGLの表示関連を利用する場合、Shaderプログラムも必要となりますが、これらは<MGCLtoRelease>\binに含まれています。
実行時に必要となりますので、利用プログラムのパスが通るところに複写する必要があります。
