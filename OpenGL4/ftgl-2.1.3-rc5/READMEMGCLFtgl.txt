MGCLV12のための制約
MGCLV12ではftglのRendering modesのうち、Polygon meshesのみをサポートしている。
MGCLV12ではOpenGLV4を使用するため、mgVBOのクラスを使用して従来のglBegin, glEndによる
display listの代替としている。
Rendering modesの下記は未対応：
- Bit maps
- Antialiased Pix maps
- Outlines
- Extruded polygon meshes
- Texture maps
- Buffer maps

