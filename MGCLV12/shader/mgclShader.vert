#version 420

// functionID definitions. see MGGL mgGL/glslprogram.h
const int FUNC_DRAW   = 0; // mgGLSL::FuncID::standard(no shading functions);
const int FUNC_SELECT = 1;   // mgGLSL::FuncID::Select
const int FUNC_ZEBRA  = 5; 

const int COORDINATE_WORLD = 0;
const int COORDINATE_SCREENCOORD = 2;// mgGLSL::FuncID::ScreenDrawing
const int COORDINATE_ANCHORPOINT = 3;//mgGLSL::FuncID::AnchorPoint
const int COORDINATE_ANCHORPTSCREEN = 4;

const int DRAW_PRIMITIVE = 0;
const int DRAW_TEXTURE = 11; // mgGLSL::FuncID::Texture
	//Texture that uses texture1 sampler2D.(MGPlaneImage Texture uses this one.)
const int DRAW_GRIDTEXTURE = 12; //mgGLSL::FuncID::gridTexture,GridCursor drawing.

uniform mat4 modelViewProjMatrix;//=projMatrix*modelViewMatrix.
uniform mat4 modelViewMatrix;
uniform mat4 projMatrix;
uniform mat3 normalMatrix;
uniform mat4 ndcMatrix;
uniform mat3 ndcScaleMatrix;
uniform float dpiFactor;

uniform int drawType;		//描画種別(Draw/Texture2D/gridTexture)
uniform int funcType;		// DRAW, SELECT, ANALYSIS
uniform int coordinateType;	// coordinate Type, WORLD,NDC,POINTS(出力座標系) 

uniform float pointSize;
uniform vec3 anchorPoint;


layout(location = 0) in vec3 vPosition;
layout(location = 1) in vec4 vColor;
layout(location = 2) in vec3 vNormal;
layout(location = 3) in vec2 vTextureCoord;

out vec4 fPosition;
out vec4 fColor;
out vec3 fNormal;
out vec2 fTextureCoord;

subroutine vec4 CoordinateModelType();

subroutine (CoordinateModelType) 
vec4 WorldCoordinateModel()
{
	// for shading. without projection(pre-perspective)
	vec4 glPos=modelViewProjMatrix*vec4(vPosition,1.0f);
	fPosition = modelViewMatrix * vec4(vPosition, 1.0f);
	fNormal = normalize(normalMatrix * vNormal);
	return glPos;
}

subroutine (CoordinateModelType) 
vec4 ScreenCodeModel()
{
	// for shading. without projection(pre-perspective)
	vec4 glPos=vec4(vPosition,1.0f);
	fPosition = vec4(vPosition, 1.0f);
	fNormal = vNormal;
	return glPos;
}

subroutine (CoordinateModelType) 
vec4 AnchorPointModel()
{
	fNormal = vNormal;

	vec4 basePos = modelViewProjMatrix*vec4(anchorPoint,1.0f);
	basePos.z = -1.0f*basePos.w;
	vec4 glPos = basePos + vec4(ndcScaleMatrix * (vPosition * dpiFactor),0.0f) * basePos.w;

	vec4 basePosNP = modelViewMatrix*vec4(anchorPoint,1.0f);
	basePosNP.z = -1.0f*basePosNP.w;
	fPosition = basePosNP + vec4(ndcScaleMatrix * (vPosition * dpiFactor) , 0.0f) * basePosNP.w;

	return glPos;
}
 
subroutine (CoordinateModelType) 
vec4 AnchorPointScreenModel()
{
	fNormal = vNormal;

	vec4 basePos = vec4(anchorPoint, 1.0f);
	basePos.z = -1.0f*basePos.w;
	vec4 glPos = basePos + vec4(ndcScaleMatrix * (vPosition * dpiFactor), 0.0f) * basePos.w;

	vec4 basePosNP = vec4(anchorPoint,1.0f);
	basePosNP.z = -1.0f*basePosNP.w;
	fPosition = basePosNP   + vec4(ndcScaleMatrix * (vPosition * dpiFactor), 0.0f) * basePosNP.w;

	return glPos;
}
 
 
subroutine (CoordinateModelType) 
vec4 AnchorPointPtsModel()
{
	fNormal = vNormal;

	vec4 basePos = vec4(ndcScaleMatrix * (anchorPoint * dpiFactor), 1.0f);
	basePos.z = -1.0f*basePos.w;
	vec4 glPos = basePos + vec4(ndcScaleMatrix * (vPosition * dpiFactor), 0.0f) * basePos.w;

	vec4 basePosNP = vec4(ndcScaleMatrix * (anchorPoint * dpiFactor), 1.0f);
	basePosNP.z = -1.0f*basePosNP.w;
	fPosition = basePosNP   + vec4(ndcScaleMatrix * (vPosition * dpiFactor), 0.0f) * basePosNP.w;

	return glPos;
}

void main(){
	gl_PointSize=pointSize;
	fColor=vColor;
	fTextureCoord=vTextureCoord;

	if(coordinateType==COORDINATE_SCREENCOORD){
		gl_Position=ScreenCodeModel();
	} else if(coordinateType==COORDINATE_ANCHORPOINT){
		gl_Position=AnchorPointModel();
	} else if (coordinateType==COORDINATE_ANCHORPTSCREEN){
		gl_Position=AnchorPointScreenModel();
	}else{
		gl_Position=WorldCoordinateModel();
	}
}
