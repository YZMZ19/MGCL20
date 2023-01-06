#version 420

// functionID definitions. see MGGL mgGL/glslprogram.h
const int FUNC_DRAW    = 0; // mgGLSL::FuncID::standard(no shading functions);
const int FUNC_SELECT = 1;   // mgGLSL::FuncID::Select
const int FUNC_ZEBRA = 5; 

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

uniform int drawType;//描画種別(Draw/Texture2D/gridTexture)
uniform int funcType;		// DRAW, SELECT, ANALYSIS
uniform int coordinateType; // coordinate Type, WORLD,NDC,POINTS(出力座標系) 

uniform sampler2D texture2D;


uniform int ShaderMode;
// ShaderMode Values See MGGL mgGL/glslprogram.h
const int NOSHADING = 0; // mgGLSL::ShadeMode::NoShading
const int SHADING = 1;   // mgGLSL::ShadeMode::Shading


// Force Lighing mode 
// if if value is false and no enable Lights, use GL_LIGHT0 Default vavlue
uniform bool ForceLight;

///--- LightModel Grobals.
uniform bool LightTwoSides;

/// ゼブラ処理パラメータ
uniform int ZebraAxis;
uniform float ZebraSize;

// Lighting props. compatible OpenGL 1.1
struct LightProperties {
	bool isEnabled;

	vec4 ambientColor;
	vec4 diffuseColor;
	vec4 specularColor;

	vec4 position;
	vec3 spotDirection;
	float spotCutoff;
	float spotExponent;
    float constantAttenuation;
	float linearAttenuation;
    float quadraticAttenuation;
};

// the set of lights to apply, per invocation of this shader
const int MaxLights = 10; // See. MGCL mgGL/glslprogram.h  LIGHT_NUM
uniform LightProperties Lights[MaxLights];

/// LightMode Constants.
const vec4 LightModelAmbientColor = vec4(0.2f, 0.2f, 0.2f, 1.0f);
const bool LightModelLocalViewer = false;
const vec3 EyeDirection = vec3(0.0f, 0.0f, 1.0f);

//// Material Defaults.
// Materialは両面で同じプロパティ。
// おそらくこれはデフォルトのまま使ってるので
const vec4 MaterialAmbient  = vec4(0.2f, 0.2f, 0.2f, 1.0f);
const vec4 MaterialDiffulse = vec4(0.8f, 0.8f, 0.8f, 1.0f);
const vec4 MaterialSpecular = vec4(0.0f, 0.0f, 0.0f, 1.0f);
const vec4 MaterialEmission = vec4(0.0f, 0.0f, 0.0f, 1.0f);
const float MaterialShiness = 0.0f;

const vec4 globalAmbientColor = MaterialAmbient * LightModelAmbientColor;

/// attributes.

in vec4 fColor;
in vec3 fNormal;
in vec2 fTextureCoord;
in vec4 fPosition;

out vec4 FragColor;

vec4 ApplyLightEffects(vec4 Position, vec3 Normal, vec4 Color)
{
	vec4 lightEffectColor = vec4(0.0); // or, to a global ambient light

    // loop over all the lights
	bool bLightEnabled = false;
	for(int i = 0; i <MaxLights; ++i){
		if(! Lights[i].isEnabled) continue;
		
		bLightEnabled = true;

		// Each light calculates.
//		vec4 lightPos = modelViewMatrix *  Lights[i].position;
		vec4 lightPos = Lights[i].position;
		vec3 lightDir = - normalize(lightPos.xyz);
        
		// Light.
		float attenuation = 1.0f;

		// Spot.effect.
		float spotCoef = 1.0;

		if(Lights[i].position.w>0.0f){
//			vec3 dir = Position.xyz - lightPos.xyz;
			vec3 dir = (Position -  (lightPos + modelViewMatrix *  vec4(0.0f, 0.0f, 0.0f, 1.0f))).xyz;

			float dist = length(dir);
			lightDir = dir/dist;

			attenuation = 1.0f /
				( Lights[i].constantAttenuation
				+ Lights[i].linearAttenuation    * dist
				+ Lights[i].quadraticAttenuation * dist * dist);

			if(Lights[i].spotCutoff < 180.0f){
//				float spotCos = dot(lightDir, normalize(normalMatrix * Lights[i].spotDirection));
				float spotCos = dot(lightDir, normalize(Lights[i].spotDirection));

				if(spotCos < cos(radians(Lights[i].spotCutoff)) ){
					spotCoef = 0.0f;
				} else {
					spotCoef = pow(max(spotCos, 0.0f), Lights[i].spotExponent);
				}
			}
		}

		/// Ambient 
		vec4 colorAmbient = Lights[i].ambientColor * MaterialAmbient;

		// Diffuse. 
		float diffuseRow = dot(-lightDir, Normal);
		vec4 colorDiffuse = max(diffuseRow, 0.0f) * Lights[i].diffuseColor * MaterialDiffulse;

		// Specular.
		float specularRow = 0.0f;
		vec4 colorSpecular = vec4(0.0f);

		if(diffuseRow > 0.0f){
//			vec3 halfVector = normalize(lightDir + normalMatrix * EyeDirection);
			vec3 halfVector = normalize(lightDir + EyeDirection);

			specularRow = dot(halfVector, Normal);
			colorSpecular = pow(specularRow, MaterialShiness) * Lights[i].specularColor * MaterialSpecular;
		}

        // Accumulate all the lights’ effects
		lightEffectColor +=  attenuation *  spotCoef * (colorAmbient + colorDiffuse + colorSpecular);
	}

	// GL_LIGHT0のデフォルト値で計算
	if( !bLightEnabled && ForceLight){

		// Light.
		float attenuation = 1.0f;

		// Spot.effect.
		float spotCoef = 1.0;
		
		vec4 defaultAmbient =  vec4(0.0f, 0.0f, 0.0f, 1.0f); 
		vec4 defaultDiffuse =  vec4(1.0f, 1.0f, 1.0f, 1.0f);
		vec4 defaultSpecular = vec4(1.0f, 1.0f, 1.0f, 1.0f);
		vec4 defaultPosition = vec4(0.0f, 0.0f, 1.0f, 0.0f);

		vec4 colorAmbient = defaultAmbient * MaterialAmbient;

//		vec3 lightDir = -normalize(normalMatrix * defaultPosition.xyz);
		vec3 lightDir = -normalize(defaultPosition.xyz);
		float diffuseRow = dot(-lightDir, Normal);
		vec4 colorDiffuse = max(diffuseRow, 0.0f) * defaultDiffuse * MaterialDiffulse;

		float specularRow = 0.0f;
		vec4 colorSpecular = vec4(0.0f);

		if(diffuseRow > 0.0f){
//			vec3 halfVector = normalize(lightDir + normalMatrix * EyeDirection);
			vec3 halfVector = normalize(lightDir + EyeDirection);
			specularRow = dot(halfVector, Normal);
			colorSpecular = pow(specularRow, MaterialShiness) * defaultSpecular * MaterialSpecular;
		}
	
		lightEffectColor =  attenuation *  spotCoef * (colorAmbient + colorDiffuse + colorSpecular);
	} 

	if(!bLightEnabled && !ForceLight){
		return Color;
	}

	vec3 ResultColor = min( vec3(Color * lightEffectColor + MaterialAmbient * LightModelAmbientColor + MaterialEmission), vec3(1.0f) );
	return vec4(ResultColor.rgb , Color.a);

}

vec4 LighitingProc(vec4 Position, vec3 Normal, vec4 Color){
	vec3 prjNormal = normalMatrix * Normal;	
	if(LightTwoSides){
		if(gl_FrontFacing){
			return ApplyLightEffects(Position, Normal, Color);
		} else {
			return ApplyLightEffects(Position, -Normal, Color);
		}
	} else {
		return ApplyLightEffects(Position, Normal, Color);
	}
}

subroutine vec4 Draw();
subroutine vec4 BaseColor();

//subroutine uniform Draw DrawFunc;
//subroutine uniform BaseColor GetBaseColor;

subroutine (Draw)
vec4 SelectionFunc()
{
	return fColor;
}


subroutine (Draw)
vec4 ZebraAnalysisFunc()
{
	vec3 v = fPosition.xyz / fPosition.w;
	vec3 fReflect = reflect(v, fNormal);
	vec3 r = fReflect;

	float coef = 0.0f;
	if(ZebraAxis == 0){
		coef = step(0.5f, fract(1.0f / ZebraSize * atan(r.z, r.y)));
	} else {
		coef = step(0.5f, fract(1.0f / ZebraSize * atan(r.y, r.x)));
	}

	return (fColor* coef);
}

subroutine (Draw)
vec4 DisplayFunc(){
	vec4 BaseColor=fColor;
	// Texture processing
	if(drawType == DRAW_TEXTURE ){
		// When using texture1.
		BaseColor=texture(texture2D,fTextureCoord);
	}else if(drawType == DRAW_GRIDTEXTURE ){
		BaseColor=texture(texture2D,gl_PointCoord);
	}

	if(BaseColor.a == 0.0f){
		discard;
	}

	if( ShaderMode == NOSHADING){
		return BaseColor;
	} else if (ShaderMode == SHADING){
		return LighitingProc(fPosition, fNormal, BaseColor);
	}
}

void main(){
	if(funcType == FUNC_SELECT){
		FragColor = SelectionFunc();
	} else if (funcType == FUNC_ZEBRA){
		FragColor = ZebraAnalysisFunc();
	} else if (funcType == FUNC_DRAW){
		FragColor = DisplayFunc();
	} else {
		discard;
	}
}

/** 廃止
void main1(){

	if(funcType == FUNC_SELECT){
		FragColor = fColor;
	}else if(funcType == 5 ){
		vec3 v = fPosition.xyz / fPosition.w;
		vec3 fReflect = reflect(v, fNormal);
		vec3 r = fReflect;
//		FragColor = vec4(step(0.5, fract(2.0 * fReflect.z / fReflect.y)));
		FragColor = vec4(vec3(step(0.5f, fract(10.0 * atan(r.z, r.y)))),1.0f);
	} else {

		vec4 Color=fColor;
		// Texture processing
		if(functionID == FUNC_TEXTURE ){
			// When using texture1.
			Color=texture(texture2D,fTextureCoord);
		}else if(functionID == FUNC_GRIDTEXTURE ){
			Color=texture(texture2D,gl_PointCoord);
		}

		if(Color.a == 0.0f){
			discard;
		}

		if( ShaderMode == NOSHADING){
			FragColor = Color;
		} else if (ShaderMode == SHADING){
			FragColor = LighitingProc(fPosition, fNormal, Color);
		}
	}
}
**/