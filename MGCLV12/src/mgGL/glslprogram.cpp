#include "stdafx.h"
#include <sys/stat.h>
#include "mg/Position.h"
#include "mgGL/Color.h"
#include "mgGL/StaticGLAttrib.h"
#include "mgGL/glslprogram.h"

namespace{GLenum glErr;};

mgGLSLProgram::mgGLSLProgram():m_handle(0), m_linked(false){;}
mgGLSLProgram::~mgGLSLProgram(){
	freeProgram();
}

///Delete this program and initialize this.
void mgGLSLProgram::freeProgram(){
	if(m_handle)
		glDeleteProgram(m_handle);
	m_handle=0;
	m_linked=false;
	m_logString="";
	if (getCurrentGLSLProgram() == this)
		setCurrentGLSLProgram(nullptr);
}

bool mgGLSLProgram::compileShaderFromFile(
	const char* fileName,
	GLSLShaderType type
){
    if( !fileExists(fileName) ){
        m_logString = "File not found.";
        return false;
    }

    if( m_handle <= 0 ){
        m_handle = glCreateProgram();
        if( m_handle == 0) {
            m_logString = "Unable to create shader program.";
            return false;
        }
    }

    std::ifstream inFile( fileName, std::ios::in );
    if( !inFile ) {
        return false;
    }

    std::ostringstream code;
    while( inFile.good() ) {
        int c = inFile.get();
        if( ! inFile.eof() ) code << (char) c;
    }
    inFile.close();
    return compileShaderFromString(code.str(), type);
}

bool mgGLSLProgram::compileShaderFromString(
	const std::string& source, GLSLShaderType type
){
    GLuint shaderHandle = 0;
    switch( type ) {
    case VERTEX:
        shaderHandle = glCreateShader(GL_VERTEX_SHADER);
        break;
    case FRAGMENT:
        shaderHandle = glCreateShader(GL_FRAGMENT_SHADER);
        break;
    case GEOMETRY:
        shaderHandle = glCreateShader(GL_GEOMETRY_SHADER);
        break;
    case TESS_CONTROL:
        shaderHandle = glCreateShader(GL_TESS_CONTROL_SHADER);
        break;
    case TESS_EVALUATION:
        shaderHandle = glCreateShader(GL_TESS_EVALUATION_SHADER);
        break;
    default:
        return false;
    }

    const char* c_code = source.c_str();
    glShaderSource( shaderHandle, 1, &c_code, NULL );

    // Compile the shader
    glCompileShader(shaderHandle );

    // Check for errors
    int result;
    glGetShaderiv( shaderHandle, GL_COMPILE_STATUS, &result );
    if( GL_FALSE == result ) {
        // Compile failed, store log and return false
        int length = 0;
        m_logString = "";
        glGetShaderiv(shaderHandle, GL_INFO_LOG_LENGTH, &length );
        if( length > 0 ) {
            char* c_log = new char[length];
            int written = 0;
            glGetShaderInfoLog(shaderHandle, length, &written, c_log);
            m_logString = c_log;
            delete[] c_log;
        }

        return false;
    }else{
        // Compile succeeded, attach shader and return true
        glAttachShader(m_handle, shaderHandle);
        return true;
    }
}

bool mgGLSLProgram::link(){
    if( m_linked ) return true;
    if( m_handle <= 0 ) return false;

    glLinkProgram(m_handle);

    int status = 0;
    glGetProgramiv( m_handle, GL_LINK_STATUS, &status);
    if( GL_FALSE == status ) {
        // Store log and return false
        int length = 0;
        m_logString = "";

        glGetProgramiv(m_handle, GL_INFO_LOG_LENGTH, &length );

        if( length > 0 ) {
            char* c_log = new char[length];
            int written = 0;
            glGetProgramInfoLog(m_handle, length, &written, c_log);
            m_logString = c_log;
            delete[] c_log;
        }

        return false;
    }else{
		build_attribLocations();
		build_attribLocationsLight();
        m_linked = true;
        return m_linked;
    }
}

void DefaultLightsSetup();
void mgGLSLProgram::use(){
    if( m_handle <= 0 || (! m_linked) ) return;
    glUseProgram( m_handle );
	setCurrentGLSLProgram(this);

	// Set default lighting.
	DefaultLightsSetup();
}

void mgGLSLProgram::setCurrentGLSLProgram(mgGLSLProgram* glsl){
	m_CurrrentGLSL=glsl;
}
mgGLSLProgram* mgGLSLProgram::getCurrentGLSLProgram(){
	return m_CurrrentGLSL;
}

std::string& mgGLSLProgram::log(){return m_logString;}

int mgGLSLProgram::getHandle(){return m_handle;}

bool mgGLSLProgram::isLinked(){return m_linked;}

void mgGLSLProgram::bindAttribLocation( GLuint location, const char* name){
    glBindAttribLocation(m_handle, location, name);
}

void mgGLSLProgram::bindFragDataLocation( GLuint location, const char* name ){
    glBindFragDataLocation(m_handle, location, name);
}

///in vec3 vPositionのlocationを求める。
int mgGLSLProgram::getvPositionLocation()const{
	const static char name[] ={"vPosition"};
    int loc = getAttribLocation(name);
	assert(loc>=0);
	return loc;
}

///in vec4 vColorのlocationを求める。
int mgGLSLProgram::getvColorLocation()const{
	const static char name[] ={"vColor"};
    int loc = getAttribLocation(name);
	assert(loc>=0);
	return loc;
}

///in vec3 vNormalのlocationを求める。
int mgGLSLProgram::getvNormalLocation()const{
	const static char name[] ={"vNormal"};
    int loc = getAttribLocation(name);
	//assert(loc>=0);
	return loc;
}

///in vec2 vTextureのlocationを求める。
int mgGLSLProgram::getvTextureCoordLocation()const{
	const static char name[] ={"vTextureCoord"};
    int loc = getAttribLocation(name);
	//assert(loc>=0);
	return loc;
}

void mgGLSLProgram::setUniform( GLint loc, float x, float y, float z){
	//ASSERT(loc>=0);
	glUniform3f(loc,x,y,z);
}

void mgGLSLProgram::setUniform( GLint loc, const glm::vec3 & v){
	//ASSERT(loc>=0);
	glUniform3f(loc,v.x,v.y,v.z);
}

void mgGLSLProgram::setUniform( GLint loc, const glm::vec4 & v){
	//ASSERT(loc>=0);
	glUniform4f(loc,v.x,v.y,v.z,v.w);
}

void mgGLSLProgram::setUniform( GLint loc, const glm::mat4 & m){
	//ASSERT(loc>=0);
	glUniformMatrix4fv(loc, 1, GL_FALSE, &m[0][0]);
}

void mgGLSLProgram::setUniform( GLint loc, const glm::mat3 & m){
	//ASSERT(loc>=0);
	glUniformMatrix3fv(loc, 1, GL_FALSE, &m[0][0]);
}

void mgGLSLProgram::setUniform( GLint loc, float val ){
	//ASSERT(loc>=0);
	glUniform1f(loc, val);
}

void mgGLSLProgram::setUniform( GLint loc, int val ){
	//ASSERT(loc>=0);
	glUniform1i(loc, val);
}

void mgGLSLProgram::setUniform( GLint loc, bool val ){
	//ASSERT(loc>=0);
	glUniform1i(loc, val);
}

void mgGLSLProgram::printActiveUniforms(){
    GLint nUniforms, size, location, maxLen;
    GLchar* name;
    GLsizei written;
    GLenum type;

    glGetProgramiv( m_handle, GL_ACTIVE_UNIFORM_MAX_LENGTH, &maxLen);
    glGetProgramiv( m_handle, GL_ACTIVE_UNIFORMS, &nUniforms);

    name = (GLchar*) malloc( maxLen );

    COUT<<"Active uniforms, Name:Index"<<std::endl;
    for( int i = 0; i < nUniforms; ++i ) {
        glGetActiveUniform( m_handle, i, maxLen, &written, &size, &type, name );
        location = glGetUniformLocation(m_handle, name);
		COUT<<name<<":"<<location<<std::endl;
    }
    free(name);
}

void mgGLSLProgram::printActiveAttribs(){
    GLint written, size, location, maxLength, nAttribs;
    GLenum type;
    GLchar* name;

    glGetProgramiv(m_handle, GL_ACTIVE_ATTRIBUTE_MAX_LENGTH, &maxLength);
    glGetProgramiv(m_handle, GL_ACTIVE_ATTRIBUTES, &nAttribs);

    name = (GLchar *) malloc( maxLength );
    COUT<<"Active attribs, Name:Index"<<std::endl;
    for( int i = 0; i < nAttribs; i++ ) {
        glGetActiveAttrib( m_handle, i, maxLength, &written, &size, &type, name );
        location = glGetAttribLocation(m_handle, name);
		COUT<<name<<":"<<location<<std::endl;
    }

    free(name);
}

int mgGLSLProgram::getUniformLocation(const char* name)const{
#if _DEBUG
	GLint loc = glGetUniformLocation(m_handle, name);
	if(loc<0){
		std::cerr << "GLSL::getUniformLocation: ERROR!! " << name << std::endl;
	}
    return loc;
#else
    return glGetUniformLocation(m_handle, name);
#endif
}

int mgGLSLProgram::getAttribLocation(const char* name)const{
#if _DEBUG
	GLint loc = glGetAttribLocation(m_handle, name);
	if(loc<0){
		std::cerr << "GLSL::getAttribLocation: ERROR!! " << name << std::endl;
	}
    return loc;
#else
    return glGetAttribLocation(m_handle, name);
#endif
}

bool mgGLSLProgram::fileExists( const std::string& fileName)const{
    struct stat info;
    int ret = -1;

    ret = stat(fileName.c_str(), &info);
    return 0 == ret;
}

void mgGLSL::printOpenGLError(int errorCode){
	CString msg(gluErrorString(errorCode));
	COUT<<"glError:"<<(TCAST)msg<<std::endl;
}
int mgGLSL::checkForOpenGLError(const char* file, int line) {
    //
    // Returns 1 if an OpenGL error occurred, 0 otherwise.
    //
    GLenum glErr;
    int    retCode = 0;
	glErr = glGetError();
    while (glErr != GL_NO_ERROR){
        glErr = glGetError();
		CString msg(gluErrorString(glErr));
		COUT<<"glError in file:"<<file<<" "<<line<<" "<<(TCAST)msg<<std::endl;
        retCode = 1;
    }
    return retCode;
}

void mgGLSL::dumpGLInfo(bool dumpExtensions) {
    CString renderer(glGetString( GL_RENDERER ));
    CString vendor(glGetString( GL_VENDOR ));
    CString version(glGetString( GL_VERSION ));
    CString glslVersion(glGetString( GL_SHADING_LANGUAGE_VERSION ));

    GLint major, minor;
    glGetIntegerv(GL_MAJOR_VERSION, &major);
    glGetIntegerv(GL_MINOR_VERSION, &minor);

	COUT<<"GL Vendor    :"<<(TCAST)vendor<<std::endl;
	COUT<<"GL Renderer  :"<<(TCAST)renderer<<std::endl;
	COUT<<"GL Version   :"<<(TCAST)version<<","<<major<<"."<<minor<<std::endl;
	COUT<<"GLSL Version :"<<(TCAST)glslVersion<<std::endl;

    if( dumpExtensions ) {
        GLint nExtensions;
        glGetIntegerv(GL_NUM_EXTENSIONS, &nExtensions);
        for( int i = 0; i < nExtensions; i++ ) {
			CString msg(glGetStringi(GL_EXTENSIONS, i));
			COUT<<(TCAST)msg<<std::endl;
        }
    }
}

///Set static color as selection name.
void mgGLSL::setColorAsSelectionName(unsigned name){
	if(!name)
		return;

	mgGLSLProgram::SELECT_NAME nub;
	nub.uiName=name;

	GLuint vColorLoc=mgGLSLProgram::vColor;
	glVertexAttribPointer(vColorLoc,4,GL_FLOAT,GL_FALSE,0,0);
	glErr=glGetError();
	glVertexAttrib4Nub(vColorLoc,
		nub.ucName[0],nub.ucName[1],nub.ucName[2],nub.ucName[3]);

	glDisableVertexAttribArray(vColorLoc);
}

/// UniformLocationのIndexをせっと
/// Shaderがロードされてから呼ばれる

void mgGLSLProgram::build_attribLocations(){
	memset(m_uniform_locations, -1,MAX_ATTRIB_LOC); 
	m_uniform_locations[modelViewProjMatrix]= getUniformLocation("modelViewProjMatrix");
	m_uniform_locations[modelViewMatrix]	= getUniformLocation("modelViewMatrix");
	m_uniform_locations[projMatrix]			= getUniformLocation("projMatrix");
	m_uniform_locations[normalMatrix]		= getUniformLocation("normalMatrix");

	m_uniform_locations[ndcMarix]			= getUniformLocation("ndcMarix");
	m_uniform_locations[ndcScaleMatrix]		= getUniformLocation("ndcScaleMatrix");
	m_uniform_locations[dpiFactor]			= getUniformLocation("dpiFactor");

	m_uniform_locations[DrawType]			= getUniformLocation("drawType");
	m_uniform_locations[FuncType]			= getUniformLocation("funcType");
	m_uniform_locations[ShaderMode]			= getUniformLocation("ShaderMode");

	m_uniform_locations[CoordinateType]		= getUniformLocation("coordinateType");
	
	m_uniform_locations[texture2D]			= getUniformLocation("texture2D");
	m_uniform_locations[pointSize]			= getUniformLocation("pointSize");
	m_uniform_locations[anchorPoint]		= getUniformLocation("anchorPoint");

	m_uniform_locations[LightTwoSides]	= getUniformLocation("LightTwoSides");
	m_uniform_locations[ForceLight]		= getUniformLocation("ForceLight");

	m_uniform_locations[ZebraAxis]		= getUniformLocation("ZebraAxis");
	m_uniform_locations[ZebraSize]		= getUniformLocation("ZebraSize");
}

void mgGLSLProgram::build_attribLocationsLight(){

	const size_t NAMESIZE = 80;
	char unifornName[NAMESIZE];
	for(int lightNo = 0; lightNo<LIGHT_NUM; ++lightNo){
		sprintf_s(unifornName, NAMESIZE, "Lights[%u].isEnabled", lightNo);
		m_LightPropsIds[lightNo][isEnabled]	= getUniformLocation(unifornName);

		sprintf_s(unifornName, NAMESIZE, "Lights[%u].ambientColor", lightNo);
		m_LightPropsIds[lightNo][ambientColor]	= getUniformLocation(unifornName);

		sprintf_s(unifornName, NAMESIZE, "Lights[%u].diffuseColor", lightNo);
		m_LightPropsIds[lightNo][diffuseColor]	= getUniformLocation(unifornName);

		sprintf_s(unifornName, NAMESIZE, "Lights[%u].specularColor", lightNo);
		m_LightPropsIds[lightNo][specularColor]	= getUniformLocation(unifornName);

		sprintf_s(unifornName, NAMESIZE, "Lights[%u].position", lightNo);
		m_LightPropsIds[lightNo][position]	= getUniformLocation(unifornName);

		sprintf_s(unifornName, NAMESIZE, "Lights[%u].spotDirection", lightNo);
		m_LightPropsIds[lightNo][spotDirection]	= getUniformLocation(unifornName);

		sprintf_s(unifornName, NAMESIZE, "Lights[%u].spotExponent", lightNo);
		m_LightPropsIds[lightNo][spotExponent]	= getUniformLocation(unifornName);

		sprintf_s(unifornName, NAMESIZE, "Lights[%u].spotCutoff", lightNo);
		m_LightPropsIds[lightNo][spotCutoff]	= getUniformLocation(unifornName);

		sprintf_s(unifornName, NAMESIZE, "Lights[%u].constantAttenuation", lightNo);
		m_LightPropsIds[lightNo][constantAttenuation]	= getUniformLocation(unifornName);

		sprintf_s(unifornName, NAMESIZE, "Lights[%u].linearAttenuation", lightNo);
		m_LightPropsIds[lightNo][linearAttenuation]	= getUniformLocation(unifornName);

		sprintf_s(unifornName, NAMESIZE, "Lights[%u].quadraticAttenuation", lightNo);
		m_LightPropsIds[lightNo][quadraticAttenuation]	= getUniformLocation(unifornName);

	}
}

void mgGLSLProgram::setUniform( mgGLSLProgram::UniformName name, float x, float y, float z){
	GLint loc = m_uniform_locations[name];
	setUniform(loc, x, y, z);
}

void mgGLSLProgram::setUniform( mgGLSLProgram::UniformName name, const glm::vec3 & v){
	GLint loc = m_uniform_locations[name];
	setUniform(loc, v);
}

void mgGLSLProgram::setUniform( mgGLSLProgram::UniformName name, const glm::vec4 & v){
	GLint loc = m_uniform_locations[name];
	setUniform(loc, v);
}

void mgGLSLProgram::setUniform( mgGLSLProgram::UniformName name, const glm::mat4 & m){
	GLint loc = m_uniform_locations[name];
	setUniform(loc, m);
}

void mgGLSLProgram::setUniform( mgGLSLProgram::UniformName name, const glm::mat3 & m){
	GLint loc = m_uniform_locations[name];
	setUniform(loc, m);
}

void mgGLSLProgram::setUniform( mgGLSLProgram::UniformName name, float val ){
	GLint loc = m_uniform_locations[name];
	setUniform(loc, val);
}

void mgGLSLProgram::setUniform( mgGLSLProgram::UniformName name, int val ){
	GLint loc = m_uniform_locations[name];
	setUniform(loc, val);
}

void mgGLSLProgram::setUniform( mgGLSLProgram::UniformName name, bool val ){
	GLint loc = m_uniform_locations[name];
	setUniform(loc, val);
}

//// LightProps
void mgGLSLProgram::setUniformLights( GLint lightNo, mgGLSLProgram::LightProps name, float x, float y, float z){
	ASSERT(lightNo>=0 && lightNo<LIGHT_NUM);
	GLint loc = m_LightPropsIds[lightNo][name];
	setUniform(loc, x, y, z);
}

void mgGLSLProgram::setUniformLights( GLint lightNo, mgGLSLProgram::LightProps name, const glm::vec3 & v){
	ASSERT(lightNo>=0 && lightNo<LIGHT_NUM);
	GLint loc = m_LightPropsIds[lightNo][name];
	setUniform(loc, v);
}

void mgGLSLProgram::setUniformLights( GLint lightNo, mgGLSLProgram::LightProps name, const glm::vec4 & v){
	ASSERT(lightNo>=0 && lightNo<LIGHT_NUM);
	GLint loc = m_LightPropsIds[lightNo][name];
	setUniform(loc, v);
}

void mgGLSLProgram::setUniformLights( GLint lightNo, mgGLSLProgram::LightProps name, const glm::mat4 & m){
	ASSERT(lightNo>=0 && lightNo<LIGHT_NUM);
	GLint loc = m_LightPropsIds[lightNo][name];
	setUniform(loc, m);
}

void mgGLSLProgram::setUniformLights( GLint lightNo, mgGLSLProgram::LightProps name, const glm::mat3 & m){
	ASSERT(lightNo>=0 && lightNo<LIGHT_NUM);
	GLint loc = m_LightPropsIds[lightNo][name];
	setUniform(loc, m);
}

void mgGLSLProgram::setUniformLights( GLint lightNo, mgGLSLProgram::LightProps name, float val ){
	ASSERT(lightNo>=0 && lightNo<LIGHT_NUM);
	GLint loc = m_LightPropsIds[lightNo][name];
	setUniform(loc, val);
}

void mgGLSLProgram::setUniformLights( GLint lightNo, mgGLSLProgram::LightProps name, int val ){
	ASSERT(lightNo>=0 && lightNo<LIGHT_NUM);
	GLint loc = m_LightPropsIds[lightNo][name];
	setUniform(loc, val);
}

void mgGLSLProgram::setUniformLights( GLint lightNo, mgGLSLProgram::LightProps name, bool val ){
	ASSERT(lightNo>=0 && lightNo<LIGHT_NUM);
	GLint loc = m_LightPropsIds[lightNo][name];
	setUniform(loc, val);
}

///Get version by glGetInteger(), and set the info
///in this program.
void mgGLSLProgram::setOpenGLVersion(){
    glGetIntegerv(GL_MAJOR_VERSION, &m_VersionMajor);
    glGetIntegerv(GL_MINOR_VERSION, &m_VersionMinor);
}
void mgGLSLProgram::getOpenGLVerion(GLint& major, GLint& minor)const{
    major=m_VersionMajor;
    minor=m_VersionMinor;
}

/// Debugまわり
void CALLBACK mgGLSL::debugCallback(GLenum source,
	GLenum type,
	GLuint id,
	GLenum severity,
	GLsizei length,
	const GLchar* message,
	void* userParam)
{
	std::ostream* stream = static_cast<std::ostream*>(userParam);
	if(stream!=NULL){
		std::ostream& out = *stream;

		out<<"Severity:";
		if(severity == GL_DEBUG_SEVERITY_HIGH_ARB)
			out<<"High";
		else if(severity == GL_DEBUG_SEVERITY_MEDIUM_ARB)
			out<<"Medium";
		else if(severity == GL_DEBUG_SEVERITY_LOW_ARB)
			out<<"Low";

		out<<"\tType (DEBUG_TYPE):";

		if(type == GL_DEBUG_TYPE_ERROR_ARB)
			out<<"ERROR";
		else if(type == GL_DEBUG_TYPE_DEPRECATED_BEHAVIOR_ARB)
			out <<"DEPRECATED_BEHAVIOR";
		else if(type == GL_DEBUG_TYPE_UNDEFINED_BEHAVIOR_ARB)
			out <<"UNDEFINED_BEHAVIOR";
		else if(type == GL_DEBUG_TYPE_PORTABILITY_ARB)
			out <<"PORTABILITY ";
		else if(type == GL_DEBUG_TYPE_PERFORMANCE_ARB)
			out <<"PERFORMANCE";
		else if(type == GL_DEBUG_TYPE_OTHER_ARB)
			out <<"OTHER";

		out<< "\tSource:";
		if(source == GL_DEBUG_SOURCE_API_ARB)
			out<<"GL_DEBUG_SOURCE_API_ARB";
		else if(source == GL_DEBUG_SOURCE_WINDOW_SYSTEM_ARB)
			out<<"GL_DEBUG_SOURCE_WINDOW_SYSTEM_ARB";
		else if(source == GL_DEBUG_SOURCE_SHADER_COMPILER_ARB)
			out<<"GL_DEBUG_SOURCE_SHADER_COMPILER_ARB";
		else if(source == GL_DEBUG_SOURCE_THIRD_PARTY_ARB)
			out<<"GL_DEBUG_SOURCE_THIRD_PARTY_ARB";
		else if(source == GL_DEBUG_SOURCE_APPLICATION_ARB)
			out<<"GL_DEBUG_SOURCE_APPLICATION_ARB";
		else if(source == GL_DEBUG_SOURCE_OTHER_ARB)
			out<<"GL_DEBUG_SOURCE_OTHER_ARB";
		else
			out<< source;

		out<<std::endl
			<< "Message: "<< message << std::endl;

	}
}

// ShaderモードのOn/Offを設定する
void mgGLSLProgram::EnableLights(bool bEnabled){
	mgGLSL::ShadeMode mode = (bEnabled)?mgGLSL::Shading:mgGLSL::NoShading;
	setUniform(ShaderMode, mode);
}

// ShaderモードのOn/Offを取得する
bool mgGLSLProgram::LightEnabled(){
	ASSERT(isLinked());
	GLint loc = m_uniform_locations[mgGLSLProgram::ShaderMode];
	GLint lightEnabled;

	glGetUniformiv(m_handle, loc, &lightEnabled);

	return (lightEnabled == mgGLSL::Shading) ? true : false;
}

void mgGLSL::execLightMode(
	int mode/// <0: undefined, =0:Light is disabled, >0:Light is enabled.
){
	if(mode<0)
		return;

	ASSERT(mgGLSLProgram::getCurrentGLSLProgram()!=NULL);
	mgGLSLProgram* pGLSL = mgGLSLProgram::getCurrentGLSLProgram();
	ASSERT(pGLSL->isLinked());

	mgGLSL::ShadeMode shadeMode = (mode) ? mgGLSL::Shading:mgGLSL::NoShading;
	pGLSL->setUniform(mgGLSLProgram::ShaderMode, shadeMode);

	mgStaticGLAttrib& cattr=mgGLSL::getCurrentStaticGLAttrib();
	cattr.setLightMode(mode);
}

bool mgGLSL::LightEnabled(){
	ASSERT(mgGLSLProgram::getCurrentGLSLProgram()!=NULL);
	mgGLSLProgram* pGLSL = mgGLSLProgram::getCurrentGLSLProgram();
	ASSERT(pGLSL->isLinked());

	return pGLSL->LightEnabled();
}

void DefaultLightsSetup(){

	mgGLSLProgram* pGLSL = mgGLSLProgram::getCurrentGLSLProgram();

	// 両面描画
	pGLSL->setUniform(mgGLSLProgram::LightTwoSides, true);

	// ForceLightモードをFALSEに。有効なライトがない場合は、
	// Light-OFF状態で描画するように変更。Fugenの場合はTRUEとするところ。
	pGLSL->setUniform(mgGLSLProgram::ForceLight, true);
}

void mgGLSL::execStaticColorAttrib(
	const float colr[4]
){
	int vColorLoc=mgGLSLProgram::vColor;
	glDisableVertexAttribArray(vColorLoc);	
	glErr=glGetError();

	glVertexAttribPointer(vColorLoc,4,GL_FLOAT,GL_FALSE,0,0);
	glVertexAttrib4fv(vColorLoc,colr);
	mgStaticGLAttrib& cattr=mgGLSL::getCurrentStaticGLAttrib();
	cattr.setColor(colr);
}

void mgGLSL::execStaticColorAttrib(
	const MGColor& color
){
	if(color.defined()){
		mgGLSL::execStaticColorAttrib(color.color());
	}
}

void mgGLSL::execStaticLineWidth(float lineWidth){
	if(lineWidth<=0.)
		return;

	glLineWidth(lineWidth);
	mgStaticGLAttrib& cattr=mgGLSL::getCurrentStaticGLAttrib();
	cattr.setLineWidth(lineWidth);
}

///Line stipple属性をセットする。
///When factor=0 is input, line pattern is disabled. This means the line is solid.
///When factor<0, the stipple attribute is undefined. This means the attribute
///is defined by the environment.
///When factor<=0, pattern is unnecessary.
void mgGLSL::execStaticLineStipple(short int factor, GLuint pattern){
	if(factor<0)
		return;

	if(factor>0){
		glEnable(GL_LINE_STIPPLE);
		glLineStipple((GLint)factor,pattern);
	}else
		glDisable(GL_LINE_STIPPLE);
	mgStaticGLAttrib& cattr=mgGLSL::getCurrentStaticGLAttrib();
	cattr.setLineStipple(factor,pattern);
}

void mgGLSL::initializeStaticGLAttribStack(){
	std::stack<mgStaticGLAttrib>& attribs=mgGLSL::getStaticGLAttribStack();
	if(attribs.empty()){
		attribs.push(mgStaticGLAttrib());
	}
}
std::stack<mgStaticGLAttrib>& mgGLSL::getStaticGLAttribStack(){
	static std::stack<mgStaticGLAttrib> m_attribStack;
	return m_attribStack;
}

mgStaticGLAttrib& mgGLSL::getCurrentStaticGLAttrib(){
	std::stack<mgStaticGLAttrib>& attribs=mgGLSL::getStaticGLAttribStack();
	return attribs.top();
}
void mgGLSL::pushStaticGLAttrib(){
	mgStaticGLAttrib& cattriv=mgGLSL::getCurrentStaticGLAttrib();
	std::stack<mgStaticGLAttrib>& attribs=mgGLSL::getStaticGLAttribStack();
	attribs.push(cattriv);
}
void mgGLSL::popStaticGLAttrib(){
	std::stack<mgStaticGLAttrib>& attribs=mgGLSL::getStaticGLAttribStack();
	if(attribs.size()>1)
		attribs.pop();
	mgStaticGLAttrib& attrib=attribs.top();
	mgGLSL::execStaticGLAttrib(attrib);
}
void mgGLSL::execStaticGLAttrib(const mgStaticGLAttrib& attrib){
	mgGLSL::execStaticColorAttrib(attrib.color());
	mgGLSL::execStaticLineWidth(attrib.getLineWidth());
	short int factor; GLushort pattern;
	attrib.getLineStipple(factor, pattern);
	mgGLSL::execStaticLineStipple(factor, pattern);
	mgGLSL::execLightMode(attrib.getLightMode());
}

void mgGLSLProgram:: setFuncType(mgGLSL::FuncType type){
	GLint loc = m_uniform_locations[FuncType];
	setUniform(loc, type);
}

int mgGLSLProgram:: getFuncType(){
	GLint loc = m_uniform_locations[FuncType];
	GLint funcType;
	glGetUniformiv(m_handle, loc, &funcType);
	return funcType;
}

void mgGLSLProgram::setCoordinateType(mgGLSL::CoordinateType type){
	GLint loc = m_uniform_locations[CoordinateType];
	setUniform(loc, type);
}

mgGLSL::CoordinateType mgGLSLProgram::getCoordinateType(){
	GLint loc = m_uniform_locations[CoordinateType];
	GLint coordinateType;
	glGetUniformiv(m_handle, loc, &coordinateType);
	return (mgGLSL::CoordinateType)coordinateType;
}

void mgGLSLProgram::getAnchor(MGPosition& anchor){
	GLint loc = getUniformLocation("funcType");
	float xyz[3];
	glGetUniformfv(m_handle, loc, xyz);
	anchor=MGPosition(xyz[0],xyz[1],xyz[2]);
}

void printMatrix(glm::mat4& mat){
	for(int i = 0; i<4;++i){
	COUT<<"|";
		for(int ii = 0; ii<4;++ii){
			COUT <<mat[i][ii] <<" " ;
		}
	COUT<<std::endl;
	}
}

void printMatrix3(glm::mat3& mat){
	for(int i = 0; i<3;++i){
	COUT<<"|";
		for(int ii = 0; ii<3;++ii){
			COUT <<mat[i][ii] <<" " ;
		}
	COUT<<std::endl;
	}
}

void mgGLSLProgram::printPrjMatrix(){

	glm::mat4 mvp;
	glGetUniformfv( m_handle, m_uniform_locations[mgGLSLProgram::modelViewProjMatrix], &mvp[0][0]);

	COUT<<"modelViewProjMatrix"<<":"<<std::endl;
	printMatrix(mvp);

	glm::mat4 mv;
	glGetUniformfv( m_handle, m_uniform_locations[mgGLSLProgram::modelViewMatrix], &mv[0][0]);

	COUT<<"modelViewMatrix"<<":"<<std::endl;
	printMatrix(mv);

	glm::mat4 prj;
	glGetUniformfv( m_handle, m_uniform_locations[mgGLSLProgram::projMatrix], &prj[0][0]);

	COUT<<"projMatrix"<<":"<<std::endl;
	printMatrix(prj);

	glm::mat3 mn;
	glGetUniformfv( m_handle, m_uniform_locations[mgGLSLProgram::normalMatrix], &mn[0][0]);

	COUT<<"normalMatrix"<<":"<<std::endl;
	printMatrix3(mn);

}

void mgGLSLProgram::printLightProps()
{
	COUT<<"--Status-------"<<std::endl;
	GLint status = -1;
	glGetUniformiv( m_handle, m_uniform_locations[mgGLSLProgram::CoordinateType], &status);

	GLint drawType = -1;
	glGetUniformiv( m_handle, m_uniform_locations[mgGLSLProgram::DrawType], &drawType);

	COUT<<"DrawType:"<<drawType<<"  CoordType:"<<status<< std::endl;
	glm::vec3 pos;
	glGetUniformfv( m_handle, m_uniform_locations[mgGLSLProgram::anchorPoint], &pos[0]);
	COUT<<"AnchoPos:"<<pos.x <<", "<< pos.y<<", "<<pos.z<<std::endl;
}


///Utility class to invoke glsl's setFuncType.
///mgFuncTypeSwitcher saves the current function type and invoke input function type.
///When mgFuncTypeSwitcher is destructed, the saved original type is restored.
mgFuncTypeSwitcher::mgFuncTypeSwitcher(mgGLSL::FuncType type){
	m_pGLSL=mgGLSLProgram::getCurrentGLSLProgram();
	m_orgFuncType = m_pGLSL->getFuncType();
	m_pGLSL->setFuncType((mgGLSL::FuncType)type);
}

mgFuncTypeSwitcher::~mgFuncTypeSwitcher(){
	m_pGLSL->setFuncType((mgGLSL::FuncType)m_orgFuncType);
}

///mgCoordinateTypeSwitcher saves the current coordinate type and invoke input coordinate type.
///When mgCoordinateTypeSwitcher is destructed, the saved original type is restored.
mgCoordinateTypeSwitcher::mgCoordinateTypeSwitcher(
	mgGLSL::CoordinateType coordinateType, const MGPosition* anchorP
){
	m_pGLSL = mgGLSLProgram::getCurrentGLSLProgram();
	m_coordinateType = m_pGLSL->getCoordinateType();//save the old one.
	if(m_coordinateType!=coordinateType)
		m_pGLSL->setUniform(mgGLSLProgram::CoordinateType, coordinateType);//set new one.
	if(coordinateType==mgGLSL::AnchorPoint ||coordinateType==mgGLSL::AnchorPointScreen){
		//When new one is to update the anchor point.
		if(m_coordinateType==mgGLSL::AnchorPoint ||m_coordinateType==mgGLSL::AnchorPointScreen)			
			m_pGLSL->getAnchor(m_anchorPoint);//save only the old one that is valid.

		glm::vec3 anchor = glm::vec3((*anchorP)[0], (*anchorP)[1], (*anchorP)[2]);
		m_pGLSL->setUniform(mgGLSLProgram::anchorPoint, anchor);//set new one.
	}
}

mgCoordinateTypeSwitcher::~mgCoordinateTypeSwitcher(){
	mgGLSL::CoordinateType coordinateType = m_pGLSL->getCoordinateType();//save the current one.
	if((coordinateType==mgGLSL::AnchorPoint || coordinateType==mgGLSL::AnchorPointScreen)
		&& (m_coordinateType==mgGLSL::AnchorPoint ||m_coordinateType==mgGLSL::AnchorPointScreen)
		){//When new one is to update the anchor point.
		glm::vec3 anchor = glm::vec3(m_anchorPoint[0], m_anchorPoint[1], m_anchorPoint[2]);
		m_pGLSL->setUniform(mgGLSLProgram::anchorPoint, anchor);//restore the old one.
	}
	if(coordinateType!=m_coordinateType)
		m_pGLSL->setCoordinateType(m_coordinateType);//restore the old coordinate type.
}
