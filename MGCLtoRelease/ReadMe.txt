MGCL�@Version12 Release�t�@�C��

1. �r���h�ɂ���
�EOpenGL��V4.2�ȏオ�ғ�����K�v������܂��B
�E���p���Ă���C++�̃o�[�W������C++20�ł��B

�Q�DMGCL V12�̎�ȓ���
(1) C++ Move Semantics�Ή�
MGVLV12�ł̓��[�u�Z�}���e�B�N�X�𗘗p���Ă��܂��BIGES�֘A�v���O�����ł͂܂����[�u�Z�}���e�B�N�X�Ή���
���������K�v�ł����A���ł́A����B-Spline�֘A�N���X�iMGNDDArray, MGKnotVector, MGBPointSeq, MGSPointSeq�Ȃǁj�́A
���̋@�\�𗘗p����悤���ǂ���Ă��܂��B

(2) unique_ptr�̗��p
MGCLV12�ł͐ϋɓI��unique_ptr�𗘗p���Ă��܂��B���̂��߂�MGCLDefs.h�ŁA������UniqueXxxxx���`���Ă��܂��B
�����𗘗p����unique_ptr��vector, list�Ȃǂ𗘗p���Ă��܂��B
�]���g�p���Ă���auto_ptr�̗��p�͂��ׂ�unique_ptr�̗��p�ɕύX���Ă��܂��B

(3) Tessellation�v���O��������
Tessellation�v���O������啝�ɉ��ǂ��Ă��܂��B���܂ł��܂������Ȃ����������̗���ƂȂ��Ă��܂��B

(4) B-Spline�N���X MGLBRep, MGSBRep��ctor����
���܂ŃR���X�g���N�^�[(ctor)�Ƃ��Ē񋟂��Ă���B-Spline�̐����iCurve, Surface�Ƃ��Ɂj��
buildxxxxx�Ƃ��ă����o�[�֐������Ă킩��₷�������B

�R�D�g�p���Ă���O���\�t�g
OpenGL V4.2�̂��߁Aglm, glew, freetype, ftgl �𗘗p���Ă��܂��B
MGCL�̗��p��MIT���C�Z���X�ŁA���p�ɍۂ��Đ��񂪂���܂��񂪁A���p�ŗ��p����ꍇ��L�̃��C�Z���X�󋵂��悭���ׂė��p���Ă��������B

freetype, glew, glm�͑S����������邱�ƂȂ��I���W�i���Ȃ��̂����̂܂ܗ��p���Ă��܂����A
ftgl�ɂ�MGCL��VBO�N���X�𗘗p�ł���悤��������Ă��܂��B
freetype, glew��dll, .lib��MGCL��download�œ���ł��܂����A�\�[�X�v���O�����͊e�T�C�g������肵�Ă��������B
MGCLtoRelease�ɂ�DLL�Ƃ��̂��߂�lib�t�@�C���A����уR���p�C���ɕK�v��include�t�@�C���݂̂�����Ă���܂��B
ftgl��MGCLV12.dll�Ƃ���MGCL�̈ꕔ�ɑg�ݍ��܂�Ă��܂��B

(1) glm
glm�̓R���p�C�����ɂ̂ݎg�p���܂��Binclude�t�@�C�������𗘗p���Acpp�t�@�C���i���s�t�@�C���j�͂���܂���B
Link�̂��߂�lib�t�@�C������ю��s�̂��߂�dll�͂���܂���B

(2) glew
"The OpenGL Extension Wrangler Library (GLEW) is a cross-platform open-source C/C++ extension loading library"
�Ƃ���悤��glew�͎�X��OpenGL���W���[���̃A�h���X���������Ă���܂��B
OpenGL�̃��W���[�����p�̋��n�������Ă���ALink�̂��߂�lib�t�@�C���Ǝ��s�̂��߂�dll��K�v�Ƃ��܂��B
.lib�t�@�C����dll�������Ă��܂��B

(3) freetype
freetype��ftgl�����p���Ă��邽�ߕK�v�ŁAMGCL�͒��ڗ��p���Ă��܂���B
lib�t�@�C�����񋟂���܂��Bdll�͂���܂���B

(4) ftgl
ftgl�͕����`��̂��߂̃v���O�����ł��B
ftgl�͊e��t�H���g�𗘗p�ł��܂����AMGCL�́Aftgl�̃t�H���g���ׂĂ��T�|�[�g���Ă���킯�ł͂Ȃ��A
�ꕔ�����T�|�[�g���Ă��܂���B����̉ۑ�ł��B

�S�DOpenGL��shader�v���O����
MGCL�́A�ꕔWindows�Ɉˑ����Ă��܂��B
�摜�����i.png�Ȃǂ̉摜�t�@�C���̓ǂݍ��݁j�AOpenGL��Windows�pcontext��Windows�Ή��ƂȂ��Ă��܂��B
OpenGL�̕\���֘A�𗘗p����ꍇ�AShader�v���O�������K�v�ƂȂ�܂����A������<MGCLtoRelease>\bin�Ɋ܂܂�Ă��܂��B
���s���ɕK�v�ƂȂ�܂��̂ŁA���p�v���O�����̃p�X���ʂ�Ƃ���ɕ��ʂ���K�v������܂��B
