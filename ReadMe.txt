MGCL20�\�����[�V����

�{�\�����[�V�����́A���ʂƂ���MGCLtoRelease�t�H���_�[�z���� include, lib, dll��
�e�t�H���_�[�z���̃t�@�C�����X�V����B�r���h�ł�lib, bin�t�H���_�[�̃t�@�C���𒼐ڍX�V���A
�r���h��̃C�x���g�Ƃ��āAMGCL�̊J���ŕύX����Ă���\���̂���t�@�C�� (include�t�@�C�����ׂ�)��
MGCLV12/include����MGCLtoRelease/include�փR�s�[����(�������AGL(glew�pinclude�t�@�C��)
�����glm��include�t�@�C���̓R�s�[���Ȃ�)�B

MGCL��DLL��MGCLtoRease�t�H���_�[�S�̂��R�s�[���ė��p����BMGCLtoRelease�z���ɂ̓R���p�C���ɕK�v��
include�t�@�C���A�����N�ɕK�v��lib�t�@�C���A����ю��s���ɕK�v�Ȏ��s�t�@�C��
�iMGCL��glew��dll, OpenGL Shader�v���O����)���܂܂�Ă��邽�߁A���ꂾ���ŕK�v�\���ł���B

MGCLtoReleas�͉��L�̍\���ƂȂ�F

MGCLtoRelease
 |
 +- bin(MGCLDLL�AOpenGL�֘A��DLL�AOpenGL Shader�v���O����)
 |
 +- lib(MGCL��OpenGL�֘A��DLL�p��lib�Q)
 |
 +- include(MGCL��OpenGL��include�t�@�C���Q)

�P�jinclude�t�@�C��(MGCLtoRelease/include)
include�t�H���_�[�z���Ɏ��̃T�u�t�H���_�[������F
 GL(glew�j, glm, mg(MGCL geometry and other basic classes), 
topo(MGCL topology classes), mgGL(MGCL Graphic classes), Tl2(MGCL tesellation classes), 
mgiges(MGCL iges classes)

2) lib�t�@�C��(MGCLtoRelease/lib)
DLL�Ƃ��Ē񋟂����Aglew��MGCL�̃����N�p�t�@�C��
glew32.lib�@��MGCLVxx.lib(xx��MGCL�̃o�[�W�����ԍ�)

3) dll�t�@�C����OpenGL Shader�v���O����(MGCLtoRelease/bin)
�@1) MGCLVxx.dll, 2) glew32.dll, 3) mgclShade.frag, 4) mgclShader.vert
�@
