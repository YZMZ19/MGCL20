MGCL�𗘗p�����J���ɂ�include�t�@�C��(.h�t�@�C���j�A�����N�̂��߂�.lib�t�@�C���A
����ю��s�t�@�C��(DLL��OpenGL Shader�v���O����)���K�v�Ƃ����B

�P�jinclude�t�@�C��
include�t�H���_�[�z���Ɏ��̃T�u�t�H���_�[������F
 GL(glew�j, glm, mg(MGCL geometry and other basic classes), 
topo(MGCL topology classes), mgGL(MGCL Graphic classes), Tl2(MGCL tesellation classes), 
mgiges(MGCL iges classes)

2) lib�t�@�C��
DLL�Ƃ��Ē񋟂����Aglew��MGCL�̃����N�p�t�@�C��
glew32.lib�@��MGCLVxx.lib(xx��MGCL�̃o�[�W�����ԍ�)

3) dll�t�@�C����OpenGL Shader�v���O����
�@1) MGCLVxx.dll, 2) glew32.dll, 3) mgclShade.frag, 4) mgclShader.vert

MGCL20�̃v���W�F�N�g�ł̓r���h��̃C�x���g�Ƃ��āAMGCL�̊J���ŕύX����Ă���\���̂���t�@�C����
fugen�\�����[�V������MGCLtoRelease�փR�s�[���邽�߂̃o�b�`�R�}���h�𑖂点�Ă���F

�i�P�j�쐬���ꂽdll�Ƃ��̗��p�̂��߂�lib�t�@�C���̃R�s�[
 $(FUGEN_MGCL_DIR)bin\�@�Ɓ@$(FUGEN_MGCL_DIR)lib\�@��
�i�Q�jMGCL��dll�ŗ��p���邽�߂�include�t�@�C��
mg, mgGL, tl2, topo, mgiges�̊e�t�H���_�[��$(FUGEN_MGCL_DIR)include��

�����ŁAFUGEN_MGCL_DIR�� MGCLSolution.props�Œ�`����Ă���}�N���ŁA
���L�̃t�H���_�[�\����O��Ƃ��Ă���F

MGCLtoRelease :fugen�v���W�F�N�g�ŗ��p����MGCL��folder.
 |
 +- bin(MGCLDLL��OpenGL�֘A��DLL)
 |
 +- lib(MGCL��OpenGL�֘A��DLL�p��lib�Q)
 |
 +- include(MGCL��OpenGL��include�t�@�C���Q)

