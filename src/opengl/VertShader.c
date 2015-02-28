#version 330
//Version specified by LoadShader !!!
uniform mat4 mvp;
uniform vec4 fixed_color;	//if alpha is not 0 this will be used instead of v_color
in vec3 coord3d;
in vec3 v_color;
flat out vec4 f_color;
void main(void) {
	gl_Position = mvp * vec4(coord3d, 1.0);
	if (fixed_color.w == 0)
		f_color = vec4(v_color,1.0);
	else
		f_color = fixed_color;
}
