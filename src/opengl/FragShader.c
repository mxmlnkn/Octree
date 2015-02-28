#version 330
//Version specified by LoadShader !!!
flat in vec4 f_color;
void main(void) {
	gl_FragColor = f_color;
	//gl_FragColor = vec4(z,z,z,fade);
	//gl_FragColor[0] = gl_FragCoord.x/640.0;		//red
	//gl_FragColor[1] = gl_FragCoord.y/480.0;		//green
	//gl_FragColor[2] = 0.5;						//blue
	//gl_FragColor[3] = mod(gl_FragCoord.y/20, 1.0);	//alpha
}