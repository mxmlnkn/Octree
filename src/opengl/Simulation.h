#ifndef SIMULATION_H
#define SIMULATION_H

#define OBJ_FILE_NAME "cityandwindmill.obj"
GLuint program;
//attributes: Variables for our shaders
GLint attribute_coord3d, attribute_v_color,uniform_mvp,uniform_fixed_color;
//Virtual Buffer Object IDs
GLuint vbo_vertices, ibo_triangles, vbo_triangle_colors, vbo_normal_pts;
//Own storage arrays for Objects
float *projection_matrix=0, *view_matrix=0, *model_matrix=0, *mvp=0;
float obj_rot=0, x_transl=0, y_transl=0, z_transl=0, w_transl=1, obj_scale=1;

struct {	//Initialized by ChangeSize
	unsigned int Width,Height;
} Screen;

typedef float vec3[3];
typedef float vec4[4];
//asdf translates all vectors. arrow keys only translate Center(left,right) and Up (up,down) NOTE: Up is a direction and therefore mustn't be 0!!!
typedef struct {
	float Eye[3];					// Position of the Camera (initialize: {0,0,0})
	float Center[3];				// Point where we look at (initialize: {0,0,-1})
	float Up[3];					// specifiying what up is (relative!!!) to be able to tile the camera along the viewing-axis (KEY_LEFT, KEY_RIGHT) (initialize: {0,1,0})
	float Phi,Theta,Psi;		// spherial coordinates with theta ranging from -pi/2 (-z-axis) to +pi/2 (+z-axis) (initialize: theta=-M_PI, PSI=0)
	// the initial values let the camera look to the -z-axis
	// => Center = Eye + [Cos Phi Sin Theta, Sin Phi Sin Theta, Cos Theta]
	// Psi is the angle from the y-axis and the vector specified by Theta and Phi ... It's to tilt the Camera.Up
	// => Camera.Up = RotateMatrixAround(Camera.Center-Camera.Eye, Degree)*(0,1,0)
} CameraData;
CameraData Camera = {{0,1.7,0},{0,0,-3},{0,2.7,0},0,-0.51555,0};	//Eye at 1.70m height looks at a point on the ground 3m distant (Theta =atan(-1.7/3))

typedef struct {
	char Name[255];
	GLuint* Triangles;		//a,b,c,a,b,c,...(counterclockwise face)
	GLuint TriangleCount;
	GLuint QuadCount;		//can be neglected, because every quad also increases TriangleCount by 2
	GLfloat* NormalLines;	//Face-1-Center, Center+Normal; Face-2-Center, ...
	GLfloat* Normals;		//1st Normal, ...
} Object3D;
typedef struct {
	GLfloat* Vertices;	//x,y,z;x,y,z;...
	GLuint VertexCount;
	Object3D* Objects;	//beliebig langer Array mitObject3D-Strukturen (die konstanter Größe sind), also Object3D*
	GLuint ObjectCount;
} OBJData;

OBJData Scene;

#define SetMatrix4fToIdentity(m) axpby(1, NULL, 0, m, 4,4);
#define MovVec3f(a,b) memcpy(a,b,3*sizeof(float));
#define ScalarProduct3f(a,b,res) MultiplicateMat(&(a)[0],&(b)[0],1,3,3,1,res);
void MultiplicateMat(const float* a, const float* b, unsigned int a_rows, unsigned int a_cols, unsigned int b_rows, unsigned int b_cols, float** c);
#define ScaleMat(lambda,m,r,c) axpby(0, NULL, lambda, m, r,c);
void AddMat(float* a, float* b, unsigned int rows, unsigned int cols, float** c);
void axpby(const float sc_a, const float* x, const float sc_b, float* y, const unsigned int rows, const unsigned int cols);
float VectorNorm(const float* vec, const unsigned int dim);
void CrossProduct(float* a, float* b, float* c);
void DebugDisplayMat(const float* mat, const unsigned int rows, const unsigned int cols);
void CreateRotationMat3f(const float ux, const float uy, const float uz, const float theta, float** m);
void ProcessAsciiKeys(unsigned char key, int x, int y);
void ProcessSpecialKeys(int key, int x, int y);
OBJData LoadOBJFile(const char* filename);
char* LoadSourceFile(const char* filename);
GLuint LoadShader(const char* filename, GLenum type);
GLuint GetAttributeID(const char* attribute_name);
GLuint GetUniformID(const char* attribute_name);
int init_resources(void);
void ChangeSize(int w, int h);
void CalcView(const CameraData* Camera, float* mvp);
void CalcProjection(float fovy, float aspect, float zNear, float zFar, float* view_m);
void RenderScene(void);
void DoOnIdle(void);
void free_resources(void);

#endif