#define DEBUG_SIM
//When compiling freeglut and glew use msys not cygwin!
#include <stdio.h>

#include <time.h>
#include <math.h>
#include <omp.h>

#include <windows.h>        // https://stackoverflow.com/questions/20441581/getting-errors-when-compiling-glew-sdl-program
#include <GL/glew.h>		// Loads OpenGL2.0+ functionality for us
#include <GL/glut.h>		// Windows FreeGlut equivalent
#include <GL/freeglut_ext.h>

#include "Simulation.h"

/********************** TODO ***************
 * Diffuse Lambert Lighting                *
 * Actual Simulation                       *
 * Collision-Recognization:                *
 *  -for Sphere and Cylinder               *
 *  -General: two triangle faces (easy)    *
 * Jumping (use Cylindrical human)         *
 * Normal Mapping                          *
 * Clean Code                              *
 *******************************************/

/******************************************************************
 *************************** Definitions **************************
 ******************************************************************/

void MultiplicateMat(const float* a, const float* b, unsigned int a_rows, unsigned int a_cols, unsigned int b_rows, unsigned int b_cols, float** c) {
	if(a_cols != b_rows) {
		*c = 0;
		return;
	}
	*c = realloc(*c, sizeof(float)*a_rows*b_cols);
	memset(*c,0,sizeof(float)*a_rows*b_cols);
	unsigned int j;
	#pragma omp parallel for	//parallelizes only the j loop ... therefore won't work if result is row-vector
	for (j=0; j<b_cols; j++) {
		unsigned int i;
		for (i=0; i<a_rows; i++) {
			unsigned int k;
			for (k=0; k<a_cols; k++) {
				//fprintf(stderr,"c[%u,%u] = %f + %f \n",i,j,(*c)[b_cols*i+j],a[i*a_cols+k]*b[k*b_cols+j]);
				(*c)[b_cols*i+j] += a[i*a_cols+k]*b[k*b_cols+j];
			}
		}
	}
	return;
}
/* With the following you can do: a=b or a=a+b or a=a-b or a=l*a and many more which would be extra functions normally ...
 * You can even use this for b=I with axpby(1,NULL,0,&mat,...), because axpby will use Identity if nothing is specified for the first matrix */
void axpby(const float sc_a, const float* x, const float sc_b, float* y, const unsigned int rows, const unsigned int cols) {
	if (y == 0) return;
	unsigned int j;//axpby(1,Camera.Center,-1,(float**)&f,3,1);
	#pragma omp parallel for	//parallelizes only the j loop ... therefore won't work if result is row-vector
	for (j=0; j<cols; j++) {
		unsigned int i;
		for (i=0; i<rows; i++) {
			//fprintf(stderr, "(*y)[%u,%u] = ... (y=%p, *y=%p)\n",i,j,y,*y);
			float X;
			if (x==NULL) {
				if(i==j)
					X=1;
				else
					X=0;
			} else
				X = x[i*cols+j];	//if no matrix specified for first one, then use identity matrix
			y[i*cols+j] = sc_a *X + sc_b * y[cols*i+j];
		}
	}
	return;
}
float VectorNorm(const float* vec, const unsigned int dim) {
	int i; float sum=0;
	for (i=0;i<dim;i++) {
		//fprintf(stderr,"VectorNorm[%u]=%f, ",i,vec[i]);
		sum += vec[i]*vec[i];
	}
	return sqrt(sum);
}
void CrossProduct(float* a, float*b, float* c) {
	c[0] = a[1]*b[2]-a[2]*b[1];
	c[1] = a[2]*b[0]-a[0]*b[2];
	c[2] = a[0]*b[1]-a[1]*b[0];
	return;
}
void DebugDisplayMat(const float* mat, const unsigned int rows, const unsigned int cols) {
	unsigned int i,j;
	fprintf(stderr, "%ux%u-Matrix: {",rows,cols);
	for (i=0; i<rows; i++) {
		fprintf(stderr,"{");
		for (j=0; j<cols; j++) {
			fprintf(stderr, "%f",mat[i*cols+j]);
			if(j!=cols-1)
				fprintf(stderr,",");
		}
		fprintf(stderr, "}");
		if(i!=rows-1)
			fprintf(stderr,",");
	}
	fprintf(stderr,"}\n");
	return;
}
void CreateRotationMat3f(const float ux, const float uy, const float uz, const float theta, float** m) {
	*m = realloc(*m,9*sizeof(float));
	(*m)[0*3+0]=0;		(*m)[0*3+1]=-uz;	(*m)[0*3+2]=uy;
	(*m)[1*3+0]=uz;		(*m)[1*3+1]=0;		(*m)[1*3+2]=-ux;
	(*m)[2*3+0]=-uy;	(*m)[2*3+1]=ux;		(*m)[2*3+2]=0;
	//m=cos theta*I + sin theta * m
	axpby(cos(theta*M_PI/180.0),NULL, sin(theta*M_PI/180.0) ,*m,3,3);
	//+ uu^T*(1-cos theta)
	float* u = malloc(3*sizeof(float));
	u[0]=ux; u[1]=uy; u[2]=uz; //needed as an array for multiplicatemat
	float* m1=0; MultiplicateMat(u,u,3,1,1,3,&m1);
	axpby(1-cos(theta*M_PI/180.0),m1,1,*m,3,3);
	free(u); free(m1);
	return;
}
void CreateRotationMat4f(const float ux, const float uy, const float uz, const float theta, float** m) {
	#define ROT_DIM 4
	#if ROT_DIM >= 3
		*m = realloc(*m,ROT_DIM*ROT_DIM*sizeof(float));
		memset(*m,0,sizeof(float)*ROT_DIM*ROT_DIM);
		(*m)[0*ROT_DIM+0]=0;		(*m)[0*ROT_DIM+1]=-uz;		(*m)[0*ROT_DIM+2]=uy;
		(*m)[1*ROT_DIM+0]=uz;		(*m)[1*ROT_DIM+1]=0;		(*m)[1*ROT_DIM+2]=-ux;
		(*m)[2*ROT_DIM+0]=-uy;		(*m)[2*ROT_DIM+1]=ux;		(*m)[2*ROT_DIM+2]=0;
		//m=cos theta*I + sin theta * m
		axpby(cos(theta*M_PI/180.0),NULL, sin(theta*M_PI/180.0) ,*m,ROT_DIM,ROT_DIM);
		//+ uu^T*(1-cos theta)
		float* u = malloc(ROT_DIM*sizeof(float));
		memset(u,0,ROT_DIM*sizeof(float));
		u[0]=ux; u[1]=uy; u[2]=uz; //needed as an array for MultiplicateMat
		float* m1=0;
		MultiplicateMat(u,u,ROT_DIM,1,1,ROT_DIM,&m1);
		axpby(1-cos(theta*M_PI/180.0),m1,1,*m,ROT_DIM,ROT_DIM);
		free(u); free(m1);
		//if ROT_DIM>3 then the matrix should contain a 3x3 matrix in the upper left corner and be zero in the other elements
		#if ROT_DIM == 4
			(*m)[(ROT_DIM-1)*ROT_DIM+(ROT_DIM-1)]=1;
		#endif
	#endif
	#undef ROT_DIM
	return;
}

int CollisionSphereObject(const float* Eye, const float Radius, const Object3D* Object, const float* Vertices) {
	unsigned short int Collides = 0;
	//Go through every Face in Object
	unsigned int triangle_num;
	#pragma omp parallel for
	for (triangle_num=0; triangle_num<Object->TriangleCount; triangle_num++) {
		GLuint* cur_tri = &(Object->Triangles[3*triangle_num]);
		//Initialize some vector: AB,AC,AEye,Normal
		float AB[3],AC[3],BC[3],AEye[3],BEye[3],Normal[3];			// float* Normal = &Object->Normals[3*triangle_num];	//could also be used
		// OA->AB,AC,AEye and OB->BEye,BC
			MovVec3f(&AB[0], &Vertices[3*cur_tri[0]]);
			MovVec3f(&AC[0],&AB[0]);
			MovVec3f(&AEye[0],&AB[0]);
			MovVec3f(&BC[0], &Vertices[3*cur_tri[1]]);
			MovVec3f(&BEye[0], &Vertices[3*cur_tri[1]]);
		//subtract other location vector
			axpby(1,&Vertices[3*cur_tri[1]],-1,&AB[0], 3,1);
			axpby(1,&Vertices[3*cur_tri[2]],-1,&AC[0], 3,1);
			axpby(1,&Vertices[3*cur_tri[2]],-1,&BC[0], 3,1);
			axpby(1,&Eye[0],-1,&AEye[0], 3,1);
			axpby(1,&Eye[0],-1,&BEye[0], 3,1);
		CrossProduct(&AB[0], &AC[0], &Normal[0]);
			float scale_normal = (VectorNorm(&Normal[0],3)==0) ? (0) : (1.0/VectorNorm(&Normal[0],3));
			ScaleMat(scale_normal,&Normal[0],3,1);
		/***TEST FOR WALL************************************
		 * Define Wall as Planes slanted more than 50° from *
		 * the y=0 plane. Therefore assuming the man is a   *
		 * cylindiracalobject with axis: (Eye.x,0,Eye.z) to *
		 * Eye and radius 0.5m                              *
		 * SCP gives angle diverging from (0,1,0). It's a   *
		 * Wall if this angle is e.g. 90°+-50° => 40-140°   *
		 ****************************************************/
		#define WALL_MARGIN 50
		float y_angle = acos(Normal[1]/VectorNorm(&Normal[0],3));	// here the y-coordinate is height
		if ( y_angle > (90-WALL_MARGIN)*M_PI/180
			 && y_angle < (90+WALL_MARGIN)*M_PI/180 ) {		// Pi>acos(x)>0
			/***TEST FOR COLLISON WITH WALLFACE*******************
			 * if distance to plane is positive and smaller than *
			 * COLLISION_MARGIN and vertical distance vector     *
			 * goes through triangle-face then it collides       *
			 *****************************************************
			 * Imagine Prism instead of triangle with normals    *
			 * pointing outwards. A Point lies in this infinte   *
			 * long prism if the distance to these normals is    *
			 * negative for all (=> Normalization doesn't matter *
			 * Normal of a line which lies in the triangle-plane:*
			 *  CB_nrml = Normal x CB is approx parallel to      *
			 *  AM_{BC} => dist:=<BEye,CB_orth> must be =< 0     *
			 *  AEye doesn't work because A isn't part of the    *
			 *  line or rather the plane perpendicular to the    *
			 *  triangle and containing CB                       *
			 * For the other two we can use AEye                 *

			 *****************************************************/
			#define COLLISION_MARGIN 0.5
			float *distance=0;
			ScalarProduct3f(&AEye[0],&Normal[0],&distance);
			//AB_nrml=AB x Normal, AC_nrml=CA x Normal, BC_nrml=BC x Normal
				float AB_nrml[3], BC_nrml[3], CA_nrml[3];
				CrossProduct(&AB[0], &Normal[0], &AB_nrml[0]);
				CrossProduct(&BC[0], &Normal[0], &BC_nrml[0]);
				CrossProduct(&Normal[0], &AC[0], &CA_nrml[0]);	//Because we use AC instead of CA
				float *dist_AB=0,*dist_BC=0,*dist_CA=0;
				ScalarProduct3f(&AB_nrml[0],&AEye[0],&dist_AB);
				ScalarProduct3f(&BC_nrml[0],&BEye[0],&dist_BC);
				ScalarProduct3f(&CA_nrml[0],&AEye[0],&dist_CA);
			//has to be adjusted for real sphere collision (it's more like a decting collision with an infinitesimal thin cylinder whose axis is parallel to the normal of the plane)
			if (dist_AB[0]<0 && dist_BC[0]<0 && dist_CA[0]<0) {
				//if <AEye,AB/|AB|> < 1 (other corner point is OA+1*AB and OA+1*AC)
				if (distance[0]<COLLISION_MARGIN && distance[0]>0.001) {	//Normal dot (0,0,1) == Normal[2]
					#ifdef DEBUG_SIM
						fprintf(stderr, "Camera.Eye(%2.2f,%2.2f,%2.2f), Face: A:(%2.2f,%2.2f,%2.2f) B:(%2.2f,%2.2f,%2.2f) C:(%2.2f,%2.2f,%2.2f)\n",
							Eye[0], Eye[1], Eye[2],
							Vertices[3*cur_tri[0]+0],Vertices[3*cur_tri[0]+1],Vertices[3*cur_tri[0]+2],
							Vertices[3*cur_tri[1]+0],Vertices[3*cur_tri[1]+1],Vertices[3*cur_tri[1]+2],
							Vertices[3*cur_tri[2]+0],Vertices[3*cur_tri[2]+1],Vertices[3*cur_tri[2]+2]);
						fprintf(stderr, "\tAB_nrml(%2.2f,%2.2f,%2.2f), BC_nrml(%2.2f,%2.2f,%2.2f), CA_nrml(%2.2f,%2.2f,%2.2f) => dist_AB:%2.2f, dist_BC:%2.2f, dist_CA:%2.2f\n",
							AB_nrml[0],AB_nrml[1],AB_nrml[2],
							BC_nrml[0],BC_nrml[1],BC_nrml[2],
							CA_nrml[0],CA_nrml[1],CA_nrml[2],
							dist_AB[0],dist_BC[0],dist_CA[0]);
						//if like in the case of the city only paralelogramlike quads are used it is always in the triangle if it is in the parallogram
						fprintf(stderr, "\tAre we near this face? Distance:%f, Camera.Eye: (%f,%f,%f)\n",
							distance[0],Eye[0], Eye[1], Eye[2]);
						fprintf(stderr, "\tWe are near this face!");
					#endif
					Collides = 1;
				}
			}
			free(distance);
		}
	}
	return Collides;
}
void RecalculateCameraFromAngles(CameraData *Camera) {
	//calculate Camera.Center and Camera.Up from Phi(Yaw),Theta(Pitch) and Psi(Roll) and Camera.Eye
	//Calculate new Camera.Center (isn't actually needed if only Psi changes)
	//(Phi,Theta)=(0,0)=>Center-Eye=(0,0,-1) and (Phi,Theta)=(45,45)=>Center-Eye=(0.7,0.7,-0.7)
	Camera->Center[0] = Camera->Eye[0] +  sin(Camera->Phi)*cos(Camera->Theta);
	Camera->Center[1] = Camera->Eye[1] +  sin(Camera->Theta);
	Camera->Center[2] = Camera->Eye[2] + -cos(Camera->Phi)*cos(Camera->Theta);
	//Calculate new Camera.Up which depends on (in this order) Phi(around y-axis) and Theta(up-down) and Psi(around viewing-axis)
	//Dependance of Phi and Theta is exactly the same as Center only that Up is pitched up 90°=Pi/2 (Up is a relative vector!)
	Camera->Up[0] =  sin(Camera->Phi)*cos(Camera->Theta+M_PI/2);
	Camera->Up[1] =  sin(Camera->Theta+M_PI/2);
	Camera->Up[2] = -cos(Camera->Phi)*cos(Camera->Theta+M_PI/2);
	//Caclulate Psi by rotating Camera.Up around Viewing-Axis=Camera.Center-Camera.Eye
	float *rot_mat=0, *tmp_vec=0;
	CreateRotationMat3f(
		Camera->Center[0]-Camera->Eye[0],
		Camera->Center[1]-Camera->Eye[1],
		Camera->Center[2]-Camera->Eye[2], 180*Camera->Psi/M_PI, &rot_mat);
	MultiplicateMat(rot_mat,Camera->Up,3,3,3,1,&tmp_vec);
	memcpy(&(Camera->Up),tmp_vec,3*sizeof(float));
	free(rot_mat); free(tmp_vec);
	#ifdef DEBUG_SIM2
		fprintf(stderr,"Phi:%f, Theta:%f, Psi:%f, Camera.Up=(%f,%f,%f), Camera.Center=(%f,%f,%f), Camera.Eye=(%f,%f,%f)\n",
		Camera.Phi, Camera.Theta,Camera.Psi,Camera.Up[0],Camera.Up[1],Camera.Up[2],Camera.Center[0],Camera.Center[1],Camera.Center[2],Camera.Eye[0],Camera.Eye[1],Camera.Eye[2]);
		float *tmp=malloc(16*sizeof(float));
		CalcView(&Camera,tmp);
		fprintf(stderr,"Look-At-Matrix: ");
		DebugDisplayMat(tmp,4,4);
		free(tmp);
	#endif
	return;
}
void ProcessAsciiKeys(unsigned char key, int x, int y) {
	static unsigned char LastKey=0, IsCalculating=0;
	//Test whether this functions was called when this function was being executed...
	if (IsCalculating) {
		fprintf(stderr, "Another Key pressed, while we still processed the last one\n");
		return;
	} else
		IsCalculating = 1;
	fprintf(stderr, "%c pressed\n",key);
	float cam_transl_sign = 0;
	unsigned short int cam_transl_coord = 0xFFFF;
	float* transl_coord=0;
	switch(LastKey) {
		case 'x': transl_coord = &x_transl; break;
		case 'y': transl_coord = &y_transl; break;
		case 'z': transl_coord = &z_transl; break;
		case 'W': transl_coord = &w_transl; break;
		case 'S': transl_coord = &obj_scale; break;
		default: transl_coord = 0; break;
	}
	switch(key) {
		case 'w': cam_transl_coord = 1; cam_transl_sign = +1; break;
		case 's': cam_transl_coord = 1; cam_transl_sign = -1; break;
		case 'd': cam_transl_coord = 0; cam_transl_sign = +1; break;
		case 'a': cam_transl_coord = 0; cam_transl_sign = -1; break;
		case 'e': cam_transl_coord = 2; cam_transl_sign = +1; break;
		case 'q': cam_transl_coord = 2; cam_transl_sign = -1; break;
		case 'x': LastKey=key; break;
		case 'y': LastKey=key; break;
		case 'z': LastKey=key; break;
		case 'r': LastKey=key; break;
		case 'W': LastKey=key; break;
		case 'S': LastKey=key; break;
		case '+':
			if (LastKey == 'r') {
				obj_rot += 5;
			} else if (transl_coord!=0)
				(*transl_coord)+=0.1;
			fprintf(stderr,"Translating: X:%f, Y:%f, Z:%f, W:%f, R:%f, S:%f\n",x_transl,y_transl,z_transl,w_transl,obj_rot,obj_scale);
			break;
		case '-':
			if (LastKey == 'r') {
				obj_rot -= 5;
			} else if (transl_coord!=0)
				(*transl_coord)-=0.1;
			fprintf(stderr,"Translating: X:%f, Y:%f, Z:%f, W:%f, R:%f, S:%f\n",x_transl,y_transl,z_transl,w_transl,obj_rot,obj_scale);
			break;
		case 'R':
			obj_rot=0; x_transl=0; y_transl=0; z_transl=0; w_transl=1; obj_scale=1;
			Camera.Eye[0]=0;Camera.Eye[1]=1.7;Camera.Eye[2]=0;
			Camera.Center[0]=0;Camera.Center[1]=0;Camera.Center[2]=-3;
			Camera.Up[0]=0;Camera.Up[1]=2.7;Camera.Up[2]=0;
			Camera.Theta = -0.51555;
			Camera.Phi = 0;
			Camera.Psi = 0;
			break;
		case 27:  glutLeaveMainLoop(); break;
		default: break;
	}
	if (cam_transl_coord==1) {	//w,s
		//Remember old Camera.Center and Camera.Eye in case we collide with somethin
			float CenterSAV[3], EyeSAV[3];
			MovVec3f(&CenterSAV[0],&Camera.Center[0]);
			MovVec3f(&EyeSAV[0],&Camera.Eye[0]);
		float* tmp_vec = malloc(3*sizeof(float));
		memcpy(tmp_vec,Camera.Center,3*sizeof(float));
		axpby(-1,Camera.Eye,1,tmp_vec,3,1);				//tmp_vec = <Center-Eye,(1,0,0)> (1,0,0) + <Center-Eye,(0,0,1)> (0,0,1)
		tmp_vec[1] = 0;	//project onto x,z-plane because we are human and only can move horizontally
		float norm = 0;
		if (VectorNorm(tmp_vec,3)!=0)
			norm = cam_transl_sign*0.1f/VectorNorm(tmp_vec,3);
		axpby(norm, tmp_vec,1,Camera.Eye,3,1);
		axpby(norm, tmp_vec,1,Camera.Center,3,1);
		#ifdef DEBUG_SIM2
			fprintf(stderr,"<(1,0,1),Center-Eye> = ");
			DebugDisplayMat(tmp_vec,3,1);
		#endif
		free(tmp_vec);
		//Test for Collision and revert Camera if true
			if (CollisionSphereObject(&Camera.Eye[0], 0.3, &(Scene.Objects[0]), Scene.Vertices)) {
				fprintf(stderr, "We did collide with something. Revert Movement\n");
				MovVec3f(&Camera.Center[0], &CenterSAV[0]);
				MovVec3f(&Camera.Eye[0]   , &EyeSAV[0]   );
			}
	} else if (cam_transl_coord==0) {	//a,d
		//Remember old Camera.Center and Camera.Eye in case we collide with somethin
			float CenterSAV[3], EyeSAV[3];
			MovVec3f(&CenterSAV[0],&Camera.Center[0]);
			MovVec3f(&EyeSAV[0],&Camera.Eye[0]);
		//tv2 = <(Center-Eye) x Up,(1,0,1)> (projects onto x,z-plane because we are human and only can move horizontally)
			float tv1[3],tv2[3];
			memcpy(&tv1[0],Camera.Center,3*sizeof(float));
			axpby(-1,Camera.Eye,1,&tv1[0],3,1);
			CrossProduct(&tv1[0], Camera.Up, &tv2[0]);
			tv2[1] = 0;
		//Camera = Camera+dist*tv2 (Move along tv2 (temp vector 2))
			float norm = 0;
			if (VectorNorm(&tv2[0],3)!=0)
				norm = cam_transl_sign*0.1f/VectorNorm(&tv2[0],3);
			axpby(norm, &tv2[0],1,Camera.Eye,3,1);
			axpby(norm, &tv2[0],1,Camera.Center,3,1);
			#ifdef DEBUG_SIM
				fprintf(stderr,"<(1,0,1),(Center-Eye) x Camera.Up> = ");
				DebugDisplayMat(&tv2[0],3,1);
			#endif
		//Test for Collision and revert Camera if true
			if (CollisionSphereObject(&Camera.Eye[0], 0.3, &Scene.Objects[0], Scene.Vertices)) {
				fprintf(stderr, "We did collide with something. Revert Movement\n");
				MovVec3f(&Camera.Center[0], &CenterSAV[0]);
				MovVec3f(&Camera.Eye[0]   , &EyeSAV[0]   );
			}
	} else {
		/*Camera.Eye[cam_transl_coord]		+= cam_transl_sign * 0.1f;
		Camera.Center[cam_transl_coord]		+= cam_transl_sign * 0.1f;*/
	}
	#ifdef DEBUG_SIM2
		fprintf(stderr,"Phi:%f, Theta:%f, Psi:%f, Camera.Up=(%f,%f,%f), Camera.Center=(%f,%f,%f), Camera.Eye=(%f,%f,%f)\n",
		Camera.Phi, Camera.Theta,Camera.Psi,Camera.Up[0],Camera.Up[1],Camera.Up[2],Camera.Center[0],Camera.Center[1],Camera.Center[2],Camera.Eye[0],Camera.Eye[1],Camera.Eye[2]);
	#endif
	//glutPostRedisplay();
	RenderScene();
	#ifdef DEBUG_SIM2
		fprintf(stderr, "glutPostRedisplay.\n");
	#endif
	IsCalculating = 0;	//show that we are ready and the next key can be processed
	return;
}
/* key The key whose press triggers the callback
 * x,y The coordinates of the mouse relative to the window at the time the key is pressed */
void ProcessSpecialKeys(int key, int x, int y) {
    /*GLUT_KEY_F1, GLUT_KEY_F2, ..., GLUT_KEY_F12 - F1 through F12 keys
    GLUT_KEY_PAGE_UP, GLUT_KEY_PAGE_DOWN - Page Up and Page Down keys
    GLUT_KEY_HOME, GLUT_KEY_END - Home and End keys
    GLUT_KEY_INSERT - Insert key*/
	switch(key) {
		case GLUT_KEY_F9: fprintf(stderr, "Toggle Fullscreen\n"); glutFullScreenToggle(); break;
		case GLUT_KEY_LEFT:
			fprintf(stderr, "LEFT_KEY pressed.\n");
			Camera.Phi -= 5.0f*M_PI/180.0f;	//head rotation (around z-axis)
#if 1==0
			//Center = Eye + Rotate*(Center-Eye)
			float* vec_ec = malloc(3*sizeof(float));
			float* rot_mat=0; CreateRotationMat3f(0,1,0, 10, &rot_mat);
			memcpy(vec_ec, Camera.Eye, 3*sizeof(float));		// vec_ec = Camera.Eye
			axpby(1,Camera.Center,-1,vec_ec,3,1);				// vec_ec = Camera.Center - vec_ec
			float* tmp_vec=0;
			MultiplicateMat(rot_mat,vec_ec,3,3,3,1,&tmp_vec);	// tmp_vec = rot_mat * vec_ec
			axpby(1,Camera.Eye,1,tmp_vec,3,1);					// tmp_vec = Camera.Eye + tmp_vec
			axpby(1,tmp_vec,0,Camera.Center,3,1);				// Camera.Center = tmp_vec
			#ifdef DEBUG_SIM
				fprintf(stderr, "\tNew Camera.Center: ");
				DebugDisplayMat(Camera.Center,3,1);
				float* projection_matrix = malloc(16*sizeof(float));
				CalcView(&Camera, projection_matrix);
				fprintf(stderr, "\tProjection View: ");
				DebugDisplayMat(projection_matrix,3,1);
				fprintf(stderr,"\n");
				free(projection_matrix);
			#endif
			//Calculate new Camera.Up
			MultiplicateMat(rot_mat,Camera.Up,3,3,3,1,&tmp_vec);
			memcpy(Camera.Up,tmp_vec,3*sizeof(float));
			free(vec_ec); free(rot_mat); free(tmp_vec);
#endif
			break;
		case GLUT_KEY_RIGHT:
			fprintf(stderr, "RIGHT_KEY pressed.\n");
			Camera.Phi += 5.0f*M_PI/180.0f;			//head rotation (around z-axis)
			break;
		case GLUT_KEY_DOWN:
			fprintf(stderr, "RIGHT_DOWN pressed.\n");
			if (Camera.Theta > -40.0f)
				Camera.Theta -= 5.0f*M_PI/180.0f;			//head "nodding"
			break;
		case GLUT_KEY_UP:
			fprintf(stderr, "RIGHT_UP pressed.\n");
			if (Camera.Theta < 40.0f)
				Camera.Theta += 5.0f*M_PI/180.0f;
			break;
		/*case GLUT_KEY_HOME:
			fprintf(stderr, "POS1/HOME pressed.\n");	//head tilt ... not possible for humans
			Camera.Psi -= 5.0f*M_PI/180.0f;
			break;
		case GLUT_KEY_PAGE_UP:
			fprintf(stderr, "PAGE_UP pressed.\n");
			Camera.Psi += 5.0f*M_PI/180.0f;
			break;*/
		default: break;
	}
	//calculate Camera.Center and Camera.Up from Phi,Theta and Camera.Eye
	if(key==GLUT_KEY_LEFT || key==GLUT_KEY_RIGHT || key==GLUT_KEY_UP || key==GLUT_KEY_DOWN || key==GLUT_KEY_HOME || key==GLUT_KEY_PAGE_UP ) {
		RecalculateCameraFromAngles(&Camera);
	}
	glutPostRedisplay();
	return;
}
void MouseWhileKeyPressed(int x, int y) {
	#ifdef DEBUG_SIM
		fprintf(stderr, "Key Pressed. Mouse: X:%i ,Y:%i\n",x,y);
	#endif
	return;
}
void MouseWhileNoKeyPressed(int x, int y) {
	//GLUT LEFT BUTTON, GLUT MIDDLE BUTTON or GLUT RIGHT BUTTON. state is GLUT UP xor GLUT DOWN
	static int lastX=0xFFFF, lastY=0xFFFF;
	int deltaX=0, deltaY=0;
	if (lastX!=0xFFFF) {
		deltaX = x - lastX;
		deltaY = y - lastY;
	}
	lastX = x; lastY = y;
	if( deltaX == 0 && deltaY == 0 ) return;

	#define err_border 20
	if( x<=err_border || y<=err_border || x>=Screen.Width-err_border || y>=Screen.Height-err_border) {
		lastX = Screen.Width/2; lastY = Screen.Height/2;
		glutWarpPointer( lastX, lastY );
	}
	#undef err_border

	Camera.Phi   +=  deltaX*M_PI/Screen.Width;			//300 pixel per 180°
	Camera.Theta += -deltaY*M_PI/Screen.Height;			//300 pixel per 180°
	if(Camera.Theta<-90) Camera.Theta=-90;
	if(Camera.Theta>+90) Camera.Theta=+90;
	RecalculateCameraFromAngles(&Camera);
	glutPostRedisplay();
	#ifdef DEBUG_SIM2
		fprintf(stderr, "No Key Pressed. Mouse: X:%i ,Y:%i\n",x,y);
	#endif
	return;
}
void ProcessMouseButtons(int button, int state,int x, int y) {
	#ifdef DEBUG_SIM
		fprintf(stderr, "Button: %i, State:%i Mouse-Position: X:%i ,Y:%i\n",button,state,x,y);
	#endif
	return;
}

void CalculateNormals(OBJData Scene) {
	//Go through every Face in 'Scene'
	unsigned int obj_ind;
	for (obj_ind=0; obj_ind<Scene.ObjectCount; obj_ind++) {
		Object3D* cur_obj = &(Scene.Objects[obj_ind]);
		cur_obj->Normals = malloc(3*sizeof(float)*cur_obj->TriangleCount);
		cur_obj->NormalLines = malloc(2*3*sizeof(float)*cur_obj->TriangleCount);
		unsigned int triangle_num;
		#pragma omp parallel for
		for (triangle_num=0; triangle_num<cur_obj->TriangleCount; triangle_num++) {
			//Initialize some vector: Center=(OA+OB+OC)/3,Normal
			//and store these in cur_obj->Normals (Face-1-Center, Center+Normal, Normal; Face-2-Center, ...)
			GLuint* cur_tri = &cur_obj->Triangles[3*triangle_num];
			float* Center   = &cur_obj->NormalLines[2*3*triangle_num];
			float* NrmPoint = &cur_obj->NormalLines[2*3*triangle_num+3];
			float* Normal   = &cur_obj->Normals[3*triangle_num];
			float AB[3],AC[3];
			MovVec3f(&AB[0], &Scene.Vertices[3*cur_tri[0]]);				// AB = OA
			MovVec3f(&AC[0], &AB[0]);										// AC = OA
			MovVec3f(Center, &AB[0]);										// Center = OA
			axpby(-1,&Scene.Vertices[3*cur_tri[1]],1,&AB[0], 3,1);			// AB = AB-Vertices[B_ID]
			axpby(-1,&Scene.Vertices[3*cur_tri[2]],1,&AC[0], 3,1);			// AC = AC-Vertices[C_ID]
			CrossProduct(&AB[0], &AC[0], Normal);
			float scale = 0;
			if (VectorNorm(Normal,3)!=0)
				scale = 1/VectorNorm(Normal,3);
			ScaleMat(scale,Normal,3,1);
			axpby(1,&Scene.Vertices[3*cur_tri[1]],1,Center, 3,1);			// Center += Vertices[B_ID]
			axpby(1.0/3.0,&Scene.Vertices[3*cur_tri[2]],1.0/3.0,Center, 3,1);		// Center = Vertices[C_ID]/3 + Center/3
			MovVec3f(NrmPoint, Center);										// NrmPoint = Center
			axpby(1,Normal,1,NrmPoint, 3,1);								// NrmPoint += Normal
		}
	}
	return;
}
OBJData LoadOBJFile(const char* filename) {
	/*********************************************************
	 * converts quad 1-2-3-4 to 2 triangles 1-2-3 and 4-1-3  *
	 *********************************************************/
	FILE* obj_file = fopen(filename, "rb");
	OBJData OutputData;

	#ifdef DEBUG_SIM
		fprintf(stderr, "Counting Numbers of Objects and Vertices...");
	#endif
	//count number of objects and Vertices
	OutputData.ObjectCount = 0;
	OutputData.VertexCount = 0;
	while(!feof(obj_file)) {
		char first_char;
		do {
			first_char = fgetc(obj_file);
		} while( (first_char==13||first_char==10) && !feof(obj_file) );
		if(first_char=='o')
			OutputData.ObjectCount++;
		else if(first_char=='v')
			OutputData.VertexCount++;
		//Read until EOL and don't store read date (*)
		fscanf(obj_file,"%*[^\n]%*c");	//assuming the newline is composed of ascii: [13(CR=\r)]10(LF=\n) (%*c reads \r, because the other command stops BEFOER \r)
	}
	rewind(obj_file);
	//Allocate Object-Array
	OutputData.Objects = malloc(OutputData.ObjectCount * sizeof(Object3D));
	OutputData.Vertices = malloc(OutputData.VertexCount * 3*sizeof(float));
	#ifdef DEBUG_SIM
		fprintf(stderr, "OK\nObjects: %u, Vertices: %u\nCounting Number of Faces per Object:\n", OutputData.ObjectCount, OutputData.VertexCount);
	#endif


	//count f to allocate memory (completely fill in Object3D (not arrays, but pointers to arrays))
	signed int cur_object_index=-1, cur_vertex_index=-1;
	Object3D *cur_object = 0;	//there should be no 'v' or 'f' statements without a prior 'o' statement :S
	while(!feof(obj_file)) {
		char first_char;
		do {first_char = fgetc(obj_file);} while( (first_char==13||first_char==10) && !feof(obj_file) );	//read EOL-Characters
		switch (first_char) {
			case 'o':
				cur_object_index++;
				if (cur_object_index >= OutputData.ObjectCount)
					fprintf(stderr,"More Objects than initally thought! Can't handle that.\n");
				cur_object = &OutputData.Objects[cur_object_index];
				fscanf(obj_file, " %s", (char*)&(cur_object->Name));
				cur_object->TriangleCount = 0;
				cur_object->QuadCount = 0;
				break;
			case 'v':
				cur_vertex_index++;
				if (cur_vertex_index >= OutputData.VertexCount) {
					fprintf(stderr, "Calculated Array-Length for Vertices exceeded when loading %s\n",filename);
					break;
				}
				fscanf(obj_file, "%f %f %f%*[^\n]",
					&(OutputData.Vertices[3*cur_vertex_index]),
					&(OutputData.Vertices[3*cur_vertex_index+1]),
					&(OutputData.Vertices[3*cur_vertex_index+2])
				);	//reads until EOL
				break;
			case 'f':
				if (cur_object==0) {
					fprintf(stderr, "Error no object-statement prior to Vertex-Statement in %s\n",filename);
					break;
				}
				char tmp_sr[256];
				memset(tmp_sr,0,256*sizeof(char));
				fgets(&tmp_sr[0], 256, obj_file);
				fseek(obj_file,-1,SEEK_CUR);	//can't jump back to prior line, because at least f was read from this line. But it only should jump back "\n"
				unsigned int i = 0, numbers = 0;
				for(i=0;i<=255-1;i++) {	//-1 because we also want to analyze the next string
					if (tmp_sr[i]==' ' && tmp_sr[i+1]>='0' && tmp_sr[i+1]<='9') {
						numbers++;
						i++;		//KMP like
					}
				}
				if (numbers==3)
					cur_object->TriangleCount++;
				else if (numbers==4) {
					cur_object->TriangleCount += 2;
					cur_object->QuadCount++;
				} else
					fprintf(stderr, "Invalid Number after face-directive: %u, String:%s\n",numbers,tmp_sr);
				break;
			default:
				break;
		}
		//Read until EOL and don't store read date (*)
		fscanf(obj_file,"%*[^\n]%*c");	//assuming the newline is composed of ascii: [13(CR=\r)]10(LF=\n) (%*c reads \r, because the other command stops BEFOER \r)
	}
	rewind(obj_file);
	#ifdef DEBUG_SIM
		for (cur_object_index=0; cur_object_index<OutputData.ObjectCount; cur_object_index++) {
			fprintf(stderr, "\t%s: %u Triangles and %u Quads\n", OutputData.Objects[cur_object_index].Name, OutputData.Objects[cur_object_index].TriangleCount,OutputData.Objects[cur_object_index].QuadCount);
		}
	#endif


	unsigned int cur_triangle_index=-1;
	cur_object_index = -1;
	//Load data into allocated arrays
	#ifdef DEBUG_SIM
		fprintf(stderr, "First Chars while loading: ");
	#endif
	while(!feof(obj_file)) {
		char first_char;
		do {
			first_char = fgetc(obj_file);
		} while( (first_char==13||first_char==10) && !feof(obj_file) );
		#ifdef DEBUG_SIM
			fprintf(stderr, "%c, ",first_char);
		#endif
		switch (first_char) {
			case 'o':
				//select cur_object to write in for 'v' and 'f'
				cur_object_index++;
				if (cur_object_index >= OutputData.ObjectCount)
					fprintf(stderr, "More Objects than initally thought! Can't handle that.\n");
				//Allocate Memory for Triangle-Data
				cur_object = &OutputData.Objects[cur_object_index];
				cur_object->Triangles = malloc( 3*sizeof(GLuint)*cur_object->TriangleCount );
				cur_triangle_index = -1;
				break;
			case 'f':
				if (cur_object==0) {
					fprintf(stderr, "Error no object-statement prior to face-statements in %s\n",filename);
					break;
				}
				if (cur_triangle_index+1 >= cur_object->TriangleCount) {
					fprintf(stderr, "Calculated Array-Length for Triangles exceeded when loading %s\n",cur_object->Name);
					break;
				}
				char tmp_sr[256];
				memset(tmp_sr,0,256*sizeof(char));
				fgets(&tmp_sr[0], 256, obj_file);
				fseek(obj_file,-1,SEEK_CUR);	//can't jump back to prior line, because at least f was read from this line. But it only should jump back "\n"
				unsigned int i = 0, numbers = 0;
				for(i=0;i<=255-1;i++) {	//-1 because we also want to analyze the next string
					if (tmp_sr[i]==' ' && tmp_sr[i+1]>='0' && tmp_sr[i+1]<='9') {
						numbers++;
						i++;		//KMP like
					}
					if (tmp_sr[i]==0)
						break;
				}
				if (numbers==3) {
					cur_triangle_index++;
					sscanf(&tmp_sr[0], "%u %u %u%*[^\n]",
						&(cur_object->Triangles[3*cur_triangle_index]),
						&(cur_object->Triangles[3*cur_triangle_index+1]),
						&(cur_object->Triangles[3*cur_triangle_index+2])
					);	//reads until EOL
					//decrement vertexID, because .obj-format starts counting from 1 but OpenGL from 0
					cur_object->Triangles[3*cur_triangle_index]--;
					cur_object->Triangles[3*cur_triangle_index+1]--;
					cur_object->Triangles[3*cur_triangle_index+2]--;
				} else if (numbers==4) {
					cur_triangle_index++;
					sscanf(&tmp_sr[0], "%u %u %u %u%*[^\n]",
						&(cur_object->Triangles[3*cur_triangle_index]),
						&(cur_object->Triangles[3*cur_triangle_index+1]),
						&(cur_object->Triangles[3*cur_triangle_index+2]),
						&(cur_object->Triangles[3*(cur_triangle_index+1)])
					);	//reads until EOL
					//decrement vertexID, because .obj-format starts counting from 1 but OpenGL from 0
					//convert quad 1-2-3-4 to 2 triangles 1-2-3 and 4-1-3
					cur_object->Triangles[3*cur_triangle_index]--;
					cur_object->Triangles[3*cur_triangle_index+1]--;
					cur_object->Triangles[3*cur_triangle_index+2]--;
					cur_object->Triangles[3*(cur_triangle_index+1)]--;
					cur_object->Triangles[3*(cur_triangle_index+1)+1] = cur_object->Triangles[3*cur_triangle_index];
					cur_object->Triangles[3*(cur_triangle_index+1)+2] = cur_object->Triangles[3*cur_triangle_index+2];
					cur_triangle_index++;
				} else
					fprintf(stderr, "Invalid Number after face-directive: %u, String:%s\n",numbers,tmp_sr);
				break;
			default:
				fscanf(obj_file, "%*[^\n]%*c");
				break;
		}
	}
#ifdef DEBUG_SIM
		fprintf(stderr, "\n");
#endif

#ifdef DEBUG_SIM2
	fprintf(stderr, "Vertice Array loaded from %s, VertexCount:%u:\n",filename,OutputData.VertexCount);
	for (cur_vertex_index=0; cur_vertex_index<OutputData.VertexCount; cur_vertex_index++) {
		fprintf(stderr, "v %f %f %f\n",OutputData.Vertices[3*cur_vertex_index],OutputData.Vertices[3*cur_vertex_index+1],OutputData.Vertices[3*cur_vertex_index+2]);
	}
	for(cur_object_index=0; cur_object_index<OutputData.ObjectCount; cur_object_index++) {
		cur_object = &(OutputData.Objects[cur_object_index]);
		fprintf(stderr, "Object: %s, TriangleCount: %i\n",cur_object->Name, cur_object->TriangleCount);
		for (cur_triangle_index=0; cur_triangle_index<cur_object->TriangleCount; cur_triangle_index++) {
			fprintf(stderr, "f %u %u %u\n",cur_object->Triangles[3*cur_triangle_index+0], cur_object->Triangles[3*cur_triangle_index+1], cur_object->Triangles[3*cur_triangle_index+2]);
		}
	}
#endif
	return OutputData;
}
/* Reads complete binary file into a string allocated with malloc and returns it returns 0 pointer if there is an error */
char* LoadSourceFile(const char* filename) {
	FILE* source_file = fopen(filename,"rb");
	if (source_file == 0) {
		fprintf(stderr, "Couldn't open %s\n",filename);
		return 0;
	}
	//Determine File Lenght:
	fseek(source_file, 0, SEEK_END);
	int file_size = ftell(source_file);
	char* buffer = malloc(file_size*sizeof(char)+1);	//+1 for 0-end-byte
	fseek(source_file, 0, SEEK_SET);
	//Read File and write errors into stderr
	int bytes_read = fread( buffer, sizeof(char), file_size, source_file );
	fclose(source_file);
	buffer[file_size] = 0;
	fprintf(stderr, "FragShader.c Size: %i\nNumber of Bytes read from file: %i\n",file_size, bytes_read);
	char* buf_ptr = 0;
	do {	//replace all ascii 13 with space, because only ascii 10 is needed for newline in textmode (or else we would get output: 13,13,10)
		buf_ptr = memchr(buffer, 13, file_size);	//substituting buffer with buf_ptr would be better, but then we have to calculate the third argument by subtracting pointers from each other :S
		if (buf_ptr!=0)
			*buf_ptr = ' ';
		else
			break;
	} while(1);
	fwrite( buffer, sizeof(char), file_size, stderr );
	return buffer;
}

/***********************************************************************************
 * Loads Shader SourceFile and compiles it returning the ID to the compiled Shader.*
 * Needs to know what kind of Shader and Filename of Source                        *
 ***********************************************************************************/
GLuint LoadShader(const char* filename, GLenum type) {
	char* source = LoadSourceFile(filename);
	if (!source)
		return 0;
	GLuint shaderID = glCreateShader(type);
	#define NO_AUTOMATIC_GL_VERSION_DETECTION
	#ifdef NO_AUTOMATIC_GL_VERSION_DETECTION
	const GLchar* sources[2] = {source, ""};
	#else
	const GLchar* sources[2] = {
		#ifdef GL_ES_VERSION_2_0
			"#version 100\n"
			"#define GLES2\n",
		#else
			"#version 400\n",
		#endif
		source
	};
	#endif
	glShaderSource(shaderID, 2, (const char**)&sources, NULL);	//implicit const declaring won't "dig deep enough" therefore the explicit cast
	fprintf(stderr,"Freed Source Char Array.\n");
	free(source);
	fprintf(stderr,"Now Compiling Shader.\n");

	glCompileShader(shaderID);
#ifdef DEBUG_SIM
	#ifdef GL_ES_VERSION_2_0
		fprintf(stderr,"Not embedded System\n");
	#else
		GLint version;
		glGetIntegerv(GL_MAJOR_VERSION, &version);
		fprintf(stderr,"OpenGL Major Version: %i\n",version);
	#endif
	//Get Compiler Log
	GLint infologlength;
	glGetShaderiv(shaderID,GL_INFO_LOG_LENGTH,&infologlength);
	char* infolog = malloc( sizeof(char)*infologlength );
	glGetShaderInfoLog(shaderID, infologlength, NULL,  infolog);
	fprintf(stderr, "\n------------- LOG -------------\n-- InfoLog-Length: %i\n%s-------------------------------\n\n",infologlength, infolog);
#endif

	GLint compile_ok;
	glGetShaderiv(shaderID, GL_COMPILE_STATUS, &compile_ok);
	if (compile_ok == 0) {
		fprintf(stderr, "Error compiling shader: %s\n",filename);
		return 0;
	}
	return shaderID;
}

/* Bind Input Variables used in Shaders */
GLuint GetAttributeID(const char* attribute_name) {
	GLuint attribute_var_id = glGetAttribLocation(program, attribute_name);	//global variable GLint
	fprintf(stderr, "Binding %s ...", attribute_name);
	if (attribute_var_id == -1) {
		fprintf(stderr, "FAILED\n");
		return 0;
	} else {
		fprintf(stderr, "OK (ID=%u)\n",attribute_var_id);
		return attribute_var_id;
	}
}
/* Bind Input Variables used in Shaders */
GLuint GetUniformID(const char* attribute_name) {
	GLuint attribute_var_id = glGetUniformLocation(program, attribute_name);	//global variable GLint
	fprintf(stderr, "Binding %s ...", attribute_name);
	if (attribute_var_id == -1) {
		fprintf(stderr, "FAILED\n");
		return 0;
	} else {
		fprintf(stderr, "OK (ID=%u)\n",attribute_var_id);
		return attribute_var_id;
	}
}

/* This function creates all GLSL related stuff
Returns 1 when all is ok, 0 with a displayed error */
int init_resources(void) {
	//Initialize Camera:
	/*Camera.Eye[0]=0; Camera.Eye[1]=0; Camera.Eye[2]=0;
	Camera.Center[0]=0; Camera.Center[1]=0; Camera.Center[2]=-1;
	Camera.Up[0]=0; Camera.Up[1]=1; Camera.Up[2]=0; */
	//Initialize Matrices
	projection_matrix = malloc(16*sizeof(float));
	view_matrix = malloc(16*sizeof(float));
	model_matrix = malloc(16*sizeof(float));
	mvp = malloc(16*sizeof(float));
	//SetMatrix4fToIdentity(projection_matrix);
	CalcView(&Camera,view_matrix);
	//SetMatrix4fToIdentity(model_matrix);
	//SetMatrix4fToIdentity(mvp);

	// Enable alpha
	glEnable(GL_BLEND);
	glClearDepth(1.0f);				// Depth Buffer Setup
	glEnable(GL_DEPTH_TEST);		// Enables Depth Testing
	glDepthFunc(GL_LEQUAL);			// The Type Of Depth Test To Do
	//glEnable(GL_DEBUG_OUTPUT);	//only availbale für OpenGL 4.3+
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	//Compile Shaders
	GLuint vs = LoadShader("VertShader.c",GL_VERTEX_SHADER);	//Vertex Shader will be called for every vertex
	GLuint fs = LoadShader("FragShader.c",GL_FRAGMENT_SHADER);	//Fragment Shader will be called for every Pixel
	//Link Program
	GLint link_ok=0;
	program = glCreateProgram();	//Gloabl Variable program
	glAttachShader(program, vs);
	glAttachShader(program, fs);
	glLinkProgram(program);
#ifdef DEBUG_SIM
	{
		GLint infologlength;
		glGetProgramiv(program,GL_INFO_LOG_LENGTH,&infologlength);
		char* infolog = malloc( sizeof(char)*infologlength );
		glGetProgramInfoLog(program, infologlength, NULL,  infolog);
		fprintf(stderr, "Program Linker InfoLog-Length: %i\n------------- LOG -------------\n%s-------------------------------\n\n",infologlength, infolog);
	}
#endif
	glGetProgramiv(program, GL_LINK_STATUS, &link_ok);	//iv = integer value
	if (!link_ok) {
		fprintf(stderr, "glLinkProgram: not successfull");
		return 0;
	}

	/********************************************************
	 * Bind Input Variables used in Shaders                 *
	 ********************************************************/
	glGetError();	//clear error to not get errors from prior commands
	attribute_coord3d	= GetAttributeID("coord3d");
	attribute_v_color	= GetAttributeID("v_color");
	glUseProgram(program);
	uniform_fixed_color	= GetUniformID("fixed_color");
	uniform_mvp	= GetUniformID("mvp");
	glUniformMatrix4fv(uniform_mvp, 1, 1, mvp);
	#ifdef DEBUG_SIM
		fprintf(stderr,"glGetError after projection/modelview matrix: %i (GL_INVALID_ENUM=%i, GL_INVALID_VALUE=%i)\n",glGetError(),GL_INVALID_ENUM,GL_INVALID_VALUE);
	#endif

	//Load Icosaeder into an array of triangles (note that vertices will get duplicated many times most of the time)
	Scene = LoadOBJFile(OBJ_FILE_NAME);
	glGenBuffers(1, &vbo_vertices);
	glBindBuffer(GL_ARRAY_BUFFER, vbo_vertices);
	#ifdef DEBUG_SIM
		fprintf(stderr, "glGetError after BindBuffer(vbo_vertices): %i\n", glGetError());
	#endif
	//  glBufferData(GL_ARRAY_BUFFER, sizeof(triangle_vertices), triangle_vertices, GL_STATIC_DRAW);
	glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)*3*Scene.VertexCount, Scene.Vertices, GL_STATIC_DRAW);
	#ifdef DEBUG_SIM
		fprintf(stderr, "glGetError after BufferData Vertices: %i (GL_ARRAY_BUFFER=%i, Scene.Vertices=%p)\n", glGetError(),GL_ARRAY_BUFFER,Scene.Vertices);
	#endif
	glBindBuffer(GL_ARRAY_BUFFER,0);
	glGenBuffers(1, &ibo_triangles);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo_triangles);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(GLuint)*3*Scene.Objects[0].TriangleCount, Scene.Objects[0].Triangles, GL_STATIC_DRAW);
	#ifdef DEBUG_SIM
		fprintf(stderr, "glGetError after BufferData IBO Trianges: %i\n", glGetError());
	#endif
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,0);
	/****************************************
	 * Load Normal-Points in Virtual Buffer *
	 ****************************************/
	CalculateNormals(Scene);
	glGenBuffers(1, &vbo_normal_pts);
	glBindBuffer(GL_ARRAY_BUFFER, vbo_normal_pts);
	glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)*3*2*Scene.Objects[0].TriangleCount, Scene.Objects[0].NormalLines, GL_STATIC_DRAW);
	#ifdef DEBUG_SIM
		fprintf(stderr, "glGetError after BufferData NormalLines: %i\n", glGetError());
	#endif
	glBindBuffer(GL_ARRAY_BUFFER,0);

	//Load random colors for every fucking vertex for every fucking face (20*3 rgb values for the icosaeder)
	float* triangle_colors = malloc(sizeof(float)*3*Scene.VertexCount);
	int i;
	srand(time(NULL));
	for(i=0;i<3*Scene.VertexCount;i++) {
		triangle_colors[i] = ((float)(rand()%10000))/10000.0;
	}
	glGenBuffers(1, &vbo_triangle_colors);
	glBindBuffer(GL_ARRAY_BUFFER, vbo_triangle_colors);		//sets already loaded buffer (in init resources) for further adjustments
	glBufferData(GL_ARRAY_BUFFER, sizeof(float)*3*Scene.VertexCount, triangle_colors, GL_STATIC_DRAW);
	#ifdef DEBUG_SIM
		fprintf(stderr, "glGetError after BufferData triangle Colors: %i\n", glGetError());
	#endif
	glBindBuffer(GL_ARRAY_BUFFER, 0);		//we are finished setting the options of vbo_triangle_colors. Shouldn't doing this one time be enough?
	free(triangle_colors);

	return 1;
}

void CalcProjection(float fovy, float aspect, float zNear, float zFar, float* view_m) {	//zNear and zFar always positive! aspect=w/h
	//Note: OpenGL normally views orthogonal from infinity and only shows -w<[x,y,z]<+w but depth buffer from 0 to 1 oO?
	#ifdef DEBUG_SIM
		if(zNear<0 || zFar<0)
			fprintf(stderr, "zNear and zFar aren't positive for this CalcProjecton-Call! zNear:%f, zFar:%f\n",zNear,zFar);
	#endif
	if (view_m==0)
		return;
	/*****************************************************
	 * Assumptions: -Camera is at(0,0,0) and looks to -z *
	 *  zRar > zNear > 0, w=1                            *
	 * This Matrix does the following to a 4x1 Vector:   *
	 ****** FOVY *****************************************
	 * Wanted: Scale x,y                                 *
	 * f = cos(fovy/2)/sin(fovy/2) = A/O = span_Z/span_Y *
	 * y = f*y; x=f*x; z=z; w=w;                         *
	 * => y is scaled to match the fovy-angle at z=1     *
	 *    (and x is scaled to keep the aspect ratio)     *
	 ****** ASPECT-RATIO: ********************************
	 * (aspect = w/h = span_X/span_Y)                    *
	 * x = x/aspect; y=y; z=z; w=w;                      *
	 * => if w>h then the x-axis will be "elongated"     *
	 ****** Depth Scale **********************************
	 * Wanted: -a 2D object takes the same percent area  *
	 *   on the zNear-Plane as on the zFar-Plane         *
	 *  -for z=pyramide-tip=0 => w=0(scale=inf)          *
	 * height at zF|zN shall be tan theta*zF|zN          *
	 * => h(z)=tan theta*z*h0 => y=y/w; w(z)=1/tan th*z  *
	 * Increase/z_Span = [tan th*zF-tan th*zN]/(zF-zN) = *
	 *  tan theta = f; => w = -z*f;                      *
	 * <=> view_m[3][2]=-f; (other in this row is 0)     *
	 ****** INTERIM RESULTS ******************************
	 * |y|<tan(fovy/2), |x|<f*h/w => -1<(x,y)<+1         *
	 *****************************************************/
	memset(view_m,0,16*sizeof(float));
	float f = 1.0f/tanf((fovy/180*M_PI)/2);
	view_m[0*4+0] = f/aspect;	//x-scale
	view_m[1*4+1] = f;			//y-scale
	view_m[3*4+2] = -1;			//w=-z
	/****** CLIPPING: ************************************
	 * Unfortunately we can't really think in scaling    *
	 *  and translating anymore, because of W which      *
	 *  scales our vector by 1/z and is a fucking pain   *
	 *  to calculate. Therefore we simply look at it in  *
	 *  an abstract way, calculating a linear fit with   *
	 *  two points: z'(-zNear)=-1 and z'(-zFar)=+1       *
	 * We have two variables: a:=view_m[2*4+2],          *
	 *  b:=view_m[3*4+3] which gives us z'=1/w*(a*z+b)   *
	 * => z'(-zN)=-1/zN*(-a*zN+b)=a+b/zN=-1  |           *
	 * => z'(-zF)=               =a+b/zF=+1  v -         *
	 *  => 2=b/zF-b/zN <=> b=2*zF*zN/(zN-zF)             *
	 *  => input b into equation 2 => a=(zN-zF-2*zN)/... *
	 *      a=(zF+zN)/(zN-zF)                            *
	 *****************************************************/
	view_m[2*4+2] = (zFar+zNear)	/(zNear-zFar);	//a: zscale
	view_m[2*4+3] = (2*zFar*zNear)	/(zNear-zFar);	//b: z-translate (heavily interdependent on w=-z!)
	return;
}
///////////////////////////////////////////////////////////////////////////////
// Window has changed size, or has just been created. In either case, we need
// to use the window dimensions to set the viewport and the projection matrix.
void ChangeSize(int w, int h) {
	glViewport(0, 0, w, h);
	Screen.Width = w;
	Screen.Height = h;
	#ifdef DEBUG_SIM
		fprintf(stderr, "Changed Window Size to: Width: %u, Height: %u\nCalcProjection...",w,h);
	#endif
	CalcProjection(45.0f, 1.0f*w/h, 0.1f, 50.0f, projection_matrix);
	#ifdef DEBUG_SIM
		fprintf(stderr, "OK\n");
	#endif
	return;
}

void CalcView(const CameraData* Camera, float* mvp) {
	/*****************************************************
	 * This function rotates and translates the world to *
	 * fit our camera position                           *
	 *****************************************************
	 * new z-Axis = -(Camera.Center-Camera.Eye)          *
	 *  no minus because for Camera.Center=-10,          *
	 *  Camera.Eye=0 => z-axis = (0,0,-1) but in this    *
	 *  case z should be (0,0,1) because we look to its  *
	 *  negative end                                     *
	 * we can't use: new y-Axis = Camera.Up , because    *
	 *  then the user could specify Camera.Up to get a   *
	 *  non-orthonormal system => use only a projection  *
	 *  of Camera.Up. Start by:                          *
	 * new x-Axis = Camera.Up(kinda y-axis) x z_new      *
	 * new y-Axis = z_new cross x_new                    *
	 * => write these into the rows of a matrix to get   *
	 *  the wanted transfomration                        *
	 *     ( x-axis, 0 )   This is the case, because when*
	 *     ( y-axis, 0 )   multiplying it is like SCP    *
	 *  m= ( z-axis, 0 )   for every row in the vector   *
	 *     ( 0,0,0,  1 )                                 *
	 * Now we have the correct rotation, but we need to  *
	 * translate it using Eye. If Eye=(1,1,1) we will    *
	 * translate everything including the camera by      *
	 * (-1,-1,-1) to get our Eye to 0 like anticipated   *
	 * when calculating the perspective and Rotation!    *
	 * Multiplying our matrix to a translation-matrix    *
	 *  does the trick. Only the last column may change, *
	 *  but not easily e.g.: mvp[4*0+3]=-Camera->Eye[0]* *
	 *   *x[0]-Camera->Eye[1]*x[1]-Camera->Eye[2]*x[2]   *
	 *          ( 1,0,0, -Eye[0] )     ( x-axis, T[0] )  *
	 *          ( 0,1,0, -Eye[1] )     ( y-axis, T[1] )  *
	 *  mvp = m*( 0,0,1, -Eye[2] )  =  ( z-axis, T[2] )  *
	 *          ( 0,0,0,     1   )     ( 0,0,0,   1   )  *
	 * T[0] = <x-Axis,-Eye>; T[1] = <y-Axis,-Eye>;       *
	 *****************************************************/
	#ifdef DEBUG_SIM2
		fprintf(stderr, "Calc Projection Mat: \n");
		fprintf(stderr, "\tCamera.Eye: ");
		DebugDisplayMat(Camera->Eye, 3,1);
		fprintf(stderr, "\n\tCamera.Center: ");
		DebugDisplayMat(Camera->Center, 3,1);
		fprintf(stderr, "\n\tCamera.Up: ");
		DebugDisplayMat(Camera->Up, 3,1);
	#endif
	float *z	= malloc(3*sizeof(float));
	float *upn	= malloc(3*sizeof(float));
	float *x	= malloc(3*sizeof(float));
	float *y	= malloc(3*sizeof(float));
	//z=-(Center-Eye)/|Center-Eye|
	memcpy(z, Camera->Eye, 3*sizeof(float) );
	axpby(-1,Camera->Center,1,z,3,1);
	ScaleMat(1.0/VectorNorm(z,3), z, 3, 1);
	//x=UP/|UP| \times z
	memcpy(upn, &(Camera->Up), 3*sizeof(float) );
	ScaleMat(1/VectorNorm(upn,3), upn,3,1);
	CrossProduct(upn,z,x);
	//y=z\times x
	CrossProduct(z,x,y);
	//mvp = {{x,0},{y,0},{z,0},{0,0,0,1}}
	memset(mvp,0,16*sizeof(float));
	float T[3];
	MultiplicateMat(x,Camera->Eye,1,3,3,1,&upn); T[0]=-upn[0];
	MultiplicateMat(y,Camera->Eye,1,3,3,1,&upn); T[1]=-upn[0];
	MultiplicateMat(z,Camera->Eye,1,3,3,1,&upn); T[2]=-upn[0];
	mvp[4*0+0]=x[0];	mvp[4*0+1]=x[1];	mvp[4*0+2]=x[2];	mvp[4*0+3]=T[0];
	mvp[4*1+0]=y[0];	mvp[4*1+1]=y[1];	mvp[4*1+2]=y[2];	mvp[4*1+3]=T[1];
	mvp[4*2+0]=z[0];	mvp[4*2+1]=z[1];	mvp[4*2+2]=z[2];	mvp[4*2+3]=T[2];
	mvp[4*3+0]=0;		mvp[4*3+1]=0;		mvp[4*3+2]=0;		mvp[4*3+3]=1;
	free(upn); free(x); free(y); free(z);
	return;
}
void RenderScene(void) {
	// Clear the window with black
	glClearColor(0.0, 0.0, 0.0, 1.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
	//Installs a program object as part of current rendering state(Nothing is rendered here)
	glUseProgram(program);

	/***********************************
	 * mvp = projection * view * model *
	 ***********************************/
	float *tmp_m1=0, *tmp_m2=malloc(16*sizeof(float));
	//Translation-Matrix
	axpby(1,NULL,0,tmp_m2,4,4);			// tmp_m2=I
	tmp_m2[0*4+3] = x_transl;			// x-translation
	tmp_m2[1*4+3] = y_transl;			// y-translation
	tmp_m2[2*4+3] = z_transl;			// z-translation
	tmp_m2[3*4+3] = w_transl;			// 1/Scale

	CalcView(&Camera, view_matrix);
	CreateRotationMat4f(0,0,1, obj_rot, &tmp_m1);
	MultiplicateMat(tmp_m2,tmp_m1,4,4,4,4,&model_matrix);	//model = tmp_m2(transl)*tmp_m1(rotation)
	axpby(obj_scale,NULL,0,tmp_m1,4,4);						//tmp_m1 = scale*I
	tmp_m1[3*4+3] = 1;										//except w
	MultiplicateMat(model_matrix,tmp_m1,4,4,4,4,&tmp_m2);	//tmp_m2 = model_matrix*tmp_m1
	memcpy(model_matrix,tmp_m2,16*sizeof(float));			//model_matrix = tmp_m2

	MultiplicateMat(projection_matrix,view_matrix,4,4,4,4,&tmp_m1);	// projection * view *
	MultiplicateMat(tmp_m1,model_matrix,4,4,4,4,&mvp);				// * model = mvp;
	glUniformMatrix4fv(uniform_mvp, 1, 1, mvp);
	free(tmp_m1); free(tmp_m2);

	/******************
	 * Draw Triangles *
	 ******************/
	glUniform4f(uniform_fixed_color, 0,0,0,0);		//alpha=0 => don't use fixed_color
	//Enable all necessary attributes needed by the shaders and only disable them after calls like DrawArray!
	glEnableVertexAttribArray(attribute_v_color);
	glBindBuffer(GL_ARRAY_BUFFER, vbo_triangle_colors);
	glVertexAttribPointer(
		attribute_v_color, // attribute
		3,                 // number of elements per vertex, here (r,g,b)
		GL_FLOAT,          // the type of each element
		GL_FALSE,          // take our values as-is
		0,                 // no extra data between each position
		0                  // offset of first element
	);

	glEnableVertexAttribArray(attribute_coord3d);
	glBindBuffer(GL_ARRAY_BUFFER, vbo_vertices);
	glVertexAttribPointer(	// Describe our vertices array to OpenGL (it can't guess its format automatically)
		attribute_coord3d,	// attribute
		3,					// number of elements per vertex, here (x,y,z)
		GL_FLOAT,			// the type of each element
		GL_FALSE,			// take our values as-is
		0,					// no extra data between each position
		0					// offset of first element
	);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo_triangles);
	int size; glGetBufferParameteriv(GL_ELEMENT_ARRAY_BUFFER, GL_BUFFER_SIZE, &size);
	glDrawElements(GL_TRIANGLES, size/sizeof(GLuint), GL_UNSIGNED_INT, 0);	//Push 3*triangle_count elements starting with index 0 from the enabled array (see glEnableVertexAttribArray) from GL_ARRAY_BUFFER with ID vbo_triangle

	glDisableVertexAttribArray(attribute_coord3d);
	glDisableVertexAttribArray(attribute_v_color);

	/****************
	 * Draw Normals *
	 ****************/
	glUniform4f(uniform_fixed_color, 1,0,0,1);		//alpha!=0 => use fixed_color
	glEnableVertexAttribArray(attribute_coord3d);	//if enabled we also need to initialize it or else unpredictable errors will happen
	glBindBuffer(GL_ARRAY_BUFFER, vbo_normal_pts);	//vbo initialized in init_ressources
	glVertexAttribPointer(attribute_coord3d,3,GL_FLOAT,GL_FALSE,0,0);
	int normal_size; glGetBufferParameteriv(GL_ARRAY_BUFFER, GL_BUFFER_SIZE, &normal_size);
	glDrawArrays(GL_LINES,0,normal_size/(3*sizeof(GLfloat)));
	glDisableVertexAttribArray(attribute_coord3d);

	/* Display the result */
	glutSwapBuffers();
	#ifdef DEBUG_SIM
		fprintf(stderr, "glutSwapBuffers.\n");
	#endif
}

void DoOnIdle(void) {
//	ico_rot = 90*(sinf(glutGet(GLUT_ELAPSED_TIME)/1000.0 *(2*M_PI)/5)/ 2 + 0.5); // 0->1->0 every 5 seconds;
	/*float cur_fade = sinf(glutGet(GLUT_ELAPSED_TIME)/1000.0 *(2*M_PI)/5)/ 2 + 0.5; // 0->1->0 every 5 seconds
	glUseProgram(program);	//specifiy current programm for next OpenGL command
	glUniform1f(uniform_fade, cur_fade);*/
//	glutPostRedisplay();
}

void free_resources() {
	glDeleteProgram(program);
	glDeleteBuffers(1,&vbo_vertices);
	glDeleteBuffers(1,&ibo_triangles);
	glDeleteBuffers(1,&vbo_triangle_colors);
	fprintf(stderr, "Deleted Programm and Buffers\n");

	GLint logged_msg, max_length, max_logged_msg;
	glGetIntegerv(GL_DEBUG_LOGGED_MESSAGES, &logged_msg);
	glGetIntegerv(GL_MAX_DEBUG_MESSAGE_LENGTH, &max_length);
	glGetIntegerv(GL_MAX_DEBUG_LOGGED_MESSAGES, &max_logged_msg);
	fprintf(stderr, "Logged Debug Messages: %i, Max Debug Msg Length: %i, Max Logged Messages: %i\n",logged_msg,max_length,max_logged_msg);

	free(projection_matrix);
	free(view_matrix);
	free(model_matrix);
	free(mvp);
}

int main(int argc, char* argv[]) {
	//Prepare Log File (append - created if it doesn't exist)
	freopen("Simulation.log","a", stderr);	//returns zero if not successfull, else it returns argument file ptr
	setvbuf(stderr,0,_IONBF,0);	//with no buffer the programm can crash without executing fclose() and we still will be able to read the logs
	time_t cur_time = time(NULL);
	struct tm* ptm = gmtime(&cur_time);
	fprintf(stderr, "\n============ %04i-%02i-%02i %02i-%02i ============\n",ptm->tm_year,ptm->tm_mon,ptm->tm_mday, ptm->tm_hour,ptm->tm_min);

	//Matrix-Tests
	/*float mat_a[] = {1,2,3,5,7,11};	//3x2
	float mat_b[] = {0,4,6,8,10,12};	//2x3
	float* mat_c = 0;	//must be NULL or initialized with malloc or else realloc will crash
	MultiplicateMat(mat_a,mat_b,2,3,3,2, &mat_c);
	DebugDisplayMat(mat_c,2,2);
	MultiplicateMat(mat_a,mat_b,3,2,2,3, &mat_c);
	DebugDisplayMat(mat_c,3,3);
	#define SETVEC3(vec,x,y,z) vec[0]=x; vec[1]=y; vec[2]=z;
	SETVEC3(Camera.Eye,0,1,2);
	SETVEC3(Camera.Center,3,4,5);
	SETVEC3(Camera.Up,6,7,8);
	float* rot_mat = 0;
	CreateRotationMat3f(0,1,0, 45, &rot_mat);
	DebugDisplayMat(rot_mat,3,3);
	return 1;*/

	//gltSetWorkingDirectory(argv[0]);
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGBA|GLUT_ALPHA|GLUT_DOUBLE|GLUT_DEPTH|GLUT_STENCIL);
	glutInitWindowSize(800, 600);
	glutCreateWindow("Simulation");
	glutReshapeFunc(&ChangeSize);
	glutDisplayFunc(&RenderScene);
	glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE, GLUT_ACTION_GLUTMAINLOOP_RETURNS);
	glutSpecialFunc(&ProcessSpecialKeys);
	glutKeyboardFunc(&ProcessAsciiKeys);
	glutIdleFunc(&DoOnIdle);
	glutSetCursor(GLUT_CURSOR_NONE);	//Hide Cursor
	glutMotionFunc(&MouseWhileKeyPressed);
	glutMouseFunc(&ProcessMouseButtons);
	glutPassiveMotionFunc(&MouseWhileNoKeyPressed);

	//Just a short feel for the things GLEW will do for us :3
	fprintf(stderr, "OpenGL-Version: %s\n", glGetString(GL_VERSION));
	fprintf(stderr, "Extensions: ");
	const char *supported = NULL;
	void* wglGetExtString = wglGetProcAddress("wglGetExtensionsStringARB");	// Try To Use wglGetExtensionStringARB On Current DC, If Possible
	if (wglGetExtString)
		supported = ((char*(__stdcall*)(HDC))wglGetExtString)(wglGetCurrentDC());
	fprintf(stderr,"%s\n",supported);

	GLenum glew_status = glewInit();
	if (glew_status!=GLEW_OK) {
		fprintf(stderr, "Glew Error: %s\n", glewGetErrorString(glew_status));
		return 1;
	}

	/* When all init functions runs without errors, the program can initialise the resources */
	if (init_resources() == 1) {
		fprintf(stderr,"Entered Glut Main Loop\n");
		glutMainLoop();	//doesn't return :S
	}

	free_resources();
	return 0;
}