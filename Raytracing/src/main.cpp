
#include <fstream>
#include <string>
#include <iostream>
#include <fstream>
#include <GL/glut.h>
#include <GL/gl.h>
#include <stdio.h>
#include <cstdio>
#include <ctime>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <omp.h>

using namespace std;

const int wsize = 512; //window size (height amd width)
const int depth_global=2; //ray depth

int mycount = 0;

struct myColor //color structure (r,g,b)
{
	float r,g,b;
	myColor()
	{
		r=0;
		g=0;
		b=0;
	}
	myColor(float _r,float _g,float _b)
	{
		r=_r;
		g=_g;
		b=_b;
	}
	myColor & operator =(const myColor & c1)
	{
		r=c1.r;
		g=c1.g;
		b=c1.b;
		return *this;
	}
	myColor operator +(myColor & c1)
	{
		myColor res;
		res.r = r + c1.r;
		res.g = g + c1.g;
		res.b = b + c1.b;
		return res;
	}
	void scalar(float s)
	{
		r=r*s;
		g=g*s;
		b=b*s;
	}
	void setColor(float _r,float _g,float _b)
	{
		r=_r;
		g=_g;
		b=_b;
	}
	bool isZero()
	{
		return r==0 && g==0 && b==0;
	}
};

struct myVector //vector structure with some operations (x,y,z)
{
	double x,y,z;
	myVector(){};
	myVector (double _x,double _y,double _z)
	{
		x=_x;
		y=_y;
		z=_z;
	}
	myVector (const myVector & p1)
	{
		x=p1.x;
		y=p1.y;
		z=p1.z;
	}
	myVector & operator =(const myVector & p1)
	{
		x=p1.x;
		y=p1.y;
		z=p1.z;
		return *this;
	}
	myVector operator +(myVector & p1)
	{
		myVector res;
		res.x = x + p1.x;
		res.y = y + p1.y;
		res.z = z + p1.z;
		return res;
	}
	myVector operator -(myVector & p1)
	{
		myVector res;
		res.x = x - p1.x;
		res.y = y - p1.y;
		res.z = z - p1.z;
		return res;
	}
	myVector operator /(const double & c) //scalar div
	{
		myVector res;
		res.x = x/c;
		res.y = y/c;
		res.z = z/c;
		return res;
	}
	bool operator ==(const myVector & v)
	{
		double tol = 0.01;
		if (x<=(v.x + tol) && x>=(v.x - tol) && y<=(v.y + tol) && y>=(v.y - tol) && z<=(v.z + tol) && z>=(v.z - tol))
			return true;
		else 
			return false; 
	}
	bool operator !=(const myVector & v)
	{
		return !((*this)==v);
	}
	myVector scalar(const double & c) // scalar product
	{
		myVector res;
		res.x = x*c;
		res.y = y*c;
		res.z = z*c;
		return res;
	}
	double operator*(const myVector & p1) // dot product
	{
		return x*p1.x+y*p1.y+z*p1.z;
	}
	myVector crossProduct(const myVector & p) // cross product
	{
		myVector res;
		res.x = y*p.z - z*p.y;
		res.y = z*p.x - x*p.z;
		res.z = x*p.y - y*p.x;
		return res;
	}
	double distance(const myVector & p1) //distance between two vectors
	{
		return sqrt(pow(x-p1.x,2.0)+pow(y-p1.y,2.0)+pow(z-p1.z,2.0));
	}
	double norm() //norm of vector
	{
		return sqrt(pow(x,2.0)+pow(y,2.0)+pow(z,2.0));
	}
	void toUnit() //transform to unit vector
	{
		double _norm = norm();
		x=x/_norm;
		y=y/_norm;
		z=z/_norm;
	}
	void print()
	{
		cout<<x<<" "<<y<<" "<<z<<" ";
	}
};

struct BarycentricCoordinate //braycentric coordinate structure (v,w,u)
{
	double v,w,u;
	void print()
	{
		cout<<v<<" "<<w<<" "<<u<<" ";
	}
	bool liesOutside() //return if barycentric coordinate is inside de triangle
	{
		return v>1 || v<0 || w>1 || w<0 || u>1 || u<0 ; 
	}
};

struct myMesh //mesh structure
{
	int numVertices,numFaces,numEdges;
	double * vertices; // x1,y1,z1,...,xn,yn,zn
	int * faces; // f11 f12 f13 f21 f22 f23 ... fn1 fn2 fn3 
	myVector * normals; // normals for each face
	void readOFF(const char * filename) // read OFF file
	{
		string OFFWord;
		int d;
		ifstream file(filename);
		file>>OFFWord>>numVertices>>numFaces>>numEdges;
		vertices = new double[numVertices*3];
		faces = new int[numFaces*3];
		normals = new myVector[numFaces];
		for(int i=0;i<numVertices*3;i++)
			file>>vertices[i];
		int j=0;
		for(int i=0;i<numFaces*4;i+=4,j+=3)
		{
			file>>d>>faces[j]>>faces[j+1]>>faces[j+2];
		}
	}
	void computeNormals() //normal computation
	{
		for(int i=0;i<numFaces;i++)
		{
			myVector p1(vertices[3*faces[3*i]],vertices[3*faces[3*i]+1],vertices[3*faces[3*i]+2]);
			myVector p2(vertices[3*faces[3*i+1]],vertices[3*faces[3*i+1]+1],vertices[3*faces[3*i+1]+2]);
			myVector p3(vertices[3*faces[3*i+2]],vertices[3*faces[3*i+2]+1],vertices[3*faces[3*i+2]+2]);
			myVector normal = ((p2-p1).crossProduct(p3-p2))/((p2-p1).crossProduct(p3-p2)).norm();
			normals[i]=normal;
			//cout<<"p1: ";p1.print();cout<<endl;
			//cout<<"p2: ";p2.print();cout<<endl;
			//cout<<"p3: ";p3.print();cout<<endl;
			//normal.print();cout<<endl;
		}
	}
	void getPoints(int t,myVector & a,myVector & b,myVector & c) //get points of a triangle
	{
		a.x = vertices[faces[t*3]*3];
		a.y = vertices[faces[t*3]*3+1];
		a.z = vertices[faces[t*3]*3+2];
		b.x = vertices[faces[t*3+1]*3];
		b.y = vertices[faces[t*3+1]*3+1];
		b.z = vertices[faces[t*3+1]*3+2];
		c.x = vertices[faces[t*3+2]*3];
		c.y = vertices[faces[t*3+2]*3+1];
		c.z = vertices[faces[t*3+2]*3+2];
	}
	void printMesh()
	{
		int j=0;
		cout<<endl<<"Vertices:"<<endl;
		for(int i=0;i<numVertices*3;i+=3)
			cout<<vertices[i]<<" "<<vertices[i+1]<<" "<<vertices[i+2]<<endl;
		cout<<endl<<"Faces:"<<endl;
		for(int i=0;i<numFaces*3;i+=3)
			cout<<faces[i]<<" "<<faces[i+1]<<" "<<faces[i+2]<<endl;
	}
	myMesh(const char * filename)
	{
		readOFF(filename);
	}
	void normalize() //normalize mesh, all values between 0 and 1
	{
		double maxvalue = -999999;
		for(int i=0;i<numVertices*3;i++)
		{
			if (abs(vertices[i])>maxvalue)
				maxvalue = abs(vertices[i]);
		}
		for(int i=0;i<numVertices*3;i++)
		{
			vertices[i]=vertices[i]/maxvalue;
		}
	}
	~myMesh()
	{
		delete[] vertices;
		delete[] faces;
	}
};

myMesh * pMesh; // global mesh value

float * pBuffer; // global buffer of pixels 

						//light color
						double lr=1.0; 
						double lg=1.0;
						double lb=1.0;
						//specular light reduction factor
						double k_sr=1.0;
						double k_sg=1.0;
						double k_sb=1.0;
						//ambient light color
						double I_ar=0.1;
						double I_ag=0.1;
						double I_ab=0.1;

//typedef myVector myVector;

struct myCamera // camera structure, name of variables equal to slides
{
	myVector eye;
	myVector center;
	myVector up;
	double fovy;
	double n,f;
	double wp,hp;
	
	myVector xe;
	myVector ye;
	myVector ze;
	double df;
	double h;
	double w;
	
	
	myCamera(myVector _eye,myVector _center,myVector _up,double _fovy,double _n,double _f,double _wp,double _hp)
	{
		eye = _eye;
		center = _center;
		up = _up;
		fovy = _fovy;
		n = _n;
		f = _f;
		wp = _wp;
		hp = _hp;
	}
	
	void update() //compute intrinsic variables
	{
		df = n;
		h = 2.0*df*tan(fovy/2.0);
		w = (wp/hp)*h;
		
		ze = (eye-center)/((eye-center).norm());
		xe = (up.crossProduct(ze))/((up.crossProduct(ze)).norm());
		ye = ze.crossProduct(xe);
		//xe.toUnit();
		//ye.toUnit();
		//ze.toUnit();
	}
	
	void print()
	{
		cout<<endl<<endl<<"*******Camera*******"<<endl<<endl;
		cout<<"eye: "; eye.print(); cout<<endl;
		cout<<"center: "; center.print(); cout<<endl;
		cout<<"up: "; up.print(); cout<<endl;
		cout<<"fovy: "; cout<<fovy; cout<<endl;
		cout<<"n: "; cout<<n; cout<<endl;
		cout<<"f: "; cout<<f; cout<<endl;
		cout<<"wp: "; cout<<wp; cout<<endl;
		cout<<"hp: "; cout<<hp; cout<<endl;
		cout<<endl; 
		cout<<"df: "; cout<<df; cout<<endl;
		cout<<"h: "; cout<<h; cout<<endl;
		cout<<"w: "; cout<<w; cout<<endl;
		cout<<"xe: "; xe.print(); cout<<endl;
		cout<<"ye: "; ye.print(); cout<<endl;
		cout<<"ze: "; ze.print(); cout<<endl;
		cout<<"*******************"<<endl;
		
	}
};

struct myRender // render structure
{
	float * buffer; //buffer of pixels
	myCamera * cam; //render camera
	vector<myMesh*> data; //objects (triangular meshes)
	vector<vector<myColor> > colors; //color of each object
	myVector lightSource; //light source position
	myColor getColor(int mesh_id,int vertex_id) //get color of a mesh vertex
	{
		int _mesh_id  = mesh_id%colors.size();
		int _vertex_id = vertex_id%colors[_mesh_id].size();
		return colors[_mesh_id][_vertex_id];
	}
	myVector getRayDirection(int x,int y) //get ray direction of a pixel
	{
		myVector v1 = ((cam->ze).scalar(-(cam->df)));
		myVector v2 = ((cam->ye).scalar((cam->h)*((double)y/cam->hp-0.5)));
		myVector v3 = (cam->xe.scalar((cam->w)*((double)x/cam->wp-0.5)));
		myVector dir = v1 + v2 + v3;
		return dir;
	}
	
	
	bool intersect(myVector & o,myVector & dir,myVector & normal,myVector & a,myVector & b,myVector & c,
			myVector & pi,BarycentricCoordinate & res) //detect intersection, returning some useful values
	{
		double ti = ((a-o)*normal)/(dir*normal);
		myVector temp = (dir.scalar(ti));
		pi = (o)+temp;
   		myVector v0 = b - a, v1 = c - a, v2 = pi - a;
    	double d00 = v0*v0;
    	double d01 = v0*v1;
    	double d11 = v1*v1;
    	double d20 = v2*v0;
    	double d21 = v2*v1;
    	double denom = d00 * d11 - d01 * d01;
    	res.v = (d11 * d20 - d01 * d21) / denom;
    	res.w = (d00 * d21 - d01 * d20) / denom;
    	res.u = 1.0 - res.v - res.w;
    	//res.print();
		//cout<<endl;
    	if (res.liesOutside())
    		return 0;
    	else
    	{
    		//cout<<"a";
			return 1;
		}
	}
	
	bool rayShadow(myVector & o,myVector & dir,int mi,int ti) // ray for shadows, given the intial point, the direction, and the triangle
	{
		double dist = lightSource.distance(o);
		for(int i=0;i<data.size();i++) //for each mesh
		{
			for(int j=0;j<data[i]->numFaces;j++) //for each face
			{
				myVector a,b,c,pi;
				data[i]->getPoints(j,a,b,c);
				int v1id = data[i]->faces[j*3];
				int v2id = data[i]->faces[j*3+1];
				int v3id = data[i]->faces[j*3+2];
				myVector normal = data[i]->normals[j];
				BarycentricCoordinate res;
				if(intersect(o,dir,normal,a,b,c,pi,res))
				{
					if(i!=mi&&j!=ti&&lightSource.distance(pi)<dist) //ignore same point
						return 1;
				}
			}
		}
		return 0;
	}
	
	myColor ray(myVector & o,myVector & dir,int depth) //trace ray returning the corresponding color
	{
		myColor color;
		double max_dist = 999999;
		if(depth<=0) return color;
		for(int i=0;i<data.size();i++) //for each mesh
		{
			for(int j=0;j<data[i]->numFaces;j++) //for each face
			{
				mycount++;
				myVector a,b,c,pi;
				data[i]->getPoints(j,a,b,c);
				int v1id = data[i]->faces[j*3];
				int v2id = data[i]->faces[j*3+1];
				int v3id = data[i]->faces[j*3+2];
				myVector normal = data[i]->normals[j];
				BarycentricCoordinate res;
				if(intersect(o,dir,normal,a,b,c,pi,res)) //for each intersection
				{
					if(cam->eye.distance(pi)<max_dist) //for each nearest intersection
					{
								
						max_dist = cam->eye.distance(pi);
					
						double k_dr = getColor(i,v1id).r*(res.u) + getColor(i,v2id).r*(res.v) + getColor(i,v3id).r*(res.w);
						double k_dg = getColor(i,v1id).g*(res.u) + getColor(i,v2id).g*(res.v) + getColor(i,v3id).g*(res.w);
						double k_db = getColor(i,v1id).b*(res.u) + getColor(i,v2id).b*(res.v) + getColor(i,v3id).b*(res.w);
						double n_specular = 50.0;

						myVector vL = lightSource-pi;
						vL.toUnit();
						double normal_vL = normal*vL;
						
						myVector vv = cam->eye - pi;
						vv.toUnit();
						double vL_normal = vL*normal;
						myVector vr = normal.scalar(2.0*vL_normal)-vL;
						double vr_vv = vr*vv;
						
						double vv_normal = vv*normal;
						myVector rr = normal.scalar(2.0*vv_normal)-vv;
								
						double ambient_r = I_ar*k_dr;
						double ambient_g = I_ag*k_dg;
						double ambient_b = I_ab*k_db;
								
						double diffuse_r = lr*k_dr*normal_vL;
						double diffuse_g = lg*k_dg*normal_vL;
						double diffuse_b = lb*k_db*normal_vL;
								
						double specular_r = lr*k_sr*pow(vr_vv,n_specular);
						double specular_g = lg*k_sg*pow(vr_vv,n_specular);
						double specular_b = lb*k_sb*pow(vr_vv,n_specular);
						
						//myVector pit = pi.scalar(0.1);
						myVector vrs = lightSource-pi;
						vrs.toUnit();
						if (depth==depth_global) //compute shadow only on the first level of trace
						{
							if(rayShadow(pi,vrs,i,j))
							{
								//cout<<"entra "<<endl;
								color.setColor(ambient_r,ambient_g,ambient_b);
							}
							else
							{
								color.setColor(ambient_r+diffuse_r+specular_r,ambient_g+diffuse_g+specular_g,ambient_b+diffuse_b+specular_b);
							}
							myColor reflection = ray(pi,rr,depth-1);
							reflection.scalar(0.2);
							color = color + reflection;
						}
						else
						{
							color.setColor(ambient_r+diffuse_r+specular_r,ambient_g+diffuse_g+specular_g,ambient_b+diffuse_b+specular_b);
							myColor reflection = ray(pi,rr,depth-1);
							reflection.scalar(0.2);
							color = color + reflection;
						}
						/*buffer[(int)(y*(cam->wp)+x)*3] = ambient_r+diffuse_r+specular_r;
						buffer[(int)(y*(cam->wp)+x)*3+1] = ambient_g+diffuse_g+specular_g;
						buffer[(int)(y*(cam->wp)+x)*3+2] = ambient_b+diffuse_b+specular_b;*/
					}
							
				}
			}
		}
		return color;
	}
	
	
	
	void run() //trace all rays (for each pixel)
	{
		int nx = (cam->wp);
		#pragma omp parallel for
		for (int x=0;x<nx;x+=1)  //for each pixel
		{
			for (int y=0;y<(cam->hp);y+=1)  //for each pixel
			{
				myVector dir = getRayDirection(x,y);
				myColor color = ray(cam->eye,dir,depth_global);
				//dir.toUnit();
				//cout<<x<<" "<<y<<" : ";dir.print(); cout<<endl;
				buffer[(int)(y*(cam->wp)+x)*3] = color.r;
				buffer[(int)(y*(cam->wp)+x)*3+1] = color.g;
				buffer[(int)(y*(cam->wp)+x)*3+2] = color.b;
			}
		}
	}
};


myRender * pRender; //ray tracing render

void renderObjects(void) //Z-Buffer objects
{
	glBegin(GL_TRIANGLES);
	for(int i=0;i<pRender->data.size();i++) //for each mesh
	{
		float temp[] = {pRender->colors[i][0].r,pRender->colors[i][0].g,pRender->colors[i][0].b};
		float white[] = {1.0,1.0,1.0};
		float black[] = {0.1,0.1,0.1};
		glMaterialfv(GL_FRONT, GL_AMBIENT, black);
		glMaterialfv(GL_FRONT, GL_DIFFUSE, temp);
  		glMaterialfv(GL_FRONT, GL_SPECULAR, white);
  		glMaterialf(GL_FRONT, GL_SHININESS, 50.0);
  		glColor3fv(temp);
  		glBindTexture(GL_TEXTURE_2D, 0);
		for(int j=0;j<pRender->data[i]->numFaces;j++) //for each face
		{
			myVector a,b,c;
			pRender->data[i]->getPoints(j,a,b,c);
			glNormal3f(pRender->data[i]->normals[j].x,pRender->data[i]->normals[j].y,pRender->data[i]->normals[j].z);
			glVertex3f(a.x,a.y,a.z);
			glVertex3f(b.x,b.y,b.z);
			glVertex3f(c.x,c.y,c.z);
		}
	}
	glEnd();
}

void displayZBuffer(void)  //Z-Buffer display
{
	glClearColor(0.0, 0.0, 0.0, 0.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glLoadIdentity();
	gluLookAt(pRender->cam->eye.x, pRender->cam->eye.y, pRender->cam->eye.z, pRender->cam->center.x, pRender->cam->center.y, pRender->cam->center.z, pRender->cam->up.x, pRender->cam->up.y, pRender->cam->up.z);
	
	float g_lightPos[3] = {pRender->lightSource.x,pRender->lightSource.y,pRender->lightSource.z};
	glLightfv(GL_LIGHT0, GL_POSITION, g_lightPos);
	
	renderObjects();
	glFlush();
	
	/*glMatrixMode(GL_MODELVIEW);
	float luzamb[3] = {0.1,0.1,0.1};
	float luzdif[3] = { 1, 1, 1 };
	float luzspec[3] = {1, 1, 1 };
	
	
	float vlightSource[3] = {-3.0,8.0,10.0};

	//glLightfv(GL_LIGHT0, GL_POSITION, vlightSource);
	//glLightfv(GL_LIGHT0, GL_AMBIENT, luzamb);
	//glLightfv(GL_LIGHT0, GL_DIFFUSE, luzdif);
	//glLightfv(GL_LIGHT0, GL_SPECULAR, luzspec);
	
	glLightfv(GL_LIGHT0, GL_DIFFUSE, luzdif);
    glLightfv(GL_LIGHT0, GL_POSITION, vlightSource);
 
    glLightfv(GL_LIGHT0, GL_SPECULAR, luzspec);
	
	
	glEnable(GL_LIGHT0);
	//glEnable(GL_COLOR_MATERIAL);
	
	float cor[] = {0.9,0.1,0.1} ;
	
	glColorMaterial(GL_FRONT, GL_DIFFUSE);
	glEnable(GL_COLOR_MATERIAL);*/
	
	
	//glColor3f(0.9, 0.1, 0.1);
	//glutWireSphere(0.5, 64, 64);
	
	//glutSwapBuffers();
	//glFlush();
	

	//glutSwapBuffers();
}

void reshape(GLint width, GLint height) //Z-Buffer reshape
{
	glViewport(0, 0, width, height);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective((pRender->cam->fovy)*(180.0/3.1416), pRender->cam->wp/pRender->cam->hp, pRender->cam->n, pRender->cam->f);
	glMatrixMode(GL_MODELVIEW);
}

void saveImage() //save current image in ppm format
{
	float * pixels = new float[ 3 * wsize * wsize];
	glReadPixels(0, 0, wsize, wsize, GL_RGB, GL_FLOAT, pixels);
	ofstream file("output.ppm");
	file<<"P3"<<endl<<wsize<<" "<<wsize<<endl<<255<<endl;
	for(int i=wsize-1;i>=0;i--)
	{
		for(int j=0;j<3*wsize;j++)
			file<<(int)(pixels[i*(3*wsize)+j]*255)<<" ";
		file<<endl;
	}
}

bool flag = 0; //flag for Z-Buffer and raytracing control

void keyboard(unsigned char key, int x, int y) {
	switch (key) {
		case 'q': flag=0; 
		break;
		case 'w': flag=1; cout<<"w";
		break;
		case 's': saveImage(); 
		break;
		case 27:
		exit(0);
		break;
	}
	glutPostRedisplay();
}


void display(void) //display function regarding flag
{
	if (flag==0)
		displayZBuffer();
	else
	{
		glClearColor(1.0, 1.0, 1.0, 0.0);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		glDrawPixels(wsize,wsize,GL_RGB,GL_FLOAT,pBuffer);
		glutSwapBuffers();
	}
}

void init(void) //Z-Buffer initialization
{
	glEnable(GL_DEPTH_TEST);
   	glDepthFunc(GL_LESS);
   	glShadeModel(GL_FLAT);
   	glEnable(GL_LIGHTING);
   	glEnable(GL_LIGHT0);
   	float amb[3] = {0.1,0.1,0.1};
	float dif[3] = { 1, 1, 1 };
	float spe[3] = {1, 1, 1 };

	glLightfv(GL_LIGHT0, GL_AMBIENT, amb);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, dif);
	glLightfv(GL_LIGHT0, GL_SPECULAR, spe);
}

 
int main(int argc, char **argv)
{
		pBuffer = new float [wsize*wsize*3];
		myVector eye(1.0,0.5,1);//view1
		//myVector eye(-0.5,1.3,1.5);//view2
		//myVector eye(0,0.5,1);//view3
		//myVector eye(0.5,2.5,1.5);//view4
		
		
		//myVector eye(0.0,0.5,2.0);
		myVector center(0,0,0);
		myVector up(0,1,0);
		//myVector up(0.7071,0.7071,-5);
		//(eye+center+up).print();
		double fovy = 3.1416/2.0;
		double n = 0.1;
		double f = 100;
		double wp=wsize;
		double hp=wsize;
		//load camera
		myCamera cam(eye,center,up,fovy,n,f,wp,hp);
		cam.update();
		cam.print();
		//load meshes
		pMesh = new myMesh("200.off");
		pMesh->computeNormals();
		myMesh * pMesh2 = new myMesh("plane2.off");
		pMesh2->computeNormals();
		myMesh * pMesh3 = new myMesh("plane.off");
		pMesh3->computeNormals();
		//initialize mesh colors
		vector<myColor> col;
		col.push_back(myColor(0.8,0.4,0.1));
		vector<myColor> col2;
		col2.push_back(myColor(0.2,0.2,0.7));
		vector<myColor> col3;
		col3.push_back(myColor(0.2,0.2,0.7));
		//myVector light(6.0,12.0,4.0);
		//initialize light position
		myVector light(-3,8,10);
		//set render variables
		myRender * render = new myRender;
		render->cam = &cam;
		render->data.push_back(pMesh);
		render->data.push_back(pMesh2);
		render->data.push_back(pMesh3);
		render->buffer = pBuffer;
		render->colors.push_back(col);
		render->colors.push_back(col2);
		render->colors.push_back(col3);
		render->lightSource = light;
		pRender = render;
		//run render
		render->run();
		cout<<endl<<"Esto: "<<mycount<<endl;
		glutInit(&argc, argv);
		glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB | GLUT_DEPTH);
		glutInitWindowSize(wsize, wsize);
		glutInitWindowPosition(20, 20);
		glutCreateWindow("Mesh");
		init();
		glutDisplayFunc(display);
		glutReshapeFunc(reshape);
		glutKeyboardFunc(keyboard);
		glutMainLoop();
}
