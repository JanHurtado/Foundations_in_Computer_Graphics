
#include <fstream>
#include <string>
#include <iostream>
#include <GL/glut.h>
#include <GL/gl.h>
#include <stdio.h>
#include <cstdio>
#include <ctime>
#include <cstdlib>

#include "CornerTable.h"

using namespace std;

void pointGenerate(double &x, double &y)
{
    double a1 = ((rand() % 1000001) / 1000000.0) / 5.0;
    double a2 = ((rand() % 1000001) / 1000000.0) / 5.0;
    double a3 = ((rand() % 1000001) / 1000000.0) / 5.0;
    double a4 = ((rand() % 1000001) / 1000000.0) / 5.0;
    double a5 = ((rand() % 1000001) / 1000000.0) / 5.0;
    double a6 = 1.0 - (a1 + a2 + a3 + a4 + a5);

    x = a1 * (2) + a2 * (1) + a3 * (-1) + a4 * (-2) + a5 * (-1) + a6 * (1);
    y = a1 * (0) + a2 * (1.73205) + a3 * (1.73205) + a4 * (0) + a5 * (-1.73205) + a6 * (-1.73205);    
}

struct myMesh
{
	int numVertices,numFaces,numEdges;
	double * vertices;
	int * faces;
	void readOFF(const char * filename)
	{
		string OFFWord;
		int d;
		ifstream file(filename);
		file>>OFFWord>>numVertices>>numFaces>>numEdges;
		vertices = new double[numVertices*3];
		faces = new int[numFaces*3];
		for(int i=0;i<numVertices*3;i++)
			file>>vertices[i];
		int j=0;
		for(int i=0;i<numFaces*4;i+=4,j+=3)
		{
			file>>d>>faces[j]>>faces[j+1]>>faces[j+2];
		}
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
	~myMesh()
	{
		delete[] vertices;
		delete[] faces;
	}
};

myMesh * pMesh;
CornerTable * pCT;
double scaleFactor = 0.35;
double tx = 0.0;
double ty = 0.0;
double tz = 0.0;
vector<Point> * pPath;
int current_file=0;
int num_files = 6;
string filenames[] = {"Meshes/mesh2.mesh","Meshes/mesh5.mesh","Meshes/mesh10.mesh","Meshes/mesh15.mesh","Meshes/mesh20.mesh","Meshes/mesh26.mesh"};


void newMesh()
{
	delete pMesh;
	delete pCT;
	pMesh = new myMesh();
	pMesh->readOFF(filenames[current_file].c_str());
	pCT = new CornerTable(pMesh->faces,pMesh->vertices,pMesh->numFaces,pMesh->numVertices,3);
}

void newPoint()
{
	delete pPath;
	Point p; p.z=0; 
	pointGenerate(p.x,p.y);
	vector<CornerType> path;
	pCT->findPointDFS(p,path);
	pCT->getPathPoints(p,path);
	
	pPath = new vector<Point>(pCT->getPathPoints(p,path));
	cout<<"POINT: "<<p.x<<" "<<p.y<<" "<<p.z<<endl;
}

void keyboard(unsigned char key, int x, int y) {
	switch (key) {
		case '+': scaleFactor=scaleFactor*1.1; 
		break;
		case '-': scaleFactor=scaleFactor*0.9; 
		break;
		case 'w': ty=ty-0.1/scaleFactor; 
		break;
		case 's': ty=ty+0.1/scaleFactor; 
		break;
		case 'a': tx=tx+0.1/scaleFactor; 
		break;
		case 'd': tx=tx-0.1/scaleFactor; 
		break;
		case 'n': newMesh(); current_file=(current_file+1)%num_files; newPoint();
		break;
		case 'm': newPoint();
		break;
		case 27:
		exit(0);
		break;
	}
	glutPostRedisplay();
}


void display(void) 
{
	glPushMatrix();
	glScalef(scaleFactor,scaleFactor,scaleFactor);
	glTranslatef(tx,ty,tz);
	glClearColor(1.0, 1.0, 1.0, 0.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
	glBegin(GL_TRIANGLES);
	for(int i=0;i<pMesh->numFaces;i++)
	{
	    glColor3f(0.1, 0.2, 0.3);
		glVertex3f(pMesh->vertices[pMesh->faces[i*3]*3],
					pMesh->vertices[pMesh->faces[i*3]*3+1],
					pMesh->vertices[pMesh->faces[i*3]*3+2]);
		glVertex3f(pMesh->vertices[pMesh->faces[i*3+1]*3],
					pMesh->vertices[pMesh->faces[i*3+1]*3+1],
					pMesh->vertices[pMesh->faces[i*3+1]*3+2]);
		glVertex3f(pMesh->vertices[pMesh->faces[i*3+2]*3],
					pMesh->vertices[pMesh->faces[i*3+2]*3+1],
					pMesh->vertices[pMesh->faces[i*3+2]*3+2]);
	}
	glEnd();
	glColor3f(0.1, 0.1, 0.8);
	glBegin(GL_LINE_STRIP);
	for(int i=0;i<pPath->size();i++)
	{
		glVertex3f((*pPath)[i].x,(*pPath)[i].y,(*pPath)[i].z);
	}
	glEnd();
	glColor3f(0.8, 0.1, 0.1);
	glPushMatrix();
   	glTranslatef ((*pPath)[pPath->size()-1].x,(*pPath)[pPath->size()-1].y,(*pPath)[pPath->size()-1].z);
	glutWireSphere(0.01/scaleFactor,20, 20);
	glPopMatrix();
   	glFlush ();
   	glPopMatrix();
   	glFlush ();
	glutSwapBuffers();
}



int main(int argc, char **argv)
{
	
		srand(time(NULL));
		newMesh();
		newPoint();
		glutInit(&argc, argv);
		glutInitDisplayMode(GLUT_RGB | GLUT_DEPTH | GLUT_DOUBLE);
		glutInitWindowSize(512, 512);
		glutInitWindowPosition(20, 20);
		glutCreateWindow("Mesh");
		glutDisplayFunc(display);
		glutKeyboardFunc(keyboard);
		glutMainLoop();		
}
