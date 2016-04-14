#include <pcl/visualization/cloud_viewer.h>
#include <pcl/io/io.h>
#include <pcl/io/pcd_io.h>
#include <pcl/io/ifs_io.h>
#include <pcl/point_types.h>
#include <vector>
#include "myColor.h"

vector<float> xs;
vector<float> ys;
vector<float> zs;
vector<float> rs;
vector<float> gs;
vector<float> bs;
		
pcl::visualization::PCLVisualizer::Ptr viewer;

enum nodes{e_rgb=0,e_xyz,e_xyy,e_lab,e_srgb};

int currentNode = 0;

int transitions[5][5] {{0,1,0,0,0},{1,0,1,1,1},{0,1,0,0,0},{0,1,0,0,0},{0,1,0,0,0}};
	
void initialize(vector<float> & _xs,vector<float> & _ys,vector<float> & _zs,
				vector<float> & _rs,vector<float> & _gs,vector<float> & _bs)
{
	xs.resize(_xs.size());
	ys.resize(_xs.size());
	zs.resize(_xs.size());
	rs.resize(_xs.size());
	gs.resize(_xs.size());
	bs.resize(_xs.size());
	for(int i=0;i<_xs.size();i++)
	{
		xs[i]=_xs[i];
		ys[i]=_ys[i];
		zs[i]=_zs[i];
		rs[i]=_rs[i];
		gs[i]=_gs[i];
		bs[i]=_bs[i];
	}
}
							
pcl::PointXYZRGB genWhitePoint(float _x,float _y,float _z) //Function for PCL white point generation
{
	pcl::PointXYZRGB p(255,255,255);
	p.x = _x;
	p.y = _y;
	p.z = _z;
	return p;
}

void setCube()
{
	viewer->addLine(genWhitePoint(0,0,0),genWhitePoint(1,0,0),"l1");
	viewer->addLine(genWhitePoint(0,0,0),genWhitePoint(0,1,0),"l2");
	viewer->addLine(genWhitePoint(0,0,0),genWhitePoint(0,0,1),"l3");
	viewer->addLine(genWhitePoint(1,0,0),genWhitePoint(1,1,0),"l4");
	viewer->addLine(genWhitePoint(1,0,0),genWhitePoint(1,0,1),"l5");
	viewer->addLine(genWhitePoint(0,1,0),genWhitePoint(0,1,1),"l6");
	viewer->addLine(genWhitePoint(0,1,0),genWhitePoint(1,1,0),"l7");
	viewer->addLine(genWhitePoint(0,0,1),genWhitePoint(0,1,1),"l8");
	viewer->addLine(genWhitePoint(0,0,1),genWhitePoint(1,0,1),"l9");
	viewer->addLine(genWhitePoint(1,0,1),genWhitePoint(1,1,1),"l10");
	viewer->addLine(genWhitePoint(1,1,0),genWhitePoint(1,1,1),"l11");
	viewer->addLine(genWhitePoint(0,1,1),genWhitePoint(1,1,1),"l12"); 
}

void updatePointCloud()
{
	viewer->removePointCloud("main");
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud (new pcl::PointCloud<pcl::PointXYZRGB>);
	for (int i=0;i<xs.size();i++)
	{
		pcl::PointXYZRGB o(rs[i],gs[i],bs[i]); //point RGB colors
		//point coordinates...
		o.x = xs[i];
		o.y = ys[i];
		o.z = zs[i];
		cloud->push_back(o);
	}
	viewer->addPointCloud(cloud,"main");
	viewer->resetCamera();
}

void updateLabel(string s)
{
	//viewer->removeText("label");
	viewer->updateText(s,10,50,"label");
}

void keyboardEventOccurred (const pcl::visualization::KeyboardEvent &event,
	                        void* viewer_void)
{
	boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer = *static_cast<boost::shared_ptr<pcl::visualization::PCLVisualizer> *> (viewer_void);
	if (event.getKeySym () == "t" && event.keyDown ())
	{
		if(transitions[currentNode][e_rgb])
		{
			myXYZ2RGB(xs,ys,zs);
			updatePointCloud();
			currentNode = e_rgb;
			updateLabel("RGB");
		}
		else
			updateLabel("Try another color space");
	}
	if (event.getKeySym () == "y" && event.keyDown ())
	{
		if(transitions[currentNode][e_xyz])
		{
			if(currentNode==e_rgb)
			{
				myRGB2XYZ(xs,ys,zs);
				updatePointCloud();
				currentNode = e_xyz;
			}
			if(currentNode==e_xyy)
			{
				myxyY2XYZ(xs,ys,zs);
				updatePointCloud();
				currentNode = e_xyz;
			}
			if(currentNode==e_lab)
			{
				myLab2XYZ(xs,ys,zs);
				updatePointCloud();
				currentNode = e_xyz;
			}
			if(currentNode==e_srgb)
			{
				mysRGB2XYZ(xs,ys,zs);
				updatePointCloud();
				currentNode = e_xyz;
			}
			updateLabel("XYZ");
		}
		else
			updateLabel("Try another color space");
	}
	if (event.getKeySym () == "v" && event.keyDown ())
	{
		if(transitions[currentNode][e_xyy])
		{
			myXYZ2xyY(xs,ys,zs);
			updatePointCloud();
			currentNode = e_xyy;
			updateLabel("xyY");
		}
		else
			updateLabel("Try another color space");
	}
	if (event.getKeySym () == "b" && event.keyDown ())
	{
		if(transitions[currentNode][e_lab])
		{
			myXYZ2Lab(xs,ys,zs);
			updatePointCloud();
			currentNode = e_lab;
			updateLabel("Lab");
		}
		else
			updateLabel("Try another color space");
	}
	if (event.getKeySym () == "n" && event.keyDown ())
	{
		if(transitions[currentNode][e_srgb])
		{
			myXYZ2sRGB(xs,ys,zs);
			updatePointCloud();
			currentNode = e_srgb;
			updateLabel("sRGB");
		}
		else
			updateLabel("Try another color space");
	}
}

void runVis()
{
	//PCL visualizer
	viewer.reset(new pcl::visualization::PCLVisualizer);
	updatePointCloud();
	viewer->addText("t: RGB     y: XYZ     v: xyY     b: Lab     n: sRGB",10,20);
	viewer->addText("RGB",10,50,"label");
 	setCube(); // Draw cube
 	viewer->registerKeyboardCallback (keyboardEventOccurred, (void*)&viewer);
	viewer->addCoordinateSystem(); //Draw axis
	viewer->setBackgroundColor (0, 0, 0); //set background color
	viewer->spin(); //run visualizer
}

