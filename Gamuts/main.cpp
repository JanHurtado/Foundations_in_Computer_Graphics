
//#include "myColor.h"
#include "myVisualization.h"

int main()
{
	spectrumSet ss;
	ss.loadFromFile("result");
	vector<float> v1;
	vector<float> v2;
	vector<float> v3;
	vector<float> cr;
	vector<float> cg;
	vector<float> cb;
	gamutGeneration(ss,1,1,"R","G","B","1",v1,v2,v3,cr,cg,cb);
	myRGB2XYZ(v1,v2,v3);
	myXYZ2RGB(v1,v2,v3);
	//mysRGB2XYZ(v1,v2,v3);
	//myXYZ2RGB(v1,v2,v3);
	//myXYZ2Lab(v1,v2,v3);
	//myLab2XYZ(v1,v2,v3);
	//myXYZ2xyY(v1,v2,v3);
	//myxyY2XYZ(v1,v2,v3);
	initialize(v1,v2,v3,cr,cg,cb);
	runVis();
	//showPointCloud(v1,v2,v3,cr,cg,cb);
}
