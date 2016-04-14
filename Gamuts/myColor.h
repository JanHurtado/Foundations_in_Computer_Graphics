#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <unordered_map>
#include <cmath>

using namespace std;


typedef float TValue;
typedef float TWavelength;

const TValue defaultValue = 0.0;

struct spectrumSet
{
	size_t nSpectrums;
	size_t nSampledValues;
	unordered_map<string,size_t> names;
	vector<TWavelength> wavelengths;
	vector<vector<TValue> > data;
	TWavelength minWavelength;
	TWavelength maxWavelength;
	void setMinWavelength(){minWavelength = wavelengths[0];}
	void setMaxWavelength(){maxWavelength = wavelengths[wavelengths.size()-1];}
	void update()
	{	
		nSampledValues = wavelengths.size();
		setMinWavelength();
		setMaxWavelength();
	}
	vector<TValue> * getVectorPointer(string & spectrum_name)
	{
		return &(data[names[spectrum_name]]);
	}
	TValue getValue(size_t specPos,TWavelength wl)
	{
		if(wl<minWavelength || wl>maxWavelength)
		{
			return defaultValue;
		}
		else
		{
			return data[specPos][(int) wl-minWavelength];
		}
	}
	TValue getValue(string & spectrum_name,TWavelength wl)
	{
		return getValue(names[spectrum_name],wl);
	}
	TValue getValue(vector<TValue> & spec,TWavelength wl)
	{
		if(wl<minWavelength || wl>maxWavelength)
		{
			return defaultValue;
		}
		else
		{
			return spec[(int) wl-minWavelength];
		}
	}
	void loadFromFile(const char * filename)
	{
		ifstream file(filename);
		file>>nSpectrums>>nSampledValues;
		wavelengths.resize(nSampledValues);
		data.resize(nSpectrums);
		for(size_t i=0;i<nSpectrums;i++)
		{
			data[i].resize(nSampledValues);
			string temp;
			file>>temp;
			names[temp]=i;
		}
		for(size_t i=0;i<nSampledValues;i++)
		{
			file>>wavelengths[i];
			//cout<<wavelengths[i]<<" ";
			for(size_t j=0;j<nSpectrums;j++)
			{
				file>>data[j][i];
				//cout<<data[j][i]<<" ";
			}
			//cout<<endl;
		}
		update();
	}
	void printValues(string _name)
	{
		for(int i=0;i<nSampledValues;i++)
		{
			//cout<<data[0][i]<<" ";
			//cout<<minWavelength<<endl;
			cout<<getValue(_name,(float)i+minWavelength)<<" ";
		}
	}
};

void gamutGeneration	(spectrumSet & ss, TWavelength delta_lambda, TWavelength delta_lambda_increment, 
						string nS1, string nS2, string nS3, string nWS,
						vector<TValue> & xs,vector<TValue> & ys,vector<TValue> & zs,
						vector<TValue> & rs,vector<TValue> & gs,vector<TValue> & bs)
{
	TValue k = 0;
    TValue sum_1 = 0;
    TValue sum_2 = 0;
    TValue sum_3 = 0;   
    for (TWavelength i=ss.minWavelength;i<ss.maxWavelength;i+=delta_lambda)
    {
    	TValue white_value = ss.getValue(nWS,i);
    	//TValue white_value = 1.0; //White spectrum value
    	sum_1 += white_value*ss.getValue(nS1,i)*delta_lambda;
    	sum_2 += white_value*ss.getValue(nS2,i)*delta_lambda;
    	sum_3 += white_value*ss.getValue(nS3,i)*delta_lambda;
    } 
    TValue current_max = max(sum_1,sum_2);
    current_max = max(current_max,sum_3);
    k = (TValue)1.0/current_max;

    TValue minx = 99999; //min R
    TValue miny = 99999; //min G
    TValue minz = 99999; //min B
    TValue maxx = -99999; //max R
    TValue maxy = -99999; //max G
    TValue maxz = -99999; //max B
    
    vector<TValue> * v1 = ss.getVectorPointer(nS1);
    vector<TValue> * v2 = ss.getVectorPointer(nS2);
    vector<TValue> * v3 = ss.getVectorPointer(nS3);
    vector<TValue> * lum = ss.getVectorPointer(nWS);
    int num = 0;
    for (TWavelength i = 0;i<ss.maxWavelength-ss.minWavelength;i+=delta_lambda_increment)
    {
    	for (TWavelength ii=ss.minWavelength;ii<ss.maxWavelength;ii+=delta_lambda)
    	{
    		vector<TValue> tSpectrum(ss.nSampledValues);
    		for (TWavelength iii = ii; iii < ii + i ; iii += delta_lambda)
    		{
    			tSpectrum[((int) (iii-ss.minWavelength))%((int)ss.nSampledValues)]=1.0;
    			//cout<<((int) (iii-ss.minWavelength))%((int)ss.nSampledValues)<<" ";
    		}
    		TValue sum_1 = 0.0;
			TValue sum_2 = 0.0;
			TValue sum_3 = 0.0;
			for (TWavelength iii=ss.minWavelength;iii<ss.maxWavelength;iii+=delta_lambda) //spectrum summation
			{
				num++;
				sum_1 += ss.getValue(tSpectrum,iii)*(ss.getValue(*v1,iii))*ss.getValue(*lum,iii);
				sum_2 += ss.getValue(tSpectrum,iii)*(ss.getValue(*v2,iii))*ss.getValue(*lum,iii);
				sum_3 += ss.getValue(tSpectrum,iii)*(ss.getValue(*v3,iii))*ss.getValue(*lum,iii);
			} 
			TValue tristimulus1 = k*delta_lambda*sum_1; //R tristimulus value
			TValue tristimulus2 = k*delta_lambda*sum_2; //G tristimulus value
			TValue tristimulus3 = k*delta_lambda*sum_3; //B tristimulus value
			xs.push_back(tristimulus1);
			ys.push_back(tristimulus2);
			zs.push_back(tristimulus3);
			//maximum and minimum computations
			if (tristimulus1 > maxx)
				maxx = tristimulus1;
			if (tristimulus2 > maxy)
				maxy = tristimulus2;
			if (tristimulus3 > maxz)
				maxz = tristimulus3;	
			if (tristimulus1 < minx)
				minx = tristimulus1;
			if (tristimulus2 < miny)
				miny = tristimulus2;
			if (tristimulus3 < minz)
				minz = tristimulus3;
    	}
    }
    for (int i = 0;i<xs.size();i++)
    {
    	rs.push_back(((xs[i]-minx)/(maxx-minx))*255.0);
    	gs.push_back(((ys[i]-miny)/(maxy-miny))*255.0);
    	bs.push_back(((zs[i]-minz)/(maxz-minz))*255.0);
    }
    cout<<num<<endl;
}


void myRGB2XYZ(vector<TValue> & xs,vector<TValue> & ys,vector<TValue> & zs)
{
	for (int i=0;i<xs.size();i++)
	{
		TValue r_bar = xs[i];
		TValue g_bar = ys[i];
		TValue b_bar = zs[i];
		TValue m[3*3] = { 0.490, 0.310, 0.200,  0.177, 0.813, 0.011, 0.000, 0.010, 0.990 }; 
		//TValue m[3*3] = {2.76883,1.75171,1.13014,1.00017,4.59400,0.06216,0.00000,0.05651,5.59417}; 
		xs[i] = (TValue) (m[0]*r_bar+m[1]*g_bar+m[2]*b_bar);
		ys[i] = (TValue) (m[3]*r_bar+m[4]*g_bar+m[5]*b_bar);
		zs[i] = (TValue) (m[6]*r_bar+m[7]*g_bar+m[8]*b_bar);
	}
}

void myXYZ2RGB(vector<TValue> & xs,vector<TValue> & ys,vector<TValue> & zs)
{
	for (int i=0;i<xs.size();i++)
	{
		TValue X = xs[i];
		TValue Y = ys[i];
		TValue Z = zs[i];
		TValue m[3*3] = { 2.36440, -0.89580, -0.46770, -0.51483, 1.42523, 0.08817, 0.00520, -0.01440, 1.00921 };
		//TValue m[3*3] = {0.41847,-0.15866,-0.082835,-0.091169,0.25243,0.015708,0.0009209,-0.0025498,0.1786};
   		xs[i] = (TValue) (m[0]*X+m[1]*Y+m[2]*Z);
    	ys[i] = (TValue) (m[3]*X+m[4]*Y+m[5]*Z);
    	zs[i] = (TValue) (m[6]*X+m[7]*Y+m[8]*Z);
	}
}

void myXYZ2xyY(vector<TValue> & xs,vector<TValue> & ys,vector<TValue> & zs)
{
	for(int i=0;i<xs.size();i++)
	{
		TValue sum = (TValue)xs[i]+ys[i]+zs[i];
		zs[i]=ys[i];
    	if (sum>=1.e-6)
    	{
   			xs[i] = (TValue)(xs[i]/sum);
    		ys[i] = (TValue)(ys[i]/sum);
    	}
	}
}

void myxyY2XYZ(vector<TValue> & xs,vector<TValue> & ys,vector<TValue> & zs)
{
	for(int i=0;i<xs.size();i++)
	{
		TValue scale = (TValue)zs[i]/ys[i];
  		TValue z = 1.0-xs[i]-ys[i];
  		ys[i]=zs[i];
  		if (ys[i]>=1.e-6f)
  		{
  			xs[i] = (TValue)(xs[i]*scale);
  			zs[i] = (TValue)(z*scale);
  		}	
	}
}

TValue gamma_Lab(TValue t) {
	TValue ft, e=216.0/24389, k= 24389.0/27;
	if (t>e) 
		ft=pow(t,1.0/3); 
	else
		ft=(k*t+16)/116;
	return ft;
}

void myXYZ2Lab(vector<TValue> & xs,vector<TValue> & ys,vector<TValue> & zs)
{
	TValue ref_white1 = 1;
	TValue ref_white2 = 1;
	TValue ref_white3 = 1;
	for(int i=0;i<xs.size();i++)
	{
		TValue xr = (TValue)xs[i]/ref_white1;
	  	TValue yr = (TValue)ys[i]/ref_white2;
	  	TValue zr = (TValue)zs[i]/ref_white3;
	  	TValue fx = gamma_Lab(xr);
  		TValue fy = gamma_Lab(yr);
  		TValue fz = gamma_Lab(zr);
  		xs[i] = (TValue)(116*fy-16);
  		ys[i] = (TValue)(500*(fx-fy));
  		zs[i] = (TValue)(200*(fy-fz));
	}
}

void myLab2XYZ(vector<TValue> & xs,vector<TValue> & ys,vector<TValue> & zs)
{
	TValue ref_white1 = 1;
	TValue ref_white2 = 1;
	TValue ref_white3 = 1;
	for(int i=0;i<xs.size();i++)
	{
		TValue fx,fy,fz, fx3,fy3,fz3, xr,yr,zr;
	  	TValue e=216.0/24389;
		TValue k=24389.0/27;

		fy = (xs[i] + 16.0)/116;
		fx = ys[i]/500.0 + fy;
		fz = fy - zs[i]/200;

		fx3 = fx*fx*fx;
		fy3 = fy*fy*fy;
		fz3 = fz*fz*fz;

		xr=(fx3>e)?fx3:(116*fx-16)/k;
		yr=((TValue)xs[i]>k*e)?fy3:xs[i]/k;
		zr=(fz3>e)?fz3:(116*fz-16)/k;

		xs[i] = (TValue) (ref_white1 * xr);
		ys[i] = (TValue) (ref_white2 * yr);
		zs[i] = (TValue) (ref_white3 * zr);
	}
}

TValue gamma_sRGB(TValue x){
	TValue ft,t = (x>0)?x:-x;
	if (t>0.0031308) 
		ft = 1.055*pow(t,1.0/2.4)-0.055;
	else                  
		ft = 12.92*t;
	return (x>0)?ft:-ft;
}

void myXYZ2sRGB(vector<TValue> & xs,vector<TValue> & ys,vector<TValue> & zs)
{
	for(int i=0;i<xs.size();i++)
	{	
		TValue r,g,b;
		if (false) 
		{  
			r = 3.3921940*xs[i] - 1.6168667*ys[i] - 0.4906146*zs[i];
			g =-0.9787684*xs[i] + 1.9161415*ys[i] + 0.0334540*zs[i];
			b = 0.0719453*xs[i] - 0.2289914*ys[i] + 1.4052427*zs[i];
	  	} 
	  	else 
	  	{
			r = 3.2404542*xs[i] - 1.5371385*ys[i] - 0.4985314*zs[i];
			g =-0.9692660*xs[i] + 1.8760108*ys[i] + 0.0415560*zs[i];
			b = 0.0556434*xs[i] - 0.2040259*ys[i] + 1.0572252*zs[i];
	  	}
	  	r = gamma_sRGB(r);
	  	g = gamma_sRGB(g);
	  	b = gamma_sRGB(b);

	  	xs[i] = (TValue) r;
	  	ys[i] = (TValue) g;
	  	zs[i] = (TValue) b;
  	}
}

TValue inv_gamma_sRGB(TValue x) {
	TValue ft,t=(TValue) (x>0)?x:-x;
	if ( t > 0.04045 ) 
	   ft = pow((t+0.055)/1.055,2.4);
	else
       ft =  t/12.92;
           
    return (x>0)?ft:-ft;
}

void mysRGB2XYZ(vector<TValue> & xs,vector<TValue> & ys,vector<TValue> & zs)
{
	for(int i=0;i<xs.size();i++)
	{	
		double Rc = inv_gamma_sRGB(xs[i]);
	   	double Gc = inv_gamma_sRGB(ys[i]);
	   	double Bc = inv_gamma_sRGB(zs[i]);

	   	xs[i] = (float) (Rc*0.4124564 + Gc*0.3575761 + Bc*0.1804375);
	   	ys[i] = (float) (Rc*0.2126729 + Gc*0.7151522 + Bc*0.0721750);
	   	zs[i] = (float) (Rc*0.0193339 + Gc*0.1191920 + Bc*0.9503041);
	}
}




