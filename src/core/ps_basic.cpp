#include "../pntwrks.h"

//----------------------------------------------------------------------
// FUNCTION DEFINITIONS FOR POINT CLASS
//----------------------------------------------------------------------
POINT::POINT()
{
}
POINT::POINT(double X , double Y , double Z)
{
    x = X;
    y = Y;
    z = Z;
}
POINT::POINT(double X , double Y)
{
    x = X;
    y = Y;
    z = 0;
}
POINT::POINT(double X)
{
    x = X;
    y = 0;
    z = 0;
}
double POINT::distance(POINT point)
{
    double dx = x-point.x;
    double dy = y-point.y;
    double dz = z-point.z;
    double distance = sqrt(dx*dx + dy*dy + dz*dz);
    return distance;
}
POINT::~POINT()
{
    //destructor
}
RESOLUTION::RESOLUTION()
{
}
RESOLUTION::RESOLUTION(int nx, int ny, int nz)
{
	if (nx < 1 || ny < 1 || nz < 1)
	{
		cout << "ERROR: resolution can not be less than 1 point per axis" << endl;
		exit(0);
	}
	else
	{
		x = nx;
		y = ny;
		z = nz;
	}
}
RESOLUTION::RESOLUTION(int nx, int ny)
{
	if (nx < 1 || ny < 1)
	{
		cout << "ERROR: resolution can not be less than 1 point per axis" << endl;
		exit(0);
	}
	else
	{
		x = nx;
		y = ny;
		z = 1;
	}
}
RESOLUTION::RESOLUTION(int nx)
{
	if (nx < 1)
	{
		cout << "ERROR: resolution can not be less than 1 point per axis" << endl;
		exit(0);
	}
	else
	{
		x = nx;
		y = 1;
		z = 1;
	}
}
RESOLUTION::~RESOLUTION()
{
    //destructor
}
//----------------------------------------------------------------------
// FUNCTION DEFINITIONS FOR DOMAIN CLASS
//----------------------------------------------------------------------
DOMAIN::DOMAIN()
{

}
DOMAIN::DOMAIN(double lx, double ly, double lz)
{
	if (lx < 0 || ly < 0 || lz < 0)
	{
		cout << "ERROR: domain size can not be less than 0" << endl;
		exit(0);
	}
	else
	{
		x = lx;
		y = ly;
		z = lz;
	}
}
DOMAIN::DOMAIN(double lx, double ly)
{
	if (lx < 0 || ly < 0)
	{
		cout << "ERROR: domain size can not be less than 0" << endl;
		exit(0);
	}
	else
	{
		x = lx;
		y = ly;
		z = 0;
	}
}
DOMAIN::DOMAIN(double lx)
{
	if (lx < 0)
	{
		cout << "ERROR: domain size can not be less than 0" << endl;
		exit(0);
	}
	else
	{
		x = lx;
		y = 0;
		z = 0;
	}
}
DOMAIN::~DOMAIN()
{
    //destructor
}

//----------------------------------------------------------------------
// FUNCTION DEFINITIONS FOR POINTSET CLASS
//----------------------------------------------------------------------
POINTSET::POINTSET() //constructor
{
	TotalPoints = 0;
    TotalLines= 0;
    TotalSurfaces = 0;
    TotalVolumes = 0;
    dim = 0;
}
POINTSET::~POINTSET()
{
    //destructor
}
void POINTSET::computeAvgPntSpacing()
{
    double sum_d = 0;
    int TotSamples = 50;
    for (int i = 1; i <= TotSamples; i++)
    {
        double d_min = 1e6;
        for (int j = 1; j <= TotSamples; j++)
        {
            double d = points[i].distance(points[j]);
            if (d <= d_min && d > 0.00001)
            {d_min = d;}
        }
        sum_d = sum_d + d_min;
    }
    AvgPntSpacing = sum_d/TotSamples;
    cout << "AvgPntSpacing = " << AvgPntSpacing << endl;
}
