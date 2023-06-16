#include "../pntwrks.h"

//********************************************************************************
// FUNCTIONS FOR COMPUTING THE EXTENDED INTERFACE VELOCITY
//********************************************************************************
vector<double> assignConstantExtendedVelocity(int TotalPoints, double b)
{
    vector<double> V_ext(TotalPoints+1);
    for (int i = 1; i <= TotalPoints; i++)
    {V_ext[i] = b;}
    return V_ext;
}
vector<double> assign4BranchExtendedVelocity(vector<POINT> points,int TotalPoints, double b)
{
    vector<double> V_ext(TotalPoints+1);
    for (int i = 1; i <= TotalPoints; i++)
    {
        double theta = atan(points[i].y/points[i].x);
        V_ext[i] = b*cos(4*theta);
    }
    return V_ext;
}
vector<double> assign6BranchExtendedVelocity(vector<POINT> points,int TotalPoints, double b)
{
    vector<double> V_ext(TotalPoints+1);
    for (int i = 1; i <= TotalPoints; i++)
    {
        double theta = atan(points[i].y/points[i].x);
        V_ext[i] = b*cos(6*theta);
    }
    return V_ext;
}

//********************************************************************************
// FUNCTIONS FOR COMPUTING THE VELOCITY FIELD IN THE DOMAIN
//********************************************************************************
tuple<vector<double>,vector<double>,vector<double> > assignSingleVortexVelocity(vector<POINT> points, int TotalPoints)
{
	vector<double> Vx(TotalPoints+1);
	vector<double> Vy(TotalPoints+1);
	vector<double> Vz(TotalPoints+1);
    for (int i = 1; i <= TotalPoints; i++)
    {
        Vx[i] = -pow(sin(3.14*points[i].x),2)*sin(2*3.14*points[i].y);
        Vy[i] =  pow(sin(3.14*points[i].y),2)*sin(2*3.14*points[i].x);
        Vz[i] =  0;
    }
    return make_tuple(Vx,Vy,Vz);
}
tuple<vector<double>,vector<double>,vector<double> > assignSingleVortexVelocity(vector<POINT> points, int TotalPoints, double u, double v, double w)
{
    vector<double> Vx(TotalPoints+1);
	vector<double> Vy(TotalPoints+1);
	vector<double> Vz(TotalPoints+1);
    double c =  0.45;
    double pi = 3.14;
    for (int i = 1; i <= TotalPoints; i++)
    {
        //swirl profile
        Vx[i] = u*(1.0-cos(c*pi*points[i].x))*pow((1.0-points[i].x),2)*(c*pi*sin(c*pi*points[i].y)*pow((1.0-points[i].y),2)-(1.0-cos(c*pi*points[i].y))*2.0*(1.0-points[i].y));
        Vy[i] =-v*(1.0-cos(c*pi*points[i].y))*pow((1.0-points[i].y),2)*(c*pi*sin(c*pi*points[i].x)*pow((1.0-points[i].x),2)-(1.0-cos(c*pi*points[i].x))*2.0*(1.0-points[i].x));
        Vz[i] = w;
    }
    return make_tuple(Vx,Vy,Vz);
}
tuple<vector<double>,vector<double>,vector<double> > assignSingleVortexVelocity2(vector<POINT> points, int TotalPoints, double u, double v, double w)
{
    vector<double> Vx(TotalPoints+1);
	vector<double> Vy(TotalPoints+1);
	vector<double> Vz(TotalPoints+1);
    for (int i = 1; i <= TotalPoints; i++)
    {
        //swirl profile
        Vx[i] = u*cos(points[i].x)*sin(points[i].y);
        Vy[i] =-v*cos(points[i].x)*sin(points[i].y);
        Vz[i] = w;
    }
    return make_tuple(Vx,Vy,Vz);
}
tuple<vector<double>,vector<double>,vector<double> > assignExtremeDeformationVelocity(vector<POINT> points, int TotalPoints)
{
    vector<double> Vx(TotalPoints+1);
	vector<double> Vy(TotalPoints+1);
	vector<double> Vz(TotalPoints+1);
    for (int i = 1; i <= TotalPoints; i++)
    {
        Vx[i] = sin(4*3.14*(points[i].x+0.5))*sin(4*3.14*(points[i].y+0.5));
        Vy[i] = cos(4*3.14*(points[i].x+0.5))*cos(4*3.14*(points[i].y+0.5));
        Vz[i] = 0;
    }
    return make_tuple(Vx,Vy,Vz);
}
tuple<vector<double>,vector<double>,vector<double> > assignRotatingAdvectiveVelocity(vector<POINT> points, int TotalPoints, double u, double v, double w)
{
    vector<double> Vx(TotalPoints+1);
	vector<double> Vy(TotalPoints+1);
	vector<double> Vz(TotalPoints+1);
    for (int i = 1; i <= TotalPoints; i++)
    {
        //simple rotation around center
        Vx[i] = u*(points[i].y);
        Vy[i] = -v*(points[i].x);
        Vz[i] = 0;
    }
    return make_tuple(Vx,Vy,Vz);
}
tuple<vector<double>,vector<double>,vector<double> > assignConstantAdvectiveVelocity(int TotalPoints, double u, double v, double w)
{
    vector<double> Vx(TotalPoints+1);
	vector<double> Vy(TotalPoints+1);
	vector<double> Vz(TotalPoints+1);
    for (int i = 1; i <= TotalPoints; i++)
    {
        Vx[i] = u;
        Vy[i] = v;
        Vz[i] = w;
    }
    return make_tuple(Vx,Vy,Vz);
}
tuple<vector<double>,vector<double>,vector<double> > assign3DSingleVortexFlow(vector<POINT> points, int TotalPoints)
{
    vector<double> Vx(TotalPoints+1);
	vector<double> Vy(TotalPoints+1);
	vector<double> Vz(TotalPoints+1);
    for (int i = 1; i <= TotalPoints; i++)
    {
        Vx[i] = 2*pow(sin(3.14*points[i].x),2)*sin(2*3.14*points[i].y)*sin(2*3.14*points[i].z);
        Vy[i] = -sin(2*3.14*points[i].x)*pow(sin(3.14*points[i].y),2)*sin(2*3.14*points[i].z);
        Vz[i] = -sin(2*3.14*points[i].x)*sin(2*3.14*points[i].y)*pow(sin(3.14*points[i].z),2);
    }
    return make_tuple(Vx,Vy,Vz);
}
