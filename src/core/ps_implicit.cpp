#include "../pntwrks.h"

void POINTSET::addCube(POINT base_pnt, DOMAIN size)
{
    double x = base_pnt.x+size.x/2;
    double y = base_pnt.y+size.y/2;
    double z = base_pnt.z+size.z/2;
    double width = size.x/2;
    double height = size.y/2;
    double depth=1;
    if(size.z == 0)
    {depth = 1;}
    else
    {depth = size.z/2;}
    for (int i = 1; i <= TotalPoints; i++)
    {
        phi[i] = minnum(phi[i],pow((points[i].x-x)/width,100)+pow((points[i].y-y)/height,100)+pow((points[i].z-z)/depth,6)-pow(1,100));
    }
}

void POINTSET::addSphere(POINT cntr_pnt, double a)
{
    for (int i = 1; i <= TotalPoints; i++)
    {
        phi[i] = minnum(phi[i],pow(points[i].x-cntr_pnt.x,2)+pow(points[i].y-cntr_pnt.y,2)+pow(points[i].z-cntr_pnt.z,2)-pow(a,2));
    }
}
void POINTSET::addCylinderZ(POINT cntr_pnt, double a)
{
    for (int i = 1; i <= TotalPoints; i++)
    {
        phi[i] = minnum(phi[i],pow(points[i].x-cntr_pnt.x,2)+pow(points[i].y-cntr_pnt.y,2)-pow(a,2));
    }
}
void POINTSET::addCylinderY(POINT cntr_pnt, double a)
{
    for (int i = 1; i <= TotalPoints; i++)
    {
        phi[i] = minnum(phi[i],pow(points[i].x-cntr_pnt.x,2)+pow(points[i].z-cntr_pnt.z,2)-pow(a,2));
    }
}
void POINTSET::addCylinderX(POINT cntr_pnt, double a)
{
    for (int i = 1; i <= TotalPoints; i++)
    {
        phi[i] = minnum(phi[i],pow(points[i].y-cntr_pnt.y,2)+pow(points[i].z-cntr_pnt.z,2)-pow(a,2));
    }
}
void POINTSET::addStar(POINT cntr, double a, double b, double d, double theta0, int branches)
{
    double pi = 3.14;
    int iter = (branches/2) - 1;
    double dtheta = 360/branches;
    for (int i = 1; i <= TotalPoints; i++)
    {
        double theta = theta0;
        double rad = theta0*pi/180;
        double x_p = ((points[i].x)*cos(rad)-(points[i].y)*sin(rad))-((cntr.x)*cos(rad)-(cntr.y)*sin(rad));
        double y_p = ((points[i].x)*sin(rad)+(points[i].y)*cos(rad))-((cntr.x)*sin(rad)+cntr.y*cos(rad));
        phi[i] = minnum(phi[i],(pow(x_p,2)/(a*a)+pow(y_p,2)/(b*b))-d);
        for (int j = 1; j <= iter; j++)
        {
            theta = theta+dtheta;
            rad = theta*pi/180;
            x_p = (points[i].x*cos(rad)-points[i].y*sin(rad))-(cntr.x*cos(rad)-cntr.y*sin(rad));
            y_p = (points[i].x*sin(rad)+points[i].y*cos(rad))-(cntr.x*sin(rad)+cntr.y*cos(rad));
            double phi1 = (pow(x_p,2)/(a*a)+pow(y_p,2)/(b*b))-d;
            phi[i] = minnum(phi[i],(minnum(phi1, phi[i])));
        }
    }
}

void POINTSET::subCube(POINT base_pnt, DOMAIN size)
{
    double x = base_pnt.x+size.x/2;
    double y = base_pnt.y+size.y/2;
    double z = base_pnt.z+size.z/2;
    double width = size.x/2;
    double height = size.y/2;
    double depth;
    if(size.z == 0)
    {depth = 1;}
    else
    {depth = size.z/2;}
    for (int i = 1; i <= TotalPoints; i++)
    {
        phi[i] = maxnum(phi[i],-(pow((points[i].x-x)/width,10)+pow((points[i].y-y)/height,10)+pow((points[i].z-z)/depth,10)-pow(1,10)));
    }
}
void POINTSET::subSphere(POINT cntr_pnt, double a)
{
    for (int i = 1; i <= TotalPoints; i++)
    {
        phi[i] = maxnum(phi[i],-(pow(points[i].x-cntr_pnt.x,2)+pow(points[i].y-cntr_pnt.y,2)+pow(points[i].z-cntr_pnt.z,2)-pow(a,2)));
    }
}
void POINTSET::subCylinderZ(POINT cntr_pnt, double a)
{
    for (int i = 1; i <= TotalPoints; i++)
    {
        phi[i] = maxnum(phi[i],-(pow(points[i].x-cntr_pnt.x,2)+pow(points[i].y-cntr_pnt.y,2)-pow(a,2)));
    }
}
void POINTSET::subCylinderY(POINT cntr_pnt, double a)
{
    for (int i = 1; i <= TotalPoints; i++)
    {
        phi[i] = maxnum(phi[i],-(pow(points[i].x-cntr_pnt.x,2)+pow(points[i].z-cntr_pnt.z,2)-pow(a,2)));
    }
}
void POINTSET::subCylinderX(POINT cntr_pnt, double a)
{
    for (int i = 1; i <= TotalPoints; i++)
    {
        phi[i] = maxnum(phi[i],-(pow(points[i].z-cntr_pnt.z,2)+pow(points[i].y-cntr_pnt.y,2)-pow(a,2)));
    }
}
void POINTSET::computeHyperbolicTangent(double W)
{
    for (int i = 1; i <= TotalPoints; i++)
    {
        phi[i] = tanh(phi[i]/W);
        if (phi[i] > 1)
        {phi[i] = 1;}
        if (phi[i] < -1)
        {phi[i] = -1;}
    }
}
void POINTSET::computeVOF(double W)
{
    for (int i = 1; i <= TotalPoints; i++)
    {
        phi[i] = 0.5*(1-tanh(phi[i]/W));
        if (phi[i] > 1)
        {phi[i] = 1;}
        if (phi[i] < 0)
        {phi[i] = 0;}
    }
}

