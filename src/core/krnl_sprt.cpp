#include "../pntwrks.h"

void KERNEL::collectNbrsBruteForce(POINT g, vector<POINT> points0, int TotalPoints0, double radius)
{
    vector<int> c(1);
    int c_nds = 0;
    double tol = 1e-6;
    for (int iter = 1; iter <= 5; iter++)
    {
        for (int j = 1; j <= TotalPoints0; j++)
        {
            double distance = g.distance(points0[j]);
            if (distance <= radius && distance >= tol)
            {
                c_nds++;
                c.push_back(1);
                c[c_nds]=j;
            }
        }
        if (c_nds >= 3)
        {break;}
        else
        {radius = 1.5*radius;}
    }
    //ENSURE UNIQUE POINTS ARE CHOSEN
    tie(nbrs,TotalNbrs) = UniqueNums(c, c_nds);
}
void KERNEL::collectNbrsSamePhaseBruteForce(POINT g, double gphi, vector<POINT> points0, vector<double> phi, int TotalPoints0, double radius)
{
    vector<int> c(1);
    int c_nds = 0;
    double tol = 1e-6;
    for (int iter = 1; iter <= 5; iter++)
    {
        for (int j = 1; j <= TotalPoints0; j++)
        {
            double distance = g.distance(points0[j]);
            if (distance <= radius && distance >= tol)
            {
                if (gphi <= 0.5 && phi[j] <= 0.5)
                {
                    c_nds++;
                    c.push_back(1);
                    c[c_nds]=j;
                }
                        if (gphi >= 0.5 && phi[j] >= 0.5)
                {
                    c_nds++;
                    c.push_back(1);
                    c[c_nds]=j;
                }
            }
        }
        if (c_nds >= 3)
        {break;}
        else
        {radius = 1.5*radius;}
    }
    //ENSURE UNIQUE POINTS ARE CHOSEN    
    tie(nbrs, TotalNbrs) = UniqueNums(c, c_nds);
}
