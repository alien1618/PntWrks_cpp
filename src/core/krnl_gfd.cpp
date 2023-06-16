#include "../pntwrks.h"

vector<double> computeGFD(double pnt1x, double pnt1y, double U1, vector<POINT> points, vector<double> U, vector<int> supp_pnt_nums, int TotSuppPnts, double dm)
{
    vector<double> b1(6);
    vector<double> b2(6);
    vector<double> b(6);
    vector<double> D(6);
    vector<vector<double> > A(6, vector<double> (6));
    for (int m = 1; m <= 5; m++)
    {
        for (int n = 1; n <= 5; n++)
        {A[m][n] = 0;}
        b[m] = 0;
        b1[m] = 0;
        b2[m] = 0;
    }
    for (int s = 1; s <= TotSuppPnts; s++)
    {
        int nbr = supp_pnt_nums[s];
        double dx = pnt1x-points[nbr].x;
        double dy = pnt1y-points[nbr].y;
        double d = sqrt(dx*dx+dy*dy);
        double w = 0;
        if (d <= 0.5*dm)
        {w = 0.66667-4.0*pow((d/dm),2)-4.0*pow((d/dm),3);}
        else if (d > 0.5*dm && d <= dm)
        {w = 1.333333-4.0*(d/dm)+4.0*pow((d/dm),2)-1.33333*pow((d/dm),3);}
        double h = points[nbr].x-pnt1x;
        double k = points[nbr].y-pnt1y;

        A[1][1] = A[1][1] + pow(w,2)*pow(h,2);
        A[2][2] = A[2][2] + pow(w,2)*pow(k,2);
        A[3][3] = A[3][3] + pow(w,2)*0.25*pow(h,4);
        A[4][4] = A[4][4] + pow(w,2)*0.25*pow(k,4);
        A[5][5] = A[5][5] + pow(w,2)*0.25*pow(h,2)*pow(k,2);

        A[2][1] = A[2][1] + pow(w,2)*h*k;
        A[1][2] = A[1][2] + pow(w,2)*h*k;

        A[3][1] = A[3][1] + pow(w,2)*0.5*pow(h,3);
        A[1][3] = A[1][3] + pow(w,2)*0.5*pow(h,3);

        A[4][1] = A[4][1] + pow(w,2)*0.5*pow(k,2)*h;
        A[1][4] = A[1][4] + pow(w,2)*0.5*pow(k,2)*h;

        A[5][1] = A[5][1] + pow(w,2)*pow(h,2)*k;
        A[1][5] = A[1][5] + pow(w,2)*pow(h,2)*k;

        A[3][2] = A[3][2] + pow(w,2)*0.5*pow(h,2)*k;
        A[2][3] = A[2][3] + pow(w,2)*0.5*pow(h,2)*k;

        A[4][2] = A[4][2] + pow(w,2)*0.5*pow(k,3);
        A[2][4] = A[2][4] + pow(w,2)*0.5*pow(k,3);

        A[5][2] = A[5][2] + pow(w,2)*pow(k,2)*h;
        A[2][5] = A[2][5] + pow(w,2)*pow(k,2)*h;

        A[4][3] = A[4][3] + pow(w,2)*0.25*pow(h,2)*pow(k,2);
        A[3][4] = A[3][4] + pow(w,2)*0.25*pow(h,2)*pow(k,2);

        A[5][3] = A[5][3] + pow(w,2)*0.5*pow(h,3)*k;
        A[3][5] = A[3][5] + pow(w,2)*0.5*pow(h,3)*k;

        A[5][4] = A[5][4] + pow(w,2)*0.5*pow(k,3)*h;
        A[4][5] = A[4][5] + pow(w,2)*0.5*pow(k,3)*h;

        b1[1] = b1[1] + pow(w,2)*h;
        b1[2] = b1[2] + pow(w,2)*k;
        b1[3] = b1[3] + pow(w,2)*0.5*pow(h,2);
        b1[4] = b1[4] + pow(w,2)*0.5*pow(k,2);
        b1[5] = b1[5] + pow(w,2)*h*k;

        b2[1] = b2[1] + U[nbr]*pow(w,2)*h;
        b2[2] = b2[2] + U[nbr]*pow(w,2)*k;
        b2[3] = b2[3] + U[nbr]*pow(w,2)*0.5*pow(h,2);
        b2[4] = b2[4] + U[nbr]*pow(w,2)*0.5*pow(k,2);
        b2[5] = b2[5] + U[nbr]*pow(w,2)*h*k;
    }
    for (int p = 1; p <= 5; p++)
    {b[p] = -U1*b1[p] + b2[p];}

    D = solve(A,b,5); // D stores dUdx, dUdy, dUdxx, dUdyy
    return D;
}
