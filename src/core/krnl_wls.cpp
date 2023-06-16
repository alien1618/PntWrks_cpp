#include "../pntwrks.h"

void  KERNEL::computeWLS(POINT pnt, vector<POINT> points, vector<int> supp_pnt_nums, int TotSuppPnts, int weight_type, int order, double dmax)
{
    //--------------------------------------------------------------------------
    //cout << "CALCULATE BASES FOR POINT" << endl;
    //--------------------------------------------------------------------------
    KERNEL g_base;
    g_base.bases(order, pnt.x, pnt.y, pnt.z);
    int monomials = g_base.monomials;

    //--------------------------------------------------------------------------
    //cout << "ALLOCATE MEMORY FOR VCETORS" << endl;
    //--------------------------------------------------------------------------
    vector<vector<double> > A(monomials+1, vector<double>(monomials+1));
    vector<vector<double> > B(TotSuppPnts+1, vector<double>(monomials+1));

    N.resize(TotSuppPnts+1);
    dNdx.resize(TotSuppPnts+1);
    dNdy.resize(TotSuppPnts+1);
    dNdz.resize(TotSuppPnts+1);
    dNdxx.resize(TotSuppPnts+1);
    dNdyy.resize(TotSuppPnts+1);
    dNdzz.resize(TotSuppPnts+1);
    dNdxy.resize(TotSuppPnts+1);
    dNdxz.resize(TotSuppPnts+1);
    dNdyz.resize(TotSuppPnts+1);

    vector<vector<double> > P(monomials+1, vector<double>(TotSuppPnts+1));
    vector<vector<double> > P_trans(TotSuppPnts+1, vector<double>(monomials+1));
    vector<vector<double> > W(TotSuppPnts+1, vector<double>(TotSuppPnts+1));
    KERNEL base;
    for (int j = 1; j <= TotSuppPnts; j++)
    {
        base.bases(order, points[supp_pnt_nums[j]].x, points[supp_pnt_nums[j]].y, points[supp_pnt_nums[j]].z);
        for (int i = 1; i <= monomials; i++)
        {P[i][j] = base.N[i];}
    }
    for (int j = 1; j <= TotSuppPnts; j++)
    {
        for (int i = 1; i <= monomials; i++)
        {
            P_trans[j][i] = P[i][j];
        } //L X monomials
    }
    double c = 0.1;
    double temp = exp(-(1.1*dmax/c)*(1.1*dmax/c));
    for (int i = 1; i <= TotSuppPnts; i++)
    {
        if (weight_type == 0)
        {W[i][i] = 1;} // just typical least squares method
        else if (weight_type == 1) // weighted least squares with gaussian exponential weights
        {
            double r = pnt.distance(points[supp_pnt_nums[i]]);
            double w = exp(-(r/c)*(r/c)) - temp;
            w = w / (1 - temp);
            W[i][i] = w;
        }
        else if (weight_type == 2) //regularized spline
        {
            double eps = 0.0001;
            double dx = pnt.x-points[supp_pnt_nums[i]].x;
            double dy = pnt.y-points[supp_pnt_nums[i]].y;
            double dz = pnt.z-points[supp_pnt_nums[i]].z;

            double rx = abs(dx)/dmax;
            double ry = abs(dy)/dmax;
            double rz = abs(dz)/dmax;

            double term = (pow(eps,-2) - pow((1+eps),-2));
            double wrx = pow((pow(rx,2)+eps),-2)-pow((1+eps),-2);
            wrx = wrx/term;

            double wry = pow((pow(ry,2)+eps),-2)-pow((1+eps),-2);
            wry = wry/term;

            double wrz = pow((pow(rz,2)+eps),-2)-pow((1+eps),-2);
            wrz = wrz/term;

            W[i][i] = wrx*wry*wrz;   //w
        }
        else if (weight_type == 3) //cubic spline weights
        {
            double wx = 1;
            double wy = 1;
            double wz = 1;
            double dx = pnt.x-points[supp_pnt_nums[i]].x;
            double dy = pnt.y-points[supp_pnt_nums[i]].y;
            double dz = pnt.z-points[supp_pnt_nums[i]].z;
            double rx = abs(dx)/dmax;
            double ry = abs(dy)/dmax;
            double rz = abs(dz)/dmax;
            if (rx>0.5)
            {wx = 1.33333333-4*rx+4*rx*rx -1.333333*pow(rx,3);}
            else if (rx<=0.5)
            {wx = 0.6666667 - 4*rx*rx + 4*pow(rx,3);}

            if (ry>0.5)
            {wy = 1.3333333-4*ry+4*ry*ry -1.3333333*pow(ry,3);}
            else if (ry<=0.5)
            {wy = 0.6666667 - 4*ry*ry + 4*pow(ry,3);}

            if (rz>0.5)
            {wz = 1.3333333-4*rz+4*rz*rz -1.3333333*pow(rz,3);}
            else if (rz<=0.5)
            {wz = 0.6666667 - 4*rz*rz + 4*pow(rz,3);}

            W[i][i] = wx*wy*wz;   //w
        }
        else if(weight_type == 4) //quartic spline weights
        {
            double wx = 1;
            double wy = 1;
            double wz = 1;
            double dx = pnt.x-points[supp_pnt_nums[i]].x;
            double dy = pnt.y-points[supp_pnt_nums[i]].y;
            double dz = pnt.z-points[supp_pnt_nums[i]].z;

            double rx = abs(dx)/dmax;
            double ry = abs(dy)/dmax;
            double rz = abs(dz)/dmax;

            wx=1-6*rx*rx+8*rx*rx*rx-3*rx*rx*rx*rx;
            wy=1-6*ry*ry+8*ry*ry*ry-3*ry*ry*ry*ry;
            wz=1-6*rz*rz+8*rz*rz*rz-3*rz*rz*rz*rz;

            W[i][i] = wx*wy*wz;   //w
        }
        else if (weight_type == 5) //quintic spline weights
        {
            double wx = 1;
            double wy = 1;
            double wz = 1;
            double dx = pnt.x-points[supp_pnt_nums[i]].x;
            double dy = pnt.y-points[supp_pnt_nums[i]].y;
            double dz = pnt.z-points[supp_pnt_nums[i]].z;

            double rx = abs(dx)/dmax;
            double ry = abs(dy)/dmax;
            double rz = abs(dz)/dmax;

            wx= 1-10*rx*rx+20*rx*rx*rx-15*rx*rx*rx*rx+4*rx*rx*rx*rx*rx;
            wy= 1-10*ry*ry+20*ry*ry*ry-15*ry*ry*ry*ry+4*ry*ry*ry*ry*ry;
            wz= 1-10*rz*rz+20*rz*rz*rz-15*rz*rz*rz*rz+4*rz*rz*rz*rz*rz;
            W[i][i] = wx*wy*wz;   //w
        }
    }

#pragma omp parallel for schedule (static)
    for (int j = 1; j <= monomials; j++)
    {
        for (int k = 1; k <= TotSuppPnts; k++)
        {
            double sum = 0;
            for (int i = 1; i <= TotSuppPnts; i++)
            {sum = sum + P_trans[i][j] * W[k][i];}
            B[k][j] = sum;
        }
    }
#pragma omp parallel for schedule (static)
    for (int j = 1; j <= monomials; j++)
    {
        for (int k = 1; k <= monomials; k++)
        {
            double sum = 0;
            for (int i = 1; i <= TotSuppPnts; i++)
            {sum = sum + B[i][j] * P[k][i];}
            A[k][j] = sum;
        }
    }

    //--------------------------------------------------------------------------
    //cout <<"CALCULATE THE SHAPE FUNCTIONS AND THEIR DERIVATIVES" << endl;
    //--------------------------------------------------------------------------
    if (order >= 1)
    {
        vector<double> gamma(monomials+1);
        vector<double> dgammadx(monomials+1);
        vector<double> dgammady(monomials+1);
        vector<double> dgammadz(monomials+1);
        gamma = solve(A, g_base.N, monomials);
        dgammadx = solve(A, g_base.dNdx, monomials);
        dgammady = solve(A, g_base.dNdy, monomials);
        dgammadz = solve(A, g_base.dNdz, monomials);
        for(int j = 1; j <= TotSuppPnts; j++)
        {
            N[j] = 0;
            dNdx[j] = 0;
            dNdy[j] = 0;
            dNdz[j] = 0;
            for (int i = 1; i <= monomials; i++)
            {
                N[j] = N[j] + gamma[i]*B[j][i];
                dNdx[j] = dNdx[j] + dgammadx[i]*B[j][i];
                dNdy[j] = dNdy[j] + dgammady[i]*B[j][i];
                dNdz[j] = dNdz[j] + dgammadz[i]*B[j][i];
            }
        }
    }
    if (order >= 2)
    {
        vector<double> dgammadxx(monomials+1);
        vector<double> dgammadyy(monomials+1);
        vector<double> dgammadzz(monomials+1);
        vector<double> dgammadxy(monomials+1);
        vector<double> dgammadxz(monomials+1);
        vector<double> dgammadyz(monomials+1);
        dgammadxx = solve(A, g_base.dNdxx, monomials);
        dgammadyy = solve(A, g_base.dNdyy, monomials);
        dgammadzz = solve(A, g_base.dNdzz, monomials);
        dgammadxy = solve(A, g_base.dNdxy, monomials);
        dgammadxz = solve(A, g_base.dNdxz, monomials);
        dgammadyz = solve(A, g_base.dNdyz, monomials);
        for(int j = 1; j <= TotSuppPnts; j++)
        {
            dNdxx[j] = 0;
            dNdyy[j] = 0;
            dNdzz[j] = 0;
            dNdxy[j] = 0;
            dNdxz[j] = 0;
            dNdyz[j] = 0;
            for (int i = 1; i <= monomials; i++)
            {
                dNdxx[j] = dNdxx[j] + dgammadxx[i]*B[j][i];
                dNdyy[j] = dNdyy[j] + dgammadyy[i]*B[j][i];
                dNdzz[j] = dNdzz[j] + dgammadzz[i]*B[j][i];
                dNdxy[j] = dNdxy[j] + dgammadxy[i]*B[j][i];
                dNdxz[j] = dNdxz[j] + dgammadxz[i]*B[j][i];
                dNdyz[j] = dNdyz[j] + dgammadyz[i]*B[j][i];
            }
        }
    }
    double sum = 0;
    double sumx = 0;
    double sumy = 0;
    double sumz = 0;
    for (int i = 1; i <= TotSuppPnts; i++)
    {
        sum =sum + N[i];
        sumx = sumx + dNdx[i];
        sumy = sumy + dNdy[i];
        sumz = sumz + dNdz[i];
    }
    if (sum <= 0.9 || sum >= 1.1)
    {
        cout << "ERROR in constructed WLS interpolants. Partition of Unity is Not satisfied. sum = " << sum << endl;
        exit(0);
    }
    double eps = 1e-5;
    if (abs(sumx) >= eps || abs(sumy) >= eps || abs(sumz) >= eps)
    {
        cout << "ERROR in constructed WLS interpolants. Sum of derivatives not equal to zero. " << sumx << "\t" << sumy << "\t" << sumz << endl; exit(0);
    }

}
