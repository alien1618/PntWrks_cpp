#include "../pntwrks.h"

void  KERNEL::computeMLS(POINT pnt, vector<POINT> points, vector<int> supp_pnt_nums, int TotSuppPnts, int weight_type, int order, double dmax)
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
    vector<double> p(monomials+1);
    vector<vector<double> > pp(monomials+1, vector<double>(monomials+1));

    vector<vector<double> > A(monomials+1, vector<double>(monomials+1));
    vector<vector<double> > dAdx(monomials+1, vector<double>(monomials+1));
    vector<vector<double> > dAdy(monomials+1, vector<double>(monomials+1));
    vector<vector<double> > dAdz(monomials+1, vector<double>(monomials+1));
    vector<vector<double> > dAdxx(monomials+1, vector<double>(monomials+1));
    vector<vector<double> > dAdyy(monomials+1, vector<double>(monomials+1));
    vector<vector<double> > dAdzz(monomials+1, vector<double>(monomials+1));
    vector<vector<double> > dAdxy(monomials+1, vector<double>(monomials+1));
    vector<vector<double> > dAdxz(monomials+1, vector<double>(monomials+1));
    vector<vector<double> > dAdyz(monomials+1, vector<double>(monomials+1));

    vector<vector<double> > B(TotSuppPnts+1, vector<double>(monomials+1));
    vector<vector<double> > dBdx(TotSuppPnts+1, vector<double>(monomials+1));
    vector<vector<double> > dBdy(TotSuppPnts+1, vector<double>(monomials+1));
    vector<vector<double> > dBdz(TotSuppPnts+1, vector<double>(monomials+1));
    vector<vector<double> > dBdxx(TotSuppPnts+1, vector<double>(monomials+1));
    vector<vector<double> > dBdyy(TotSuppPnts+1, vector<double>(monomials+1));
    vector<vector<double> > dBdzz(TotSuppPnts+1, vector<double>(monomials+1));
    vector<vector<double> > dBdxy(TotSuppPnts+1, vector<double>(monomials+1));
    vector<vector<double> > dBdxz(TotSuppPnts+1, vector<double>(monomials+1));
    vector<vector<double> > dBdyz(TotSuppPnts+1, vector<double>(monomials+1));

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

    vector<vector<double> > weights(TotSuppPnts+1, vector<double>(monomials+1));

    vector<double> p_dAdx_gamma(monomials+1);
    vector<double> p_dAdy_gamma(monomials+1);
    vector<double> p_dAdz_gamma(monomials+1);

    vector<double> term_xx(monomials+1);
    vector<double> term_yy(monomials+1);
    vector<double> term_zz(monomials+1);
    vector<double> term_xy(monomials+1);
    vector<double> term_xz(monomials+1);
    vector<double> term_yz(monomials+1);

    vector<double> gamma(monomials+1);
    vector<double> dgammadx(monomials+1);
    vector<double> dgammady(monomials+1);
    vector<double> dgammadz(monomials+1);
    vector<double> dgammadxx(monomials+1);
    vector<double> dgammadyy(monomials+1);
    vector<double> dgammadzz(monomials+1);
    vector<double> dgammadxy(monomials+1);
    vector<double> dgammadxz(monomials+1);
    vector<double> dgammadyz(monomials+1);


    //--------------------------------------------------------------------------
    //cout << "CALCULATE DISTANCE BETWEEN GAUSS POINT AND SUPPORT POINTS" << endl;
    //--------------------------------------------------------------------------
    //cout << "weight_type " << weight_type << endl;
    if (weight_type == 1)
    {weights = MLS_Spline3(pnt, points, supp_pnt_nums, TotSuppPnts, dmax);}
    else if (weight_type == 2)
    {weights = MLS_Spline4(pnt, points, supp_pnt_nums, TotSuppPnts, dmax);}
    else if (weight_type == 3)
    {weights = MLS_Spline5(pnt, points, supp_pnt_nums, TotSuppPnts, dmax);}
    else if (weight_type == 4)
    {weights = MLS_RegSpline4(pnt, points, supp_pnt_nums, TotSuppPnts, dmax);}
    else
    {
        cout << "ERROR in INTERPOLATION_FUNCTIONS::computeMLS_interp" << endl;
        cout << "MLS Weight_type " << weight_type << " is undefined." << endl; exit(0);
    }
    //--------------------------------------------------------------------------
    //cout <<  "CALCULATE THE A AND B MATRICES" << endl;
    //--------------------------------------------------------------------------
    for (int j = 1; j <= monomials; j++)
    {
        for (int i = 1; i <= monomials; i++)
        {
            A[i][j] = 0;
            dAdx[i][j] = 0;
            dAdy[i][j] = 0;
            dAdz[i][j] = 0;
            dAdxx[i][j] = 0;
            dAdyy[i][j] = 0;
            dAdzz[i][j] = 0;
            dAdxy[i][j] = 0;
            dAdxz[i][j] = 0;
            dAdyz[i][j] = 0;
        }
    }

    KERNEL base;
    for (int i = 1; i <= TotSuppPnts; i++)
    {
        base.bases(order, points[supp_pnt_nums[i]].x, points[supp_pnt_nums[i]].y, points[supp_pnt_nums[i]].z);
        p = base.N;
        //pp = multiply_transp(p,p,monomials);
        for (int n = 1; n <= monomials; n++)
        {
            for (int m = 1; m <= monomials; m++)
            {pp[m][n] = p[m] * p[n];}
        }

        if (order >= 1)
        {
            double w = weights[i][1];
            double dwdx = weights[i][2];
            double dwdy = weights[i][3];
            double dwdz = weights[i][4];

#pragma omp parallel for
            for (int n = 1; n <= monomials; n++)
            {
#pragma omp parallel for
                for (int m = 1; m <= monomials; m++)
                {
                    A[m][n] = A[m][n] + pp[m][n] * w;
                    dAdx[m][n] = dAdx[m][n] + pp[m][n] * dwdx;
                    dAdy[m][n] = dAdy[m][n] + pp[m][n] * dwdy;
                    dAdz[m][n] = dAdz[m][n] + pp[m][n] * dwdz;
                }
            }
            for (int k = 1; k <= monomials; k++)
            {
                B[i][k] = p[k]* w;
                dBdx[i][k] = p[k]* dwdx;
                dBdy[i][k] = p[k]* dwdy;
                dBdz[i][k] = p[k]* dwdz;
            }

        }
        if (order >= 2)
        {
            double dwdxx = weights[i][5];
            double dwdyy = weights[i][6];
            double dwdzz = weights[i][7];
            double dwdxy = weights[i][8];
            double dwdxz = weights[i][9];
            double dwdyz = weights[i][10];

#pragma omp parallel for
            for (int n = 1; n <= monomials; n++)
            {
                for (int m = 1; m <= monomials; m++)
                {
                    dAdxx[m][n] = dAdxx[m][n] + pp[m][n] * dwdxx;
                    dAdyy[m][n] = dAdyy[m][n] + pp[m][n] * dwdyy;
                    dAdzz[m][n] = dAdzz[m][n] + pp[m][n] * dwdzz;
                    dAdxy[m][n] = dAdxy[m][n] + pp[m][n] * dwdxy;
                    dAdxz[m][n] = dAdxz[m][n] + pp[m][n] * dwdxz;
                    dAdyz[m][n] = dAdyz[m][n] + pp[m][n] * dwdyz;
                }
            }
#pragma omp parallel for
            for (int k = 1; k <= monomials; k++)
            {

                dBdxx[i][k] = p[k] * dwdxx;
                dBdyy[i][k] = p[k] * dwdyy;
                dBdzz[i][k] = p[k] * dwdzz;
                dBdxy[i][k] = p[k] * dwdxy;
                dBdxz[i][k] = p[k] * dwdxz;
                dBdyz[i][k] = p[k] * dwdyz;
            }

        }
    }

    //--------------------------------------------------------------------------
    //cout <<"CALCULATE THE SHAPE FUNCTIONS AND THEIR DERIVATIVES" << endl;
    //--------------------------------------------------------------------------
    gamma = solve(A, g_base.N, monomials);
    if (order >= 1)
    {
        for(int n = 1; n<= monomials ; n++)
        {
            double dAdx_gamma = 0;
            double dAdy_gamma = 0;
            double dAdz_gamma = 0;
            for(int m = 1; m<= monomials ; m++)
            {
                dAdx_gamma = dAdx_gamma + (dAdx[m][n] * gamma[m]);
                dAdy_gamma = dAdy_gamma + (dAdy[m][n] * gamma[m]);
                dAdz_gamma = dAdz_gamma + (dAdz[m][n] * gamma[m]);
            }
            p_dAdx_gamma[n] = g_base.dNdx[n]-dAdx_gamma;
            p_dAdy_gamma[n] = g_base.dNdy[n]-dAdy_gamma;
            p_dAdz_gamma[n] = g_base.dNdz[n]-dAdz_gamma;
        }

        dgammadx = solve(A, p_dAdx_gamma, monomials);
        dgammady = solve(A, p_dAdy_gamma, monomials);
        dgammadz = solve(A, p_dAdz_gamma, monomials);

#pragma omp parallel for
        for(int n = 1; n <= TotSuppPnts ; n++)
        {
            double gamma_B = 0;
            double dgammadx_B = 0;
            double dgammady_B = 0;
            double dgammadz_B = 0;
            double gamma_dBdx = 0;
            double gamma_dBdy = 0;
            double gamma_dBdz = 0;
            for(int m = 1; m <= monomials ; m++)
            {
                gamma_B = gamma_B + (gamma[m] * B[n][m]);
                dgammadx_B = dgammadx_B + (dgammadx[m] * B[n][m]);
                dgammady_B = dgammady_B + (dgammady[m] * B[n][m]);
                dgammadz_B = dgammadz_B + (dgammadz[m] * B[n][m]);
                gamma_dBdx = gamma_dBdx + (gamma[m] * dBdx[n][m]);
                gamma_dBdy = gamma_dBdy + (gamma[m] * dBdy[n][m]);
                gamma_dBdz = gamma_dBdz + (gamma[m] * dBdz[n][m]);
            }
            N[n] = gamma_B;
            dNdx[n] = dgammadx_B+gamma_dBdx;
            dNdy[n] = dgammady_B+gamma_dBdy;
            dNdz[n] = dgammadz_B+gamma_dBdz;
        }
    }

    if (order >= 2)
    {
#pragma omp parallel for
        for (int n = 1; n <= monomials; n++)
        {
            double dAdxx_gamma = 0;
            double dAdyy_gamma = 0;
            double dAdzz_gamma = 0;
            double dAdxy_gamma = 0;
            double dAdxz_gamma = 0;
            double dAdyz_gamma = 0;
            double dAdx_dgammadx = 0;
            double dAdx_dgammady = 0;
            double dAdx_dgammadz = 0;
            double dAdy_dgammadx = 0;
            double dAdy_dgammady = 0;
            double dAdy_dgammadz = 0;
            double dAdz_dgammadx = 0;
            double dAdz_dgammady = 0;
            double dAdz_dgammadz = 0;
            for (int m = 1; m <= monomials; m++)
            {
                dAdxx_gamma = dAdxx_gamma + dAdxx[m][n]*gamma[m];
                dAdyy_gamma = dAdyy_gamma + dAdyy[m][n]*gamma[m];
                dAdzz_gamma = dAdzz_gamma + dAdzz[m][n]*gamma[m];
                dAdxy_gamma = dAdxy_gamma + dAdxy[m][n]*gamma[m];
                dAdxz_gamma = dAdxz_gamma + dAdxz[m][n]*gamma[m];
                dAdyz_gamma = dAdyz_gamma + dAdyz[m][n]*gamma[m];

                dAdx_dgammadx = dAdx_dgammadx + dAdx[m][n]*dgammadx[m];
                dAdx_dgammady = dAdx_dgammady + dAdx[m][n]*dgammady[m];
                dAdx_dgammadz = dAdx_dgammadz + dAdx[m][n]*dgammadz[m];

                dAdy_dgammadx = dAdy_dgammadx + dAdy[m][n]*dgammadx[m];
                dAdy_dgammady = dAdy_dgammady + dAdy[m][n]*dgammady[m];
                dAdy_dgammadz = dAdy_dgammadz + dAdy[m][n]*dgammadz[m];

                dAdz_dgammadx = dAdz_dgammadx + dAdz[m][n]*dgammadx[m];
                dAdz_dgammady = dAdz_dgammady + dAdz[m][n]*dgammady[m];
                dAdz_dgammadz = dAdz_dgammadz + dAdz[m][n]*dgammadz[m];
            }
            term_xx[n] = g_base.dNdxx[n] - (dAdx_dgammadx + dAdx_dgammadx + dAdxx_gamma);
            term_yy[n] = g_base.dNdyy[n] - (dAdy_dgammady + dAdy_dgammady + dAdyy_gamma);
            term_zz[n] = g_base.dNdzz[n] - (dAdz_dgammadz + dAdz_dgammadz + dAdzz_gamma);
            term_xy[n] = g_base.dNdxy[n] - (dAdx_dgammady + dAdy_dgammadx + dAdxy_gamma);
            term_xz[n] = g_base.dNdxz[n] - (dAdx_dgammadz + dAdz_dgammadx + dAdxz_gamma);
            term_yz[n] = g_base.dNdyz[n] - (dAdy_dgammadz + dAdz_dgammady + dAdyz_gamma);
        }
        dgammadxx = solve(A, term_xx, monomials);
        dgammadyy = solve(A, term_yy, monomials);
        dgammadzz = solve(A, term_zz, monomials);
        dgammadxy = solve(A, term_xy, monomials);
        dgammadxz = solve(A, term_xz, monomials);
        dgammadyz = solve(A, term_yz, monomials);

#pragma omp parallel for
        for(int n = 1; n <= TotSuppPnts ; n++)
        {
            double dgammadxx_B = 0;
            double dgammadyy_B = 0;
            double dgammadzz_B = 0;
            double dgammadxy_B = 0;
            double dgammadxz_B = 0;
            double dgammadyz_B = 0;
            double gamma_dBdxx = 0;
            double gamma_dBdyy = 0;
            double gamma_dBdzz = 0;
            double gamma_dBdxy = 0;
            double gamma_dBdxz = 0;
            double gamma_dBdyz = 0;
            double dgammadx_dBdx = 0;
            double dgammadx_dBdy = 0;
            double dgammadx_dBdz = 0;
            double dgammady_dBdx = 0;
            double dgammady_dBdy = 0;
            double dgammady_dBdz = 0;
            double dgammadz_dBdx = 0;
            double dgammadz_dBdy = 0;
            double dgammadz_dBdz = 0;
            for(int m = 1; m <= monomials ; m++)
            {
                dgammadxx_B = dgammadxx_B + (dgammadxx[m] * B[n][m]);
                dgammadyy_B = dgammadyy_B + (dgammadyy[m] * B[n][m]);
                dgammadzz_B = dgammadzz_B + (dgammadzz[m] * B[n][m]);
                dgammadxy_B = dgammadxy_B + (dgammadxy[m] * B[n][m]);
                dgammadxz_B = dgammadxz_B + (dgammadxz[m] * B[n][m]);
                dgammadyz_B = dgammadyz_B + (dgammadyz[m] * B[n][m]);
                gamma_dBdxx = gamma_dBdxx + (gamma[m] * dBdxx[n][m]);
                gamma_dBdyy = gamma_dBdyy + (gamma[m] * dBdyy[n][m]);
                gamma_dBdzz = gamma_dBdzz + (gamma[m] * dBdzz[n][m]);
                gamma_dBdxy = gamma_dBdxy + (gamma[m] * dBdxy[n][m]);
                gamma_dBdxz = gamma_dBdxz + (gamma[m] * dBdxz[n][m]);
                gamma_dBdyz = gamma_dBdyz + (gamma[m] * dBdyz[n][m]);
                dgammadx_dBdx = dgammadx_dBdx + (dgammadx[m] * dBdx[n][m]);
                dgammadx_dBdy = dgammadx_dBdy + (dgammadx[m] * dBdy[n][m]);
                dgammadx_dBdz = dgammadx_dBdz + (dgammadx[m] * dBdz[n][m]);
                dgammady_dBdx = dgammady_dBdx + (dgammady[m] * dBdx[n][m]);
                dgammady_dBdy = dgammady_dBdy + (dgammady[m] * dBdy[n][m]);
                dgammady_dBdz = dgammady_dBdz + (dgammady[m] * dBdz[n][m]);
                dgammadz_dBdx = dgammadz_dBdx + (dgammadz[m] * dBdx[n][m]);
                dgammadz_dBdy = dgammadz_dBdy + (dgammadz[m] * dBdy[n][m]);
                dgammadz_dBdz = dgammadz_dBdz + (dgammadz[m] * dBdz[n][m]);

            }
            dNdxx[n] = dgammadxx_B + dgammadx_dBdx + dgammadx_dBdx + gamma_dBdxx;
            dNdyy[n] = dgammadyy_B + dgammady_dBdy + dgammady_dBdy + gamma_dBdyy;
            dNdzz[n] = dgammadzz_B + dgammadz_dBdz + dgammadz_dBdz + gamma_dBdzz;
            dNdxy[n] = dgammadxy_B + dgammadx_dBdy + dgammady_dBdx + gamma_dBdxy;
            dNdxz[n] = dgammadxz_B + dgammadx_dBdz + dgammadz_dBdx + gamma_dBdxz;
            dNdyz[n] = dgammadyz_B + dgammady_dBdz + dgammadz_dBdy + gamma_dBdyz;
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

    if (sum <= 0.85 || sum >= 1.2)
    {
        cout << "ERROR in constructed MLS interpolants. Partition of Unity is Not satisfied. sum = " << sum << endl;
        exit(0);
    }
    double eps = 1e-5;
    if (abs(sumx) >= eps || abs(sumy) >= eps || abs(sumz) >= eps)
    {
        cout << "ERROR in constructed MLS interpolants. Sum of derivatives not equal to zero. " << sumx << "\t" << sumy << "\t" << sumz << endl; exit(0);
    }
}
vector<vector<double> > KERNEL::MLS_Spline3(POINT pnt, vector<POINT> points, vector<int> supp_pnt_nums, int TotSuppPnts, double dmax)
{
    //---------------------------------------------------------------
    // CUBIC SPLINE WEIGHT FUNCTION
    //---------------------------------------------------------------
    double wx=0, dwx=0, dwxx = 0, wy = 0, dwy=0, dwyy=0, wz=0, dwz=0, dwzz=0;
    vector<vector<double> > w(TotSuppPnts+1, vector<double>(11));
    for (int i=1; i <= TotSuppPnts; i++)
    {
        double dx = pnt.x-points[supp_pnt_nums[i]].x;
        double dy = pnt.y-points[supp_pnt_nums[i]].y;
        double dz = pnt.z-points[supp_pnt_nums[i]].z;
        double drdx = double(sign <double> (dx))/dmax;
        double drdy = double(sign <double> (dy))/dmax;
        double drdz = double(sign <double> (dz))/dmax;
        double rx = abs(dx)/dmax;
        double ry = abs(dy)/dmax;
        double rz = abs(dz)/dmax;
        if (rx>0.5)
        {
            wx = 1.33333333-4*rx+4*rx*rx -1.333333*pow(rx,3);
            dwx = (-4 + 8*rx-4*pow(rx,2))*drdx;
            dwxx = (8-8*rx)*drdx*drdx;
        }
        else if (rx<=0.5)
        {
            wx = 0.6666667 - 4*rx*rx + 4*pow(rx,3);
            dwx = (-8*rx + 12*pow(rx,2))*drdx;
            dwxx=(-8+24*rx)*drdx*drdx;
        }

        if (ry>0.5)
        {
            wy = 1.3333333-4*ry+4*ry*ry -1.3333333*pow(ry,3);
            dwy = (-4 + 8*ry-4*pow(ry,2))*drdy;
            dwyy=(8-8*ry)*drdy*drdy;
        }
        else if (ry<=0.5)
        {
            wy = 0.6666667 - 4*ry*ry + 4*pow(ry,3);
            dwy = (-8*ry + 12*pow(ry,2))*drdy;
            dwyy=(-8+24*ry)*drdy*drdy;
        }
        if (rz>0.5)
        {
            wz = 1.3333333-4*rz+4*rz*rz -1.3333333*pow(rz,3);
            dwz = (-4 + 8*rz-4*pow(rz,2))*drdz;
            dwzz=(8-8*rz)*drdz*drdz;
        }
        else if (rz<=0.5)
        {
            wz = 0.6666667 - 4*rz*rz + 4*pow(rz,3);
            dwz = (-8*rz + 12*pow(rz,2))*drdz;
            dwzz=(-8+24*rz)*drdz*drdz;
        }

        w[i][1] = wx*wy*wz;   //w
        w[i][2] = wy*wz*dwx;  //dwdx
        w[i][3] = wx*wz*dwy;  //dwdy
        w[i][4] = wx*wy*dwz;  //dwdz
        w[i][5] = wy*wz*dwxx;   //dwdxx
        w[i][6] = wx*wz*dwyy;   //dwdyy
        w[i][7] = wx*wy*dwzz;   //dwdzz
        w[i][8] = dwx*dwy*wz;   //dwdxdwdy
        w[i][9] = dwx*wy*dwz;   //dwdxdwdz
        w[i][10]=wx*dwy*dwz;   //dwdydwdz
    }
    return w;
}

vector<vector<double> > KERNEL::MLS_Spline4(POINT pnt, vector<POINT> points, vector<int> supp_pnt_nums, int TotSuppPnts, double dmax)
{
    //---------------------------------------------------------------
    //QUARTIC SPLINES WEIGHT FUNCTIONS
    //---------------------------------------------------------------
    double wx, dwx, dwxx, wy, dwy, dwyy, wz, dwz, dwzz;
    vector<vector<double> > w(TotSuppPnts+1, vector<double>(11));
    for (int i=1; i <= TotSuppPnts; i++)
    {
        double dx = pnt.x-points[supp_pnt_nums[i]].x;
        double dy = pnt.y-points[supp_pnt_nums[i]].y;
        double dz = pnt.z-points[supp_pnt_nums[i]].z;

        double drdx = double(sign <double> (dx))/dmax;
        double drdy = double(sign <double> (dy))/dmax;
        double drdz = double(sign <double> (dz))/dmax;

        double rx = abs(dx)/dmax;
        double ry = abs(dy)/dmax;
        double rz = abs(dz)/dmax;

        wx=1-6*rx*rx+8*rx*rx*rx-3*rx*rx*rx*rx;
        dwx = (-12*rx+24*rx*rx-12*rx*rx*rx)*drdx;
        dwxx=(-12+48*rx-36*rx*rx)*drdx*drdx;

        wy=1-6*ry*ry+8*ry*ry*ry-3*ry*ry*ry*ry;
        dwy = (-12*ry+24*ry*ry-12*ry*ry*ry)*drdy;
        dwyy=(-12+48*ry-36*ry*ry)*drdy*drdy;

        wz=1-6*rz*rz+8*rz*rz*rz-3*rz*rz*rz*rz;
        dwz = (-12*rz+24*rz*rz-12*rz*rz*rz)*drdz;
        dwzz=(-12+48*rz-36*rz*rz)*drdz*drdz;

        w[i][1] = wx*wy*wz;   //w
        w[i][2] = wy*wz*dwx;  //dwdx
        w[i][3] = wx*wz*dwy;  //dwdy
        w[i][4] = wx*wy*dwz;  //dwdz
        w[i][5]=wy*wz*dwxx;   //dwdxx
        w[i][6]=wx*wz*dwyy;   //dwdyy
        w[i][7]=wx*wy*dwzz;   //dwdzz
        w[i][8]=dwx*dwy*wz;   //dwdxdwdy
        w[i][9]=dwx*wy*dwz;   //dwdxdwdz
        w[i][10]=wx*dwy*dwz;   //dwdydwdz

    }
    return w;
}
vector<vector<double> > KERNEL::MLS_Spline5(POINT pnt, vector<POINT> points, vector<int> supp_pnt_nums, int TotSuppPnts, double dmax)
{
    //---------------------------------------------------------------
    //FIFTH ORDER SPLINES WEIGHT FUNCTIONS
    //---------------------------------------------------------------
    vector<vector<double> > w(TotSuppPnts+1, vector<double>(11));
    for (int i=1; i <= TotSuppPnts; i++)
    {
        double dx = pnt.x-points[supp_pnt_nums[i]].x;
        double dy = pnt.y-points[supp_pnt_nums[i]].y;
        double dz = pnt.z-points[supp_pnt_nums[i]].z;

        double drdx = double(sign <double> (dx))/dmax;
        double drdy = double(sign <double> (dy))/dmax;
        double drdz = double(sign <double> (dz))/dmax;

        double rx = abs(dx)/dmax;
        double ry = abs(dy)/dmax;
        double rz = abs(dz)/dmax;

        double wx= 1-10*rx*rx+20*rx*rx*rx-15*rx*rx*rx*rx+4*rx*rx*rx*rx*rx;
        double dwx = (-20*rx+60*rx*rx-60*rx*rx*rx+20*rx*rx*rx*rx)*drdx;
        double dwxx =(-20+120*rx-180*rx*rx+80*rx*rx*rx)*drdx*drdx;

        double wy= 1-10*ry*ry+20*ry*ry*ry-15*ry*ry*ry*ry+4*ry*ry*ry*ry*ry;
        double dwy = (-20*ry+60*ry*ry-60*ry*ry*ry+20*ry*ry*ry*ry)*drdy;
        double dwyy =(-20+120*ry-180*ry*ry+80*ry*ry*ry)*drdy*drdy;

        double wz= 1-10*rz*rz+20*rz*rz*rz-15*rz*rz*rz*rz+4*rz*rz*rz*rz*rz;
        double dwz = (-20*rz+60*rz*rz-60*rz*rz*rz+20*rz*rz*rz*rz)*drdz;
        double dwzz =(-20+120*rz-180*rz*rz+80*rz*rz*rz)*drdz*drdz;

        w[i][1] = wx*wy*wz;   //w
        w[i][2] = wy*wz*dwx;  //dwdx
        w[i][3] = wx*wz*dwy;  //dwdy
        w[i][4] = wx*wy*dwz;  //dwdz
        w[i][5]=wy*wz*dwxx;   //dwdxx
        w[i][6]=wx*wz*dwyy;   //dwdyy
        w[i][7]=wx*wy*dwzz;   //dwdzz
        w[i][8]=dwx*dwy*wz;   //dwdxdwdy
        w[i][9]=dwx*wy*dwz;   //dwdxdwdz
        w[i][10]=wx*dwy*dwz;   //dwdydwdz
    }
    return w;
}
vector<vector<double> > KERNEL::MLS_RegSpline4(POINT pnt, vector<POINT> points, vector<int> supp_pnt_nums, int TotSuppPnts, double dmax)
{
    //REGULARIZED QUARTIC SPLINES WEIGHT FUNCTIONS
    vector<vector<double> > w(TotSuppPnts+1, vector<double>(11));
    double eps = 1e-3;

    for (int i=1; i <= TotSuppPnts; i++)
    {
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

        double dwrx, dwry, dwrz, dwrxx, dwryy, dwrzz;

        double term1_x = pow((pow(rx,2) + eps),-3);
        double term1_y = pow((pow(ry,2) + eps),-3);
        double term1_z = pow((pow(rz,2) + eps),-3);

        double term2_x = pow((pow(rx,2) + eps),-4);
        double term2_y = pow((pow(ry,2) + eps),-4);
        double term2_z = pow((pow(rz,2) + eps),-4);

        double term3 = pow(eps,-2)-pow((1+eps),-2);

        if (dx >= 1e-8)
        {
            dwrx = 4*term1_x*rx/(dmax*term3);
            dwrxx = -24*(pow(rx,2)/pow(dmax,2))*term2_x - (4/pow(dmax,2))*term1_x;
            dwrxx = dwrxx/term3;
        }
        else
        {
            dwrx = 4*pow(eps,2)*(1/dmax);
            dwrxx = 20*pow(eps,2)*(1/pow(dmax,2));
        }

        if (dy >= 1e-8)
        {
            dwry = 4*term1_y*ry/(dmax*term3);
            dwryy = -24*(pow(ry,2)/pow(dmax,2))*term2_y + (4/pow(dmax,2))*term1_y;
            dwryy = dwryy/term3;
        }
        else
        {
            dwry = 4*pow(eps,2)*(1/dmax);
            dwryy = 20*pow(eps,2)*(1/pow(dmax,2));
        }

        if (dz >= 1e-8)
        {
            dwrz = 4*term1_z*rz/(dmax*term3);
            dwrzz = -24*(pow(rz,2)/pow(dmax,2))*term2_z + (4/pow(dmax,2))*term1_z;
            dwrzz = dwrzz/term3;
        }
        else
        {
            dwrz = 4*pow(eps,2)*(1/dmax);
            dwrzz = 20*pow(eps,2)*(1/pow(dmax,2));
        }

        w[i][1] = wrx*wry*wrz;   //w
        w[i][2] = wry*wrz*dwrx;  //dwdx
        w[i][3] = wrx*wrz*dwry;  //dwdy
        w[i][4] = wrx*wry*dwrz;  //dwdz
        w[i][5]=wry*wrz*dwrxx;   //dwdxx
        w[i][6]=wrx*wrz*dwryy;   //dwdyy
        w[i][7]=wrx*wry*dwrzz;   //dwdzz
        w[i][8]=dwrx*dwry*wrz;   //dwdxdwdy
        w[i][9]=dwrx*wry*dwrz;   //dwdxdwdz
        w[i][10]=wrx*dwry*dwrz;   //dwdydwdz
    }
    return w;
}
