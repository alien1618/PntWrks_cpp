#include "../pntwrks.h"

void KERNEL::computeRBF(POINT pnt, vector<POINT> points, vector<int> supp_pnt_nums, int TotSuppPnts, int RBF_type, int RBF_order, double dmax, double alfc)
{
    if (TotSuppPnts <= 0)
    {
        cout << "ERROR: insufficient support TotalPoints for interpolation" << endl;
        cout << "sum_TotalPoints = " << TotSuppPnts << endl;
        exit(0);
    }

    KERNEL g_base;
    g_base.bases(RBF_order, pnt.x, pnt.y, pnt.z);
    int monomials = g_base.monomials;

    int mg = TotSuppPnts+monomials;
    vector<vector<double> > A(mg+1, vector<double>(mg+1));
    vector<vector<double> > RBF_mat(mg+1, vector<double>(11));
    vector<vector<double> > weights(mg+1, vector<double>(11));

    vector<double> rk(mg+1);
    vector<double> drkdx(mg+1);
    vector<double> drkdy(mg+1);
    vector<double> drkdz(mg+1);
    vector<double> drkdxx(mg+1);
    vector<double> drkdyy(mg+1);
    vector<double> drkdzz(mg+1);
    vector<double> drkdxy(mg+1);
    vector<double> drkdxz(mg+1);
    vector<double> drkdyz(mg+1);

    for (int j = 1; j <= mg; j++)
    {
        for (int i = 1; i <= mg; i++)
        {A[i][j] = 0;}
    }

    //cout << " assembling matrix A " << endl;
    for (int i = 1; i <= TotSuppPnts; i++)
    {
        for (int k=1; k <= TotSuppPnts; k++)
        {
            double dx = points[supp_pnt_nums[i]].x-points[supp_pnt_nums[k]].x;
            double dy = points[supp_pnt_nums[i]].y-points[supp_pnt_nums[k]].y;
            double dz = points[supp_pnt_nums[i]].z-points[supp_pnt_nums[k]].z;
            double r= dx*dx+dy*dy+dz*dz;
            if (RBF_type == 1)
            {
                //Multi-quartics
                double q = -0.5;
                double rc = alfc*dmax;
                A[i][k]=pow((r+rc*rc),q);
            }
            else   if (RBF_type == 2)
            {
                double qc = alfc/dmax/dmax;
                //Gaussian exponential
                A[i][k] = exp(-qc*r);
            }
            else if (RBF_type == 3)
            {
                // Thin plate spline
                double q = alfc;
                A[i][k] =pow(r,(0.5*q));
            }
        }
        if (monomials > 0)
        {
            g_base.bases(RBF_order, points[supp_pnt_nums[i]].x, points[supp_pnt_nums[i]].y, points[supp_pnt_nums[i]].z);
            for (int j = 1; j <= monomials; j++)
            {
                A[i][TotSuppPnts+j] = g_base.N[j];
                A[TotSuppPnts+j][i] = g_base.N[j];
            }
        }
    }

    //cout << "solve the linear equations to get shape functions" << endl;
    if (RBF_type == 1)
    {RBF_mat = RBF_MQ(pnt, points, supp_pnt_nums, TotSuppPnts, dmax, alfc, RBF_order);}
    else if (RBF_type == 2)
    {RBF_mat = RBF_GE(pnt, points, supp_pnt_nums, TotSuppPnts, dmax, alfc, RBF_order);}
    else if (RBF_type == 3)
    {RBF_mat = RBF_TPS(pnt, points, supp_pnt_nums, TotSuppPnts, dmax, alfc, RBF_order);}
    for (int j = 1; j <= mg; j++)
    {
        rk[j] = RBF_mat[j][1];
        drkdx[j] = RBF_mat[j][2];
        drkdy[j] = RBF_mat[j][3];
        drkdz[j] = RBF_mat[j][4];
        drkdxx[j] = RBF_mat[j][5];
        drkdyy[j] = RBF_mat[j][6];
        drkdzz[j] = RBF_mat[j][7];
        drkdxy[j] = RBF_mat[j][8];
        drkdxz[j] = RBF_mat[j][9];
        drkdyz[j] = RBF_mat[j][10];
    }
 
    if (RBF_order >= 1)
    {
        N = solve(A,rk,mg);
        dNdx = solve(A,drkdx,mg);
        dNdy = solve(A,drkdy,mg);
        dNdz = solve(A,drkdz,mg);
        dNdxx = rk;
        dNdyy = rk;
        dNdzz = rk;
        dNdxy = rk;
        dNdxz = rk;
        dNdyz = rk;
    }
    if (RBF_order >= 2)
    {
        dNdxx = solve(A,drkdxx,mg);
        dNdyy = solve(A,drkdyy,mg);
        dNdzz = solve(A,drkdzz,mg);
        dNdxy = solve(A,drkdxy,mg);
        dNdxz = solve(A,drkdxz,mg);
        dNdyz = solve(A,drkdyz,mg);
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
        cout << "ERROR in constructed RBF interpolants. Partition of Unity is Not satisfied. sum = " << sum << endl;
        exit(0);
    }
    double eps = 1e-5;
    if (abs(sumx) >= eps || abs(sumy) >= eps || abs(sumz) >= eps)
    {
        cout << "ERROR in constructed RBF interpolants. Sum of derivatives not equal to zero. " << sumx << "\t" << sumy << "\t" << sumz << endl; exit(0);
    }
}
vector<vector<double> > KERNEL::RBF_MQ(POINT pnt, vector<POINT> points, vector<int> supp_pnt_nums, int TotSuppPnts, double dmax, double alfc, int RBF_order)
{
    int monomials;
    KERNEL g_base;
    g_base.bases(RBF_order, pnt.x, pnt.y, pnt.z);
    monomials = g_base.monomials;
    vector<vector<double> > RBF(TotSuppPnts+monomials+1, vector<double>(11));
    for(int j = 1; j <= 10; j++)
    {
        for (int i = 1; i <= TotSuppPnts+monomials; i++)
        {RBF[i][j] = 0;}
    }
    for (int i=1; i <= TotSuppPnts; i++)
    {
        double dx = pnt.x-points[supp_pnt_nums[i]].x;
        double dy = pnt.y-points[supp_pnt_nums[i]].y;
        double dz = pnt.z-points[supp_pnt_nums[i]].z;
        double rr= dx*dx+dy*dy+dz*dz;

        //Multi-quartics
        double q = -0.5;
        double rcrc = alfc*dmax*alfc*dmax;
        RBF[i][1]=pow((rr+rcrc),q);                                                            //W
        RBF[i][2]=2*q*pow((rr+rcrc),(q-1.))*dx;                                                //dWdx
        RBF[i][3]=2*q*pow((rr+rcrc),(q-1.))*dy;                                                //dWdy
        RBF[i][4]=2*q*pow((rr+rcrc),(q-1.))*dz;                                                //dWdz
        RBF[i][5]=2*q*pow((rr+rcrc),(q-1.))+4.*(q-1.)*q*dx*dx*pow((rr+rcrc),(q-2.));     //dWdxx
        RBF[i][6]=2*q*pow((rr+rcrc),(q-1.))+4.*(q-1.)*q*dy*dy*pow((rr+rcrc),(q-2.));     //dWdyy
        RBF[i][7]=2*q*pow((rr+rcrc),(q-1.))+4.*(q-1.)*q*dz*dz*pow((rr+rcrc),(q-2.));     //dWdzz
        RBF[i][8]=4*q*(q-1)*pow((rr+rcrc),(q-2.))*dx*dy;     //dWdxy
        RBF[i][9]=4*q*(q-1)*pow((rr+rcrc),(q-2.))*dx*dz;     //dWdxz
        RBF[i][10]=4*q*(q-1)*pow((rr+rcrc),(q-2.))*dy*dz;    //dWdyz
    }

    if (monomials > 0)
    {
        for (int k=1; k <= monomials; k++)
        {
            RBF[TotSuppPnts+k][1]=g_base.N[k];
            RBF[TotSuppPnts+k][2]=g_base.dNdx[k];
            RBF[TotSuppPnts+k][3]=g_base.dNdy[k];
            RBF[TotSuppPnts+k][4]=g_base.dNdz[k];
            RBF[TotSuppPnts+k][5]=g_base.dNdxx[k];
            RBF[TotSuppPnts+k][6]=g_base.dNdyy[k];
            RBF[TotSuppPnts+k][7]=g_base.dNdzz[k];
            RBF[TotSuppPnts+k][8]=g_base.dNdxy[k];
            RBF[TotSuppPnts+k][9]=g_base.dNdxz[k];
            RBF[TotSuppPnts+k][10]=g_base.dNdyz[k];
        }
    }
    return RBF;
}
vector<vector<double> > KERNEL::RBF_GE(POINT pnt, vector<POINT> points, vector<int> supp_pnt_nums, int TotSuppPnts, double dmax, double alfc, int RBF_order)
{
    int monomials;
    KERNEL g_base;
    g_base.bases(RBF_order, pnt.x, pnt.y, pnt.z);
    monomials = g_base.monomials;
    vector<vector<double> > RBF(TotSuppPnts+monomials+1, vector<double>(11));

    for(int j = 1; j <= 10; j++)
    {
        for (int i = 1; i <= TotSuppPnts+monomials; i++)
        {RBF[i][j] = 0;}
    }

    for (int i=1; i <= TotSuppPnts; i++)
    {
        double dx = pnt.x-points[supp_pnt_nums[i]].x;
        double dy = pnt.y-points[supp_pnt_nums[i]].y;
        double dz = pnt.z-points[supp_pnt_nums[i]].z;
        double rr= dx*dx+dy*dy+dz*dz;

        // Gaussian exp
        double qc = alfc/dmax/dmax;
        RBF[i][1]=exp(-qc*rr);                                                            //W
        RBF[i][2]=-2*qc*exp(-qc*rr)*dx;                                                   //dWdx
        RBF[i][3]=-2*qc*exp(-qc*rr)*dy;                                                   //dWdy
        RBF[i][4]=-2*qc*exp(-qc*rr)*dz;                                                   //dWdz
        RBF[i][5]=-2*qc*exp(-qc*rr)+4*qc*qc*exp(-qc*rr)*dx*dx;                       //dWdxx
        RBF[i][6]=-2*qc*exp(-qc*rr)+4*qc*qc*exp(-qc*rr)*dy*dy;                       //dWdyy
        RBF[i][7]=-2*qc*exp(-qc*rr)+4*qc*qc*exp(-qc*rr)*dz*dz;                       //dWdzz
        RBF[i][8]=4*qc*qc*exp(-qc*rr)*dx*dy;                       //dWdxy
        RBF[i][9]=4*qc*qc*exp(-qc*rr)*dx*dz;                       //dWdxz
        RBF[i][10]=4*qc*qc*exp(-qc*rr)*dy*dz;                      //dWdyz
    }

    if (monomials > 0)
    {
        for (int k=1; k <= monomials; k++)
        {
            RBF[TotSuppPnts+k][1]=g_base.N[k];
            RBF[TotSuppPnts+k][2]=g_base.dNdx[k];
            RBF[TotSuppPnts+k][3]=g_base.dNdy[k];
            RBF[TotSuppPnts+k][4]=g_base.dNdz[k];
            RBF[TotSuppPnts+k][5]=g_base.dNdxx[k];
            RBF[TotSuppPnts+k][6]=g_base.dNdyy[k];
            RBF[TotSuppPnts+k][7]=g_base.dNdzz[k];
            RBF[TotSuppPnts+k][8]=g_base.dNdxy[k];
            RBF[TotSuppPnts+k][9]=g_base.dNdxz[k];
            RBF[TotSuppPnts+k][10]=g_base.dNdyz[k];
        }
    }
    return RBF;
}
vector<vector<double> > KERNEL::RBF_TPS(POINT pnt, vector<POINT> points, vector<int> supp_pnt_nums, int TotSuppPnts, double dmax, double alfc, int RBF_order)
{
    int monomials;
    KERNEL g_base;
    g_base.bases(RBF_order, pnt.x, pnt.y, pnt.z);
    monomials = g_base.monomials;
    vector<vector<double> > RBF(TotSuppPnts+monomials+1, vector<double>(11));

    for(int j = 1; j <= 10; j++)
    {
        for (int i = 1; i <= TotSuppPnts+monomials; i++)
        {RBF[i][j] = 0;}
    }

    for (int i=1; i <= TotSuppPnts; i++)
    {
        double dx = pnt.x-points[supp_pnt_nums[i]].x;
        double dy = pnt.y-points[supp_pnt_nums[i]].y;
        double dz = pnt.z-points[supp_pnt_nums[i]].z;
        double num = dx*dx+dy*dy+dz*dz;
        double r= sqrt(num);
        double r1 = pow(num,-0.5);
        double r2 = 1/pow(num,-1.5);
        double drdx = dx*r1;
        double drdy = dy*r1;
        double drdz = dz*r1;
        double drdxx = r1-(pow(dx,2)*r2);
        double drdyy = r1-(pow(dy,2)*r2);
        double drdzz = r1-(pow(dz,2)*r2);
        double drdxy = -dx*dy*r2;
        double drdxz = -dx*dz*r2;
        double drdyz = -dy*dz*r2;

        // Thin plate spline
        double q = alfc;
        RBF[i][1]=pow(r,(0.5*q));                                                               //W
        RBF[i][2]=0.5*q*pow(r,(0.5*q-1))*drdx;                                                  //dWdx
        RBF[i][3]=0.5*q*pow(r,(0.5*q-1))*drdy;                                                  //dWdy
        RBF[i][4]=0.5*q*pow(r,(0.5*q-1))*drdz;                                                  //dWdz
        RBF[i][5]=0.5*q*pow(r,(0.5*q-1))*drdxx+0.5*q*(0.5*q-1)*pow(r,(0.5*q-2))*drdx*drdx;      //dWdxx
        RBF[i][6]=0.5*q*pow(r,(0.5*q-1))*drdyy+0.5*q*(0.5*q-1)*pow(r,(0.5*q-2))*drdy*drdy;      //dWdyy
        RBF[i][7]=0.5*q*pow(r,(0.5*q-1))*drdzz+0.5*q*(0.5*q-1)*pow(r,(0.5*q-2))*drdz*drdz;      //dWdzz
        RBF[i][8]=0.5*q*pow(r,(0.5*q-1))*drdxy+0.5*q*(0.5*q-1)*pow(r,(0.5*q-2))*drdx*drdy;      //dWdxy
        RBF[i][9]=0.5*q*pow(r,(0.5*q-1))*drdxz+0.5*q*(0.5*q-1)*pow(r,(0.5*q-2))*drdx*drdz;      //dWdxz
        RBF[i][10]=0.5*q*pow(r,(0.5*q-1))*drdyz+0.5*q*(0.5*q-1)*pow(r,(0.5*q-2))*drdy*drdz;     //dWdyz
    }
    if (monomials > 0)
    {
        for (int k=1; k <= monomials; k++)
        {
            RBF[TotSuppPnts+k][1]=g_base.N[k];
            RBF[TotSuppPnts+k][2]=g_base.dNdx[k];
            RBF[TotSuppPnts+k][3]=g_base.dNdy[k];
            RBF[TotSuppPnts+k][4]=g_base.dNdz[k];
            RBF[TotSuppPnts+k][5]=g_base.dNdxx[k];
            RBF[TotSuppPnts+k][6]=g_base.dNdyy[k];
            RBF[TotSuppPnts+k][7]=g_base.dNdzz[k];
            RBF[TotSuppPnts+k][8]=g_base.dNdxy[k];
            RBF[TotSuppPnts+k][9]=g_base.dNdxz[k];
            RBF[TotSuppPnts+k][10]=g_base.dNdyz[k];
        }
    }
    return RBF;
}
