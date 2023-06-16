#include "../pntwrks.h"

void KERNEL::computeSPH(POINT g, vector<POINT> points, vector<int> supp_pnt_nums, int TotSuppPnts, int type, double h, int dim)
{
    vector<double> w(TotSuppPnts+1);
    vector<double> dwdx(TotSuppPnts+1);
    vector<double> dwdy(TotSuppPnts+1);
    vector<double> dwdz(TotSuppPnts+1);
    double eps = 1.0e-6;
    if (type == 1 || type == 6) //Gaussian Kernel
    {
        double pi = 3.14;
        double alpha = 1.0/(pow(pi,0.5*dim)*pow(h,dim));
        for (int s = 1; s <= TotSuppPnts; s++)
        {
            double dx = g.x-points[supp_pnt_nums[s]].x;
            double dy = g.y-points[supp_pnt_nums[s]].y;
            double dz = g.z-points[supp_pnt_nums[s]].z;
            double r = sqrt(dx*dx+dy*dy+dz*dz)+eps;
            double q = r/h;
            double dqdx = dx/(r*h);
            double dqdy = dy/(r*h);
            double dqdz = dz/(r*h);
            if (q >= 0 && q <= 3)
            {
                double val =  (alpha * exp(-q*q));
                w[s] = val;
                dwdx[s] = -2.0*dqdx*val;
                dwdy[s] = -2.0*dqdy*val;
                dwdz[s] = -2.0*dqdz*val;
            }
            else
            {
                w[s] = 0;
                dwdx[s] = 0;
                dwdy[s] = 0;
                dwdz[s] = 0;
            }
        }
    }
    else if (type == 2 ) //cubic spline
    {
        double pi = 3.14;
        double alpha = 1;
        if (dim == 1)
        {alpha = 1/(6.0*h);}
        else if (dim == 2)
        {alpha = 15.0/(14*pi*h*h);}
        else if (dim == 3)
        {alpha = 1.0/(4*pi*h*h*h);}
        for (int s = 1; s <= TotSuppPnts; s++)
        {
            double dx = g.x-points[supp_pnt_nums[s]].x;
            double dy = g.y-points[supp_pnt_nums[s]].y;
            double dz = g.z-points[supp_pnt_nums[s]].z;
            double r = sqrt(dx*dx+dy*dy+dz*dz);
            double q = r/h;
            double dqdx = dx/(r*h);
            double dqdy = dy/(r*h);
            double dqdz = dz/(r*h);
            if (q >= 0 && q < 1)
            {
                w[s] = alpha*(pow((2-q),3)-4*pow((1.0-q),3));
                dwdx[s] = alpha*(-3*dqdx*pow((2-q),2)+12*dqdx*pow((1-q),2));
                dwdy[s] = alpha*(-3*dqdy*pow((2-q),2)+12*dqdy*pow((1-q),2));
                dwdz[s] = alpha*(-3*dqdz*pow((2-q),2)+12*dqdz*pow((1-q),2));
            }
            else if(q >= 1 && q < 2)
            {
                w[s] = alpha*(pow((2-q),3));
                dwdx[s] = alpha*(-3*dqdx*pow((2-q),2));
                dwdy[s] = alpha*(-3*dqdy*pow((2-q),2));
                dwdz[s] = alpha*(-3*dqdz*pow((2-q),2));
            }
            else if(q >= 2)
            {
                w[s] = 0;
                dwdx[s] = 0;
                dwdy[s] = 0;
                dwdz[s] = 0;
            }
        }
    }
    else if (type == 3) //quintic spline
    {
        double pi = 3.14;
        double alpha = 1;
        if (dim == 1)
        {alpha = 1.0/(120.0*h);}
        else if (dim == 2)
        {alpha = 7.0/(478.0*pi*h*h);}
        else if (dim == 3)
        {alpha = 1.0/(120.0*pi*h*h*h);}
        for (int s = 1; s <= TotSuppPnts; s++)
        {
            double dx = g.x-points[supp_pnt_nums[s]].x;
            double dy = g.y-points[supp_pnt_nums[s]].y;
            double dz = g.z-points[supp_pnt_nums[s]].z;
            double r = sqrt(dx*dx+dy*dy+dz*dz);
            double q = r/h;
            double drdx = -dx;
            double drdy = -dy;
            double drdz = -dz;
            if (q >= 0 && q < 1)
            {
                w[s] = alpha*(pow((3-q),5)-6*pow((2.0-q),5)+15*pow((1.0-q),5));
                dwdx[s] = alpha*((-120+120*q-50*q*q)/(h*h*drdx));
                dwdy[s] = alpha*((-120+120*q-50*q*q)/(h*h*drdy));
                dwdz[s] = alpha*((-120+120*q-50*q*q)/(h*h*drdz));
            }
            else if(q >= 1 && q < 2)
            {
                w[s] = alpha*(pow((3-q),5)-6*pow((2-q),5));
                dwdx[s] = alpha*((-5*pow((3-q),4)+30*pow((2-q),4))/(h*(drdx/r)));
                dwdy[s] = alpha*((-5*pow((3-q),4)+30*pow((2-q),4))/(h*(drdy/r)));
                dwdz[s] = alpha*((-5*pow((3-q),4)+30*pow((2-q),4))/(h*(drdz/r)));
            }
            else if(q >= 2 && q < 3)
            {
                w[s] = alpha*(pow((3-q),5));
                dwdx[s] = alpha*((-5*pow((3-q),4))/(h*(drdx/r)));
                dwdy[s] = alpha*((-5*pow((3-q),4))/(h*(drdy/r)));
                dwdz[s] = alpha*((-5*pow((3-q),4))/(h*(drdz/r)));
            }
            else
            {
                w[s] = 0;
                dwdx[s] = 0;
                dwdy[s] = 0;
                dwdz[s] = 0;
            }
        }
    }
    if (type == 4 )   //Wendland Quintic C4
    {
        /*
        double pi=3.14;
        double alpha=0;
        if (dim == 2)
        {alpha = 7.0/(4.0*pi*h*h);}
        if (dim == 3)
        {alpha = 21.0/(16.0*pi*h*h*h);}
        */
        double pi = 3.14;
        double alpha=0;
        if (dim == 1)
        {alpha = 3.0/(64.0*h);}
        else if (dim == 2)
        {alpha = 7.0/(64.0*pi*h*h);}
        else if (dim == 3)
        {alpha = 21.0/(16.0*pi*h*h*h);} // alpha = 21.0/(256.0*pi*h*h*h);}

        for (int s = 1; s <= TotSuppPnts; s++)
        {
            double dx = g.x-points[supp_pnt_nums[s]].x;
            double dy = g.y-points[supp_pnt_nums[s]].y;
            double dz = g.z-points[supp_pnt_nums[s]].z;
            double r = sqrt(dx*dx+dy*dy+dz*dz);
            double q = r/h;
            double dqdx = dx/(r*h);
            double dqdy = dy/(r*h);
            double dqdz = dz/(r*h);
            if (q >= 0 && q <= 2)
            {
                w[s] = alpha*pow((2-q),4)*(2*q+1);
                dwdx[s] = -alpha*4*pow((2-q),3)*dqdx*(2*q+1)+alpha*pow((2-q),4)*(2*dqdx);
                dwdy[s] = -alpha*4*pow((2-q),3)*dqdy*(2*q+1)+alpha*pow((2-q),4)*(2*dqdy);
                dwdz[s] = -alpha*4*pow((2-q),3)*dqdz*(2*q+1)+alpha*pow((2-q),4)*(2*dqdz);
            }
            else
            {
                w[s] = 0;
                dwdx[s] = 0;
                dwdy[s] = 0;
                dwdz[s] = 0;
            }
        }
    }
    if (type == 5 ) //Wendland Quintic C6
    {
        double pi=3.14;
        double alpha=0;
        if (dim == 2)
        {alpha = 9.0/(4.0*pi*h*h);}
        if (dim == 3)
        {alpha = 495.0/(256.0*pi*h*h*h);}
        for (int s = 1; s <= TotSuppPnts; s++)
        {
            double dx = g.x-points[supp_pnt_nums[s]].x;
            double dy = g.y-points[supp_pnt_nums[s]].y;
            double dz = g.z-points[supp_pnt_nums[s]].z;
            double r = sqrt(dx*dx+dy*dy+dz*dz);
            double q = r/h;
            double dqdx = dx/(r*h);
            double dqdy = dy/(r*h);
            double dqdz = dz/(r*h);
            if (q >= 0 && q <= 2)
            {
                w[s] = alpha*pow((2-q),6)*((35/12)*q*q+3*q+1);
                dwdx[s] = alpha*6*pow((2-q),5)*(-dqdx)*((35/12)*q*q+3*q+1)+alpha*pow((2-q),6)*((35/12)*2*dqdx+3*dqdx);
                dwdy[s] = alpha*6*pow((2-q),5)*(-dqdy)*((35/12)*q*q+3*q+1)+alpha*pow((2-q),6)*((35/12)*2*dqdy+3*dqdy);
                dwdz[s] = alpha*6*pow((2-q),5)*(-dqdz)*((35/12)*q*q+3*q+1)+alpha*pow((2-q),6)*((35/12)*2*dqdz+3*dqdz);
            }
            else
            {
                w[s] = 0;
                dwdx[s] = 0;
                dwdy[s] = 0;
                dwdz[s] = 0;
            }
        }
    }
    if (type == 7 || type == 8) //inverse distance
    {
        for (int s = 1; s <= TotSuppPnts; s++)
        {
            int j = supp_pnt_nums[s];
            double d = points[j].distance(g)+1e-6;
            w[s] = ((h/(d)) - 1);
            dwdx[s] = (1/(d*d))*w[s]*(points[j].x-g.x);
            dwdy[s] = (1/(d*d))*w[s]*(points[j].y-g.y);
            dwdz[s] = (1/(d*d))*w[s]*(points[j].z-g.z);
        }
    }
    N = w;
    dNdx = dwdx;
    dNdy = dwdy;
    dNdz = dwdz;

    if (type != 6 && type != 8)
    {
       
        // Normalizing the interpolants
        double sum = 0;
        for (int s = 1; s <= TotSuppPnts; s++)
        {
            double dx = g.x-points[supp_pnt_nums[s]].x;
            double dy = g.y-points[supp_pnt_nums[s]].y;
            double dz = g.z-points[supp_pnt_nums[s]].z;
            double r = sqrt(dx*dx+dy*dy+dz*dz);

            double m=0;
            if (dim == 1)
            {m = r;}
            if (dim == 2)
            {m = r*r;}
            if (dim == 3)
            {m = r*r*r;}

            N[s] = m*N[s];
            dNdx[s] = m*dNdx[s];
            dNdy[s] = m*dNdy[s];
            dNdz[s] = m*dNdz[s];
            sum = sum + N[s];
        }
       
/*
        double sum = 0;
        for (int s = 1; s <= TotSuppPnts; s++)
        {sum = sum + N[s];}
  */    
        for (int s = 1; s <= TotSuppPnts; s++)
        {
            N[s] = N[s]/sum;
            dNdx[s] = dNdx[s]/sum;
            dNdy[s] = dNdy[s]/sum;
            dNdz[s] = dNdz[s]/sum;
        }
        sum = 0;
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
            cout << "ERROR in constructed SPH interpolants. Partition of Unity is Not satisfied. sum = " << sum << endl;
            exit(0);
        }
        /*
        double eps = 1e-5;
        if (abs(sumx) >= eps || abs(sumy) >= eps || abs(sumz) >= eps)
        {
            cout << "ERROR in constructed SPH interpolants. Sum of derivatives not equal to zero. " << endl;
            cout << sumx << "\t" << sumy << "\t" << sumz << endl;
            exit(0);
        }
        */
    }
}
