#include "../pntwrks.h"

tuple<vector<double>, vector<double> > computeTransport(int dim, vector<double> T, vector<double> Vx, vector<double> Vy, vector<double> Vz, vector<double> Q, vector<double> D, double dt, int t, vector<double> RHS, vector<POINT> points, vector<double> vof, int TotalPoints, vector<KERNEL> kernels, BOUNDARY_CONDITION BC_U)
{
    int time_step_scheme = 2;
    vector<double> T_new(TotalPoints+1);
    vector<double> RHS_old(TotalPoints+1);
    RHS_old = RHS;
    
#pragma omp parallel for
        for (int i = 1; i <= TotalPoints; i++)
        {
            double Tx = 0;
            double Ty = 0;
            double Tz = 0;
			double nabla2_T = 0;
            for (int j = 1; j <= kernels[i].TotalNbrs; j++)
            {
                int pntnum = kernels[i].nbrs[j];
				double T_diff = T[pntnum]-T[i];
                Tx = Tx + kernels[i].dNdx[j]*T_diff;
                Ty = Ty + kernels[i].dNdy[j]*T_diff;
                Tz = Tz + kernels[i].dNdz[j]*T_diff;
				nabla2_T = nabla2_T + kernels[i].nabla2[j]*T_diff;
            }
            double K = vof[i]*D[1]+(1-vof[i])*D[2];//assumes vof = 1 
			RHS[i] = (K*nabla2_T-(Vx[i]*Tx+Vy[i]*Ty+Vz[i]*Tz))+Q[i];
        }
    
    if (time_step_scheme == 1) //Adam-Bashforth O1
    {
#pragma omp parallel for
        for (int i = 1; i <= TotalPoints; i++)
        {T_new[i] = T[i] + dt*RHS[i];}
    }
    if (time_step_scheme == 2) //Adam-Bashforth O2
    {
        if (t < 2)
        {
#pragma omp parallel for
            for (int i = 1; i <= TotalPoints; i++)
            {T_new[i] = T[i] + dt*RHS[i];}
        }
        else
        {
#pragma omp parallel for
            for (int i = 1; i <= TotalPoints; i++)
            {T_new[i] = T[i] + (dt/2)*(3*RHS[i] - RHS_old[i]);}
        }
    }
    for (int r = 1; r <= BC_U.total; r++)
    {
        T_new[BC_U.points[r]] = BC_U.values[r];
    }
    return make_tuple(T_new, RHS);
    
}
vector<double> computeTransportUpwind(int dim, vector<double> T, vector<double> Vx, vector<double> Vy, vector<double> Vz, vector<double> Q, vector<double> vof, vector<double> diffusivity, int TotalPhases, double AV_factor, double dt, vector<POINT> points, int TotalPoints, vector<KERNEL> kernels, BOUNDARY_CONDITION BC_U, double upwind_ratio)
{
    double D = diffusivity[1];
    vector<double> RHS(TotalPoints+1);

    #pragma omp parallel for
    for (int i = 1; i <= TotalPoints; i++)
    {
        double Tx = 0;
        double Ty = 0;
        double Tz = 0;
        double Tx_up = 0;
        double Ty_up = 0;
        double Tz_up = 0;
        double Tx_dn = 0;
        double Ty_dn = 0;
        double Tz_dn = 0;
        double nabla2_T = 0;
        //double sum_w_up = 0;
        //double sum_w_dn = 0;
        for (int j = 1; j <= kernels[i].TotalNbrs; j++)
        {
            int nbr = kernels[i].nbrs[j];
            double xj_xi = points[nbr].x-points[i].x;
            double yj_yi = points[nbr].y-points[i].y;         
            double zj_zi = points[nbr].z-points[i].z;

            //double d = computeDistance(px[i], py[i], px[nbr], py[nbr]);
            //double d_sq = d*d;
            double check = Vx[i]*xj_xi+Vy[i]*yj_yi+Vz[i]*zj_zi;
            if (check < 0) //upwind
            {
                Tx_up = Tx_up + ((T[nbr]-T[i])*kernels[i].dNdx[j]);
                Ty_up = Ty_up + ((T[nbr]-T[i])*kernels[i].dNdy[j]);
                Tz_up = Tz_up + ((T[nbr]-T[i])*kernels[i].dNdz[j]);

                //Tx_up = Tx_up + (T[nbr]-T[i])*(xj_xi)*kernels[i].N[j];
                //Ty_up = Ty_up + (T[nbr]-T[i])*(yj_yi)*kernels[i].N[j];

                //double w = exp(1/d);
                //Tx_up = Tx_up + ((T[nbr]-T[i])*(xj_xi)*w)/d_sq;
                //Ty_up = Ty_up + ((T[nbr]-T[i])*(yj_yi)*w)/d_sq;
                //sum_w_up = sum_w_up + w;
            }

            if (check >= 0) // downwind
            {
                Tx_dn = Tx_dn + ((T[nbr]-T[i])*kernels[i].dNdx[j]);
                Ty_dn = Ty_dn + ((T[nbr]-T[i])*kernels[i].dNdy[j]);
                Tz_dn = Tz_dn + ((T[nbr]-T[i])*kernels[i].dNdz[j]);

                //Tx_dn = Tx_dn + (T[nbr]-T[i])*(xj_xi)*kernels[i].N[j];
                //Ty_dn = Ty_dn + (T[nbr]-T[i])*(yj_yi)*kernels[i].N[j];

                //double w = exp(1/d);
                //Tx_dn = Tx_dn + ((T[nbr]-T[i])*(xj_xi)*w)/d_sq;
                //Ty_dn = Ty_dn + ((T[nbr]-T[i])*(yj_yi)*w)/d_sq;
                //sum_w_dn = sum_w_dn + w;
            }
            nabla2_T = nabla2_T + kernels[i].nabla2[j]*(T[nbr]-T[i]);
        }
        /*
        if (sum_w_up > 0)
        {
        Tx_up = Tx_up/sum_w;
        Ty_up = Ty_up/sum_w;
        }
        if (sum_w_dn > 0)
        {
        Tx_dn = Tx_dn/sum_w_dn;
        Ty_dn = Ty_dn/sum_w_dn;
        }
        */
        Tx = dim*(upwind_ratio*Tx_up+(1-upwind_ratio)*Tx_dn);
        Ty = dim*(upwind_ratio*Ty_up+(1-upwind_ratio)*Ty_dn);
        Tz = dim*(upwind_ratio*Tz_up+(1-upwind_ratio)*Tz_dn);
        //nabla_T = dim*nabla_T;
        if (TotalPhases > 1)
        {D = diffusivity[2]*vof[i] + diffusivity[1]*(1-vof[i]);} //assumes phi = 1 solid, phi = 0 is liquid
        //RHS[i] = (D*nabla_T-(Vx[i]*Tx+Vy[i]*Ty))+Q[i];
        RHS[i] = (D*nabla2_T-(Vx[i]*Tx+Vy[i]*Ty+Vz[i]*Tz))+Q[i];
    }

#pragma omp parallel for
    for (int i = 1; i <= TotalPoints; i++)
    {T[i] = T[i] + 10*dt*RHS[i];}
    for (int r = 1; r <= BC_U.total; r++)
    {
        T[BC_U.points[r]] = BC_U.values[r];
    }
    return T;
}
vector<double> computeTransportIterative(int dim, vector<double> T, vector<double> T_iter, vector<double> Vx, vector<double> Vy, vector<double> Vz, vector<double> Q, vector<double> vof, vector<double> diffusivity, int TotalPhases, double AV_factor, double dt, vector<POINT> points, int TotalPoints, vector<KERNEL> kernels, BOUNDARY_CONDITION BC_U)
{
    double D = diffusivity[1];
    vector<double> RHS(TotalPoints+1);
#pragma omp parallel for
    for (int i = 1; i <= TotalPoints; i++)
    {
        double Tx = 0;
        double Ty = 0;
        double Tz = 0;
        double nabla2_T = 0;
        //double art_visc_x = 0;
        //double art_visc_y = 0;
        for (int j = 1; j <= kernels[i].TotalNbrs; j++)
        {
            int nbr = kernels[i].nbrs[j];
            Tx = Tx + kernels[i].dNdx[j]*(T_iter[nbr]-T_iter[i]);
            Ty = Ty + kernels[i].dNdy[j]*(T_iter[nbr]-T_iter[i]);
            Tz = Tz + kernels[i].dNdz[j]*(T_iter[nbr]-T_iter[i]);
            nabla2_T = nabla2_T + kernels[i].nabla2[j]*(T_iter[nbr]-T_iter[i]);
        }
        if (TotalPhases > 1)
        {D = diffusivity[2]*vof[i] + diffusivity[1]*(1-vof[i]);} //assumes phi = 1 solid, phi = 0 is liquid
        RHS[i] = (D*nabla2_T)-(Vx[i]*Tx+Vy[i]*Ty+Vz[i]*Tz)+Q[i];
    }

#pragma omp parallel for
    for (int i = 1; i <= TotalPoints; i++)
    {T[i] = T[i] + dt*RHS[i];}
    for (int r = 1; r <= BC_U.total; r++)
    {
        T[BC_U.points[r]] = BC_U.values[r];
    }
    return T;
}
vector<double> computeTransportGaussSeidel(int dim, vector<double> T, vector<double> Vx, vector<double> Vy, vector<double> Vz, vector<double> Q, vector<double> D, double dt, vector<POINT> points, vector<double> vof, int TotalPoints, vector<KERNEL> kernels, BOUNDARY_CONDITION BC_U)
{
    double beta = 3;
#pragma omp parallel for
    for (int i = 1; i <= TotalPoints; i++)
    {
        double Tx = 0;
        double Ty = 0;
        double Tz = 0;
		double nabla2_T = 0;
        for (int j = 1; j <= kernels[i].TotalNbrs; j++)
        {
            int pntnum = kernels[i].nbrs[j];
			double T_diff = T[pntnum]-T[i];
            Tx = Tx + kernels[i].dNdx[j]*T_diff;
            Ty = Ty + kernels[i].dNdy[j]*T_diff;
            Tz = Tz + kernels[i].dNdz[j]*T_diff;
			nabla2_T = nabla2_T + kernels[i].nabla2[j]*T_diff;
        }
        double K = vof[i]*D[1]+(1-vof[i])*D[2];
        //assumes vof = 1
        double RHS = (K*nabla2_T-(Vx[i]*Tx+Vy[i]*Ty+Vz[i]*Tz))+Q[i];
        T[i] = beta*(T[i] + dt*RHS) + (1-beta)*T[i];
    }
    for (int r = 1; r <= BC_U.total; r++)
    {
        T[BC_U.points[r]] = BC_U.values[r];
    }
    return T;
}
vector<double> computeTransportGFD(double radius, vector<double> T, vector<double> Vx, vector<double> Vy, vector<double> Q, vector<double> vof, vector<double> diffusivity, int TotalPhases, double dt, vector<POINT> points, int TotalPoints, vector<KERNEL> kernels, BOUNDARY_CONDITION BC_U)
{   //only works for 2D problems
    double D = diffusivity[1];
    for (int i = 1; i <= TotalPoints; i++)
    {
        vector<double> G;
        G = computeGFD(points[i].x, points[i].y, T[i], points, T, kernels[i].nbrs, kernels[i].TotalNbrs, radius);
        if (TotalPhases > 1)
        {D = diffusivity[2]*vof[i] + diffusivity[1]*(1-vof[i]);} //assumes phi = 1 solid, phi = 0 is liquid
        T[i] = T[i] + dt*(D*(G[3]+G[4])-(Vx[i]*G[1]+Vy[i]*G[2])+Q[i]);
    }
    for (int r = 1; r <= BC_U.total; r++)
    {
        T[BC_U.points[r]] = BC_U.values[r];
    }
    return T;
}

tuple<vector<double>,vector<double>,vector<double>,vector<double>,vector<double>>  computePeaksFunction(vector<POINT> points, int TotalPoints)
{
    vector<double> f(TotalPoints+1);
     vector<double> dfdx(TotalPoints+1);
     vector<double> dfdy(TotalPoints+1);
     vector<double> dfdxx(TotalPoints+1);
     vector<double> dfdyy(TotalPoints+1);
    
    for (int i = 1; i <= TotalPoints; i++)
    {
        double x = points[i].x;
        double y = points[i].y;
        f[i] = 3*pow((1-x),2)*exp(-pow(x,2) - pow((y+1),2))-10*(x/5 - pow(x,3) - pow(y,5))*exp(-pow(x,2)-pow(y,2))- 1/3*exp(-pow((x+1),2) - pow(y,2));
        double term = exp(-pow((y+1),2)-x*x);
        dfdx[i] = -6*pow((1-x),2)*x*term-6*(1-x)*term +(2/3)*(x+1)*exp(-y*y-pow((x+1),2))+20*x*(-pow(y,5)-pow(x,3)+x/5)*exp(-y*y-x*x)-10*(0.2-3*x*x)*exp(-y*y-x*x);
        dfdy[i] = -6*pow((1-x),2)*(y+1)*term+(2/3)*y*exp(-y*y-pow((x+1),2))+20*y*(-pow(y,5)-pow(x,3)+x/5)*exp(-y*y-x*x)+50*pow(y,4)*exp(-y*y-x*x);
        dfdxx[i] = 12*pow((1-x),2)*x*x*term+24*(1-x)*x*term-6*pow((1-x),2)*term+6*term-(4/3)*pow((x+1),2)*exp(-y*y-pow((x+1),2))+(2/3)*exp(-y*y-pow((x+1),2))-40*x*x*(-pow(y,5)-pow(x,3)+(x/5))*exp(-y*y-x*x)+20*(-pow(y,5)-pow(x,3)+(x/5))*exp(-y*y-x*x)+40*x*(0.2-3*x*x)*exp(-y*y-x*x)+60*x*exp(-y*y-x*x);
        dfdyy[i] = 12*pow((1-x),2)*pow((y+1),2)*term-6*pow((1-x),2)*term-(4/3)*y*y*exp(-y*y-pow((x+1),2))+(2/3)*exp(-y*y-pow((x+1),2))-200*pow(y,5)*exp(-y*y-x*x)-40*y*y*(-pow(y,5)-pow(x,3)+(x/5))*exp(-y*y-x*x)+20*(-pow(y,5)-pow(x,3)+(x/5))*exp(-y*y-x*x)+200*pow(y,3)*exp(-y*y-x*x);
    }

    return make_tuple(f,dfdx,dfdy,dfdxx,dfdyy);
}
