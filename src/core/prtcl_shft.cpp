#include "../pntwrks.h"

vector<vector<double> > computeParticleShiftingVelocity(double factor, double srf_mrkr, double h, double dt, vector<POINT> pnts, int TotalPoints, vector<KERNEL> kernels)
{
    //===========================================================
    // computing particle shifting
    //===========================================================
    vector<vector<double> > VEL(TotalPoints+1, vector<double> (3));
    vector<double> prtcl_conc(TotalPoints+1);

#pragma omp parallel for
    for (int i = 1; i<= TotalPoints; i++)
    {
        prtcl_conc[i] = 0;
        for (int s = 1; s <= kernels[i].TotalNbrs; s++)
        {
            int j = kernels[i].nbrs[s];
            double d = pnts[i].distance(pnts[j])+1e-8;
            prtcl_conc[i] = prtcl_conc[i] + pow(((h/d) - 1),2);
            //prtcl_conc[i] = prtcl_conc[i] + ((h/d) - 1);
        }
    }

    for (int i = 1;  i<= TotalPoints; i++)
    {
        double dC_x = 0;
        double dC_y = 0;
        //if (prtcl_conc[i] >= srf_mrkr)
        //{
        dC_x = 0;
        dC_y = 0;
        for (int s = 1; s <= kernels[i].TotalNbrs; s++)
        {
            int j = kernels[i].nbrs[s];
            double d = pnts[i].distance(pnts[j])+1e-8;
            //double w = ((h/d) - 1)*((h/d) - 1);
            //dC_x = dC_x + (prtcl_conc[j])*((2*h*w)/(d*d))*(points[j].x-px[i]);
            //dC_y = dC_y + (prtcl_conc[j])*((2*h*w)/(d*d))*(points[j].y-py[i]);
            double w = ((h/d) - 1);
            dC_x = dC_x + (prtcl_conc[j])*(w/(d*d))*(pnts[j].x-pnts[i].x);
            dC_y = dC_y + (prtcl_conc[j])*(w/(d*d))*(pnts[j].y-pnts[i].y);
        }
        //}
        double D = factor*h;
        VEL[i][1] = - D*dC_x;
        VEL[i][2] = - D*dC_y;
    }
    return VEL;
}
vector<vector<double> > computeParticleShiftingVelocity_v2(double factor, double srf_mrkr, double h, double dt, vector<POINT> pnts,  int TotalPoints, vector<KERNEL> kernels)
{
    //===========================================================
    // computing particle shifting
    //===========================================================
    vector<vector<double> > VEL(TotalPoints+1, vector<double> (3));
    vector<double> prtcl_conc(TotalPoints+1);
    vector<double> P(TotalPoints+1);
    double gamma = 7;
    double k = 0.05;
    double rho0 = 1;
    double m = 1;
#pragma omp parallel for
    for (int i = 1; i<= TotalPoints; i++)
    {
        prtcl_conc[i] = 0;
        for (int s = 1; s <= kernels[i].TotalNbrs; s++)
        {
            int j = kernels[i].nbrs[s];
            double d = pnts[i].distance(pnts[j])+1e-8;
            //prtcl_conc[i] = prtcl_conc[i] + pow(((h/d) - 1),2);
            prtcl_conc[i] = prtcl_conc[i] + ((h/d) - 1);
        }
        P[i] = k*(pow((prtcl_conc[i]/rho0),gamma)-1);
        //P[i] = k*(prtcl_conc[i]-rho0);
    }

    for (int i = 1;  i<= TotalPoints; i++)
    {
        double A_pressure_x = 0;
        double A_pressure_y = 0;
        double rho_rho_i = prtcl_conc[i]*prtcl_conc[i];
        double p1 = P[i]/rho_rho_i;
        for (int s = 1; s <= kernels[i].TotalNbrs; s++)
        {
            int j = kernels[i].nbrs[s];
            double rho_rho_j = prtcl_conc[j]*prtcl_conc[j];
            double p2 = P[j]/rho_rho_j;
            //A_pressure_x = A_pressure_x + (m)*(p1+p2)*kernels[i].dNdx[s];
            //A_pressure_y = A_pressure_y + (m)*(p1+p2)*kernels[i].dNdy[s];

            double d = pnts[i].distance(pnts[j])+1e-8;
            double w = pow(((h/d) - 1),1);
            A_pressure_x = A_pressure_x + (m)*(p1+p2)*(w/(d*d))*(pnts[j].x-pnts[i].x);
            A_pressure_y = A_pressure_y + (m)*(p1+p2)*(w/(d*d))*(pnts[j].y-pnts[i].y);

        }
        VEL[i][1] = -0.1*factor*h*A_pressure_x;
        VEL[i][2] = -0.1*factor*h*A_pressure_y;
    }
    return VEL;
}
/*
vector<vector<double> > computeParticleShiftingPFM(double factor, double srf_mrkr, double h, double dt, vector<double> px, vector<double> py, int TotalPoints, vector<KERNEL> kernels)
{
    //===========================================================
    // computing particle shifting
    //===========================================================
    vector<vector<double> > VEL(TotalPoints+1, vector<double> (3));
    vector<double> prtcl_conc(TotalPoints+1);
    //double k = 3*h;
    double k = 0.05;
    double eps = 1e-6;
#pragma omp parallel for
    for (int i = 1; i<= TotalPoints; i++)
    {
        prtcl_conc[i] = 0;
        for (int s = 1; s <= kernels[i].TotalNbrs; s++)
        {
            int j = kernels[i].nbrs[s];
            double d = computeDistance(px[i],py[i], px[j], py[j])+1e-8;
            prtcl_conc[i] = prtcl_conc[i] + pow(((h/d) - 1),2);
        }
    }

    vector<double> P(TotalPoints+1);
    double gamma = 7;

    double rho0 = 1;
    double m = 1;
#pragma omp parallel for
    for (int i = 1; i<= TotalPoints; i++)
    {
        prtcl_conc[i] = 0;
        for (int s = 1; s <= kernels[i].TotalNbrs; s++)
        {
            int j = kernels[i].nbrs[s];
            double d = computeDistance(px[i],py[i], px[j], py[j])+1e-8;
            //prtcl_conc[i] = prtcl_conc[i] + pow(((h/d) - 1),2);
            prtcl_conc[i] = prtcl_conc[i] + ((h/d) - 1);
        }
        P[i] = k*(pow((prtcl_conc[i]/rho0),gamma)-1);
        //P[i] = k*(prtcl_conc[i]-rho0);
    }

    vector<double> drodx(TotalPoints+1);
    vector<double> drody(TotalPoints+1);
    vector<double> term(TotalPoints+1);
    vector<double> termx(TotalPoints+1);
    vector<double> termy(TotalPoints+1);
    for (int i = 1; i <= TotalPoints; i++)
    {
        drodx[i] = 0;
        drody[i] = 0;
        for (int j = 1; j <= kernels[i].TotalNbrs; j++)
        {
            int pntnum = kernels[i].nbrs[j];
            double d = computeDistance(px[i],py[i], px[pntnum], py[pntnum])+1e-8;
            double w = pow(((h/d) - 1),1);
            drodx[i] = drodx[i] + (prtcl_conc[pntnum]-prtcl_conc[i])*(w/(d*d))*(px[pntnum]-px[i]);
            drody[i] = drody[i] + (prtcl_conc[pntnum]-prtcl_conc[i])*(w/(d*d))*(py[pntnum]-py[i]);

            //  drodx[i] = drodx[i] + kernels[i].dNdx[j]*(prtcl_conc[pntnum]-prtcl_conc[i]);
            //  drody[i] = drody[i] + kernels[i].dNdy[j]*(prtcl_conc[pntnum]-prtcl_conc[i]);
        }
        double abs = sqrt(drodx[i]*drodx[i]+drody[i]*drody[i]);
        term[i] = (k*abs*abs)/(prtcl_conc[i]*prtcl_conc[i]+1e-6);
        termx[i] = 2*k*drodx[i]/prtcl_conc[i];
        termy[i] = 2*k*drody[i]/prtcl_conc[i];
    }


    for (int i = 1;  i<= TotalPoints; i++)
    {
        double A_pressure_x = 0;
        double A_pressure_y = 0;
        double rho_rho_i = prtcl_conc[i]*prtcl_conc[i];
        double p1 = P[i]/rho_rho_i;
        for (int s = 1; s <= kernels[i].TotalNbrs; s++)
        {
            int j = kernels[i].nbrs[s];
            double rho_rho_j = prtcl_conc[j]*prtcl_conc[j];
            double p2 = P[j]/rho_rho_j;
            //A_pressure_x = A_pressure_x + (m)*(p1+p2)*kernels[i].dNdx[s];
            //A_pressure_y = A_pressure_y + (m)*(p1+p2)*kernels[i].dNdy[s];

            double d = computeDistance(px[i],py[i], px[j], py[j])+1e-8;
            double w = pow(((h/d) - 1),1);
            A_pressure_x = A_pressure_x + (m)*(p1+p2)*(w/(d*d))*(px[j]-px[i]);
            A_pressure_y = A_pressure_y + (m)*(p1+p2)*(w/(d*d))*(py[j]-py[i]);
        }

        double dtermdx = 0;
        double dtermdy = 0;
        double dtermdxx = 0;
        double dtermdyy = 0;
        for (int j = 1; j <= kernels[i].TotalNbrs; j++)
        {
            int pntnum = kernels[i].nbrs[j];
            double d = computeDistance(px[i],py[i], px[pntnum], py[pntnum])+1e-8;
            double w = pow(((h/d) - 1),1);
            //dtermdx = dtermdx + (term[pntnum]+term[i])*(w/(d*d))*(px[pntnum]-px[i]);
            //dtermdy = dtermdy + (term[pntnum]+term[i])*(w/(d*d))*(py[pntnum]-py[i]);

            dtermdx = dtermdx + (prtcl_conc[pntnum]-prtcl_conc[i])*(w/(d*d))*(px[pntnum]-px[i]);
            dtermdy = dtermdy + (prtcl_conc[pntnum]-prtcl_conc[i])*(w/(d*d))*(py[pntnum]-py[i]);

            //dtermdx = dtermdx + kernels[i].dNdx[j]*(term[pntnum]+term[i]);
            //dtermdy = dtermdy + kernels[i].dNdy[j]*(term[pntnum]+term[i]);

            //double d = computeDistance(px[i],py[i], px[pntnum], py[pntnum]);
            double nijx = (px[pntnum]-px[i])/(d+eps);
            double nijy = (py[pntnum]-py[i])/(d+eps);
            //dtermdxx = dtermdxx + (termx[pntnum]-termx[i])*(nijx/(d+eps))*kernels[i].dNdx[j];
            //dtermdyy = dtermdyy + (termy[pntnum]-termy[i])*(nijy/(d+eps))*kernels[i].dNdy[j];
            //dtermdxx = dtermdxx + (termx[pntnum]+termx[i])*kernels[i].dNdx[j];
            //dtermdyy = dtermdyy + (termy[pntnum]+termy[i])*kernels[i].dNdy[j];
            dtermdxx = dtermdxx + (prtcl_conc[pntnum]-prtcl_conc[i])*(nijx/(d+eps))*kernels[i].dNdx[j];
            dtermdyy = dtermdyy + (prtcl_conc[pntnum]-prtcl_conc[i])*(nijy/(d+eps))*kernels[i].dNdy[j];
        }
        //dtermdxx = 0;
        //dtermdyy = 0;
        //VEL[i][1] = dt*(dtermdx+dtermdxx);
        //VEL[i][2] = dt*(dtermdy+dtermdyy);
        //  VEL[i][1] = dt*(dtermdx);
        //  VEL[i][2] = dt*(dtermdy);

        //VEL[i][1] = -0.1*factor*h*(A_pressure_x + (dtermdx));
        //VEL[i][2] = -0.1*factor*h*(A_pressure_y + (dtermdy));
        VEL[i][1] = -0.5*factor*h*(A_pressure_x) + 10*dt*(dtermdx);
        VEL[i][2] = -0.5*factor*h*(A_pressure_y) + 10*dt*(dtermdy);

        //  VEL[i][1] = -0.5*dt*(A_pressure_x) + 10*dt*(dtermdx);
        //VEL[i][2] = -0.5*dt*(A_pressure_y) + 10*dt*(dtermdy);

        //  VEL[i][1] = dt*(dtermdx);
        //  VEL[i][2] = dt*(dtermdy);
    }
    return VEL;
}
*/
