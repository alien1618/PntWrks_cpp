#include "../pntwrks.h"

//********************************************************************************
// FUNCTIONS FOR SOLVING THE VOLUME OF FLUID FUNCTION
//********************************************************************************
vector<double> computeVolumeOfFluid(vector<KERNEL> kernels, vector<double> phi_now, int TotalPoints,  vector<double> Vx, vector<double> Vy, vector<double> Vz, vector<double> Vn, int iterations, double dt)
{
    vector<double> U_store(iterations+1);
    vector<double> phi_iter(TotalPoints+1);
    vector<double> phi_iter_old(TotalPoints+1);
    int method = 1; //method = 1 non-conservative form, method = 2 conservative form

    phi_iter = phi_now;
    phi_iter_old = phi_iter;
    int cntr = 0;
    if (method == 1)
    {
        for (int iternum = 1; iternum <= iterations; iternum++)
        {
            cntr++;
            // cout << "Computing iteration " << iternum << " of " << iterations << endl;
            phi_iter_old = phi_iter;
            for (int i = 1; i <= TotalPoints; i++)
            {
                double dphidx = 0;
                double dphidy = 0;
                double dphidz = 0;
                for (int s = 1; s <= kernels[i].TotalNbrs; s++)
                {
                    int nbr = kernels[i].nbrs[s];
                    double dphi =  phi_iter[nbr]-phi_iter[i];
                    dphidx = dphidx + kernels[i].dNdx[s]*dphi;
                    dphidy = dphidy + kernels[i].dNdy[s]*dphi;
                    dphidz = dphidz + kernels[i].dNdz[s]*dphi;
                }
                double absgradphi = sqrt(dphidx*dphidx + dphidy*dphidy + dphidz*dphidz);
                double adv = Vx[i]*dphidx + Vy[i]*dphidy + Vz[i]*dphidz;
                phi_iter[i] = phi_now[i] - dt*(adv + Vn[i]*absgradphi);
            }
            U_store[iternum] = checkConvergence(cntr, phi_iter, phi_iter_old, TotalPoints);
        }
    }
    else if (method == 2)
    {
        for (int iternum = 1; iternum <= iterations; iternum++)
        {
            cntr++;
            // cout << "Computing iteration " << iternum << " of " << iterations << endl;
            phi_iter_old = phi_iter;
            for (int i = 1; i <= TotalPoints; i++)
            {
                double dphidx = 0;
                double dphidy = 0;
                double dphidz = 0;
                double dphidx1 = 0;
                double dphidy1 = 0;
                double dphidz1 = 0;
                for (int s = 1; s <= kernels[i].TotalNbrs; s++)
                {
                    int nbr = kernels[i].nbrs[s];
                    dphidx = dphidx + kernels[i].dNdx[s]*(phi_iter[nbr]*Vx[nbr]-phi_iter[i]*Vx[i]);
                    dphidy = dphidy + kernels[i].dNdy[s]*(phi_iter[nbr]*Vy[nbr]-phi_iter[i]*Vy[i]);
                    dphidz = dphidz + kernels[i].dNdz[s]*(phi_iter[nbr]*Vz[nbr]-phi_iter[i]*Vz[i]);
                    dphidx1 = dphidx1 + kernels[i].dNdx[s]*(phi_iter[nbr]-phi_iter[i]);
                    dphidy1 = dphidy1 + kernels[i].dNdy[s]*(phi_iter[nbr]-phi_iter[i]);
                    dphidz1 = dphidz1 + kernels[i].dNdz[s]*(phi_iter[nbr]-phi_iter[i]);
                }
                double absgradphi = sqrt(dphidx1*dphidx1 + dphidy1*dphidy1 + dphidz1*dphidz1);
                double adv = dphidx + dphidy + dphidz;
                phi_iter[i] = phi_now[i] - dt*(adv + Vn[i]*absgradphi);
            }
            U_store[iternum] = checkConvergence(cntr, phi_iter, phi_iter_old, TotalPoints);
        }
    }
    /*
    for (int i = 1; i <= TotalPoints; i++)
    {
        if (phi_iter[i] <= 0.0001)
        {phi_iter[i] = 0;}
        else if (phi_iter[i] >= 0.999)
        {phi_iter[i] = 1;}
    }
    */
    //for (int m = 2; m <= iterations; m++)
    //{cout << "Convergence difference = " << m << " = " << U_store[m] << endl;}
    return phi_iter;
}
vector<double> computeAllenCahn(double W, double M0, vector<KERNEL> kernels, vector<double> phi_now, double dim, vector<POINT> particles, int TotalPoints,  vector<double> Vx, vector<double> Vy, vector<double> Vz, vector<double> Vn, int iterations, double dt)
{
    //sharp-interface phase field of sun and beckermann
    //Phi needs to be 0 to 1
    double eps = 0.000001;
    int cntr = 0;

    vector<double> dphi_abs(TotalPoints+1);
    vector<double> termx(TotalPoints+1);
    vector<double> termy(TotalPoints+1);
    vector<double> termz(TotalPoints+1);
    vector<double> phi_iter(TotalPoints+1);
    vector<double> phi_iter_old(TotalPoints+1);

    vector<double> U_store(iterations+1);
    phi_iter = phi_now;
    phi_iter_old = phi_iter;

    for (int iternum = 1; iternum <= iterations; iternum++)
    {
        cntr++;
        // cout << "Computing iteration " << iternum << " of " << iterations << endl;
        phi_iter_old = phi_iter;
        for (int i = 1;  i<= TotalPoints; i++)
        {
            double dphidx = 0;
            double dphidy = 0;
            double dphidz = 0;
            for (int s = 1; s <= kernels[i].TotalNbrs; s++)
            {
                int nbr = kernels[i].nbrs[s];
                double dphi = phi_iter[nbr]-phi_iter[i];
                dphidx = dphidx + dphi*kernels[i].dNdx[s];
                dphidy = dphidy + dphi*kernels[i].dNdy[s];
                dphidz = dphidz + dphi*kernels[i].dNdz[s];
            }
            dphi_abs[i] = dphidx*dphidx + dphidy*dphidy + dphidz*dphidz+eps;
            termx[i] = W*dphidx - phi_iter[i]*(1-phi_iter[i])*(dphidx/dphi_abs[i]);
            termy[i] = W*dphidy - phi_iter[i]*(1-phi_iter[i])*(dphidy/dphi_abs[i]);
            termz[i] = W*dphidz - phi_iter[i]*(1-phi_iter[i])*(dphidz/dphi_abs[i]);
        }
        for (int i = 1;  i<= TotalPoints; i++)
        {
            double dtermdx = 0;
            double dtermdy = 0;
            double dtermdz = 0;
            for (int s = 1; s <= kernels[i].TotalNbrs; s++)
            {
                int nbr = kernels[i].nbrs[s];
                dtermdx = dtermdx + (termx[nbr]-termx[i])*kernels[i].dNdx[s];
                dtermdy = dtermdy + (termy[nbr]-termy[i])*kernels[i].dNdy[s];
                dtermdz = dtermdz + (termz[nbr]-termz[i])*kernels[i].dNdz[s];
            }
            double RHS = dtermdx+dtermdy+dtermdz;
            double dterm2dx = 0;
            double dterm2dy = 0;
            double dterm2dz = 0;
            for (int s = 1; s <= kernels[i].TotalNbrs; s++)
            {
                int nbr = kernels[i].nbrs[s];
                dterm2dx = dterm2dx + ((phi_iter[nbr]*Vx[nbr])-(phi_iter[i]*Vx[i]))*kernels[i].dNdx[s];
                dterm2dy = dterm2dy + ((phi_iter[nbr]*Vy[nbr])-(phi_iter[i]*Vy[i]))*kernels[i].dNdy[s];
                dterm2dz = dterm2dz + ((phi_iter[nbr]*Vz[nbr])-(phi_iter[i]*Vz[i]))*kernels[i].dNdz[s];
            }
            double adv = dterm2dx+dterm2dy+dterm2dz;
            phi_iter[i] = phi_now[i] + dt*(M0*RHS - adv + dphi_abs[i]*Vn[i]);
            if (phi_iter[i] < 0)
            {phi_iter[i] = 0;}
            if (phi_iter[i] > 1)
            {phi_iter[i] = 1;}
        }
        U_store[iternum] = checkConvergence(cntr, phi_iter, phi_iter_old, TotalPoints);
    }
    /*
    cout << "Checking convergence..." << endl;
    for (int i = 2; i <= iterations; i++)
    {cout<< "U_store["<<i<<"] = " << U_store[i] << endl;}
    */
    return phi_iter;
}

vector<double> computeAllenCahn2(double W, double M0, vector<KERNEL> kernels, vector<double> phi_now, double dim, vector<POINT> particles, int TotalPoints,  vector<double> Vx, vector<double> Vy, vector<double> Vz, vector<double> Vn, int iterations, double dt)
{
    //Conservative allen-cahn equation of brassel and bretin, also see kim, lee and choi
    //Phi needs to be -1 to 1
    double eps = 0.000001;
    int cntr = 0;

    vector<double> phi_iter(TotalPoints+1);
    vector<double> phi_iter_old(TotalPoints+1);
    vector<double> U_store(iterations+1);
    phi_iter = phi_now;
    phi_iter_old = phi_iter;

    for (int iternum = 1; iternum <= iterations; iternum++)
    {
        cntr++;
        // cout << "Computing iteration " << iternum << " of " << iterations << endl;
        phi_iter_old = phi_iter;
        double sum_F_prime = 0;
        double sum_F = 0;
        for (int i = 1; i <= TotalPoints; i++)
        {
            sum_F_prime = sum_F_prime + phi_now[i]*((phi_now[i]*phi_now[i])-1);
            sum_F = sum_F + 0.25*pow(((phi_now[i]*phi_now[i])-1),2);
        }
        double beta = sum_F_prime/(sum_F+eps);

        for (int i = 1; i <= TotalPoints; i++)
        {
            double dphidx = 0;
            double dphidy = 0;
            double dphidz = 0;
            double nabla2_phi = 0;
            //double dphiVxdx = 0;
            //double dphiVydy = 0;
            //double dphiVzdz = 0;
            //double V_max = 0;
            for (int s = 1; s <= kernels[i].TotalNbrs; s++)
            {
                int nbr = kernels[i].nbrs[s];
                dphidx = dphidx + (phi_iter[nbr]-phi_iter[i])*kernels[i].dNdx[s];
                dphidy = dphidy + (phi_iter[nbr]-phi_iter[i])*kernels[i].dNdy[s];
                dphidz = dphidz + (phi_iter[nbr]-phi_iter[i])*kernels[i].dNdz[s];
                //dphiVxdx = dphiVxdx + (phi_iter[nbr]*Vx[nbr]-phi_iter[i]*Vx[i])*kernels[i].dNdx[s];
                //dphiVydy = dphiVydy + (phi_iter[nbr]*Vy[nbr]-phi_iter[i]*Vy[i])*kernels[i].dNdy[s];
                //dphiVzdz = dphiVzdz + (phi_iter[nbr]*Vz[nbr]-phi_iter[i]*Vz[i])*kernels[i].dNdz[s];
                nabla2_phi = nabla2_phi + (phi_iter[nbr]-phi_iter[i])*kernels[i].nabla2[s];
                //double V = sqrt(Vx[nbr]*Vx[nbr] + Vy[nbr]*Vy[nbr] + Vz[nbr]*Vz[nbr] + Vn[nbr]*Vn[nbr]);
                //if (V >= V_max)
                //{V_max = V;}
            }
            //double Vi = sqrt(Vx[i]*Vx[i] + Vy[i]*Vy[i] + Vz[i]*Vz[i] + Vn[i]*Vn[i]);
           // if (Vi >= V_max)
            //{V_max = Vi;}
            //double M = M0*V_max;
            double dphi_abs = sqrt(dphidx*dphidx + dphidy*dphidy + dphidz*dphidz);
            double ADV = dphidx*Vx[i]+dphidy*Vy[i]+dphidz*Vz[i];
            //double ADV = dphiVxdx+dphiVydy+dphiVzdz;
            double DIFF = nabla2_phi;
            double k = 2*sqrt(0.25*pow(((phi_now[i]*phi_now[i])-1),2));
            DIFF =  M0*(DIFF) - phi_now[i]*((phi_now[i]*phi_now[i])-1);
            phi_iter[i] = phi_now[i] + dt*(DIFF - ADV + dphi_abs*Vn[i] + beta*k);

            if(phi_iter[i] > 1.0)
            {phi_iter[i]= 1.0;}
            if(phi_iter[i] < -1.0)
            {phi_iter[i] = -1.0;}
        }
        U_store[iternum] = checkConvergence(cntr, phi_iter, phi_iter_old, TotalPoints);
    }
    /*
   cout << "Checking convergence..." << endl;
   for (int i = 2; i <= iterations; i++)
   {cout<< "U_store["<<i<<"] = " << U_store[i] << endl;}
   */

    return phi_iter;
}
vector<double> computeCahnHilliard(double W, double M0, vector<KERNEL> kernels, vector<double> phi_now, double dim, vector<POINT> particles, int TotalPoints,  vector<double> Vx, vector<double> Vy, vector<double> Vz, vector<double> Vn, int iterations, double dt)
{
    //Phi needs to be 0 to 1
    int cntr = 0;
    vector<double> mu(TotalPoints+1);
    vector<double> phi_iter(TotalPoints+1);
    vector<double> phi_iter_old(TotalPoints+1);
    vector<double> U_store(iterations+1);
    phi_iter = phi_now;
    phi_iter_old = phi_iter;

    for (int iternum = 1; iternum <= iterations; iternum++)
    {
        cntr++;
        // cout << "Computing iteration " << iternum << " of " << iterations << endl;
        phi_iter_old = phi_iter;

#pragma omp parallel for
        for (int i = 1; i <= TotalPoints; i++)
        {
            double nabla2_phi = 0;
            for (int j = 1; j <= kernels[i].TotalNbrs; j++)
            {
                int nbr = kernels[i].nbrs[j];
                nabla2_phi = nabla2_phi + (phi_iter[nbr]-phi_iter[i])*kernels[i].nabla2[j];
            }
            double dfdcon =2.0*phi_iter[i]*pow((1-phi_iter[i]),2)-2.0*pow(phi_iter[i],2)*(1.0-phi_iter[i]);
            mu[i] = dfdcon - W*(nabla2_phi);
        }
#pragma omp parallel for
        for (int i = 1; i <= TotalPoints; i++)
        {
			double dphidx = 0;
            double dphidy = 0;
            double dphidz = 0;
            double nabla2_mu = 0;
            //double V_max = 0;
            for (int j = 1; j <= kernels[i].TotalNbrs; j++)
            {
                int nbr = kernels[i].nbrs[j];
                dphidx = dphidx + (phi_iter[nbr]-phi_iter[i])*kernels[i].dNdx[j];
                dphidy = dphidy + (phi_iter[nbr]-phi_iter[i])*kernels[i].dNdy[j];
                dphidz = dphidz + (phi_iter[nbr]-phi_iter[i])*kernels[i].dNdz[j];
                nabla2_mu = nabla2_mu + (mu[nbr]-mu[i])*kernels[i].nabla2[j];
                //double V = sqrt(Vx[nbr]*Vx[nbr] + Vy[nbr]*Vy[nbr] + Vz[nbr]*Vz[nbr] + Vn[nbr]*Vn[nbr]);
                //if (V >= V_max)
                //{V_max = V;}
            }
            /*
            double Vi = sqrt(Vx[i]*Vx[i] + Vy[i]*Vy[i] + Vz[i]*Vz[i] + Vn[i]*Vn[i]);
            if (Vi >= V_max)
            {V_max = Vi;}
            if (V_max == 0)
            {V_max = 1;}
            double M = M0*V_max;
            */
            double dphi_abs = sqrt(dphidx*dphidx + dphidy*dphidy + dphidz*dphidz);
            double ADV = dphidx*Vx[i]+dphidy*Vy[i]+dphidz*Vz[i];
            phi_iter[i]  = phi_now[i] + dt*(M0*nabla2_mu - ADV + dphi_abs*Vn[i]);
            if(phi_iter[i] > 1.0)
            {phi_iter[i]= 1.0;}
            if(phi_iter[i] < 0)
            {phi_iter[i] = 0;}
        }
        U_store[iternum] = checkConvergence(cntr, phi_iter, phi_iter_old, TotalPoints);
    }
    /*
    cout << "Checking convergence..." << endl;
    for (int i = 2; i <= iterations; i++)
    {cout<< "U_store["<<i<<"] = " << U_store[i] << endl;}
    */
    return phi_iter;
}

vector<double> computeCahnHilliard2(double W, double M0, vector<KERNEL> kernels, vector<double> phi_now, double dim, vector<POINT> particles, int TotalPoints,  vector<double> Vx, vector<double> Vy, vector<double> Vz, vector<double> Vn, int iterations, double dt)
{
    //Phi needs to be -1 to 1
    int cntr = 0;

    vector<double> mu(TotalPoints+1);
    vector<double> phi_iter(TotalPoints+1);
    vector<double> phi_iter_old(TotalPoints+1);
    vector<double> U_store(iterations+1);
    phi_iter = phi_now;
    phi_iter_old = phi_iter;

    for (int iternum = 1; iternum <= iterations; iternum++)
    {
        cntr++;
        // cout << "Computing iteration " << iternum << " of " << iterations << endl;
        phi_iter_old = phi_iter;

        for (int i = 1; i <= TotalPoints; i++)
        {
            double nabla2_phi = 0;
          //  double V_max = 0;
            for (int s = 1; s <= kernels[i].TotalNbrs; s++)
            {
                int nbr = kernels[i].nbrs[s];             
		nabla2_phi = nabla2_phi + (phi_iter[nbr]-phi_iter[i])*kernels[i].nabla2[s];
               //double V = sqrt(Vx[nbr]*Vx[nbr] + Vy[nbr]*Vy[nbr] + Vz[nbr]*Vz[nbr] + Vn[nbr]*Vn[nbr]);
                //if (V >= V_max)
                //{V_max = V;}
            }
           // double Vi = sqrt(Vx[i]*Vx[i] + Vy[i]*Vy[i] + Vz[i]*Vz[i] + Vn[i]*Vn[i]);
           // if (Vi >= V_max)
           // {V_max = Vi;}
            mu[i] =  (phi_iter[i]*((phi_iter[i]*phi_iter[i])-1) - W*W*nabla2_phi); // be careful the sign of diff
        }
        for (int i = 1; i <= TotalPoints; i++)
        {
            double dphidx = 0;
            double dphidy = 0;
            double dphidz = 0;
            double nabla2_mu = 0;
            //double dphiVxdx = 0;
            //double dphiVydy = 0;
            //double dphiVzdz = 0;
           // double V_max = 0;
            for (int s = 1; s <= kernels[i].TotalNbrs; s++)
            {
                int nbr = kernels[i].nbrs[s];
                double dphi = phi_iter[nbr]-phi_iter[i]; 
                dphidx = dphidx + dphi*kernels[i].dNdx[s];
                dphidy = dphidy + dphi*kernels[i].dNdy[s];
                dphidz = dphidz + dphi*kernels[i].dNdz[s];
                //dphiVxdx = dphiVxdx + (phi_iter[nbr]*Vx[nbr]-phi_iter[i]*Vx[i])*kernels[i].dNdx[s];
                //dphiVydy = dphiVydy + (phi_iter[nbr]*Vy[nbr]-phi_iter[i]*Vy[i])*kernels[i].dNdy[s];
                //dphiVzdz = dphiVzdz + (phi_iter[nbr]*Vz[nbr]-phi_iter[i]*Vz[i])*kernels[i].dNdz[s];
                nabla2_mu = nabla2_mu + (mu[nbr]-mu[i])*kernels[i].nabla2[s];
               // double V = sqrt(Vx[nbr]*Vx[nbr] + Vy[nbr]*Vy[nbr] + Vz[nbr]*Vz[nbr] + Vn[nbr]*Vn[nbr]);
               // if (V >= V_max)
               // {V_max = V;}
            }
          //  double Vi = sqrt(Vx[i]*Vx[i] + Vy[i]*Vy[i] + Vz[i]*Vz[i] + Vn[i]*Vn[i]);
          //  if (Vi >= V_max)
         //   {V_max = Vi;}
         //   double M = M0*V_max;

            double dphi_abs = sqrt(dphidx*dphidx + dphidy*dphidy + dphidz*dphidz);
            double ADV = dphidx*Vx[i]+dphidy*Vy[i]+dphidz*Vz[i];
            //double ADV = dphiVxdx+dphiVydy+dphiVzdz;
            phi_iter[i] = phi_now[i] + dt*(M0*nabla2_mu - ADV + dphi_abs*Vn[i]);
            if(phi_iter[i] > 1.0)
            {phi_iter[i]= 1.0;}
            if(phi_iter[i] < -1.0)
            {phi_iter[i] = -1.0;}
        }
        U_store[iternum] = checkConvergence(cntr, phi_iter, phi_iter_old, TotalPoints);
    }
    /*
    cout << "Checking convergence..." << endl;
    for (int i = 2; i <= iterations; i++)
    {cout<< "U_store["<<i<<"] = " << U_store[i] << endl;}
    */
    return phi_iter;
}
vector<double> computeAllenCahnPhaseTransformation2D(int dim, vector<double> phi, vector<double> T1, vector<double> D, INTERFACE_PRMTRS inter_param, double dt, vector<POINT> particles, int TotalParticles, vector<KERNEL> kernels)
{
    double a_2 = inter_param.surface_tension.b;
    double eps_4 = inter_param.surface_tension.d0;
    double delta = inter_param.surface_tension.A;
    double theta0 = inter_param.surface_tension.theta;
    double alpha = inter_param.PFM.alpha;
    double gamma = inter_param.PFM.gamma;
    double tau0 = inter_param.PFM.tau0;

    vector<double> phi_fut(TotalParticles+1);
    vector<double> epsilon(TotalParticles+1);
    vector<double> epsilon_termx(TotalParticles+1);
    vector<double> epsilon_termy(TotalParticles+1);

#pragma omp parallel for
    for (int i = 1; i <= TotalParticles; i++)
    {
        double dphidx = 0;
        double dphidy = 0;
        for (int s = 1; s <= kernels[i].TotalNbrs; s++)
        {
            int nbr = kernels[i].nbrs[s];
            dphidx = dphidx + kernels[i].dNdx[s]*(phi[nbr]-phi[i]);
            dphidy = dphidy + kernels[i].dNdy[s]*(phi[nbr]-phi[i]);
        }
        double theta = atan2(dphidy,dphidx);
        epsilon[i] = eps_4*(1.0+delta*cos(a_2*(theta-theta0)));
        double epsilon_deriv = -eps_4*a_2*delta*sin(a_2*(theta-theta0));
        epsilon_termx[i] = epsilon[i]*epsilon_deriv*dphidx;
        epsilon_termy[i] = epsilon[i]*epsilon_deriv*dphidy;
    }
#pragma omp parallel for
    for (int i = 1; i <= TotalParticles; i++)
    {
        double epsilon_termx_dy = 0;
        double epsilon_termy_dx = 0;
        for (int s = 1; s <= kernels[i].TotalNbrs; s++)
        {
            int nbr =  kernels[i].nbrs[s];
            epsilon_termx_dy = epsilon_termx_dy + kernels[i].dNdy[s]*(epsilon_termx[nbr]-epsilon_termx[i]);
            epsilon_termy_dx = epsilon_termy_dx + kernels[i].dNdx[s]*(epsilon_termy[nbr]-epsilon_termy[i]);
        }
        double nabla2_phi = 0;
        for (int j = 1; j <= kernels[i].TotalNbrs; j++)
        {
            int nbr = kernels[i].nbrs[j];
            nabla2_phi = nabla2_phi + (phi[nbr]-phi[i])*kernels[i].nabla2[j];
        }
        double m = alpha/3.14 * atan(gamma*(inter_param.thermal.T_eq-T1[i]));
        double Q = epsilon_termx_dy - epsilon_termy_dx + phi[i]*(1.0-phi[i])*(phi[i] -0.5 + m);
        phi_fut[i] = phi[i] + (dt/tau0)*pow(epsilon[i],2)*nabla2_phi + (dt/tau0)*Q;
        if(phi_fut[i] > 1.0)
        {phi_fut[i] = 1.0;}
        if(phi_fut[i] < 0)
        {phi_fut[i] = 0;}
    }
    phi = phi_fut;
    return phi;
}
vector<double> computeSharpenInterface(bool vof, double beta, double epsilon, vector<double> phi, int TotalPoints, vector<KERNEL> kernels, double dx,  double dt, int nt)
{
    vector<double> phi_fut(TotalPoints+1);
    vector<double> term_x(TotalPoints+1);
    vector<double> term_y(TotalPoints+1);
    vector<double> term_z(TotalPoints+1);
    double dxdx_inv = 1.0/(dx*dx);
    for (int k = 1; k <= nt; k++)
    {
        //#pragma omp parallel for
        for (int i = 1; i <= TotalPoints; i++)
        {
            double dphidx = 0;
            double dphidy = 0;
            double dphidz = 0;
            for (int s = 1; s <= kernels[i].TotalNbrs; s++)
            {
                int nbr = kernels[i].nbrs[s];
                dphidx = dphidx + kernels[i].dNdx[s]*(phi[nbr]-phi[i]);
                dphidy = dphidy + kernels[i].dNdy[s]*(phi[nbr]-phi[i]);
                dphidz = dphidz + kernels[i].dNdz[s]*(phi[nbr]-phi[i]);
            }
            double abs_phi = (sqrt(dphidx*dphidx+dphidy*dphidy+dphidz*dphidz)+1e-6);
            double nx= dphidx/abs_phi;
            double ny= dphidy/abs_phi;
            double nz= dphidz/abs_phi;
            if (vof == true)
	        {   //for phi varying between 0 to 1
                term_x[i] = phi[i]*(1-phi[i])*nx;
                term_y[i] = phi[i]*(1-phi[i])*ny;
                term_z[i] = phi[i]*(1-phi[i])*nz;
            }
            else
            {   //for phi varying between -1 to 1
                term_x[i] = 0.25*(phi[i]+1)*(1-phi[i])*nx;
                term_y[i] = 0.25*(phi[i]+1)*(1-phi[i])*ny;
                term_z[i] = 0.25*(phi[i]+1)*(1-phi[i])*nz;
            }
        }
        for (int i = 1; i <= TotalPoints; i++)
        {
            double dtermxdx = 0;
            double dtermydy = 0;
            double dtermzdz = 0;
            double nabla_phi = 0;
            for (int s = 1; s <= kernels[i].TotalNbrs; s++)
            {
                int nbr = kernels[i].nbrs[s];
                dtermxdx = dtermxdx + kernels[i].dNdx[s]*(term_x[nbr]-term_x[i]);
                dtermydy = dtermydy + kernels[i].dNdy[s]*(term_y[nbr]-term_y[i]);
                dtermzdz = dtermzdz + kernels[i].dNdz[s]*(term_z[nbr]-term_z[i]);
                nabla_phi = nabla_phi + dxdx_inv*kernels[i].N[s]*(phi[nbr]-phi[i]);
            }
            phi_fut[i] = phi[i] + dt*beta*(epsilon*nabla_phi - (dtermxdx+dtermydy+dtermzdz));
        }
        phi = phi_fut;
    }
    return phi;
}

vector<double> computeCurvature(vector<double> phi, int TotalPoints, vector<KERNEL> kernels, int dim)
{
    //computes curvature using phi of the pointset not
    double eps = 0.000001;
    vector<double> curvature(TotalPoints+1);
    for (int i = 1; i <= TotalPoints; i++)
    {
        double dphidx=0;
        double dphidy=0;
        double dphidz=0;
        double dphidxx=0;
        double dphidyy=0;
        double dphidzz=0;
        double dphidxy=0;
        double dphidxz=0;
        double dphidyz=0;
        for (int j = 1; j <=kernels[i].TotalNbrs; j++)
        {
            int nbr = kernels[i].nbrs[j];
            double phi_diff = phi[i] - phi[nbr];
            dphidx = dphidx + phi_diff*kernels[i].dNdx[j];
            dphidy = dphidy + phi_diff*kernels[i].dNdy[j];
            dphidz = dphidz + phi_diff*kernels[i].dNdz[j];
            dphidxx = dphidxx + phi_diff*kernels[i].dNdxx[j];
            dphidyy = dphidyy + phi_diff*kernels[i].dNdyy[j];
            dphidzz = dphidzz + phi_diff*kernels[i].dNdzz[j];
            dphidxy = dphidxy + phi_diff*kernels[i].dNdxy[j];
            dphidxz = dphidxz + phi_diff*kernels[i].dNdxz[j];
            dphidyz = dphidyz + phi_diff*kernels[i].dNdyz[j];
        }
        double absgradphi = sqrt(dphidx*dphidx + dphidy*dphidy + dphidz*dphidz)+eps;
        double term1 = dphidx*dphidx*dphidyy;
        double term2 = dphidx*dphidy*dphidxy;
        double term3 = dphidy*dphidy*dphidxx;
        double term4 = dphidx*dphidx*dphidzz;
        double term5 = dphidx*dphidz*dphidxz;
        double term6 = dphidz*dphidz*dphidxx;
        double term7 = dphidy*dphidy*dphidzz;
        double term8 = dphidy*dphidz*dphidyz;
        double term9 = dphidz*dphidz*dphidyy;
        curvature[i] = (term1-2*term2+term3+term4-2*term5+term6+term7-2*term8+term9)/(pow(absgradphi,3));
    }
    return curvature;
}
