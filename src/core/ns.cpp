#include "../pntwrks.h"

tuple<vector<double>, vector<double>,vector<double>, vector<double>,vector<double>, vector<double> > solveNavierStokes(vector<double> P0, vector<double> Vx1, vector<double> Vy1, vector<double> Vz1, vector<double> Fx, vector<double> Fy, vector<double> Fz, vector<double> phi, vector<POINT> points, int TotalPoints, vector<KERNEL> kernels, BOUNDARY_CONDITION BCs_Vx, BOUNDARY_CONDITION BCs_Vy, BOUNDARY_CONDITION BCs_Vz, BOUNDARY_CONDITION BCs_P, vector<double> viscosity, vector<double> density, int TotalPhases, SOLVER_SETTINGS settings, double AvgPointSpacing, double segma)
{
    //SOLVER BASED ON THE ARTIFICIAL COMPRESSIBILITY AND PROJECTION METHOD
    int iter = 1;
    double eps = 0.00001;
    double roi, Nui;   

	//do not touch
    vector<double> Vx = Vx1;
    vector<double> Vy = Vy1;
    vector<double> Vz = Vz1;
    vector<double> P = P0;
    //do not touch
    
    vector<double> V(TotalPoints+1);
    vector<double> P_fut(TotalPoints+1);
    vector<double> ro(TotalPoints+1);
    vector<double> Vx_inter(TotalPoints+1);
    vector<double> Vy_inter(TotalPoints+1);
    vector<double> Vz_inter(TotalPoints+1);
    vector<double> dphidx;
    vector<double> dphidy;
    vector<double> dphidz;
    vector<double> curv;
    double adv_x = 0;
    double adv_y = 0;
    double adv_z = 0;
    //============================================================================
    // Computing curvature for two-phase flows with surface tension
    //============================================================================
    if (segma > 0)
    {
        curv = computeCurvature(phi, TotalPoints, kernels, 2*AvgPointSpacing);
        dphidx.resize(TotalPoints+1);
        dphidy.resize(TotalPoints+1);
        dphidz.resize(TotalPoints+1);
    }//make sure interpolants have xy, yx, xz

    //============================================================================
    // Computing intermediate velocities
    //============================================================================
#pragma omp parallel for
    for (int i = 1; i <= TotalPoints; i++)
    {
        double dVxdx = 0;
        double dVxdy = 0;
        double dVxdz = 0;
        double dVydx = 0;
        double dVydy = 0;
        double dVydz = 0;
        double dVzdx = 0;
        double dVzdy = 0;
        double dVzdz = 0;
        double tension_x = 0;
        double tension_y = 0;
        double tension_z = 0;

        for (int j = 1; j <= kernels[i].TotalNbrs; j++)
        {
            int pntnum = kernels[i].nbrs[j];
            double Vx_diff = Vx[pntnum]-Vx[i];
            dVxdx = dVxdx + kernels[i].dNdx[j]*Vx_diff;
            dVxdy = dVxdy + kernels[i].dNdy[j]*Vx_diff;
            dVxdz = dVxdz + kernels[i].dNdz[j]*Vx_diff;

			double Vy_diff = Vy[pntnum]-Vy[i];
            dVydx = dVydx + kernels[i].dNdx[j]*Vy_diff;
            dVydy = dVydy + kernels[i].dNdy[j]*Vy_diff;
            dVydz = dVydz + kernels[i].dNdz[j]*Vy_diff;

			double Vz_diff = Vz[pntnum]-Vz[i];
            dVzdx = dVzdx + kernels[i].dNdx[j]*Vz_diff;
            dVzdy = dVzdy + kernels[i].dNdy[j]*Vz_diff;
            dVzdz = dVzdz + kernels[i].dNdz[j]*Vz_diff;

            if (segma > 0)
            {
                double phi_diff = phi[pntnum]-phi[i];
                dphidx[i] = dphidx[i] + kernels[i].dNdx[j]*phi_diff;
                dphidy[i] = dphidy[i] + kernels[i].dNdy[j]*phi_diff;
                dphidz[i] = dphidz[i] + kernels[i].dNdz[j]*phi_diff;
            }
        }


        double nabla2_Vx = 0;
        double nabla2_Vy = 0;
        double nabla2_Vz = 0;
        for (int j = 1; j <= kernels[i].TotalNbrs; j++)
        {
            int pntnum = kernels[i].nbrs[j];
            nabla2_Vx = nabla2_Vx + (Vx[pntnum]-Vx[i])*kernels[i].nabla2[j];
            nabla2_Vy = nabla2_Vy + (Vy[pntnum]-Vy[i])*kernels[i].nabla2[j];
            nabla2_Vz = nabla2_Vz + (Vz[pntnum]-Vz[i])*kernels[i].nabla2[j];
        }
        if (settings.lagrangian == false)
        {
			//===========================================================
			// computing acceleration due to advection
			//===========================================================
			adv_x = Vx[i]*dVxdx + Vy[i]*dVxdy + Vz[i]*dVxdz;
			adv_y = Vx[i]*dVydx + Vy[i]*dVydy + Vz[i]*dVydz;
			adv_z = Vx[i]*dVzdx + Vy[i]*dVzdy + Vz[i]*dVzdz;
		}
        //===========================================================
        // computing acceleration due to viscosity
        //===========================================================
        if (TotalPhases == 1)
        {Nui = viscosity[1];}
        else
        {Nui = phi[i]*viscosity[2] + (1-phi[i])*viscosity[1];}

        double dif_x = Nui*nabla2_Vx;
        double dif_y = Nui*nabla2_Vy;
        double dif_z = Nui*nabla2_Vz;

        if (segma > 0)
        {
            //===========================================================
            // computing acceleration due to surface tension
            //===========================================================
		    double abs_dphi = sqrt(dphidx[i]*dphidx[i]+dphidy[i]*dphidy[i]+dphidz[i]*dphidz[i]);
            tension_x = segma*curv[i]*dphidx[i]/(abs_dphi+eps);
            tension_y = segma*curv[i]*dphidy[i]/(abs_dphi+eps);
            tension_z = segma*curv[i]*dphidz[i]/(abs_dphi+eps);
            //===========================================================
        }

        //===========================================================
        // computing intermediate velocities
        //===========================================================
        Vx_inter[i] = Vx[i] + settings.dt*(-adv_x + dif_x + tension_x + Fx[i]);
        Vy_inter[i] = Vy[i] + settings.dt*(-adv_y + dif_y + tension_y + Fy[i]);
        Vz_inter[i] = Vz[i] + settings.dt*(-adv_z + dif_z + tension_z + Fz[i]);
    }

    //============================================================================
    // Computing the pressure iteratively
    //============================================================================
    double Vmax = 0;
    for (int i = 1; i <= TotalPoints; i++)
    {
		double V0 = sqrt(Vx_inter[i] * Vx_inter[i]  + Vy_inter[i] * Vy_inter[i]  + Vz_inter[i] * Vz_inter[i]);
        if (V0 >= Vmax)
        {Vmax = V0;}
    }
    double c = 10*Vmax;
        
    for (int vt = 1; vt <= iter; vt++)
    {
        // computing pressure using the artificial compressibility method
#pragma omp parallel for
        for (int i = 1; i <= TotalPoints; i++)
        {
            if (TotalPhases == 1)
            {roi = density[1];}
            else
            {roi = phi[i]*density[2] + (1-phi[i])*density[1];}
            double term = 0;
            for (int s = 1; s<= kernels[i].TotalNbrs; s++)
            {
                int nbr = kernels[i].nbrs[s];
                double Vab_x = Vx_inter[i]-Vx_inter[nbr];
                double Vab_y = Vy_inter[i]-Vy_inter[nbr];
                double Vab_z = Vz_inter[i]-Vz_inter[nbr];
                term = term + (Vab_x*kernels[i].dNdx[s]+Vab_y*kernels[i].dNdy[s]+Vab_z*kernels[i].dNdz[s]);
            }
            P_fut[i]  = P[i] + settings.P_factor*settings.dt*roi*term;
        }
        P = P_fut;

        // diffusive term acts as stabilizer
#pragma omp parallel for
        for (int i = 1; i <= TotalPoints; i++)
        {
            double nabla2_P = 0;
            /*
            double dVxdx = 0;
            double dVxdy = 0;
            double dVxdz = 0;
            double dVydx = 0;
            double dVydy = 0;
            double dVydz = 0;
            double dVzdx = 0;
            double dVzdy = 0;
            double dVzdz = 0;
            */
            for (int j = 1; j <= kernels[i].TotalNbrs; j++)
            {
                int pntnum = kernels[i].nbrs[j];
                nabla2_P = nabla2_P + (P[pntnum]-P[i])*kernels[i].nabla2[j];
				/*
                dVxdx = dVxdx + kernels[i].dNdx[j]*(Vx_inter[pntnum]-Vx_inter[i]);
                dVxdy = dVxdy + kernels[i].dNdy[j]*(Vx_inter[pntnum]-Vx_inter[i]);
                dVxdz = dVxdz + kernels[i].dNdz[j]*(Vx_inter[pntnum]-Vx_inter[i]);

                dVydx = dVydx + kernels[i].dNdx[j]*(Vy_inter[pntnum]-Vy_inter[i]);
                dVydy = dVydy + kernels[i].dNdy[j]*(Vy_inter[pntnum]-Vy_inter[i]);
                dVydz = dVydz + kernels[i].dNdz[j]*(Vy_inter[pntnum]-Vy_inter[i]);

                dVzdx = dVzdx + kernels[i].dNdx[j]*(Vz_inter[pntnum]-Vz_inter[i]);
                dVzdy = dVzdy + kernels[i].dNdy[j]*(Vz_inter[pntnum]-Vz_inter[i]);
                dVzdz = dVzdz + kernels[i].dNdz[j]*(Vz_inter[pntnum]-Vz_inter[i]);
                */
            }
            P_fut[i] = P[i] + settings.dt*(nabla2_P);// + 2*(dVydx*dVxdy+dVydz*dVzdy+dVxdz*dVzdx-dVxdx*dVydy*dVzdz));
        }
        P=P_fut;
       }
       P = setVector(P, BCs_P);

        //#pragma omp parallel for// be careful...for some reason it screwed up the artificial viscocity term
        for (int i = 1;  i<= TotalPoints; i++)
        {
            double A_pressure_x = 0;
            double A_pressure_y = 0;
            double A_pressure_z = 0;
            double art_visc_x = 0;
            double art_visc_y = 0;
            double art_visc_z = 0;

            if (TotalPhases == 1)
            {roi = density[1];}
            else
            {roi = phi[i]*density[2] + (1-phi[i])*density[1];}

            for (int s = 1; s <= kernels[i].TotalNbrs; s++)
            {
                int j = kernels[i].nbrs[s];
                //===========================================================
                // computing acceleration due to pressure
                //===========================================================
                A_pressure_x = A_pressure_x + (P[i]+P[j])*kernels[i].dNdx[s];// do not touch
                A_pressure_y = A_pressure_y + (P[i]+P[j])*kernels[i].dNdy[s];// do not touch
                A_pressure_z = A_pressure_z + (P[i]+P[j])*kernels[i].dNdz[s];// do not touch

                //===========================================================
                // computing acceleration due to artificial viscosity
                //===========================================================
                double xij = points[i].x-points[j].x;
                double yij = points[i].y-points[j].y;
               double zij = points[i].z-points[j].z;
                double Vxij = Vx_inter[i]-Vx_inter[j];
                double Vyij = Vy_inter[i]-Vy_inter[j];
                double Vzij = Vz_inter[i]-Vz_inter[j];
                double uij_rij = xij*Vxij + yij*Vyij + zij*Vzij;

                double alpha = settings.AV_factor;
                double d0 = sqrt(xij*xij + yij*yij + zij*zij);
                double h0 = 2*d0;
                double term0 = ((alpha*c*h0)/roi)/(d0*d0 + h0*1e-5);
                if (settings.lagrangian == true)
                {term0 = ((alpha*c*h0))/(d0*d0+h0*1e-5);}
                double term_f = -uij_rij*term0;
                double pii0 = term_f;
                if (term_f < 0)
                {pii0 = 0;}
                art_visc_x = art_visc_x + pii0*kernels[i].dNdx[s];
                art_visc_y = art_visc_y + pii0*kernels[i].dNdy[s];
                art_visc_z = art_visc_z + pii0*kernels[i].dNdz[s];
            }

            //===========================================================
            // computing final velocity
            //===========================================================
            Vx[i] = Vx_inter[i]+settings.dt*(-(1/(roi+eps))*A_pressure_x-art_visc_x);
            Vy[i] = Vy_inter[i]+settings.dt*(-(1/(roi+eps))*A_pressure_y-art_visc_y);
            Vz[i] = Vz_inter[i]+settings.dt*(-(1/(roi+eps))*A_pressure_z-art_visc_z);
            V[i] = sqrt(Vx[i] * Vx[i]  + Vy[i] * Vy[i] + Vz[i] * Vz[i]);
        }

        //============================================================================
        // Applying boundary conditions
        //============================================================================
        Vx = setVector(Vx, BCs_Vx);
        Vy = setVector(Vy, BCs_Vy);
        Vz = setVector(Vz, BCs_Vz);  

    for(int i = 1; i <= TotalPoints; i++)
    {
        V[i] = sqrt(Vx[i] * Vx[i]  + Vy[i] * Vy[i] + Vz[i] * Vz[i]);
        if (TotalPhases == 1)
        {ro[i] = density[1];}
        else
        {ro[i] = phi[i]*density[2] + (1-phi[i])*density[1];}
    }    
    return make_tuple(Vx, Vy, Vz, V, P, ro);
}

tuple<vector<double>, vector<double>,vector<double>, vector<double>,vector<double>, vector<double> > solveNavierStokesUpwind(vector<double> P0, vector<double> Vx1, vector<double> Vy1, vector<double> Vz1, vector<double> Fx, vector<double> Fy, vector<double> Fz, vector<double> phi, vector<POINT> points, int TotalPoints, vector<KERNEL> kernels, BOUNDARY_CONDITION BCs_Vx, BOUNDARY_CONDITION BCs_Vy, BOUNDARY_CONDITION BCs_Vz, BOUNDARY_CONDITION BCs_P, vector<double> Nu, vector<double> ro0, int TotalPhases, SOLVER_SETTINGS settings, double AvgPointSpacing, double segma, int dim)
{
    //SOLVER BASED ON THE ARTIFICIAL COMPRESSIBILITY AND PROJECTION METHOD
    int iter = 1;
    double eps = 0.00001;
    double roi, Nui;
    
    //do not touch
    vector<double> Vx = Vx1;
    vector<double> Vy = Vy1;
    vector<double> Vz = Vz1;
    vector<double> P = P0;
    //do not touch
    
    vector<double> ro(TotalPoints+1);
    vector<double> Vx_inter(TotalPoints+1);
    vector<double> Vy_inter(TotalPoints+1);
    vector<double> Vz_inter(TotalPoints+1);
    vector<double> V(TotalPoints+1);
    vector<double> P_fut(TotalPoints+1);
    vector<double> dphidx;
    vector<double> dphidy;
    vector<double> dphidz;
    vector<double> curv;
    V.resize(TotalPoints+1);
    double adv_x = 0;
    double adv_y = 0;
    double adv_z = 0;

    //============================================================================
    // Computing curvature for two-phase flows with surface tension
    //============================================================================
    if (segma > 0)
    {
        curv = computeCurvature(phi, TotalPoints, kernels, 2*AvgPointSpacing);
        dphidx.resize(TotalPoints+1);
        dphidy.resize(TotalPoints+1);
        dphidz.resize(TotalPoints+1);
    }//make sure interpolants have xy, yx, xz

    //============================================================================
    // Computing intermediate velocities
    //============================================================================
#pragma omp parallel for
    for (int i = 1; i <= TotalPoints; i++)
    {
        double dVxdx_uw = 0;
        double dVxdy_uw = 0;
        double dVxdz_uw = 0;
        double dVydx_uw = 0;
        double dVydy_uw = 0;
        double dVydz_uw = 0;
        double dVzdx_uw = 0;
        double dVzdy_uw = 0;
        double dVzdz_uw = 0;
        double dVxdx_dw = 0;
        double dVxdy_dw = 0;
        double dVxdz_dw = 0;
        double dVydx_dw = 0;
        double dVydy_dw = 0;
        double dVydz_dw = 0;
        double dVzdx_dw = 0;
        double dVzdy_dw = 0;
        double dVzdz_dw = 0;
        double tension_x = 0;
        double tension_y = 0;
        double tension_z = 0;

        for (int j = 1; j <= kernels[i].TotalNbrs; j++)
        {
            int pntnum = kernels[i].nbrs[j];
            double xj_xi = points[pntnum].x-points[i].x;
            double yj_yi = points[pntnum].y-points[i].y;
            double zj_zi = points[pntnum].z-points[i].z;
            //double d = computeDistance(px[i], py[i], px[pntnum], py[pntnum]);
            //double d_sq = d*d;
            double check = Vx[i]*xj_xi+Vy[i]*yj_yi+Vz[i]*zj_zi;
            if (check < 0) //upwind
            {
                dVxdx_uw = dVxdx_uw + kernels[i].dNdx[j]*(Vx[pntnum]-Vx[i]);
                dVxdy_uw = dVxdy_uw + kernels[i].dNdy[j]*(Vx[pntnum]-Vx[i]);
                dVxdz_uw = dVxdz_uw + kernels[i].dNdz[j]*(Vx[pntnum]-Vx[i]);

                dVydx_uw = dVydx_uw + kernels[i].dNdx[j]*(Vy[pntnum]-Vy[i]);
                dVydy_uw = dVydy_uw + kernels[i].dNdy[j]*(Vy[pntnum]-Vy[i]);
                dVydz_uw = dVydz_uw + kernels[i].dNdz[j]*(Vy[pntnum]-Vy[i]);

                dVzdx_uw = dVzdx_uw + kernels[i].dNdx[j]*(Vz[pntnum]-Vz[i]);
                dVzdy_uw = dVzdy_uw + kernels[i].dNdy[j]*(Vz[pntnum]-Vz[i]);
                dVzdz_uw = dVzdz_uw + kernels[i].dNdz[j]*(Vz[pntnum]-Vz[i]);
            }
            else //downwind
            {
                dVxdx_dw = dVxdx_dw + kernels[i].dNdx[j]*(Vx[pntnum]-Vx[i]);
                dVxdy_dw = dVxdy_dw + kernels[i].dNdy[j]*(Vx[pntnum]-Vx[i]);
                dVxdz_dw = dVxdz_dw + kernels[i].dNdz[j]*(Vx[pntnum]-Vx[i]);

                dVydx_dw = dVydx_dw + kernels[i].dNdx[j]*(Vy[pntnum]-Vy[i]);
                dVydy_dw = dVydy_dw + kernels[i].dNdy[j]*(Vy[pntnum]-Vy[i]);
                dVydz_dw = dVydz_dw + kernels[i].dNdz[j]*(Vy[pntnum]-Vy[i]);

                dVzdx_dw = dVzdx_dw + kernels[i].dNdx[j]*(Vz[pntnum]-Vz[i]);
                dVzdy_dw = dVzdy_dw + kernels[i].dNdy[j]*(Vz[pntnum]-Vz[i]);
                dVzdz_dw = dVzdz_dw + kernels[i].dNdz[j]*(Vz[pntnum]-Vz[i]);
            }

            if (segma > 0)
            {
                double phi_diff = phi[pntnum]-phi[i];
                dphidx[i] = dphidx[i] + kernels[i].dNdx[j]*phi_diff;
                dphidy[i] = dphidy[i] + kernels[i].dNdy[j]*phi_diff;
                dphidz[i] = dphidz[i] + kernels[i].dNdz[j]*phi_diff;
            }
        }
        double dVxdx = dim*(settings.kernel.upwind_ratio*(dVxdx_uw)+(1-settings.kernel.upwind_ratio)*(dVxdx_dw));
        double dVxdy = dim*(settings.kernel.upwind_ratio*(dVxdy_uw)+(1-settings.kernel.upwind_ratio)*(dVxdy_dw));
        double dVxdz = dim*(settings.kernel.upwind_ratio*(dVxdz_uw)+(1-settings.kernel.upwind_ratio)*(dVxdz_dw));
        double dVydx = dim*(settings.kernel.upwind_ratio*(dVydx_uw)+(1-settings.kernel.upwind_ratio)*(dVydx_dw));
        double dVydy = dim*(settings.kernel.upwind_ratio*(dVydy_uw)+(1-settings.kernel.upwind_ratio)*(dVydy_dw));
        double dVydz = dim*(settings.kernel.upwind_ratio*(dVydz_uw)+(1-settings.kernel.upwind_ratio)*(dVydz_dw));
		double dVzdx = dim*(settings.kernel.upwind_ratio*(dVzdx_uw)+(1-settings.kernel.upwind_ratio)*(dVzdx_dw));
        double dVzdy = dim*(settings.kernel.upwind_ratio*(dVzdy_uw)+(1-settings.kernel.upwind_ratio)*(dVzdy_dw));
        double dVzdz = dim*(settings.kernel.upwind_ratio*(dVzdz_uw)+(1-settings.kernel.upwind_ratio)*(dVzdz_dw));
        double nabla2_Vx = 0;
        double nabla2_Vy = 0;
        double nabla2_Vz = 0;
        for (int j = 1; j <= kernels[i].TotalNbrs; j++)
        {
            int pntnum = kernels[i].nbrs[j];
            nabla2_Vx = nabla2_Vx + (Vx[pntnum]-Vx[i])*kernels[i].nabla2[j];
            nabla2_Vy = nabla2_Vy + (Vy[pntnum]-Vy[i])*kernels[i].nabla2[j];
            nabla2_Vz = nabla2_Vz + (Vz[pntnum]-Vz[i])*kernels[i].nabla2[j];
        }
        if (settings.lagrangian == false)
        {
			//===========================================================
			// computing acceleration due to advection
			//===========================================================
			adv_x = Vx[i]*dVxdx + Vy[i]*dVxdy + Vz[i]*dVxdz;
			adv_y = Vx[i]*dVydx + Vy[i]*dVydy + Vz[i]*dVydz;
			adv_z = Vx[i]*dVzdx + Vy[i]*dVzdy + Vz[i]*dVzdz;
		}
        //===========================================================
        // computing acceleration due to viscosity
        //===========================================================
        if (TotalPhases == 1)
        {Nui = Nu[1];}
        else
        {Nui = phi[i]*Nu[2] + (1-phi[i])*Nu[1];}

        double dif_x = Nui*nabla2_Vx;
        double dif_y = Nui*nabla2_Vy;
        double dif_z = Nui*nabla2_Vz;

        if (segma > 0)
        {
            //===========================================================
            // computing acceleration due to surface tension
            //===========================================================
		    double abs_dphi = sqrt(dphidx[i]*dphidx[i]+dphidy[i]*dphidy[i]+dphidz[i]*dphidz[i]);
            tension_x = segma*curv[i]*dphidx[i]/(abs_dphi+eps);
            tension_y = segma*curv[i]*dphidy[i]/(abs_dphi+eps);
            tension_z = segma*curv[i]*dphidz[i]/(abs_dphi+eps);
            //===========================================================
        }

        //===========================================================
        // computing intermediate velocities
        //===========================================================
        Vx_inter[i] = Vx[i] + settings.dt*(-adv_x + dif_x + tension_x + Fx[i]);
        Vy_inter[i] = Vy[i] + settings.dt*(-adv_y + dif_y + tension_y + Fy[i]);
        Vz_inter[i] = Vz[i] + settings.dt*(-adv_z + dif_z + tension_z + Fz[i]);
    }

    //============================================================================
    // Computing the pressure iteratively
    //============================================================================
    for (int vt = 1; vt <= iter; vt++)
    {
        double Vmax = 0;
        for (int i = 1; i <= TotalPoints; i++)
        {
            double V0 = sqrt(Vx_inter[i] * Vx_inter[i]  + Vy_inter[i] * Vy_inter[i]  + Vz_inter[i] * Vz_inter[i]);
            if (V0 >= Vmax)
            {Vmax = V0;}
        }
        double c = 10*Vmax;

        // computing pressure using the artificial compressibility method
#pragma omp parallel for
        for (int i = 1; i <= TotalPoints; i++)
        {
            if (TotalPhases == 1)
            {roi = ro0[1];}
            else
            {roi = phi[i]*ro0[2] + (1-phi[i])*ro0[1];}
            double term = 0;
            for (int s = 1; s<= kernels[i].TotalNbrs; s++)
            {
                int nbr = kernels[i].nbrs[s];
                double Vab_x = Vx_inter[i]-Vx_inter[nbr];
                double Vab_y = Vy_inter[i]-Vy_inter[nbr];
                double Vab_z = Vz_inter[i]-Vz_inter[nbr];
                term = term + (Vab_x*kernels[i].dNdx[s]+Vab_y*kernels[i].dNdy[s]+Vab_z*kernels[i].dNdz[s]);
            }
            P_fut[i]  = P[i] + settings.P_factor*settings.dt*roi*term;
        }
        P = P_fut;

        // diffusive term acts as stabilizer
#pragma omp parallel for
        for (int i = 1; i <= TotalPoints; i++)
        {
            double nabla2_P = 0;
            for (int j = 1; j <= kernels[i].TotalNbrs; j++)
            {
                int pntnum = kernels[i].nbrs[j];
                nabla2_P = nabla2_P + (P[pntnum]-P[i])*kernels[i].nabla2[j];
            }
            P_fut[i] = P[i] + settings.dt*nabla2_P;
        }
        P=P_fut;
        P = setVector(P, BCs_P);

        //#pragma omp parallel for// be careful...for some reason it screwed up the artificial viscocity term
        for (int i = 1;  i<= TotalPoints; i++)
        {
            double A_pressure_x = 0;
            double A_pressure_y = 0;
            double A_pressure_z = 0;
            double art_visc_x = 0;
            double art_visc_y = 0;
            double art_visc_z = 0;

            if (TotalPhases == 1)
            {roi = ro0[1];}
            else
            {roi = phi[i]*ro0[2] + (1-phi[i])*ro0[1];}

            for (int s = 1; s <= kernels[i].TotalNbrs; s++)
            {
                int j = kernels[i].nbrs[s];
                //===========================================================
                // computing acceleration due to pressure
                //===========================================================
                A_pressure_x = A_pressure_x + (P[i]+P[j])*kernels[i].dNdx[s]; // do not touch
                A_pressure_y = A_pressure_y + (P[i]+P[j])*kernels[i].dNdy[s];// do not touch
                A_pressure_z = A_pressure_z + (P[i]+P[j])*kernels[i].dNdz[s];// do not touch

                //===========================================================
                // computing acceleration due to artificial viscosity
                //===========================================================
                double xij = points[i].x-points[j].x;
                double yij = points[i].y-points[j].y;
                double zij = points[i].z-points[j].z;
                double Vxij = Vx_inter[i]-Vx_inter[j];
                double Vyij = Vy_inter[i]-Vy_inter[j];
                double Vzij = Vz_inter[i]-Vz_inter[j];
                double uij_rij = xij*Vxij + yij*Vyij + zij*Vzij;

                double alpha = settings.AV_factor;
                double d0 = sqrt(xij*xij + yij*yij + zij*zij);
                double h0 = 2*d0;
                double term0 = ((alpha*c*h0)/roi)/(d0*d0 + h0*1e-5);
                if (settings.lagrangian == true)
                {term0 = ((alpha*c*h0))/(d0*d0+h0*1e-5);}
                double term_f = -uij_rij*term0;
                double pii0 = term_f;
                if (term_f < 0)
                {pii0 = 0;}
                art_visc_x = art_visc_x + pii0*kernels[i].dNdx[s];
                art_visc_y = art_visc_y + pii0*kernels[i].dNdy[s];
                art_visc_z = art_visc_z + pii0*kernels[i].dNdz[s];
            }

            //===========================================================
            // computing final velocity
            //===========================================================
            Vx[i] = Vx_inter[i]+settings.dt*(-(1/(roi+eps))*A_pressure_x-art_visc_x);
            Vy[i] = Vy_inter[i]+settings.dt*(-(1/(roi+eps))*A_pressure_y-art_visc_y);
            Vz[i] = Vz_inter[i]+settings.dt*(-(1/(roi+eps))*A_pressure_z-art_visc_z);
            V[i] = sqrt(Vx[i] * Vx[i]  + Vy[i] * Vy[i]  + Vz[i] * Vz[i]);
        }

        //============================================================================
        // Applying boundary conditions
        //============================================================================
        Vx = setVector(Vx, BCs_Vx);
        Vy = setVector(Vy, BCs_Vy);
        Vz = setVector(Vz, BCs_Vz);
        Vx_inter = Vx;
        Vy_inter = Vy;
        Vz_inter = Vz;
    }

    for(int i = 1; i <= TotalPoints; i++)
    {
        V[i] = sqrt(Vx[i] * Vx[i]  + Vy[i] * Vy[i]  + Vz[i] * Vz[i]);
        if (TotalPhases == 1)
        {ro[i] = ro0[1];}
        else
        {ro[i] = phi[i]*ro0[2] + (1-phi[i])*ro0[1];}
    }

    Vx = setVector(Vx, BCs_Vx);
    Vy = setVector(Vy, BCs_Vy);
    Vz = setVector(Vz, BCs_Vz);
    
    return make_tuple(Vx, Vy, Vz, V, P, ro);
}
