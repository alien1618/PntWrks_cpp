#include "../pntwrks.h"

void solveNavierStokes2Phase(POINTSET pointset, MATERIALS materials, INTERFACE_PRMTRS inter_param, BOUNDARY_CONDITIONS BC, SOLVER_SETTINGS settings, vector<double> U, vector<double> Vx, vector<double> Vy, vector<double> Vz, vector<double> Q, vector<double> Fx, vector<double> Fy, vector<double> Fz)
{
    cout << "SOLVER START..." << endl;
    cout << "Initializing variables..." << endl;
    pltctrl(pointset.points, pointset.TotalPoints, settings.nt, settings.prnt_freq);
    if (pointset.TotalSurfaces>0)
    {printElements(pointset.surfaces, pointset.SurfaceNds, pointset.TotalSurfaces, "elements");}

    // NOTE: pointset.points.phi must be at vof range (0,1)
    cout << "pointset.dim = " << pointset.dim << endl;
    cout << "pointset.AvgPntSpacing " << pointset.AvgPntSpacing << endl;
    int c = 1;
    int cc = 1;
    double W = pointset.AvgPntSpacing;
    double M = inter_param.PFM.mobility*W*W;
    double kernel_radius = settings.kernel.radius_ratio*pointset.AvgPntSpacing;
    int V_Nt = 1000;
    int phi_iter = 5;
    double sum_VOF0 = 0;
    Vx = setVector(Vx, BC.Vx);
    Vy = setVector(Vy, BC.Vy);
    Vz = setVector(Vz, BC.Vz);
    vector<double> ro(pointset.TotalPoints+1);
    vector<double> P(pointset.TotalPoints+1);
    vector<double> V(pointset.TotalPoints+1);
    vector<double> Vn_ext = setVector(pointset.TotalPoints, 0);
    vector<double> RHS = setVector(pointset.TotalPoints, 0);
    vector<double> phi = pointset.phi;
    for (int i = 1; i <= pointset.TotalPoints; i++)
    {sum_VOF0 = sum_VOF0+pointset.phi[i];}

    int method = 1;
    if (inter_param.solver == "VOF")
    {method = 1;}
    else if(inter_param.solver == "AC")
    {method = 2;}
    else if(inter_param.solver == "AC2")
    {method = 3;}
    else if(inter_param.solver == "CH")
    {method = 4;}
    else if(inter_param.solver == "CH2")
    {method = 5;}
    else
    {
        cout << "inter_param.solver = " << inter_param.solver << endl;
        cout << "ERROR: inter_param.solver is undefined..." << endl;
        exit(0);
    }

    if (inter_param.solver == "AC2" || inter_param.solver == "CH2")
    {
        //pointset.points.phi is given as vof (0,1), convert to heaveside function for interface solvers
        for (int i = 1; i <= pointset.TotalPoints; i++)
        {phi[i] = 2*pointset.phi[i] - 1;}
    }
    for (int i = 1; i <= pointset.TotalPoints; i++)
    {
        ro[i] = pointset.phi[i]*materials.ro[2] + (1-pointset.phi[i])*materials.ro[1];
        P[i] = 0;
        V[i] = sqrt(Vx[i]*Vx[i] + Vy[i]*Vy[i] + Vz[i]*Vz[i]);
    }

    if (materials.D[1] > 0)
    {
        U = setVector(U, BC.U);
        pointset.printMeshVTK(U, "U", 0);
    }
	
    cout << "Printing field variables at time 0..." << endl;
    pointset.printVectorsVTK(V,Vx,Vy,Vz,"vectors",0);
    pointset.printMeshVTK(V, "V", 0);
    pointset.printMeshVTK(ro, "ro", 0);
    pointset.printMeshVTK(P, "P", 0);
    pointset.printMeshVTK(phi, "phi", 0);

    cout << "Computing kernel approximation functions..." << endl;
    vector<KERNEL> kernels = computeKernels(pointset.points, pointset.phi, pointset.TotalPoints, pointset.dim, settings.kernel, kernel_radius);

    cout << "Starting transient solver..." << endl;
    for (int t = 1; t <= settings.nt; t++)
    {
        //cout << "Computing CFL condition..." << endl;
        double V_max = 0;
        for (int i = 1; i <= pointset.TotalPoints; i++)
        {
            double V = sqrt(Vx[i]*Vx[i]+Vy[i]*Vy[i]+Vz[i]*Vz[i]);
            if (V >= V_max)
            {V_max = V;}
            if (V >= 100)
            {
                Vx[i] = 0;
                Vy[i] = 0;
                Vz[i] = 0;
            }
        }
        double CFL = (V_max*settings.dt)/pointset.AvgPntSpacing;
        if (CFL >= 0.001)
        {settings.dt=settings.dt*0.5;}

        if (materials.D[1] > 0)
        {
            //cout << "Computing thermal field distribution..." << endl;
            double nu = sqrt(settings.Pr/settings.Ra)*settings.dt/(pointset.AvgPntSpacing*pointset.AvgPntSpacing);
            double k = sqrt(1.0/(settings.Pr*settings.Ra))*settings.dt/(pointset.AvgPntSpacing*pointset.AvgPntSpacing);
            materials.D[1] = k;
            materials.D[2] = k;
            materials.Nu[1] = nu;
            materials.Nu[2] = nu;
            
            double alpha = 0.01;
            double ro0 = materials.ro[1];
            double T0 = 0.5*(settings.Thot+settings.Tcold);
            for (int i = 1; i <= pointset.TotalPoints; i++)
            {ro[i] = ro0 - ro0*alpha*(U[i]-T0);}

            //cout << "Solving the thermal energy equation..." << endl;
            tie(U, RHS) = computeTransport(pointset.dim, U, Vx, Vy, Vz, Q, materials.D, settings.dt, t, RHS, pointset.points, pointset.phi, pointset.TotalPoints, kernels, BC.U);

            //cout << "Computing external forces terms..." << endl;
            for(int i = 1; i <= pointset.TotalPoints; i++)
            {
                Fx[i] = settings.gravity.x*settings.T_factor*(U[i]-T0)/(settings.Thot-settings.Tcold);
                Fy[i] = settings.gravity.y*settings.T_factor*(U[i]-T0)/(settings.Thot-settings.Tcold);
                Fz[i] = settings.gravity.z*settings.T_factor*(U[i]-T0)/(settings.Thot-settings.Tcold);
            }
        }

        //cout << "Solving the Navier-Stokes equations..." << endl;
        tie(Vx,Vy,Vz,V,P,ro) = solveNavierStokes(P, Vx, Vy, Vz, Fx, Fy, Fz, phi, pointset.points, pointset.TotalPoints, kernels, BC.Vx, BC.Vy, BC.Vz, BC.P, materials.Nu, materials.ro, materials.total, settings, pointset.AvgPntSpacing, inter_param.surface_tension.surface_tension_co);

        if(t/V_Nt == cc )
        {
            //cout << "Solving the phase-field/volume of fluid equations..." << endl;
            vector<double> phi_old = phi;
            //cout << "method " << method << endl;
            switch(method)
            {
                case 1:
            {phi = computeVolumeOfFluid(kernels, phi, pointset.TotalPoints, Vx, Vy, Vz, Vn_ext, phi_iter, V_Nt*settings.dt);}break;
                case 2:
            {phi = computeAllenCahn(W, M, kernels, phi, pointset.dim, pointset.points, pointset.TotalPoints, Vx, Vy, Vz, Vn_ext, phi_iter, V_Nt*settings.dt);}break;
                case 3:
            {phi = computeAllenCahn2(W, M, kernels, phi, pointset.dim, pointset.points, pointset.TotalPoints, Vx, Vy, Vz, Vn_ext, phi_iter, V_Nt*settings.dt);}break;
                case 4:
            {phi = computeCahnHilliard(W, M, kernels, phi, pointset.dim, pointset.points, pointset.TotalPoints, Vx, Vy, Vz, Vn_ext, phi_iter, V_Nt*settings.dt);}break;
                case 5:
            {
				//cout << "W " << W << " M " << M << endl;
				phi = computeCahnHilliard2(W, M, kernels, phi, pointset.dim, pointset.points, pointset.TotalPoints, Vx, Vy, Vz, Vn_ext, phi_iter, V_Nt*settings.dt);
				}break;
            }
            if (inter_param.sharpen_interface == true)
            {
                if(method == 1 || method == 2 || method == 4)
                {phi = computeSharpenInterface(true, inter_param.sharpen_interface_beta, inter_param.sharpen_interface_eps, phi, pointset.TotalPoints, kernels, pointset.AvgPntSpacing, inter_param.sharpen_interface_dt, inter_param.sharpen_interface_nt);}
                else
                {phi = computeSharpenInterface(false, inter_param.sharpen_interface_beta, inter_param.sharpen_interface_eps, phi, pointset.TotalPoints, kernels, pointset.AvgPntSpacing, inter_param.sharpen_interface_dt, inter_param.sharpen_interface_nt);}
            }

            phi = setVector(phi, BC.phi);
            if (method == 3 || method == 5)
            {
                for (int i = 1; i <= pointset.TotalPoints; i++)
                {
                    if (phi[i] < -1.0)
                    {phi[i] = -1;}
                    if (phi[i] > 1.0)
                    {phi[i] = 1;}
                }
            }
            else
            {
                for (int i = 1; i <= pointset.TotalPoints; i++)
                {
                    if (phi[i] < 0)
                    {phi[i] = 0;}
                    if (phi[i] > 1.0)
                    {phi[i] = 1.0;}
                }
            }

            //===============================================================
            if (inter_param.interp_at_bndry == true)
            {
                //interpolate interface at boundaries with zero velocity
                for (int i = 1; i <= pointset.TotalPoints; i++)
                {
                    if (abs(V[i]) <= 1e-8 && pointset.phi[i] <= 0.99999 && pointset.phi[i] >= 1e-6)
                    {
                        double sum_w = 0;
                        double sum_phi = 0;
                        for (int j = 1; j <= kernels[i].TotalNbrs; j++)
                        {
                            int nbr = kernels[i].nbrs[j];
                            double d = pointset.points[i].distance(pointset.points[nbr]);
                            if (d >= 1e-5 && d <= 1.1*pointset.AvgPntSpacing && abs(V[nbr]) >= 1e-8)
                            {
                                double w = 1/d;
                                sum_w = sum_w + w;
                                sum_phi = sum_phi + w*phi[nbr];
                            }
                        }

                        if (sum_w >= 1e-6)
                        {phi[i] = sum_phi/sum_w;}
                    }
                }
            }
            //===============================================================
            if (method == 3 || method == 5)
            {   //convert heavside function (-1,1) to VOF (0,1) for the NS solver
                for (int i = 1; i <= pointset.TotalPoints; i++)
                {pointset.phi[i] = 0.5*(phi[i]+1.0);}
            }
            else
            {pointset.phi = phi;}
            cc++;
        }

        if (t/settings.prnt_freq == c)
        {
            c++;
            cout << "Computation of time step " << t << " of " << settings.nt <<  " complete..." << endl;

            //cout << "Printing solutions to file..." << endl;
            if (pointset.TotalSurfaces == 0 && pointset.TotalVolumes == 0)
            {
                pointset.printTXT(V, "V", t);
                pointset.printVectorsVTK(V,Vx,Vy,Vz, "vectors", t);
                pointset.printTXT(P, "P", t);
                pointset.printTXT(ro, "ro", t);
                pointset.printTXT(phi, "phi", t);
                if (materials.D[1] > 0)
                {pointset.printTXT(U, "U", t);}
            }
            else
            {
                pointset.printMeshVTK(V, "V", t);
                pointset.printVectorsVTK(V, Vx, Vy, Vz, "vectors", t);
                pointset.printMeshVTK(P, "P", t);
                pointset.printMeshVTK(ro, "ro", t);
                pointset.printMeshVTK(phi, "phi", t);
                if (materials.D[1] > 0)
                {pointset.printMeshVTK(U, "U", t);}
            }

            //cout << "Computing performance monitoring data..." << endl;
            double div_sum = 0;
            double phi_sum = 0;
            double sum_VOF = 0;
            for (int i = 1; i <= pointset.TotalPoints; i++)
            {
                double term = 0;
                for (int s = 1; s<= kernels[i].TotalNbrs; s++)
                {
                    int nbr = kernels[i].nbrs[s];
                    double Vab_x = Vx[i]-Vx[nbr];
                    double Vab_y = Vy[i]-Vy[nbr];
                    double Vab_z = Vz[i]-Vz[nbr];
                    term = term + (Vab_x*kernels[i].dNdx[s]+Vab_y*kernels[i].dNdy[s]+Vab_z*kernels[i].dNdz[s]);
                }
                div_sum = div_sum+term;
                phi_sum = phi_sum+phi[i];
                sum_VOF = sum_VOF+pointset.phi[i];
            }
            cout << "Divergence sum = " << div_sum << endl;
            cout << "Percentage of volume loss = " << 100*abs(sum_VOF-sum_VOF0)/sum_VOF0 << endl;
            cout << "Maximum velocity = " << V_max << endl;
            cout << "CFL = " << CFL << endl;
        }
    }
    cout << "SOLVER COMPLETE" << endl;
}
