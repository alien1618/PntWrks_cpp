#include "../pntwrks.h"

void solveNavierStokes1Phase(POINTSET pointset, MATERIALS materials, BOUNDARY_CONDITIONS BC, SOLVER_SETTINGS settings, vector<double> U, vector<double> Vx, vector<double> Vy, vector<double> Vz, vector<double> Q, vector<double> Fx,vector<double> Fy,vector<double> Fz)
{
    cout << "SOLVER START..." << endl;
    cout << "Initializing variables..." << endl;
    pltctrl(pointset.points, pointset.TotalPoints, settings.nt, settings.prnt_freq);
    if (pointset.TotalSurfaces>0)
    {printElements(pointset.surfaces, pointset.SurfaceNds, pointset.TotalSurfaces, "elements");}

    cout << "pointset.dim = " << pointset.dim << endl;
    cout << "pointset.AvgPntSpacing " << pointset.AvgPntSpacing << endl;
    int c = 1;
    double kernel_radius = settings.kernel.radius_ratio*pointset.AvgPntSpacing;
    
	Vx = setVector(Vx, BC.Vx);
    Vy = setVector(Vy, BC.Vy);
    Vz = setVector(Vz, BC.Vz);
	vector<double> V (pointset.TotalPoints+1);
    vector<double> P = setVector(pointset.TotalPoints, 0);
    vector<double> ro =setVector(pointset.TotalPoints,materials.ro[1]);
    vector<double> RHS = setVector(pointset.TotalPoints,0);
    for (int i = 1; i <= pointset.TotalPoints; i++)
    {V[i] = sqrt(Vx[i]*Vx[i] + Vy[i]*Vy[i] + Vz[i]*Vz[i]);}
    if (materials.D[1] > 0)
    {U = setVector(U, BC.U);}

    cout << "Printing field variables at time 0..." << endl;
    if (pointset.TotalSurfaces == 0 && pointset.TotalVolumes == 0)
    {
        pointset.printTXT(V, "V", 0);
        pointset.printVectorsVTK(V, Vx,Vy,Vz, "vectors", 0);
        pointset.printTXT(P, "P", 0);
        pointset.printTXT(ro, "ro", 0);
        if (materials.D[1] > 0)
        {pointset.printTXT(U, "U", 0);}
    }
    else
    {
        pointset.printMeshVTK(V, "V", 0);
        pointset.printVectorsVTK(V,Vx,Vy,Vz,"vectors", 0);
        pointset.printMeshVTK(P, "P", 0);
        pointset.printMeshVTK(ro, "ro", 0);
        if (materials.D[1] > 0)
        {pointset.printMeshVTK(U, "U", 0);}
    }

    cout << "Computing kernel approximation functions..." << endl;
    vector<KERNEL> kernels = computeKernels(pointset.points, pointset.phi, pointset.TotalPoints, pointset.dim, settings.kernel, kernel_radius);

    cout << "Computing transient solution..." << endl;
    for (int t = 1; t <= settings.nt; t++)
    {
        if (materials.D[1] > 0)//solve the energy equation
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
            double T0 = 0.5*(settings.Thot-settings.Tcold);
			for (int i = 1; i <= pointset.TotalPoints; i++)
            {ro[i] = ro0 - ro0*alpha*(U[i]-T0);}

            //cout << "Solving the thermal energy equation..." << endl;
            tie(U, RHS) = computeTransport(pointset.dim, U, Vx, Vy, Vz, Q, materials.D, settings.dt, t, RHS, pointset.points, pointset.phi, pointset.TotalPoints, kernels, BC.U);
            for (int i = 1; i <= pointset.TotalPoints; i++)
            {
                if (pointset.phi[i] <= 0)
                {U[i] = settings.Thot+settings.Thot/5.0;}
            }

            //cout << "Computing external forces terms..." << endl;
            for(int i = 1; i <= pointset.TotalPoints; i++)
            {
                Fx[i] = settings.gravity.x*settings.T_factor*(U[i]-T0)/(settings.Thot-settings.Tcold);
                Fy[i] = settings.gravity.y*settings.T_factor*(U[i]-T0)/(settings.Thot-settings.Tcold);
                Fz[i] = settings.gravity.z*settings.T_factor*(U[i]-T0)/(settings.Thot-settings.Tcold);
            }
        }
        //cout << "Solving the Navier-Stokes equations..." << endl;
        tie(Vx,Vy,Vz,V,P,ro) = solveNavierStokes(P, Vx, Vy, Vz, Fx, Fy, Fz, pointset.phi, pointset.points, pointset.TotalPoints, kernels, BC.Vx, BC.Vy, BC.Vz, BC.P, materials.Nu, materials.ro, materials.total, settings, pointset.AvgPntSpacing, 0);
        //tie(Vx,Vy,Vz,V,P,ro) = solveNavierStokesUpwind(P, Vx, Vy, Vz, Fx, Fy, Fz, pointset.phi, pointset.points, pointset.TotalPoints, kernels, BC.Vx, BC.Vy, BC.Vz, BC.P, materials.Nu, materials.ro, materials.total, settings, pointset.AvgPntSpacing, 0, pointset.dim);
        for (int i = 1; i <= pointset.TotalPoints; i++)
        {
            if (pointset.phi[i] <= 0)
            {
                Vx[i] = 0;
                Vy[i] = 0;
                Vz[i] = 0;
            }
        }
        if (t/settings.prnt_freq == c)
        {
            c++;
            cout << "Computation of time step " << t << " of " << settings.nt <<  " complete..." << endl;
            if (pointset.TotalSurfaces == 0 && pointset.TotalVolumes == 0)
            {
                pointset.printTXT(V, "V", t);
                pointset.printVectorsVTK(V,Vx,Vy,Vz,"vectors", t);
                pointset.printTXT(P, "P", t);
                pointset.printTXT(ro, "ro", t);
                if (materials.D[1] > 0)
                {pointset.printTXT(U, "U", t);}
            }
            else
            {
                pointset.printMeshVTK(V, "V", t);
                pointset.printVectorsVTK(V,Vx,Vy,Vz, "vectors", t);
                pointset.printMeshVTK(P, "P", t);
                pointset.printMeshVTK(ro, "ro", t);
                if (materials.D[1] > 0)
                {pointset.printMeshVTK(U, "U", t);}
            }
            //cout << "Computing performance monitoring data..." << endl;
            double div_sum = 0;
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
            }
            double V_max = 0;
            for (int i = 1; i <= pointset.TotalPoints; i++)
            {
                double V = sqrt(Vx[i]*Vx[i]+Vy[i]*Vy[i]+Vz[i]*Vz[i]);
                if (V >= V_max)
                {V_max = V;}
            }
            double CFL = (V_max*settings.dt)/pointset.AvgPntSpacing;
            cout << "Divergence sum = " << div_sum << endl;
            cout << "Maximum velocity = " << V_max << endl;
            cout << "CFL = " << CFL << endl;
        }
    }
    cout << "SOLVER COMPLETE" << endl;
}

