#include "../pntwrks.h"

void solveAllenCahn(POINTSET pointset, MATERIALS materials, INTERFACE_PRMTRS inter_param, BOUNDARY_CONDITIONS BC, SOLVER_SETTINGS settings, vector<double> U, vector<double> Vx, vector<double> Vy, vector<double> Vz)
{
    cout << "SOLVER START..." << endl;
    cout << "Initializing variables..."<< endl;
    pltctrl(pointset.points, pointset.TotalPoints, settings.nt, settings.prnt_freq);
    if (pointset.TotalSurfaces>0)
    {printElements(pointset.surfaces, pointset.SurfaceNds, pointset.TotalSurfaces, "elements");}


    int c = 1;
    bool fluid_solver = false;
    double kernel_radius = pointset.AvgPntSpacing*settings.kernel.radius_ratio;
    vector<double> V0(pointset.TotalPoints+1);
    vector<double> Q(pointset.TotalPoints+1);
    vector<double> phi_old(pointset.TotalPoints+1);
	vector<double> RHS(pointset.TotalPoints+1);
	
    if (materials.Nu[1] > 0)
    {fluid_solver = true;}
    else
    {fluid_solver = false;}
    cout << "fluid solver = " << fluid_solver << endl;

	vector<double> V;
	vector<double> P;
	vector<double> ro;
    if (fluid_solver == true)
    {
        Vx = setVector(Vx, BC.Vx);
        Vy = setVector(Vy, BC.Vy);
        Vz = setVector(Vz, BC.Vz);
        ro.resize(pointset.TotalPoints+1);
        P.resize(pointset.TotalPoints+1);
        V.resize(pointset.TotalPoints+1);
        for (int i = 1; i <= pointset.TotalPoints; i++)
        {
            ro[i] = materials.ro[1];
            P[i] = 0;
            V[i] = sqrt(Vx[i]*Vx[i] + Vy[i]*Vy[i] + Vz[i]*Vz[i]);
        }
    }
    else
    {
        Vx = V0;
        Vy = V0;
        Vz = V0;
    }
    
    cout << "Printing field variables at time 0..." << endl;
    if (pointset.TotalSurfaces == 0 && pointset.TotalVolumes == 0)
    {
        pointset.printTXT(pointset.phi, "phi", 0);
        if (settings.simulation_type == "thermal")
        {pointset.printTXT(U, "U", 0);}
        if (settings.simulation_type == "hele_shaw")
        {pointset.printTXT(U, "P", 0);}
        if (fluid_solver == true)
        {
            pointset.printTXT(V, "V", 0);
            pointset.printVectorsVTK(V, Vx, Vy, Vz, "vectors", 0);
            pointset.printTXT(P, "P", 0);
            pointset.printTXT(ro, "ro", 0);
        }
    }
    else
    {
        pointset.printMeshVTK(pointset.phi, "phi", 0);
        if (settings.simulation_type == "thermal")
        {pointset.printMeshVTK(U, "U", 0);}
        if (settings.simulation_type == "hele_shaw")
        {pointset.printMeshVTK(U, "P", 0);}
        if (fluid_solver == true)
        {
            pointset.printMeshVTK(V, "V", 0);
            pointset.printVectorsVTK(V, Vx, Vy, Vz, "vectors", 0);
            pointset.printMeshVTK(P, "P", 0);
            pointset.printMeshVTK(ro, "ro", 0);
        }
    }

    cout << "Computing kernel approximation functions..." << endl;
    vector<KERNEL> kernels = computeKernels(pointset.points, pointset.phi, pointset.TotalPoints, pointset.dim, settings.kernel, kernel_radius);

    vector<KERNEL> kernels2;
    if (fluid_solver == true)
    {
        settings.kernel.approximation_method = "SPH";
        kernels2 = computeKernels(pointset.points, pointset.phi, pointset.TotalPoints, pointset.dim, settings.kernel, kernel_radius);
    }

    cout << "Computing transient solution..."<< endl;
    for (int t =1; t <= settings.nt; t++)
    {
        if (fluid_solver == true)
        {
            cout << "Computing flow solution for time step "<< t << " of " << settings.nt << endl;
            tie(Vx,Vy,Vz,V,P,ro) = solveNavierStokes(P, Vx, Vy, Vz, V0, V0, V0, pointset.phi, pointset.points, pointset.TotalPoints, kernels2, BC.Vx, BC.Vy, BC.Vz, BC.P, materials.Nu, materials.ro, materials.total, settings, pointset.AvgPntSpacing, 0);
            for (int i = 1; i <= pointset.TotalPoints; i++)
            {
                if (pointset.phi[i] >= 0.5)
                {
                    Vx[i] = 0;
                    Vy[i] = 0;
                    Vz[i] = 0;
                }
            }
        }

        phi_old = pointset.phi;
        pointset.phi = computeAllenCahnPhaseTransformation2D(pointset.dim, pointset.phi, U, materials.D, inter_param, settings.dt, pointset.points, pointset.TotalPoints, kernels);
        for (int i = 1; i <= pointset.TotalPoints; i++)
        {Q[i] = inter_param.thermal.fusion_latent_heat*(pointset.phi[i]-phi_old[i])/settings.dt;}
        tie(U,RHS) = computeTransport(pointset.dim, U, Vx, Vy, Vz, Q, materials.D, settings.dt, t, RHS, pointset.points, pointset.phi, pointset.TotalPoints, kernels, BC.U);

        if (t/settings.prnt_freq == c)
        {
            c++;
            cout << "Computation of time step " << t << " of " << settings.nt << " complete" << endl;
            if (pointset.TotalSurfaces == 0 && pointset.TotalVolumes == 0)
            {
                pointset.printTXT(pointset.phi, "phi", t);
                if (settings.simulation_type == "thermal")
                {pointset.printTXT(U, "U", t);}
                if (settings.simulation_type == "hele_shaw")
                {pointset.printTXT(U, "P", t);}
                if (fluid_solver == true)
                {
                    pointset.printTXT(V, "V", t);
                    pointset.printVectorsVTK(V,Vx,Vy,Vz, "vectors", t);
                    pointset.printTXT(P, "P", t);
                    pointset.printTXT(ro, "ro", t);
                }
            }
            else
            {
                pointset.printMeshVTK(pointset.phi, "phi", t);
                if (settings.simulation_type == "thermal")
                {pointset.printMeshVTK(U, "U", t);}
                if (settings.simulation_type == "hele_shaw")
                {pointset.printMeshVTK(U, "P", t);}
                if (fluid_solver == true)
                {
                    pointset.printMeshVTK(V, "V", t);
                    pointset.printVectorsVTK(V, Vx, Vy, Vz, "vectors", t);
                    pointset.printMeshVTK(P, "P", t);
                    pointset.printMeshVTK(ro, "ro", t);
                }
            }
        }
    }
    pointset.printTXT(pointset.phi, "phi");
    if (settings.simulation_type == "thermal")
    {pointset.printTXT(U, "U");}
    if (settings.simulation_type == "hele_shaw")
    {pointset.printTXT(U, "P");}
    if (fluid_solver == true)
    {
        pointset.printTXT(V, "V");
        pointset.printTXT(P, "P");
        pointset.printTXT(ro, "ro");
    }
    cout << "SOLVER COMPLETE" << endl;
}
void solveAllenCahnMultipleSeeds(POINTSET pointset, MATERIALS materials, INTERFACE_PRMTRS inter_param, BOUNDARY_CONDITIONS bc, SOLVER_SETTINGS settings, vector<vector<double> > phi0, vector<double> U, vector<double> Vx, vector<double> Vy, vector<double> Vz)
{
    cout << "SOLVER START..." << endl;
    cout << "Initializing variables..."<< endl;
    pltctrl(pointset.points, pointset.TotalPoints, settings.nt, settings.prnt_freq);
    if (pointset.TotalSurfaces>0)
    {printElements(pointset.surfaces, pointset.SurfaceNds, pointset.TotalSurfaces, "elements");}

    int c = 1;
    bool fluid_solver = false;
    double kernel_radius = settings.kernel.radius_ratio*pointset.AvgPntSpacing;
    vector<vector<double> > phi_seeds_old(pointset.TotalPoints+1, vector<double> (inter_param.TotalSeeds+1));
    vector<vector<double> > phi_seeds(pointset.TotalPoints+1, vector<double> (inter_param.TotalSeeds+1));
    vector<double> phi_tot(pointset.TotalPoints+1);
    vector<double> phi_tot_old(pointset.TotalPoints+1);
    vector<double> V0(pointset.TotalPoints+1);
    vector<double> Q(pointset.TotalPoints+1);
    vector<double> RHS(pointset.TotalPoints+1);

    for (int i = 1; i <= pointset.TotalPoints; i++)
    {
        phi_tot[i] = 0;
        phi_tot_old[i] = 0;
    }
    for (int i = 1; i <= pointset.TotalPoints; i++)
    {
        for(int s = 1; s <=inter_param.TotalSeeds; s++)
        {phi_seeds[i][s] = phi0[i][s];}
    }
    for (int i = 1; i <= pointset.TotalPoints; i++)
    {
        for(int s = 1;  s <=inter_param.TotalSeeds; s++)
        {
            if (phi_seeds[i][s] > 0)
            {phi_tot[i] = phi_seeds[i][s];}
        }
    }
    phi_tot_old = phi_tot;

    if (materials.Nu[1] > 0)
    {fluid_solver = true;}
    else
    {fluid_solver = false;}
    cout << "fluid solver = " << fluid_solver << endl;

	vector<double> V;
	vector<double> P;
	vector<double> ro;
    if (fluid_solver == true)
    {
        Vx = setVector(Vx, bc.Vx);
        Vy = setVector(Vy, bc.Vy);
        Vz = setVector(Vz, bc.Vz);
        ro.resize(pointset.TotalPoints+1);
        P.resize(pointset.TotalPoints+1);
        V.resize(pointset.TotalPoints+1);
        for (int i = 1; i <= pointset.TotalPoints; i++)
        {
            ro[i] = materials.ro[1];
            P[i] = 0;
            V[i] = sqrt(Vx[i]*Vx[i] + Vy[i]*Vy[i] + Vz[i]*Vz[i]);
        }
        pointset.printVectorsVTK(V,Vx,Vy,Vz, "VEL", 0);
        pointset.printMeshVTK(V, "V", 0);
        pointset.printMeshVTK(ro, "ro", 0);
    }
    else
    {
        Vx = V0;
        Vy = V0;
        Vz = V0;
    }

    cout << "Printing field variables at time 0..." << endl;
    pointset.printMeshVTK(U, "U", 0);
    pointset.printMeshVTK(phi_tot, "phi", 0);

    cout << "Computing kernel approximation functions..." << endl;
    vector<KERNEL> kernels = computeKernels(pointset.points, pointset.phi, pointset.TotalPoints, pointset.dim, settings.kernel, kernel_radius);

    vector<KERNEL> kernels2;
    if (fluid_solver == true)
    {
        settings.kernel.approximation_method = "SPH";
        kernels2 = computeKernels(pointset.points, pointset.phi, pointset.TotalPoints, pointset.dim, settings.kernel, kernel_radius);
    }

    cout << "Computing transient solution..."<< endl;
    for (int t =1; t <= settings.nt; t++)
    {
        if (fluid_solver == true)
        {
            cout << "Computing flow solution for time step "<< t << " of " << settings.nt << endl;
            tie(Vx,Vy,Vz,V,P,ro) = solveNavierStokes(P, Vx, Vy, Vz, V0, V0, V0, pointset.phi, pointset.points, pointset.TotalPoints, kernels2, bc.Vx, bc.Vy, bc.Vz, bc.P, materials.Nu, materials.ro, materials.total, settings, pointset.AvgPntSpacing, 0);
            for (int i = 1; i <= pointset.TotalPoints; i++)
            {
                if (phi_tot[i] >= 0.5)
                {
                    Vx[i] = 0;
                    Vy[i] = 0;
                    Vz[i] = 0;
                }
            }
        }

        for(int s = 1; s <=inter_param.TotalSeeds; s++)
        {
            for (int i = 1; i <= pointset.TotalPoints; i++)
            {phi_seeds_old[i][s] = phi_seeds[i][s];}
        }
        for(int s = 1; s <=inter_param.TotalSeeds; s++)
        {
            for (int i = 1; i <= pointset.TotalPoints; i++)
            {pointset.phi[i] = phi_seeds[i][s];}
            inter_param.surface_tension.theta = inter_param.seed_theta[s];
            pointset.phi = computeAllenCahnPhaseTransformation2D(pointset.dim, pointset.phi, U, materials.D, inter_param, settings.dt, pointset.points, pointset.TotalPoints, kernels);
            for (int i = 1; i <= pointset.TotalPoints; i++)
            {phi_seeds[i][s] = pointset.phi[i];}
        }
        for (int i = 1; i <= pointset.TotalPoints; i++)
        {
            for(int s = 1;  s <=inter_param.TotalSeeds; s++)
            {
                if (phi_seeds[i][s] > 0)
                {phi_tot[i] = phi_seeds[i][s];}
                if (phi_seeds_old[i][s] > 0)
                {phi_tot_old[i] = phi_seeds_old[i][s];}
            }
        }

        for (int i = 1; i <= pointset.TotalPoints; i++)
        {Q[i] = inter_param.thermal.fusion_latent_heat*(phi_tot[i]-phi_tot_old[i])/settings.dt;}
        tie(U, RHS) = computeTransport(pointset.dim, U, Vx, Vy, Vz, Q, materials.D, settings.dt, t, RHS, pointset.points, pointset.phi, pointset.TotalPoints, kernels, bc.U);

        if (t/settings.prnt_freq == c)
        {
            c++;
            cout << "Computation of Time Step " << t << " of " << settings.nt << " complete..." << endl;
            if (pointset.TotalSurfaces == 0 && pointset.TotalVolumes == 0)
            {
                pointset.printTXT(phi_tot, "phi", t);
                pointset.printTXT(U, "U", t);
                if (fluid_solver == true)
                {
                    pointset.printTXT(V, "V", t);
                    pointset.printVectorsVTK(V,Vx,Vy,Vz,"vectors", t);
                    pointset.printTXT(P, "P", t);
                    pointset.printTXT(ro, "ro", t);
                }
            }
            else
            {
                pointset.printMeshVTK(phi_tot, "phi", t);
                pointset.printMeshVTK(U, "U", t);
                if (fluid_solver == true)
                {
                    pointset.printMeshVTK(V, "V", t);
                    pointset.printVectorsVTK(V,Vx,Vy,Vz, "vectors", t);
                    pointset.printMeshVTK(P, "P", t);
                    pointset.printMeshVTK(ro, "ro", t);
                }
            }
        }
    }
    pointset.printTXT(phi_tot, "phi");
    pointset.printTXT(U, "U");
    if (fluid_solver == true)
    {
        pointset.printTXT(V, "V");
        pointset.printTXT(P, "P");
        pointset.printTXT(ro, "ro");
    }
    cout << "SOLVER COMPLETE" << endl;
}
void solveCahnHilliard(POINTSET pointset, INTERFACE_PRMTRS inter_param, SOLVER_SETTINGS settings, vector<double> Vn, vector<double> Vx, vector<double> Vy, vector<double> Vz)
{
    cout << "SOLVER START..." << endl;
    cout << "Initializing variables..."<< endl;
    pltctrl(pointset.points, pointset.TotalPoints, settings.nt, settings.prnt_freq);
    if (pointset.TotalSurfaces>0)
    {printElements(pointset.surfaces, pointset.SurfaceNds, pointset.TotalSurfaces, "elements");}

    int m = 1;
    double kernel_radius = settings.kernel.radius_ratio*pointset.AvgPntSpacing;

    cout << "Computing kernel approximation functions..." << endl;
    vector<KERNEL> kernels = computeKernels(pointset.points, pointset.phi, pointset.TotalPoints, pointset.dim, settings.kernel, kernel_radius);

    cout << "Computing transient solution..."<< endl;
    for (int t =1; t <= settings.nt; t++)
    {
        pointset.phi = computeCahnHilliard(inter_param.surface_tension.d0, inter_param.PFM.mobility, kernels, pointset.phi, pointset.dim, pointset.points, pointset.TotalPoints, Vx, Vy, Vz, Vn, 1, settings.dt);
        if (t/settings.prnt_freq == m)
        {
            m++;
            cout << "Computation of time step " << t << " of " << settings.nt << " complete" << endl;
            pointset.printMeshVTK(pointset.phi, "phi", t);
            pointset.printTXT(pointset.phi, "phi", t);
        }
    }
    pointset.printTXT(pointset.phi, "phi");
    cout << "SOLVER COMPLETE" << endl;
}


