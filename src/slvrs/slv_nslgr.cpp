#include "../pntwrks.h"

void solveNavierStokesLagrangian(POINTSET pointset, POINTSET bndry, POINTSET bg_grid, MATERIALS materials,  INTERFACE_PRMTRS inter_param, BOUNDARY_CONDITIONS bc, SOLVER_SETTINGS settings, vector<double> phi0)
{
    cout << "SOLVER START..." << endl;
    cout << "Initializing variables..."<< endl;
    pltctrl(pointset.points, pointset.TotalPoints, settings.nt, settings.prnt_freq);
    if (pointset.TotalSurfaces>0)
    {printElements(pointset.surfaces, pointset.SurfaceNds, pointset.TotalSurfaces, "elements");}

    int cntr = 1;
    double factor = 2;
    double eps = 0.000001;
    POINTSET particles;
    particles.TotalPoints = pointset.TotalPoints+bndry.TotalPoints;
    particles.points.resize(particles.TotalPoints+1);
    
    vector<double> Vx(particles.TotalPoints+1);
    vector<double> Vy(particles.TotalPoints+1);
    vector<double> Vz(particles.TotalPoints+1);
    vector<double> V(particles.TotalPoints+1);
    vector<double> P(particles.TotalPoints+1);
    vector<double> ro(particles.TotalPoints+1);
    vector<double> Fx(particles.TotalPoints+1);
    vector<double> Fy(particles.TotalPoints+1);
    vector<double> Fz(particles.TotalPoints+1);
    vector<double> phi(particles.TotalPoints+1);
	Vx = setVector(Vx, bc.Vx);
    Vy = setVector(Vy, bc.Vy);
    Vz = setVector(Vz, bc.Vz);

    vector<vector<double> > VEL(particles.TotalPoints+1, vector<double> (3));
    vector<int> bndry_pnt_num(bndry.TotalPoints+1);
    for (int i = 1; i <= pointset.TotalPoints; i++)
    {
        particles.points[i].x = pointset.points[i].x;
        particles.points[i].y = pointset.points[i].y;
        particles.points[i].z = pointset.points[i].z;
        phi[i] = phi0[i];
    }
    int s = 1;
    for (int i = pointset.TotalPoints+1; i <= pointset.TotalPoints+bndry.TotalPoints; i++)
    {
        particles.points[i].x = bndry.points[s].x;
        particles.points[i].y = bndry.points[s].y;
        particles.points[i].z = bndry.points[s].z;
        bndry_pnt_num[s] = i;
        phi[i] = 1;
        s++;
    }
    pltctrl(particles.points, particles.TotalPoints, settings.nt, settings.prnt_freq);

    double kernel_radius = settings.kernel.radius_ratio*pointset.AvgPntSpacing;
    for (int i = 1; i <= particles.TotalPoints; i++)
    {
        ro[i] = phi[i]*materials.ro[2] + (1-phi[i])*materials.ro[1];
        P[i] = 0;
        V[i] = sqrt(Vx[i]*Vx[i] + Vy[i]*Vy[i]+ Vz[i]*Vz[i]);
    }
    for (int i = 1; i <= particles.TotalPoints; i++)
    {
        Fx[i] = settings.gravity.x;
        Fy[i] = settings.gravity.y;
        Fz[i] = settings.gravity.z;
    }
    vector<double> ro0(particles.TotalPoints+1);
    particles.printTXT(V, "V", 0);
    particles.printTXT(ro, "ro", 0);
    particles.printTXT(P, "P", 0);
    particles.printTXT(phi, "phi", 0);
    particles.printVTK(V, "V", 0);
    particles.printVTK(ro,"ro", 0);
    particles.printVTK(P, "P", 0);
    particles.printVTK(phi, "phi", 0);

    int approx = 1;
    if (settings.kernel.approximation_method == "MLS")
    {approx = 1;}
    else if (settings.kernel.approximation_method == "RBF")
    {approx = 2;}
    else if (settings.kernel.approximation_method == "WLS")
    {approx = 3;}
    else if (settings.kernel.approximation_method == "SPH")
    {approx = 4;}
    else
    {cout << "ERROR: approximation method is invalid." << endl; exit(0);}

    //===========================================================
    // Main loop for transient time evolution
    //===========================================================
    for (int t = 1; t <= settings.nt; t++)
    {
        //===========================================================
        // computing kernel approximation functions
        //===========================================================
        int order = 1;
        vector<KERNEL> kernels(particles.TotalPoints+1);
#pragma omp parallel for
        for (int i = 1; i <= particles.TotalPoints; i++)
        {
            KERNEL kernel;
            kernel.collectNbrsBruteForce(particles.points[i], particles.points, particles.TotalPoints, kernel_radius);

            if (kernel.TotalNbrs <= 8)
            {kernel.computeSPH(particles.points[i], particles.points, kernel.nbrs, kernel.TotalNbrs, 1, kernel_radius, pointset.dim);}
            else
            {
                switch(approx)
                {
                    case 1:
                    {kernel.computeMLS(particles.points[i], particles.points, kernel.nbrs, kernel.TotalNbrs, 1, order, 3*kernel_radius);}break;
                    case 2:
                    {kernel.computeRBF(particles.points[i], particles.points, kernel.nbrs, kernel.TotalNbrs, 1, order, kernel_radius, 0.01);}break;
                    case 3:
                    {kernel.computeWLS(particles.points[i], particles.points, kernel.nbrs, kernel.TotalNbrs, 1, order, 3*kernel_radius);}break;
                    case 4:
                    {kernel.computeSPH(particles.points[i], particles.points,kernel.nbrs, kernel.TotalNbrs, 1, kernel_radius, pointset.dim);}break;
                }
            }
            kernels[i].TotalNbrs = kernel.TotalNbrs;
            kernels[i].nbrs = kernel.nbrs;
            kernels[i].N = kernel.N;
            kernels[i].dNdx = kernel.dNdx;
            kernels[i].dNdy = kernel.dNdy;
            kernels[i].dNdz = kernel.dNdy;
            kernels[i].dNdxx.resize(kernel.TotalNbrs+1);
            kernels[i].dNdyy.resize(kernel.TotalNbrs+1);
            kernels[i].dNdzz.resize(kernel.TotalNbrs+1);
            kernels[i].nabla2.resize(kernel.TotalNbrs+1);
            for (int j = 1; j <= kernels[i].TotalNbrs; j++)
            {
                int nbr = kernels[i].nbrs[j];
                double d = particles.points[i].distance(particles.points[nbr])+1e-6;
                double nijx = (particles.points[nbr].x-particles.points[i].x)/d;
                double nijy = (particles.points[nbr].y-particles.points[i].y)/d;
                double nijz = (particles.points[nbr].z-particles.points[i].z)/d;
                kernels[i].dNdxx[j] = factor*(nijx/(d+eps))*kernels[i].dNdx[j];
                kernels[i].dNdyy[j] = factor*(nijy/(d+eps))*kernels[i].dNdy[j];
                kernels[i].dNdzz[j] = factor*(nijz/(d+eps))*kernels[i].dNdz[j];
                kernels[i].nabla2[j] = kernels[i].dNdxx[j]+kernels[i].dNdyy[j]+kernels[i].dNdzz[j];
            }
        }
        if (t == 1)
        {
#pragma omp parallel for
            for (int i = 1; i<= particles.TotalPoints; i++)
            {
                ro0[i] = 0;
                for (int s = 1; s <= kernels[i].TotalNbrs; s++)
                {
                    int j = kernels[i].nbrs[s];
                    double d = particles.points[i].distance(particles.points[j]);
                    if (d <= kernel_radius)
                    {ro0[i] = ro0[i] + ((kernel_radius/(d+1e-8)) - 1);}
                }
            }
        }
        if (materials.total > 1 && inter_param.solver != "NONE")
        {
            //===========================================================
            // solving the Phase-Field equation
            //===========================================================
            vector<KERNEL> kernels2(pointset.TotalPoints+1);
#pragma omp parallel for
            for (int i = 1; i <= pointset.TotalPoints; i++)
            {
                KERNEL kernel;
                kernel.collectNbrsBruteForce(particles.points[i], particles.points, particles.TotalPoints, kernel_radius);
                kernel.computeSPH(particles.points[i], particles.points, kernel.nbrs, kernel.TotalNbrs, 4, kernel_radius, pointset.dim);
                kernels2[i].TotalNbrs = kernel.TotalNbrs;
                kernels2[i].nbrs = kernel.nbrs;
                kernels2[i].N = kernel.N;
                kernels2[i].dNdx = kernel.dNdx;
                kernels2[i].dNdy = kernel.dNdy;
                kernels2[i].dNdz = kernel.dNdz;
            }
            vector<double> V0(pointset.TotalPoints+1);
            phi = computeCahnHilliard(1*pointset.AvgPntSpacing, 0.001, kernels2, phi, 2, particles.points, pointset.TotalPoints, V0, V0, V0, V0, 3, settings.dt);
            //phi = computeVolumeOfFluid(kernels, phi, particles.TotalPoints, V0, V0, V0, 5, 0.001*settings.dt);
            //phi = computeSharpenInterface(inter_param.sharpen_interface, inter_param.sharpen_interface_beta, inter_param.sharpen_interface_eps, phi, pointset.TotalPoints, kernels, pointset.AvgPntSpacing, inter_param.sharpen_interface_dt, inter_param.sharpen_interface_nt);
            for (int i = 1; i <= pointset.TotalPoints; i++)
            {
                if (phi[i] <= 1e-5)
                {phi[i] = 0;}
                if (phi[i] >= 0.999)
                {phi[i] = 1;}
            }
            for (int i = pointset.TotalPoints+1;  i<= particles.TotalPoints; i++)
            {phi[i] = 0;}
        }
        //===========================================================
        // solving the Navier-Stokes equations
        //===========================================================
        tie(Vx,Vy,Vz,V,P,ro) = solveNavierStokes(P, Vx, Vy, Vz, Fx, Fy, Fz, phi, particles.points, particles.TotalPoints, kernels, bc.Vx, bc.Vy, bc.Vz, bc.P, materials.Nu, materials.ro, materials.total, settings, pointset.AvgPntSpacing, inter_param.surface_tension.surface_tension_co);
        //tie(Vx,Vy,Vz,V,P,ro) = solveNavierStokesUpwind(P, Vx, Vy, Vz, Fx, Fy, Fz, phi, particles.points, particles.TotalPoints, kernels, bc.Vx, bc.Vy, bc.Vz, bc.P, materials.Nu, materials.ro, materials.total, settings, pointset.AvgPntSpacing, inter_param.surface_tension.surface_tension_co, pointset.dim);
        //===========================================================
        // computing particle shifting velocity
        //===========================================================
        if (settings.particle_shifting == true)
        {
            VEL = computeParticleShiftingVelocity(settings.particle_shifting_factor, settings.particle_shifting_surf, kernel_radius, settings.dt, particles.points, particles.TotalPoints, kernels);
            //VEL = computeParticleShiftingVelocity_v2(settings.cfd.particle_shifting_factor, settings.cfd.particle_shifting_surf, kernel_radius, settings.dt, particles.x, particles.y, particles.TotalPoints, kernels);
            //VEL = computeParticleShiftingPFM(settings.cfd.particle_shifting_factor, settings.cfd.particle_shifting_surf, kernel_radius, settings.dt, particles.x, particles.y, particles.TotalPoints, kernels);
        }
        else
        {
#pragma omp parallel for
            for (int i = 1;  i<= pointset.TotalPoints; i++)
            {
                VEL[i][1] = 0;
                VEL[i][2] = 0;
                //VEL[i][3] = 0;
            }
        }
        //===========================================================
        // advecting particles in space
        //===========================================================
#pragma omp parallel for
        for (int i = 1;  i<= pointset.TotalPoints; i++)
        {
            particles.points[i].x = particles.points[i].x + settings.dt*(Vx[i]+VEL[i][1]);
            particles.points[i].y = particles.points[i].y + settings.dt*(Vy[i]+VEL[i][2]);
          //  particles.points[i].z = particles.points[i].z + settings.dt*(Vz[i]+VEL[i][3]);
        }
        //===========================================================
        // print data to file
        //===========================================================
        if (t/settings.prnt_freq == cntr)
        {
            cntr++;
            cout << "Processing time step " << t << " of " << settings.nt << " complete..." << endl;
            particles.printTXT(V, "V", t);
            particles.printTXT(ro, "ro", t);
            particles.printTXT(P, "P", t);
            particles.printTXT(phi, "phi", t);

            particles.printVTK(V, "V", t);
            particles.printVTK(ro, "ro", t);
            particles.printVTK(P, "P", t);
            particles.printVTK(phi, "phi", t);
        }
    }
    cout << "SOLVER COMPLETE" << endl;
    
}

