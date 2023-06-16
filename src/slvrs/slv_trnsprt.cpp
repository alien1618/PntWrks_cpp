#include "../pntwrks.h"

void solveTransport(POINTSET pointset, MATERIALS materials, BOUNDARY_CONDITION BC, SOLVER_SETTINGS settings, vector<double> U, vector<double> Vx, vector<double> Vy, vector<double> Vz, vector<double> Q)
{
    cout << "SOLVER START..." << endl;
    cout << "Initializing variables..." << endl;
    pltctrl(pointset.points, pointset.TotalPoints, settings.nt, settings.prnt_freq);
    if (pointset.TotalSurfaces>0)
    {printElements(pointset.surfaces, pointset.SurfaceNds, pointset.TotalSurfaces, "elements");}

    int s = 1;
    double kernel_radius = settings.kernel.radius_ratio*pointset.AvgPntSpacing;

    U = setVector(U, BC);
    vector<double> RHS = setVector(pointset.TotalPoints, 0);
    
    cout << "Printing field variables at time 0..." << endl;
    pointset.printTXT(U, "U", 0);
    pointset.printTXT(pointset.phi,"phi", 0);
    pointset.printMeshVTK(U, "U", 0);
    pointset.printMeshVTK(pointset.phi,"phi", 0);

    cout << "Computing kernel approximation functions..." << endl;
    vector<KERNEL> kernels = computeKernels(pointset.points, pointset.phi, pointset.TotalPoints, pointset.dim, settings.kernel, kernel_radius);

    cout << "Computing transient solution..."<< endl;
    for (int t = 1; t <= settings.nt; t++)
    {
        switch (settings.kernel.upwind)
        {
            case false:
            {tie(U, RHS) = computeTransport(pointset.dim, U, Vx, Vy, Vz, Q, materials.D, settings.dt, t, RHS, pointset.points, pointset.phi, pointset.TotalPoints, kernels, BC);}break;
            case true:
            {U = computeTransportUpwind(pointset.dim, U, Vx,Vy, Vz, Q, pointset.phi, materials.D, materials.total, settings.AV_factor, settings.dt, pointset.points, pointset.TotalPoints, kernels, BC, settings.kernel.upwind_ratio);}break;
		}
			
        if (t/settings.prnt_freq == s)
        {
            s++;
            cout << "Computation of time step " << t << " of " << settings.nt << " complete..." << endl;
            pointset.printMeshVTK(U, "U", t);
            pointset.printTXT(U, "U", t);
        }
    }
    cout << "SOLVER COMPLETE" << endl;
}
void solveTransportSteadyState(POINTSET pointset, MATERIALS materials, BOUNDARY_CONDITION BC, SOLVER_SETTINGS settings, vector<double> U, vector<double> Vx, vector<double> Vy, vector<double> Vz, vector<double> Q)
{
    cout << "SOLVER START..." << endl;
    cout << "Initializing variables..." << endl;
    pltctrl(pointset.points, pointset.TotalPoints, settings.nt, settings.prnt_freq);
    if (pointset.TotalSurfaces>0)
    {printElements(pointset.surfaces, pointset.SurfaceNds, pointset.TotalSurfaces, "elements");}
    
    double tol = 1e-3;
    int s = 1;
    double kernel_radius = settings.kernel.radius_ratio*pointset.AvgPntSpacing;
	double U_sum_old = 0;
    double U_sum = 0;
    U = setVector(U, BC);   
    for (int i = 1; i <= pointset.TotalPoints; i++)
    {U_sum = U_sum + U[i];}
    
    cout << "Printing field variables at time 0..." << endl;
    pointset.printMeshVTK(U, "U", 0);
    pointset.printMeshVTK(pointset.phi,"phi", 0);
    pointset.printTXT(U, "U", 0);
    pointset.printTXT(pointset.phi,"phi", 0);
    
    cout << "Computing kernel approximation functions..." << endl;
    vector<KERNEL> kernels = computeKernels(pointset.points, pointset.phi, pointset.TotalPoints, pointset.dim, settings.kernel, kernel_radius);

    cout << "Computing transient solution..."<< endl;
    double diff;
    diff = 0;
    U_sum_old = U_sum;
    for (int t = 1; t <= settings.nt; t++)
    {
        U = computeTransportGaussSeidel(pointset.dim, U, Vx, Vy, Vz, Q, materials.D, settings.dt, pointset.points, pointset.phi, pointset.TotalPoints, kernels, BC);
        U_sum_old = U_sum;
        U_sum = 0;
        for (int i = 1; i <= pointset.TotalPoints; i++)
        {U_sum = U_sum + U[i];}
        diff = abs(U_sum - U_sum_old);
        if (t > 10)
        {
            if (diff <= tol)
            {break;}
        }
        if (t/100 == s)
        {
            s++;
            cout << "Computation of time step " << t << " of " << settings.nt << " complete..." << endl;
            cout << "RES = " << diff << endl;
            pointset.printMeshVTK(U, "U", t);
            pointset.printTXT(U, "U", t);
        }
    }
    cout << "SOLVER COMPLETE" << endl;
}

void solveReactionDiffusion(POINTSET pointset, MATERIALS materials, BOUNDARY_CONDITION BC, SOLVER_SETTINGS settings,  double f, double k, vector<double> U_A, vector<double> U_B)
{
    cout << "SOLVER START..." << endl;
    cout << "Initializing variables..." << endl;
    pltctrl(pointset.points, pointset.TotalPoints, settings.nt, settings.prnt_freq);
    if (pointset.TotalSurfaces>0)
    {printElements(pointset.surfaces, pointset.SurfaceNds, pointset.TotalSurfaces, "elements");}

    int s = 1;
    double kernel_radius = settings.kernel.radius_ratio*pointset.AvgPntSpacing;
    MATERIALS materials_A;
    materials_A = materials;
    materials_A.D[1] = materials.D[1];
    materials_A.D[2] = materials.D[1];
    MATERIALS materials_B;
    materials_B = materials;
    materials_B.D[1] = materials.D[2];
    materials_B.D[2] = materials.D[2];
    vector<double> V0(pointset.TotalPoints+1);
    vector<double> Q_A(pointset.TotalPoints+1);
    vector<double> Q_B(pointset.TotalPoints+1);
    vector<double> RHS_A(pointset.TotalPoints,0);
    vector<double> RHS_B(pointset.TotalPoints,0);
    
#pragma omp parallel for
    for (int i = 1; i <= pointset.TotalPoints; i++)
    {
        if (pointset.phi[i] <= 0)
        {U_B[i] = U_A[i];}
    }

    cout << "Printing field variables at time 0..." << endl;
    if (pointset.TotalSurfaces > 0 || pointset.TotalVolumes > 0)
    {
        pointset.printMeshVTK(U_A, "A", 0);
        pointset.printMeshVTK(U_B, "B", 0);
    }
    else
    {
        pointset.printTXT(U_A, "A", 0);
        pointset.printTXT(U_B, "B", 0);
    }

    cout << "Computing kernel approximation functions..." << endl;
    vector<KERNEL> kernels = computeKernels(pointset.points, pointset.phi, pointset.TotalPoints, pointset.dim, settings.kernel, kernel_radius);

    cout << "Computing transient solution..."<< endl;
    for (int t = 1; t <= settings.nt; t++)
    {
#pragma omp parallel for
        for(int i = 1; i <= pointset.TotalPoints; i++)
        {
            Q_A[i] = -U_A[i]*U_B[i]*U_B[i] + f*(1-U_A[i]);
            Q_B[i] =  U_A[i]*U_B[i]*U_B[i] - (k+f)*(U_B[i]);
        }

        tie(U_A, RHS_A) = computeTransport(pointset.dim, U_A, V0, V0, V0, Q_A, materials_A.D, settings.dt, t, RHS_A, pointset.points, pointset.phi, pointset.TotalPoints, kernels, BC);
        tie(U_B, RHS_B) = computeTransport(pointset.dim, U_B, V0, V0, V0, Q_B, materials_B.D, settings.dt, t, RHS_B, pointset.points, pointset.phi, pointset.TotalPoints, kernels, BC);

#pragma omp parallel for
        for (int i = 1; i <= pointset.TotalPoints; i++)
        {
            if (U_A[i] < 0)
            {U_A[i] = 0;}
            if (U_A[i] > 1)
            {U_A[i] = 1;}
            if (U_B[i] < 0)
            {U_B[i] = 0;}
            if (U_B[i] > 1)
            {U_B[i] = 1;}
        }
        if (t/settings.prnt_freq == s)
        {
            s++;
            cout << "Computation of time step " << t << " of " << settings.nt << " complete..." << endl;
            if (pointset.TotalSurfaces > 0 || pointset.TotalVolumes > 0)
            {
                pointset.printMeshVTK(U_A, "A", t);
                pointset.printMeshVTK(U_B, "B", t);
            }
            else
            {
                pointset.printTXT(U_A, "A", t);
                pointset.printTXT(U_B, "B", t);
            }
        }
    }
    cout << "SOLVER COMPLETE" << endl;
}
void solveTransportGFD(POINTSET pointset, MATERIALS materials, BOUNDARY_CONDITION BC, SOLVER_SETTINGS settings, vector<double> U, vector<double> Vx, vector<double> Vy, vector<double> Q, vector<double> phi)
{
    cout << "SOLVER START..." << endl;
    cout << "Initializing variables..." << endl;
    pltctrl(pointset.points, pointset.TotalPoints, settings.nt, settings.prnt_freq);
    if (pointset.TotalSurfaces>0)
    {printElements(pointset.surfaces, pointset.SurfaceNds, pointset.TotalSurfaces, "elements");}
    
    int s = 1;
    double kernel_radius = settings.kernel.radius_ratio*pointset.AvgPntSpacing;

    cout << "Printing field variables at time 0..." << endl;
    pointset.printTXT(U,"U", 0);
    if (pointset.TotalSurfaces > 0)
    {pointset.printMeshVTK(U, "U", 0);}

    cout << "Collecting kernel points..." << endl;
    vector<KERNEL> kernels(pointset.TotalPoints+1);
#pragma omp parallel for
    for (int i = 1; i <= pointset.TotalPoints; i++)
    {
		KERNEL kernel;
		kernel.collectNbrsBruteForce(pointset.points[i], pointset.points, pointset.TotalPoints, kernel_radius);
		kernels[i].nbrs = kernel.nbrs;
		kernels[i].TotalNbrs = kernel.TotalNbrs;
	}

    cout << "Computing transient solution..."<< endl;
    for (int t = 1; t <= settings.nt; t++)
    {
        U = computeTransportGFD(kernel_radius, U, Vx, Vy, Q, phi, materials.D, materials.total, settings.dt, pointset.points, pointset.TotalPoints, kernels, BC);
        if (t/settings.prnt_freq == s)
        {
            s++;
            cout << "Computation of time step " << t << " of " << settings.nt << " complete..." << endl;
            pointset.printTXT(U,"U", t);
            if (pointset.TotalSurfaces > 0)
            {pointset.printMeshVTK(U,"U", t);}
        }
    }
    cout << "SOLVER COMPLETE" << endl;
}
void solveBurger(POINTSET pointset, MATERIALS materials, BOUNDARY_CONDITION BC, SOLVER_SETTINGS settings, vector<double> U)
{
    cout << "SOLVER START..." << endl;
    cout << "Initializing variables..." << endl;
    pltctrl(pointset.points, pointset.TotalPoints, settings.nt, settings.prnt_freq);
    vector<double> zero(pointset.TotalPoints+1);
    if (pointset.TotalSurfaces>0)
    {printElements(pointset.surfaces, pointset.SurfaceNds, pointset.TotalSurfaces, "elements");}
    cout << "pointset.dim = " << pointset.dim << endl;
    cout << "pointset.AvgPntSpacing " << pointset.AvgPntSpacing << endl;
    int c = 1;
    double kernel_radius = settings.kernel.radius_ratio*pointset.AvgPntSpacing;
    vector<double> phi = setVector(pointset.TotalPoints, 1);
    vector<double> RHS = setVector(pointset.TotalPoints, 0);

    cout << "Printing field variables at time 0..." << endl;
    pointset.printTXT(U,"U", 0);
    if (pointset.TotalSurfaces > 0)
    {pointset.printMeshVTK(U,"U", 0);}

    cout << "Computing kernel approximation functions..." << endl;
    vector<KERNEL> kernels = computeKernels(pointset.points, pointset.phi, pointset.TotalPoints, pointset.dim, settings.kernel, kernel_radius);

    cout << "Computing transient solution..." << endl;
    string solution = "BE";
    if (solution == "BE")//backwards Euler
    {
        for (int t = 1; t <= settings.nt; t++)
        {
			switch(settings.kernel.upwind)
			{
				case false:
				{tie(U, RHS) = computeTransport(pointset.dim, U, U, zero, zero, zero, materials.D, settings.dt, t, RHS, pointset.points, pointset.phi, pointset.TotalPoints, kernels, BC);}break;
				case true:
				{U = computeTransportUpwind(pointset.dim, U, U, zero, zero, zero, phi, materials.D, materials.total, settings.AV_factor, settings.dt, pointset.points, pointset.TotalPoints, kernels, BC, settings.kernel.upwind_ratio);}break;
			}
            if (t/settings.prnt_freq == c)
            {
                c++;
                cout << "Computation of time step " << t << " of " << settings.nt <<  " complete..." << endl;
				pointset.printTXT(U,"U", t);
				if (pointset.TotalSurfaces > 0)
				{pointset.printMeshVTK(U,"U", t);}
            }
        }
    }
    else if (solution == "implicit")//iterative implicit
    {
        vector<double> U_iter = U;
        for (int t = 1; t <= settings.nt; t++)
        {
            cout << "Computing time step " << t << " of " << settings.nt << endl;
            for (int iter = 1; iter <= 5; iter++)
            {
                vector<double> U_iter_old = U_iter;
                U_iter = computeTransportIterative(pointset.dim, U, U_iter, U_iter, zero, zero, zero, phi, materials.D, materials.total, settings.AV_factor, settings.dt, pointset.points, pointset.TotalPoints, kernels, BC);
                double res = 0;
                for (int i = 1; i <= pointset.TotalPoints; i++)
                {res = res + abs(U_iter_old[i]-U_iter[i]);}
                cout << "iteration " << iter << " res = " << res << endl;
            }
            U = U_iter;
            if (t/settings.prnt_freq == c)
            {
                c++;
                cout << "Computation of time step " << t << " of " << settings.nt <<  " complete..." << endl;
				pointset.printTXT(U,"U", t);
				if (pointset.TotalSurfaces > 0)
				{pointset.printMeshVTK(U,"U", t);}

            }
        }
    }
    cout << "SOLVER COMPLETE" << endl;
}

