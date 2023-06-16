#include "src/pntwrks.h"

int main()
{
    cout << "-----------------------------------------------------------------" << endl;
    cout << "3D poisson equation" << endl;
    cout << "-----------------------------------------------------------------" << endl;

    cout << "Constructing pointset..."<< endl;
    POINT p(0,0,0);
    DOMAIN size(1,1,0.1);
    RESOLUTION divs(20, 20,2);
    POINTSET pointset(p, size, divs, "T4");

    cout << "Assigning material properties..."<< endl;
    MATERIALS materials(1);
    materials.D[1] = 1;

    cout << "Assigning boundary conditions..."<< endl;
    BOUNDARY_CONDITION bc;
    bc.assignDBC(pointset.points, pointset.TotalPoints, "y", 0, 1);
    bc.assignDBC(pointset.points, pointset.TotalPoints, "y", 1, 1);
    bc.assignDBC(pointset.points, pointset.TotalPoints, "x", 0, 1);
    bc.assignDBC(pointset.points, pointset.TotalPoints, "x", 1, 1);

    cout << "Assigning solver settings..."<< endl;
    SOLVER_SETTINGS settings;
    settings.kernel.recompute_approximants = true;
    settings.kernel.support_method = "BF";
    settings.kernel.radius_ratio = 2.5;
    settings.kernel.approximation_method = "RBF";
    settings.kernel.approximation_order = 1;
    settings.kernel.RBF = "GE";
    settings.kernel.RBF_alpha = 0.1;
    settings.kernel.WLS = "S4";
    settings.kernel.MLS = "S4";
    settings.kernel.SPH = "WC4";
    settings.dt = 0.00001;
    settings.nt = 50000;
    settings.prnt_freq = 100;

    cout << "Initializing field variables..."<< endl;
    vector<double> U = setVector(pointset.TotalPoints,0);
	vector<double> Vx = setVector(pointset.TotalPoints,0);
	vector<double> Vy = setVector(pointset.TotalPoints,0);
	vector<double> Vz = setVector(pointset.TotalPoints,0);
	vector<double> Q = setVector(pointset.TotalPoints,0);
    tie(Vx,Vy,Vz) = assignSingleVortexVelocity(pointset.points, pointset.TotalPoints);

    vector<KERNEL> kernels(pointset.TotalPoints+1);
    kernels = computeKernels(pointset.points, pointset.phi, pointset.TotalPoints, pointset.dim, settings.kernel, settings.kernel.radius_ratio*pointset.AvgPntSpacing);

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
        Q[i] = Vx[i]*Vx[i] + Vy[i]*Vy[i] + Vz[i]*Vz[i];//(ro0/settings.dt)*(term);
    }
    
    cout << "Running solver..."<< endl;
    solveTransportSteadyState(pointset, materials, bc, settings, U, Vx, Vy, Vz, Q);
    
    return 0;
}
