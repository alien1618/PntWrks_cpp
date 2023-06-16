#include "src/pntwrks.h"

int main()
{
    cout << "-----------------------------------------------------------------" << endl;
    cout << "3D advection-diffusion equation" << endl;
    cout << "-----------------------------------------------------------------" << endl;

    cout << "Constructing pointset..."<< endl;
    POINT p(-1,-1,0);
    DOMAIN d(2,2,0.1);
    RESOLUTION r(50,50,2);
    POINTSET ps(p, d, r, "T4");

    cout << "Assigning material properties..."<< endl;
    MATERIALS materials(1);
    materials.D[1] = 0;
    
    cout << "Assigning initial conditions..."<< endl;
    vector<double> U = setVector(ps.TotalPoints,0);
    vector<double> Vx = setVector(ps.TotalPoints,0);
    vector<double> Vy = setVector(ps.TotalPoints,0);
    vector<double> Vz = setVector(ps.TotalPoints,0);
    vector<double> Q = setVector(ps.TotalPoints,0);
    POINT cntr(0,-0.5,0);
    ps.addCylinderZ(cntr, 0.25);
    tie(Vx,Vy,Vz) = assignRotatingAdvectiveVelocity(ps.points, ps.TotalPoints, 10, 10, 0);
    for (int i = 1; i <= ps.TotalPoints; i++)
    {
        if (ps.phi[i] <= 0)
        {U[i] = 1;}
        else
        {U[i] = 0;}
        ps.phi[i] = 1;
    }
    
    cout << "Assigning boundary conditions..."<< endl;
    BOUNDARY_CONDITION bc;

    cout << "Assigning solver settings..."<< endl;
    SOLVER_SETTINGS settings;
    settings.kernel.recompute_approximants = true;
    settings.kernel.support_method = "BF";
    settings.kernel.radius_ratio = 3;
    settings.kernel.approximation_method = "SPH";
    settings.kernel.approximation_order = 1;
    settings.kernel.RBF = "GE";
    settings.kernel.RBF_alpha = 0.1;
    settings.kernel.WLS = "S4";
    settings.kernel.MLS = "S4";
    settings.kernel.SPH = "WC4";
    settings.kernel.upwind = true;
    settings.kernel.upwind_ratio = 1;
    settings.dt = 0.0001;
    settings.nt = 10000;
    settings.prnt_freq = 10;
    settings.AV_factor = (settings.dt/ps.AvgPntSpacing);

    cout << "Running solver..."<< endl;
    solveTransport(ps, materials, bc, settings, U, Vx, Vy,Vz, Q);

    return 0;
}
