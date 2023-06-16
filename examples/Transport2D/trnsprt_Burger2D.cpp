#include "src/pntwrks.h"

int main()
{
    cout << "-----------------------------------------------------------------" << endl;
    cout << "2D burger equation" << endl;
    cout << "-----------------------------------------------------------------" << endl;

    cout << "Constructing pointset..."<< endl;
    POINT p(-1,-1);
    DOMAIN d(2,2);
    RESOLUTION r(100,100);
    POINTSET ps(p, d, r, "Q4");

    cout << "Assigning material properties..."<< endl;
    MATERIALS materials(1);
    materials.D[1] = 0.001;

    cout << "Assigning initial conditions..."<< endl;
    POINT cntr(0,0);
    ps.addSphere(cntr, 0.50);
    vector<double> U = setVector(ps.TotalPoints,0);
    for (int i = 1; i <= ps.TotalPoints; i++)
    {
        if (ps.phi[i] <= 0)
        {U[i] = 10;}
    }
    
    cout << "Assigning boundary conditions..."<< endl;
    BOUNDARY_CONDITION bc;

    cout << "Assigning solver settings..."<< endl;
    SOLVER_SETTINGS settings;
    settings.kernel.recompute_approximants = true;
    settings.kernel.support_method = "BF";
    settings.kernel.radius_ratio = 3;
    settings.kernel.approximation_method = "MLS";
    settings.kernel.approximation_order = 2;
    settings.kernel.RBF= "GE";
    settings.kernel.RBF_alpha= 0.1;
    settings.kernel.WLS = "LS";
    settings.kernel.MLS = "S4";
    settings.kernel.SPH = "WC4";
    settings.kernel.upwind = true;
    settings.kernel.upwind_ratio = 1;
    settings.dt = 0.00001;//use for rbf, wls, mls, nsph
    settings.nt = 20000;
    settings.AV_factor = 10*(settings.dt/ps.AvgPntSpacing);
    settings.prnt_freq = 100;

    cout << "Running solver..."<< endl;
    solveBurger(ps, materials, bc, settings, U);

    return 0;
}
