#include "src/pntwrks.h"

int main ()
{ 
    cout << "-----------------------------------------------------------------" << endl;
    cout << "3D splitting level set evolution" << endl;
    cout << "-----------------------------------------------------------------" << endl; 
    cout << "Constructing pointset..."<< endl;
    POINT p(-1,-1,-1);
    DOMAIN s(2,2,2);
    RESOLUTION divs(20,20,20);
    POINTSET pointset(p, s, divs, "H8");

    cout << "Assigning phases..."<< endl;
    POINT cntr(0.3, 0,0);
    pointset.addSphere(cntr, 0.35);
    POINT cntr1(-0.3, 0,0);
    pointset.addSphere(cntr1, 0.35);
    pointset.computeVOF(pointset.AvgPntSpacing);
    pointset.printMeshVTK(pointset.phi, "phi", 0);

    cout << "Assigning interface parameters..."<< endl;
    INTERFACE_PRMTRS interf_prmtrs;
    interf_prmtrs.solver = "VOF";

    cout << "Initializing field variables..."<< endl;
    vector<double> Vx = setVector(pointset.TotalPoints,0);
    vector<double> Vy = setVector(pointset.TotalPoints,0);
    vector<double> Vz = setVector(pointset.TotalPoints,0);
    vector<double> Vn_ext = assignConstantExtendedVelocity(pointset.TotalPoints, 0.025);

    cout << "Assigning solver parameters..."<< endl;
    SOLVER_SETTINGS settings;
    settings.kernel.recompute_approximants = true;
    settings.kernel.support_method = "BF";
    settings.kernel.radius_ratio = 1.5;
    settings.kernel.approximation_method = "SPH";
    settings.kernel.approximation_order = 1;
    settings.kernel.RBF = "GE";
    settings.kernel.RBF_alpha = 0.1;
    settings.kernel.WLS = "LS";
    settings.kernel.MLS = "S4";
    settings.kernel.SPH = "WC4";
    settings.dt = 0.01;
    settings.nt = 3000;
    settings.prnt_freq = 100;

    cout << "Running solver..."<< endl;
    solveInterfaceEvolution(pointset, interf_prmtrs, settings, Vn_ext, Vx, Vy, Vz);

    return 0;
}
