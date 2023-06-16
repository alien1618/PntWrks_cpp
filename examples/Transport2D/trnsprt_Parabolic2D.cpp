#include "src/pntwrks.h"

int main()
{
    cout << "-----------------------------------------------------------------" << endl;
    cout << "2D parabolic equation" << endl;
    cout << "-----------------------------------------------------------------" << endl;
    
    cout << "Constructing pointset..."<< endl;
    POINTSET pointset("geometry/PlateWithHole2D/mesh_400.dat");

    cout << "Assigning material properties..."<< endl;
    MATERIALS materials(1);
    materials.D[1] = 1;

    cout << "Initializing field variables..."<< endl;
	vector<double> U = setVector(pointset.TotalPoints,0);
	vector<double> Vx = setVector(pointset.TotalPoints,0);
	vector<double> Vy = setVector(pointset.TotalPoints,0);
	vector<double> Vz = setVector(pointset.TotalPoints,0);
	vector<double> Q = setVector(pointset.TotalPoints,0);
	
    cout << "Assigning boundary conditions..."<< endl;
    BOUNDARY_CONDITION bc;
    bc.assignDBC(pointset, "geometry/PlateWithHole2D/bottom_400.dat", 1);
    bc.assignDBC(pointset, "geometry/PlateWithHole2D/top_400.dat", 0);

    cout << "Assigning solver settings..."<< endl;
    SOLVER_SETTINGS settings;
    settings.kernel.recompute_approximants = true;
    settings.kernel.support_method = "BF";
    settings.kernel.radius_ratio = 2;
    settings.kernel.approximation_method = "MLS";
    settings.kernel.approximation_order = 1;
    settings.kernel.RBF = "GE";
    settings.kernel.RBF_alpha = 0.01;
    settings.kernel.WLS = "LS";
    settings.kernel.MLS = "S4";
    settings.kernel.SPH = "GE";
    settings.dt = 0.0001;
    settings.nt = 10000;
    settings.prnt_freq = 100;

    cout << "Running solver..."<< endl;
    solveTransport(pointset, materials, bc, settings, U, Vx, Vy, Vz, Q);

    return 0;
}
