#include "src/pntwrks.h"

int main()
{
    cout << "-----------------------------------------------------------------" << endl;
    cout << "3D parabolic equation" << endl;
    cout << "-----------------------------------------------------------------" << endl;

    cout << "Constructing pointset..."<< endl;
    //POINTSET pointset("geometry/PlateWithHole3D/mesh_5000.dat");
    POINT p(0,0);
    DOMAIN size(1,1,0.5);
    RESOLUTION divs(50,50,25);
    POINTSET pointset(p, size, divs, "H8");
    
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
    //bc.assignDBC(pointset, "geometry/PlateWithHole3D/bottom_5000.dat", 1);
    //bc.assignDBC(pointset, "geometry/PlateWithHole3D/top_5000.dat", 0);

    bc.assignDBC(pointset.points, pointset.TotalPoints, "y", 0, 1);
    
    cout << "Assigning solver settings..."<< endl;
    SOLVER_SETTINGS settings;
    settings.kernel.recompute_approximants = true;
    settings.kernel.support_method = "BF";
    settings.kernel.radius_ratio = 2;
    settings.kernel.approximation_method = "SPH";
    settings.kernel.approximation_order = 1;
    settings.kernel.RBF = "GE";
    settings.kernel.RBF_alpha = 0.01;
    settings.kernel.WLS = "LS";
    settings.kernel.MLS = "S4";
    settings.kernel.SPH = "GE";
    settings.dt = 0.0001;
    settings.nt = 100;
    settings.prnt_freq = 10;

    cout << "Running solver..."<< endl;
    clock_t start, end;
    start = clock();
    solveTransport(pointset, materials, bc, settings, U, Vx, Vy, Vz, Q);
    end = clock();
    double time_taken = double(end - start) / double(CLOCKS_PER_SEC);
    cout << "Time taken by program is : " << time_taken << setprecision(10) << " seconds" << endl;

    return 0;
}
