#include "src/pntwrks.h"

int main()
{
    cout << "-----------------------------------------------------------------" << endl;
    cout << "2D laplace equation" << endl;
    cout << "-----------------------------------------------------------------" << endl;

    cout << "Constructing pointset..."<< endl;
    POINTSET pointset("geometry/PlateWithHole2D/mesh_400.dat");
    //POINTSET pointset("geometry/PlateWithHole2DO2/mesh_1200.dat");
    //POINTSET pointset("geometry/SlottedRectangle2D/mesh_3000.dat");
    
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
    
    //bc.assignDBC(pointset, "geometry/PlateWithHole2DO2/bottom_1200.dat", 1);
    //bc.assignDBC(pointset, "geometry/PlateWithHole2DO2/top_1200.dat", 0);
    
    //bc.assignDBC(pointset, "geometry/SlottedRectangle2D/hot_3000.dat", 35);
    //bc.assignDBC(pointset, "geometry/SlottedRectangle2D/cold_3000.dat", 0);
    
    cout << "Assigning solver settings..."<< endl;
    SOLVER_SETTINGS settings;
    settings.kernel.recompute_approximants = true;
    settings.kernel.support_method = "BF";
    settings.kernel.radius_ratio = 2;
    settings.kernel.approximation_method = "RBF";
    settings.kernel.approximation_order = 1;
    settings.kernel.RBF = "GE";
    settings.kernel.RBF_alpha = 0.1;
    settings.kernel.WLS = "S4";
    settings.kernel.MLS = "S4";
    settings.kernel.SPH = "WC4";
    settings.dt = 0.00001;//0.01;
    settings.nt = 50000;
    settings.prnt_freq = 100;

    cout << "Running solver..."<< endl;
    solveTransportSteadyState(pointset, materials, bc, settings, U, Vx, Vy, Vz, Q);
    pointset.printMesh();

    return 0;
}
