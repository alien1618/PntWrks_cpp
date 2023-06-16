#include "src/pntwrks.h"

int main()
{
    cout << "-----------------------------------------------------------------" << endl;
    cout << "2D hour glass" << endl;
    cout << "-----------------------------------------------------------------" << endl;
    cout << "Constructing pointset..."<< endl;
    POINT p(0.1,1.3);
    DOMAIN d(0.8,0.4);
    RESOLUTION r(50,25);
    POINTSET pointset(p,d,r);

    cout << "Constructing boundaries..."<< endl;
    POINTSET bndry("geometry/SPH_HourGlass/mesh_500.dat");

    cout << "Constructing background grid..."<< endl;
    POINT p2(-0.25,-0.125);
    DOMAIN d2(1.5, 2);
    RESOLUTION r2(10, 15);
    POINTSET bg_grid(p2,d2,r2, "Q4");

    cout << "Assigning material properties..."<< endl;
    MATERIALS materials(1);
    materials.Nu[1] = 0.001;
    materials.ro[1] = 1000;

    cout << "Assigning boundary conditions..."<< endl;
    int s = 1;
    vector<int> bndry_pnt_num(pointset.TotalPoints+bndry.TotalPoints+1);
    vector<double> values(pointset.TotalPoints+bndry.TotalPoints+1);
    for (int i = pointset.TotalPoints+1; i <= pointset.TotalPoints+bndry.TotalPoints; i++)
    {
        bndry_pnt_num[s] = i;
        values[s] = 0;
        s++;
    }

    BOUNDARY_CONDITIONS bc;
    bc.Vx.total = bndry.TotalPoints;
    bc.Vx.points = bndry_pnt_num;
    bc.Vx.values = values;
    
    bc.Vy.total = bndry.TotalPoints;
    bc.Vy.points = bndry_pnt_num;
    bc.Vy.values = values;
    
    bc.Vz.total = bndry.TotalPoints;
    bc.Vz.points = bndry_pnt_num;
    bc.Vz.values = values;
    
    cout << "Assigning interface parameters..."<< endl;
    INTERFACE_PRMTRS interface_parameters;
    interface_parameters.surface_tension.surface_tension_co = 0;

    cout << "Assigning solver settings..."<< endl;
    SOLVER_SETTINGS settings;
    settings.dt = 0.000025;
    settings.nt = 300000;
    settings.kernel.radius_ratio = 1.5;
    settings.kernel.support_method = "BF";
    settings.kernel.approximation_method = "SPH";
    settings.P_factor = 1000;
    settings.P_iterations = 1;
    settings.AV_factor = 1e-3;
    settings.lagrangian = true;
    settings.particle_shifting = true;
    settings.particle_shifting_factor = 0.5;
    settings.particle_shifting_surf = 1;
    settings.gravity.x = 0;
    settings.gravity.y = -9.8;
    settings.gravity.z = 0;
    settings.prnt_freq = 100;

    cout << "Running SPH Solver..."<< endl;
    vector<double> phi = setVector(pointset.TotalPoints, 1);
    clock_t start, end;
    start = clock();
    solveNavierStokesLagrangian(pointset, bndry, bg_grid, materials, interface_parameters, bc, settings, phi);
    end = clock();
    double time_taken = double(end - start) / double(CLOCKS_PER_SEC);
    cout << "Time taken by program is : " << time_taken << setprecision(10) << endl;

    return 0;
}
