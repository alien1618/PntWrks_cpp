#include "src/pntwrks.h"

int main()
{
    cout << "-----------------------------------------------------------------" << endl;
    cout << "2D tank filling" << endl;
    cout << "-----------------------------------------------------------------" << endl;
    cout << "Constructing pointset..."<< endl;
    POINT p(0.8,2.49);
    DOMAIN d(1.125,0.725);
    RESOLUTION r(62,40);
    POINTSET pointset(p,d,r);

    cout << "Constructing boundaries..."<< endl;
    POINTSET bndry("geometry/SPH_TankFilling/mesh_1000.dat");

    cout << "Constructing background grid..."<< endl;
    POINT p2(-0.25,-0.25);
    DOMAIN d2(1.5, 2);
    RESOLUTION r2(40, 40);
    POINTSET bg_grid(p2,d2,r2, "Q4");

    cout << "Assigning material properties..."<< endl;
    MATERIALS materials(1);
    materials.Nu[1] = 0.001;
    materials.ro[1] = 1;

    cout << "Assigning interface parameters..."<< endl;
    INTERFACE_PRMTRS interface_parameters;
    interface_parameters.surface_tension.surface_tension_co = 0;

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

    cout << "Assigning solver settings..."<< endl;
    SOLVER_SETTINGS settings;
    settings.dt = 0.00005;
    settings.nt = 50000;
    settings.P_factor = 1000;
    settings.P_iterations = 1;
    settings.AV_factor = 1e-3;
    settings.kernel.radius_ratio = 1.5;
    settings.kernel.support_method = "BF";
    settings.kernel.approximation_method = "SPH";
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
    solveNavierStokesLagrangian(pointset, bndry, bg_grid, materials, interface_parameters, bc, settings, phi);
}
