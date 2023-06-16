#include "src/pntwrks.h"

int main()
{
    cout << "-----------------------------------------------------------------" << endl;
    cout << "2D bubble rise" << endl;
    cout << "-----------------------------------------------------------------" << endl;
    cout << "Constructing pointset..."<< endl;
    POINT p(0.02,0.02);
    DOMAIN d(0.96,1.5);
    RESOLUTION r(50,75);   
    POINTSET pointset(p,d,r);

    cout << "Assigning phases..."<< endl;
    POINT cntr(0.01,0.01);
    pointset.addSphere(cntr, 0.5);
    pointset.computeVOF(0.5*pointset.AvgPntSpacing);

    cout << "Constructing boundaries..."<< endl;
    POINTSET bndry("geometry/SPH_Rectangle/mesh_600.dat");

    cout << "Constructing background grid..."<< endl;
    POINT p2(-0.25,-0.25);
    DOMAIN d2(1.5, 2);
    RESOLUTION r2(40, 40);   
    POINTSET bg_grid(p2,d2,r2, "Q4");

    cout << "Assigning material properties..."<< endl;
    MATERIALS materials(2);
    materials.Nu[1] = 1e-3;
    materials.ro[1] = 1;
    materials.Nu[2] = 1e-3;
    materials.ro[2] = 100;

    cout << "Assigning interface parameters..."<< endl;
    INTERFACE_PRMTRS interface_parameters;
    interface_parameters.solver = "NONE";
    interface_parameters.surface_tension.surface_tension_co = 0.0;
    interface_parameters.sharpen_interface = true;
    interface_parameters.sharpen_interface_beta = 2*pointset.AvgPntSpacing;
    interface_parameters.sharpen_interface_eps = 1*pointset.AvgPntSpacing;
    interface_parameters.sharpen_interface_dt =  10*0.00005;
    interface_parameters.sharpen_interface_nt =  5;

    cout << "Assigning boundary conditions..."<< endl;
    int s = 1;
    vector<int> bndry_pnt_num(pointset.TotalPoints+bndry.TotalPoints+1);
    vector<double> values(pointset.TotalPoints+bndry.TotalPoints+1);
    for (int i = pointset.TotalPoints+1; i <= pointset.TotalPoints+bndry.TotalPoints; i++)
    {
        bndry_pnt_num[s] = i;
        values[i] = 0;
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
    settings.P_factor = 5000;
    settings.AV_factor = 0.02;

    settings.kernel.radius_ratio = 2;
    settings.kernel.support_method = "BF";
    settings.kernel.approximation_method = "RBF";

    settings.lagrangian = true;
    settings.particle_shifting = true;
    settings.particle_shifting_factor = 0.05;
    settings.particle_shifting_surf = 1;
    settings.gravity.x = 0;
    settings.gravity.y = -9.8;
    settings.gravity.z = 0;
    settings.prnt_freq = 100;

    cout << "Running SPH Solver..."<< endl;
    solveNavierStokesLagrangian(pointset, bndry, bg_grid, materials, interface_parameters, bc, settings, pointset.phi);
    
    return 0;
}

