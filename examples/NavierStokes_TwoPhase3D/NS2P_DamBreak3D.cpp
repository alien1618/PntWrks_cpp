#include "../fpm3d/fpm3d.h"

void testNS2P_RisingBubble3D()
{
    cout << "-----------------------------------------------------------------" << endl;
    cout << "3D dam break" << endl;
    cout << "-----------------------------------------------------------------" << endl;
    cout << "Constructing pointset..."<< endl;
    POINT p(0,0,0);
    DOMAIN s(1,1,0.1);
    RESOLUTION divs(50,50,2);
    POINTSET pointset(p, s, divs, "T4");

    cout << "Assigning phases..."<< endl;
    for (int i = 1; i <= pointset.TotalPoints; i++)
    {pointset.phi[i] = -(pointset.points[i].y -0.4);}
    POINT cntr(0.5,0.20,0);
    pointset.addCylinderZ(cntr, 0.15);
    pointset.computeVOF(0.5*pointset.AvgPntSpacing);

    cout << "Assigning material properties..."<< endl;
    MATERIALS materials(2);
    materials.Nu[1] = 1000*1e-6;
    materials.ro[1] = 1000;
    materials.Nu[2] = 1000*1.48e-5;
    materials.ro[2] = 1;

    cout << "Initializing field variables..."<< endl;
	vector<double> U = setVector(pointset.TotalPoints,0);
	vector<double> Vx = setVector(pointset.TotalPoints,0);
	vector<double> Vy = setVector(pointset.TotalPoints,0);
	vector<double> Vz = setVector(pointset.TotalPoints,0);
	vector<double> Q = setVector(pointset.TotalPoints,0);
	vector<double> Fx = setVector(pointset.TotalPoints,0);
	vector<double> Fy = setVector(pointset.TotalPoints,-9.8);
	vector<double> Fz = setVector(pointset.TotalPoints,0);
	
    cout << "Assigning interface parameters..."<< endl;
    INTERFACE_PRMTRS interf_prmtrs;
    interf_prmtrs.solver  = "VOF";
    interf_prmtrs.interp_at_bndry = true;
    interf_prmtrs.sharpen_interface = true;
    interf_prmtrs.sharpen_interface_beta = 10*pointset.AvgPntSpacing;
    interf_prmtrs.sharpen_interface_eps = 2*pointset.AvgPntSpacing;
    interf_prmtrs.sharpen_interface_dt =  1000*0.00001;
    interf_prmtrs.sharpen_interface_nt =  5;

    cout << "Assigning boundary conditions..."<< endl;
    BOUNDARY_CONDITIONS bc;
    bc.Vx.assignDBC(pointset.points, pointset.TotalPoints, "x", 0, 0);
    bc.Vy.assignDBC(pointset.points, pointset.TotalPoints, "x", 0, 0);
    bc.Vz.assignDBC(pointset.points, pointset.TotalPoints, "x", 0, 0);

    bc.Vx.assignDBC(pointset.points, pointset.TotalPoints, "x", 1, 0);
    bc.Vy.assignDBC(pointset.points, pointset.TotalPoints, "x", 1, 0);
    bc.Vz.assignDBC(pointset.points, pointset.TotalPoints, "x", 1, 0);

    bc.Vx.assignDBC(pointset.points, pointset.TotalPoints, "y", 0, 0);
    bc.Vy.assignDBC(pointset.points, pointset.TotalPoints, "y", 0, 0);
    bc.Vz.assignDBC(pointset.points, pointset.TotalPoints, "y", 0, 0);

    bc.Vx.assignDBC(pointset.points, pointset.TotalPoints, "y", 1, 0);
    bc.Vy.assignDBC(pointset.points, pointset.TotalPoints, "y", 1, 0);
    bc.Vz.assignDBC(pointset.points, pointset.TotalPoints, "y", 1, 0);

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
    settings.dt = 0.00001;
    settings.nt = 5000000;
    settings.P_factor = 1.0/settings.dt;
    settings.AV_factor = 0.1;
    settings.gravity.y = -9.8;
    settings.prnt_freq = 2000;

    cout << "Running solver..."<< endl;
    solveNavierStokes2Phase(pointset, materials, interf_prmtrs, bc, settings, U, Vx, Vy, Vz, Q, Fx, Fy, Fz);
}
void testNS2P_Droplet3D()
{
    cout << "Constructing pointset..."<< endl;
    POINT p(0,0,0);
    DOMAIN s(1,1,0.1);
    RESOLUTION divs(50,50,2);
    POINTSET pointset(p, s, divs, "T4");

    cout << "Assigning phases..."<< endl;
    for (int i = 1; i <= pointset.TotalPoints; i++)
    {pointset.phi[i] = (pointset.points[i].y -0.25);}
    POINT cntr(0.5,0.45,0);
    pointset.addCylinderZ(cntr, 0.15);
    pointset.computeVOF(0.5*pointset.AvgPntSpacing);

    cout << "Assigning material properties..."<< endl;
    MATERIALS materials(2);
    materials.Nu[1] = 1000*1.48e-5;
    materials.ro[1] = 1;
    materials.Nu[2] = 1000*1.0e-6;
    materials.ro[2] = 1000;

    cout << "Initializing field variables..."<< endl;
	vector<double> U = setVector(pointset.TotalPoints,0);
	vector<double> Vx = setVector(pointset.TotalPoints,0);
	vector<double> Vy = setVector(pointset.TotalPoints,0);
	vector<double> Vz = setVector(pointset.TotalPoints,0);
	vector<double> Q = setVector(pointset.TotalPoints,0);
	vector<double> Fx = setVector(pointset.TotalPoints,0);
	vector<double> Fy = setVector(pointset.TotalPoints,-9.8);
	vector<double> Fz = setVector(pointset.TotalPoints,0);
	
    cout << "Assigning interface parameters..."<< endl;
    INTERFACE_PRMTRS interf_prmtrs;
    interf_prmtrs.solver  = "VOF";
    interf_prmtrs.sharpen_interface = true;
    interf_prmtrs.sharpen_interface_beta = 10*pointset.AvgPntSpacing;
    interf_prmtrs.sharpen_interface_eps = 2*pointset.AvgPntSpacing;
    interf_prmtrs.sharpen_interface_dt =  1000*0.00001;
    interf_prmtrs.sharpen_interface_nt =  5;
    
    cout << "Assigning boundary conditions..."<< endl;
    BOUNDARY_CONDITIONS bc;
    bc.Vx.assignDBC(pointset.points, pointset.TotalPoints, "x", 0, 0);
    bc.Vy.assignDBC(pointset.points, pointset.TotalPoints, "x", 0, 0);
    bc.Vz.assignDBC(pointset.points, pointset.TotalPoints, "x", 0, 0);

    bc.Vx.assignDBC(pointset.points, pointset.TotalPoints, "x", 1, 0);
    bc.Vy.assignDBC(pointset.points, pointset.TotalPoints, "x", 1, 0);
    bc.Vz.assignDBC(pointset.points, pointset.TotalPoints, "x", 1, 0);

    bc.Vx.assignDBC(pointset.points, pointset.TotalPoints, "y", 0, 0);
    bc.Vy.assignDBC(pointset.points, pointset.TotalPoints, "y", 0, 0);
    bc.Vz.assignDBC(pointset.points, pointset.TotalPoints, "y", 0, 0);

    bc.Vx.assignDBC(pointset.points, pointset.TotalPoints, "y", 1, 0);
    bc.Vy.assignDBC(pointset.points, pointset.TotalPoints, "y", 1, 0);
    bc.Vz.assignDBC(pointset.points, pointset.TotalPoints, "y", 1, 0);

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
    settings.dt = 0.00001;
    settings.nt = 5000000;
    settings.P_factor = 1/(settings.dt);
    settings.AV_factor = 0.1;
    settings.gravity.y = -9.8;
    settings.prnt_freq = 2000;

    cout << "Running solver..."<< endl;
    solveNavierStokes2Phase(pointset, materials, interf_prmtrs, bc, settings, U, Vx, Vy, Vz, Q, Fx, Fy, Fz);
}
void testNS2P_RayleighTaylor3D()
{
    cout << "Constructing pointset..."<< endl;
    POINT p(0,0,0);
    DOMAIN s(0.5,1,0.1);
    RESOLUTION divs(20,40,2);
    POINTSET pointset(p, s, divs, "T4");

    cout << "Assigning phases..."<< endl;
    for (int i = 1; i <= pointset.TotalPoints; i++)
    {
        pointset.phi[i] = ((pointset.points[i].y - (0.006*cos(25*pointset.points[i].x+0.1)+0.0))-0.5);
        //pointset.phi[i] = ((pointset.points[i].y - (0.006*cos(10*pointset.points[i].x+0.6)+0.0))-0.5);
    }
    pointset.computeVOF(0.5*pointset.AvgPntSpacing);

    cout << "Assigning material properties..."<< endl;
    MATERIALS materials(2);
    materials.Nu[1] = 1e-6;
    materials.ro[1] = 10;
    materials.Nu[2] = 1.48e-5;
    materials.ro[2] = 1;

    cout << "Initializing field variables..."<< endl;
	vector<double> U = setVector(pointset.TotalPoints,0);
	vector<double> Vx = setVector(pointset.TotalPoints,0);
	vector<double> Vy = setVector(pointset.TotalPoints,0);
	vector<double> Vz = setVector(pointset.TotalPoints,0);
	vector<double> Q = setVector(pointset.TotalPoints,0);
	vector<double> Fx = setVector(pointset.TotalPoints,0);
	vector<double> Fy = setVector(pointset.TotalPoints,-9.8);
	vector<double> Fz = setVector(pointset.TotalPoints,0);
	
    cout << "Assigning interface parameters..."<< endl;
    INTERFACE_PRMTRS interf_prmtrs;
    interf_prmtrs.solver  = "VOF";
    interf_prmtrs.PFM.mobility = 100;
    
    interf_prmtrs.sharpen_interface = true;
    interf_prmtrs.sharpen_interface_beta = 10*pointset.AvgPntSpacing;
    interf_prmtrs.sharpen_interface_eps = 2*pointset.AvgPntSpacing;
    interf_prmtrs.sharpen_interface_dt =  50*0.00001;
    interf_prmtrs.sharpen_interface_nt =  5;

    cout << "Assigning boundary conditions..."<< endl;
    BOUNDARY_CONDITIONS bc;
    bc.Vx.assignDBC(pointset.points, pointset.TotalPoints, "x", 0, 0);
    bc.Vy.assignDBC(pointset.points, pointset.TotalPoints, "x", 0, 0);
    bc.Vz.assignDBC(pointset.points, pointset.TotalPoints, "x", 0, 0);

    bc.Vx.assignDBC(pointset.points, pointset.TotalPoints, "x", 0.5, 0);
    bc.Vy.assignDBC(pointset.points, pointset.TotalPoints, "x", 0.5, 0);
    bc.Vz.assignDBC(pointset.points, pointset.TotalPoints, "x", 0.5, 0);

    bc.Vx.assignDBC(pointset.points, pointset.TotalPoints, "y", 0, 0);
    bc.Vy.assignDBC(pointset.points, pointset.TotalPoints, "y", 0, 0);
    bc.Vz.assignDBC(pointset.points, pointset.TotalPoints, "y", 0, 0);

    bc.Vx.assignDBC(pointset.points, pointset.TotalPoints, "y", 1, 0);
    bc.Vy.assignDBC(pointset.points, pointset.TotalPoints, "y", 1, 0);
    bc.Vz.assignDBC(pointset.points, pointset.TotalPoints, "y", 1, 0);

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
    settings.dt = 0.000001*0.5;
    settings.nt = 3000000;
    settings.P_factor = 1/(settings.dt);
    settings.prnt_freq = 10000;
    settings.gravity.y = -9.8;

    cout << "Running solver..."<< endl;
    solveNavierStokes2Phase(pointset, materials, interf_prmtrs, bc, settings, U, Vx, Vy, Vz, Q, Fx, Fy, Fz);
}
void testNS2P_DamBreak3D()
{
    cout << "Constructing pointset..."<< endl;
    POINT p(0,0,0);
    DOMAIN s(1,1,0.1);
    RESOLUTION divs(30,30,2);
    POINTSET pointset(p, s, divs, "T4");

    cout << "Assigning phases..."<< endl;
    POINT o(-0.05,-0.05);
    DOMAIN s2(0.45, 0.45);
    pointset.addCube(o, s2);
    pointset.computeVOF(0.5*pointset.AvgPntSpacing);
    pointset.printMeshVTK(pointset.phi, "phi", 0);

    cout << "Assigning material properties..."<< endl;
    MATERIALS materials(2);
    materials.Nu[1] = 1000*1.48e-5;
    materials.ro[1] = 1;
    materials.Nu[2] = 1000*1e-6;
    materials.ro[2] = 1000;

    cout << "Initializing field variables..."<< endl;
	vector<double> U = setVector(pointset.TotalPoints,0);
	vector<double> Vx = setVector(pointset.TotalPoints,0);
	vector<double> Vy = setVector(pointset.TotalPoints,0);
	vector<double> Vz = setVector(pointset.TotalPoints,0);
	vector<double> Q = setVector(pointset.TotalPoints,0);
	vector<double> Fx = setVector(pointset.TotalPoints,0);
	vector<double> Fy = setVector(pointset.TotalPoints,-9.8);
	vector<double> Fz = setVector(pointset.TotalPoints,0);
	
    cout << "Assigning interface parameters..."<< endl;
    INTERFACE_PRMTRS interf_prmtrs;
    interf_prmtrs.solver  = "VOF";
    interf_prmtrs.PFM.mobility = 5;
    
    interf_prmtrs.sharpen_interface = true;
    interf_prmtrs.sharpen_interface_beta = 10*pointset.AvgPntSpacing;
    interf_prmtrs.sharpen_interface_eps = 2*pointset.AvgPntSpacing;
    interf_prmtrs.sharpen_interface_dt =  1000*0.00001*0.25;
    interf_prmtrs.sharpen_interface_nt =  5;
    
    cout << "Assigning boundary conditions..."<< endl;
    BOUNDARY_CONDITIONS bc; //must be all wall boundary conditions to get correct velocity profile
    bc.Vx.assignDBC(pointset.points, pointset.TotalPoints, "y", 0, 0);
    bc.Vy.assignDBC(pointset.points, pointset.TotalPoints, "y", 0, 0);
    bc.Vz.assignDBC(pointset.points, pointset.TotalPoints, "y", 0, 0);

    bc.Vx.assignDBC(pointset.points, pointset.TotalPoints, "x", 0, 0);
    bc.Vy.assignDBC(pointset.points, pointset.TotalPoints, "x", 0, 0);
    bc.Vz.assignDBC(pointset.points, pointset.TotalPoints, "x", 0, 0);

    bc.Vx.assignDBC(pointset.points, pointset.TotalPoints, "x", 1, 0);
    bc.Vy.assignDBC(pointset.points, pointset.TotalPoints, "x", 1, 0);
    bc.Vz.assignDBC(pointset.points, pointset.TotalPoints, "x", 1, 0);

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

    cout << "Smoothing out the interface corners and converting to vof..."<< endl;
    settings.dt = 0.0001;
    settings.nt = 300;
    settings.prnt_freq = 100;
    pointset = solveCurvatureDrivenInterface(pointset, interf_prmtrs, settings);

    cout << "Running solver..."<< endl;
    settings.dt =  0.00001*0.25;
    settings.nt = 5000000;
    settings.P_factor = 1000;//1.0/settings.dt;
    settings.gravity.y = -9.8;
    settings.prnt_freq = 2000;
    solveNavierStokes2Phase(pointset, materials, interf_prmtrs, bc, settings, U, Vx, Vy, Vz, Q, Fx, Fy, Fz);
}

void testNS2P_KelvinHelmholtz3D()
{
    cout << "Constructing pointset..."<< endl;
    POINT p(0,0,0);
    DOMAIN d(1,0.5,0.1);
    RESOLUTION r(100,50,2);
    POINTSET pointset(p,d,r,"H8");

    cout << "Assigning material properties..."<< endl;
    MATERIALS materials(2);
    materials.Nu[1] = 1.48e-5;
    materials.ro[1] = 2;
    materials.Nu[2] = 1e-6;
    materials.ro[2] = 1;

    cout << "Assigning phases..."<< endl;
    for (int i = 1; i <= pointset.TotalPoints; i++)
    {pointset.phi[i] = ((pointset.points[i].y - (0.01*cos(25*pointset.points[i].x+0.1)+0.0))-0.25);}
    pointset.computeVOF(pointset.AvgPntSpacing);

    cout << "Initializing field variables..."<< endl;
	vector<double> U = setVector(pointset.TotalPoints,0);
	vector<double> Vx = setVector(pointset.TotalPoints,0);
	vector<double> Vy = setVector(pointset.TotalPoints,0);
	vector<double> Vz = setVector(pointset.TotalPoints,0);
	vector<double> Q = setVector(pointset.TotalPoints,0);
	vector<double> Fx = setVector(pointset.TotalPoints,0);
	vector<double> Fy = setVector(pointset.TotalPoints,-9.8);
	vector<double> Fz = setVector(pointset.TotalPoints,0);
    for (int i = 1; i <= pointset.TotalPoints; i++)
    {
        if (pointset.phi[i] <= 0.5)
        {Vx[i] = 10;}
        else
        {Vx[i] = -10;} 
    }
    
    cout << "Assigning interface parameters..."<< endl;
    INTERFACE_PRMTRS interface_parameters;
    interface_parameters.solver  = "CH2";
    interface_parameters.PFM.mobility = 100;
    interface_parameters.sharpen_interface = false;
    interface_parameters.sharpen_interface_beta = 10*pointset.AvgPntSpacing;
    interface_parameters.sharpen_interface_eps = 2*pointset.AvgPntSpacing;
    interface_parameters.sharpen_interface_dt =  50*0.00001;
    interface_parameters.sharpen_interface_nt =  5;

    cout << "Assigning boundary conditions..."<< endl;
    BOUNDARY_CONDITIONS bc;
    bc.Vx.assignDBC(pointset.points, pointset.TotalPoints, "y", 0, -10);
    bc.Vy.assignDBC(pointset.points, pointset.TotalPoints, "y", 0, 0);
    bc.Vz.assignDBC(pointset.points, pointset.TotalPoints, "y", 0, 0);

    bc.Vx.assignDBC(pointset.points, pointset.TotalPoints, "y", 0.5, 10);
    bc.Vy.assignDBC(pointset.points, pointset.TotalPoints, "y", 0.5, 0);
    bc.Vz.assignDBC(pointset.points, pointset.TotalPoints, "y", 0.5, 0);

    bc.phi.total = 0;
    bc.phi.points.resize(1);
	bc.phi.values.resize(1);
	for (int i = 1; i<= pointset.TotalPoints; i++)
    {
        if (pointset.points[i].x <= 0.001 || pointset.points[i].x >= 0.999)
        {
            bc.phi.total++;
			bc.phi.points.push_back(1);
			bc.phi.values.push_back(1);
            bc.phi.points[bc.phi.total] = i;
            if (interface_parameters.solver == "CH2"||interface_parameters.solver == "AC2")
            {bc.phi.values[bc.phi.total] = 2*pointset.phi[i]-1;}
            else
            {bc.phi.values[bc.phi.total] = pointset.phi[i];}
        }
    }
    
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
    settings.kernel.upwind = false;
    settings.kernel.upwind_ratio = 0.75;
    settings.dt = 0.0002;
    settings.nt = 300000;
    settings.P_factor = 2000;//1/(settings.temporal.dt);
    settings.P_iterations = 1;
    settings.AV_factor = 0.1;
    settings.prnt_freq = 1000;
    settings.gravity.x = 0;
    settings.gravity.y = -9.8;
    settings.gravity.z = 0;

    cout << "Running solver..."<< endl;
    solveNavierStokes2Phase(pointset, materials, interface_parameters, bc, settings, U, Vx, Vy, Vz, Q, Fx, Fy, Fz);
}

void testNS2P_Boiling3D()
{
    cout << "Constructing pointset..."<< endl;
    POINT p(0,0,0);
    DOMAIN s(1,1,0.1);
    RESOLUTION divs(20,20,2);
    POINTSET pointset(p, s, divs, "T4");

    cout << "Assigning phases..."<< endl;
    for (int i = 1; i <= pointset.TotalPoints; i++)
    {
        pointset.phi[i] = ((pointset.points[i].y )-0.01);
    }
    pointset.computeVOF(0.5*pointset.AvgPntSpacing);
    pointset.printMeshVTK(pointset.phi,"phi", 0);

    cout << "Assigning material properties..."<< endl;
    MATERIALS materials(2);
    materials.Nu[1] = 0.1;
    materials.ro[1] = 1;
    materials.D[1] = 1;
    materials.Nu[2] = 0.1;
    materials.ro[2] = 1;
    materials.D[2] = 1;
    
    cout << "Assigning interface parameters..."<< endl;
    INTERFACE_PRMTRS interf_prmtrs;
    interf_prmtrs.solver  = "VOF";
    
    cout << "Initializing field variables..."<< endl;
	vector<double> U = setVector(pointset.TotalPoints,0);
	vector<double> Vx = setVector(pointset.TotalPoints,0);
	vector<double> Vy = setVector(pointset.TotalPoints,0);
	vector<double> Vz = setVector(pointset.TotalPoints,0);
	vector<double> Q = setVector(pointset.TotalPoints,0);
	vector<double> Fx = setVector(pointset.TotalPoints,0);
	vector<double> Fy = setVector(pointset.TotalPoints,0.001);
	vector<double> Fz = setVector(pointset.TotalPoints,0);
	
    cout << "Assigning boundary conditions..."<< endl;
    BOUNDARY_CONDITIONS bc;
    bc.Vx.assignDBC(pointset.points, pointset.TotalPoints, "x", 0, 0);
    bc.Vy.assignDBC(pointset.points, pointset.TotalPoints, "x", 0, 0);
    bc.Vz.assignDBC(pointset.points, pointset.TotalPoints, "x", 0, 0);

    bc.Vx.assignDBC(pointset.points, pointset.TotalPoints, "x", 1, 0);
    bc.Vy.assignDBC(pointset.points, pointset.TotalPoints, "x", 1, 0);
    bc.Vz.assignDBC(pointset.points, pointset.TotalPoints, "x", 1, 0);

    bc.Vx.assignDBC(pointset.points, pointset.TotalPoints, "y", 0, 0);
    bc.Vy.assignDBC(pointset.points, pointset.TotalPoints, "y", 0, 0);
    bc.Vz.assignDBC(pointset.points, pointset.TotalPoints, "y", 0, 0);

    bc.Vx.assignDBC(pointset.points, pointset.TotalPoints, "y", 1, 0);
    bc.Vy.assignDBC(pointset.points, pointset.TotalPoints, "y", 1, 0);
    bc.Vz.assignDBC(pointset.points, pointset.TotalPoints, "y", 1, 0);

    bc.U.assignDBC(pointset.points, pointset.TotalPoints, "y", 1, 0);
    bc.U.assignDBC(pointset.points, pointset.TotalPoints, "y", 0, 1);
    
    double Thot = 1;
    double Tcold = 0;
    double val = Thot+Tcold/5.0;
    
    bc.U.total++;
	bc.U.points.push_back(1);
    bc.U.points[bc.U.total] = 13;
    bc.U.values.push_back(1);
    bc.U.values[bc.U.total] = val;
    
    bc.U.total++;
	bc.U.points.push_back(1);
    bc.U.points[bc.U.total] = 25;
    bc.U.values.push_back(1);
    bc.U.values[bc.U.total] = val;
    
    bc.U.total++;
	bc.U.points.push_back(1);
    bc.U.points[bc.U.total] = 38;
    bc.U.values.push_back(1);
    bc.U.values[bc.U.total] = val;
    
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
    settings.P_factor = 1000;
    settings.T_factor = 500; //use 5000 for 200X100
    settings.Thot = 1;
    settings.Tcold = 0;
    settings.dt = 0.00005; //use 0.00001 for 200X100
    settings.nt = 1000000;
    settings.prnt_freq = 2000;
    settings.gravity.y = 0.001;

    cout << "Running solver..."<< endl;
    solveNavierStokes2Phase(pointset, materials, interf_prmtrs, bc, settings, U, Vx, Vy, Vz, Q, Fx, Fy, Fz);
}
