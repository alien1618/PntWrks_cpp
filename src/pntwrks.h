#ifndef PNTWRKS_H
#define PNTWRKS_H

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <omp.h>
#include <cstdlib>
#include <sstream>
#include <tuple>
#include <sys/stat.h>
#include <sys/types.h>
#include <iomanip>      // std::setprecision

using namespace std;

static string input_dir = "logs";
static string output_dir = "out";

struct PFM_PRMTRS
{
    double tau0; //the higher it is the more surface tension of the interface
    double alpha;     //smaller it is the smaller the solidification rate
    double gamma;   //the higher it is the more non-constant the interface thickness. smaller it is the more the surface tension
    double mobility;
};  

struct GRAVITY
{
    double x;
    double y;
    double z;
};

struct KERNEL_PRMTRS
{
    string WLS;
    string MLS;
    string SPH;
    string RBF;
    double RBF_alpha;
    string approximation_method;
	int approximation_order;
	bool recompute_approximants;
    int grad_order;
    double radius_ratio;
    string support_method;
    bool upwind;
    double upwind_ratio;
};

struct SURFACE_TENSION_PRMTRS
{
    double theta;
    double d0;
    double A;
    double b;
    double m;
    double function;
    bool curv_effect;
    double surface_tension_co;
    double kinetic_mobility_co;
    bool kinetic_mobility;
};

struct INTERFACE_THERMAL_PRMTRS
{
    double T_eq;
    double fusion_latent_heat;
};

class INTERFACE_PRMTRS
{
private:

public:
    INTERFACE_PRMTRS ();
    double delta_ratio;
    INTERFACE_THERMAL_PRMTRS thermal;
    SURFACE_TENSION_PRMTRS surface_tension;
    PFM_PRMTRS PFM;
    string solver;

    bool sharpen_interface;
    double sharpen_interface_beta;
    double sharpen_interface_eps;
    double sharpen_interface_dt;
    int sharpen_interface_nt;

    bool interp_at_bndry;

    //----------------------------------------------------------------
    // SEED PARAMETERS
    //----------------------------------------------------------------
    int TotalSeeds;
    void seeds(int num);
    vector<double> seed_theta;
    vector<double> seed_center_x;
    vector<double> seed_center_y;
    vector<double> seed_center_z;
 
    ~INTERFACE_PRMTRS();
};
class SOLVER_SETTINGS
{
private:

public:
    SOLVER_SETTINGS();

	string simulation_type;
    int prnt_freq;
    GRAVITY gravity;
    KERNEL_PRMTRS kernel;
    double dt;
    int nt;
    int P_iterations;
    double P_factor; // Artificial compressibility pressure factor
    double T_factor; //Rayleigh-Benard Temerature force factor
    double AV_factor; //Artificial viscosity factor
    double Ra; //Rayleigh number
    double Pr; //Prentle number
    double Thot;
	double Tcold;
	bool lagrangian;
    bool particle_shifting;
    double particle_shifting_factor;
    double particle_shifting_surf;
};

class POINT
{
private:

public:
    double x;
    double y;
    double z;

    POINT();
    POINT(double);
    POINT(double, double);
    POINT(double, double, double);
    double distance(POINT point);
    ~POINT();
};
class DOMAIN
{
private:

public:
    double x;
    double y;
    double z;
    
	DOMAIN();
	DOMAIN(double);
    DOMAIN(double, double);
    DOMAIN(double, double, double);
    ~DOMAIN();
};
class RESOLUTION
{
private:

public:
    int x;
    int y;
    int z;
      
    RESOLUTION();
    RESOLUTION(int);
    RESOLUTION(int, int);
    RESOLUTION(int, int, int);
    ~RESOLUTION();
};

class POINTSET
{

private:
    public:
    POINTSET ();
    POINTSET (string input_txt);
    POINTSET (POINT p, DOMAIN l, RESOLUTION size, string elemtype);
    POINTSET(POINT p, DOMAIN l, RESOLUTION divss);

    //----------------------------------------------------------------
    // GEOEMTRY DATA
    //----------------------------------------------------------------
    int dim;
    double AvgPntSpacing;
    int TotalPoints;

    int ElementOrder;

    int TotalLines;
    int LineShape;
    int LineNds;

    int TotalSurfaces;
    int SurfaceShape;
    int SurfaceNds;

    int TotalVolumes;
    int VolumeShape;
    int VolumeNds;

    vector<POINT> points;
    vector<vector<int>> lines;
	vector<vector<int>> surfaces;
    vector<vector<int>> volumes;
    vector<double> phi;

    //------------------------------------------
    //Basic Geometry Functions
    //------------------------------------------
    void computeAvgPntSpacing();
    void computeHyperbolicTangent(double W);
    void computeVOF(double W);

    //------------------------------------------
    //Generation of Structured Grids Functions
    //------------------------------------------
    void generateGridPnts1D(POINT p, DOMAIN l, RESOLUTION div);
    void generateGridPnts2D(POINT p, DOMAIN l, RESOLUTION div);
    void generateGridPnts3D(POINT p, DOMAIN l, RESOLUTION div);
    void generateGridQuad(RESOLUTION div);
    void generateGridTri(RESOLUTION div);
    void generateGridTri2(RESOLUTION div);
    void generateGridTet(POINT p, DOMAIN size, RESOLUTION divs);
    void generateGridHex(POINT p, DOMAIN size, RESOLUTION divs);

    //------------------------------------------
    //Implicit Geometry Modeling Functions
    //------------------------------------------
    void addCube(POINT base_pnt, DOMAIN size);
    void addSphere(POINT pnt, double a);
    void addCylinderZ(POINT pnt, double radius);
    void addCylinderY(POINT pnt, double radius);
    void addCylinderX(POINT pnt, double radius);
    void subCube(POINT base_pnt, DOMAIN size);
    void subSphere(POINT cntr_pnt, double a);
    void subCylinderZ(POINT cntr_pnt, double a);
    void subCylinderY(POINT cntr_pnt, double a);
    void subCylinderX(POINT cntr_pnt, double a);
	void addStar(POINT cntr, double a, double b, double d, double theta0, int branches);
    
    //------------------------------------------
    //Printing Geometric Data
    //------------------------------------------
    void printStructuredGridVTK(vector<double> U, RESOLUTION divs, string filename, int t);
    void printUnstructuredGridVTK(vector<double> U, string filename, int t);
    void printMeshVTK(vector<double> U, string filename, int t);
    void printMeshVTK(vector<int> U, string filename, int t);
    void printMeshVTK(vector<vector<double> > u, RESOLUTION g, string filename, int t);
    void printTXT(vector<double> U, string filename);
    void printTXT(vector<double> U, string filename, int t);
    void printTXT(vector<int> U, string filename, int t);
    void printVTK(vector<double> U, string filename, int t);
    void printVTK(vector<int> U, string filename, int t);
    void printVectorsVTK(vector<double> U, vector<double> Vx, vector<double> Vy, vector<double> Vz, string filename, int t);
    void printMesh();
    ~POINTSET ();
};
void pltctrl(vector<POINT> points, int TotalPoints, int n_t, int prnt_freq);
void printElements(vector<vector<int>> surfaces, int SurfaceNds, int TotalSurfaces, string filename);

class BOUNDARY_CONDITION
{
private:

public:
    int total;
	vector<int> points;
	vector<double> values;

    void assignDBC(POINTSET pointset, string salomeface_fname, double value);
    void assignDBC(vector<POINT> points, int TotalPoints, string edge, double location, double value);
   
    BOUNDARY_CONDITION();
    ~BOUNDARY_CONDITION ();
};

class BOUNDARY_CONDITIONS
{
	private:

	public:
	BOUNDARY_CONDITIONS();
	~BOUNDARY_CONDITIONS();
    
	BOUNDARY_CONDITION Vx;
    BOUNDARY_CONDITION Vy;
    BOUNDARY_CONDITION Vz;
    BOUNDARY_CONDITION P;
    BOUNDARY_CONDITION U;
    BOUNDARY_CONDITION phi;
};

class KERNEL
{
private:
    int monomials;
    void bases(int order, double x, double y, double z);

    vector<vector<double> > MLS_Spline3(POINT pnt, vector<POINT> points, vector<int> nbrpnts, int totnbrpnts, double dmax);
    vector<vector<double> > MLS_Spline4(POINT pnt, vector<POINT> points, vector<int> nbrpnts, int totnbrpnts, double dmax);
    vector<vector<double> > MLS_Spline5(POINT pnt, vector<POINT> points, vector<int> nbrpnts, int totnbrpnts, double dmax);
    vector<vector<double> > MLS_RegSpline4(POINT pnt, vector<POINT> points, vector<int> nbrpnts, int totnbrpnts, double dmax);

    vector<vector<double> > RBF_MQ(POINT pnt, vector<POINT> points, vector<int> nbrpnts, int totnbrpnts, double dmax, double alfc, int RBF_order);
    vector<vector<double> > RBF_GE(POINT pnt, vector<POINT> points, vector<int> nbrpnts, int totnbrpnts, double dmax, double alfc, int RBF_order);
    vector<vector<double> > RBF_TPS(POINT pnt, vector<POINT> points, vector<int> nbrpnts, int totnbrpnts, double dmax, double alfc, int RBF_order);


public:
    vector<int> nbrs;
    int TotalNbrs;
    vector<double> N;
    vector<double> dNdx;
    vector<double> dNdy;
    vector<double> dNdz;
    vector<double> dNdxx;
    vector<double> dNdyy;
    vector<double> dNdzz;
    vector<double> dNdxy;
    vector<double> dNdxz;
    vector<double> dNdyz;
    vector<double> nabla2;
    
    void collectNbrsBruteForce(POINT g, vector<POINT> points0, int TotalPoints0, double radius);
	void collectNbrsSamePhaseBruteForce(POINT g, double gphi, vector<POINT> points0, vector<double> phi, int TotalPoints0, double radius);
    void collectNbrsBackgroundGridO1(POINT pnt, vector<POINT> bg_nodes, vector<vector<int>> bg_elements, double radius);

    void computeRBF(POINT pnt, vector<POINT> points, vector<int> nbrpnts, int totnbrpnts, int weight, int order, double dmax, double RBF_alfc);
    void computeMLS(POINT pnt, vector<POINT> points, vector<int> nbrpnts, int totnbrpnts, int weight, int order, double dmax);
    void computeWLS(POINT pnt, vector<POINT> points, vector<int> nbrpnts, int totnbrpnts, int weight, int order, double dmax);
    void computeSPH(POINT pnt, vector<POINT> points, vector<int> nbrpnts, int totnbrpnts, int weight, double dmax, int dim);
    
    ~KERNEL ();
};

void writeKernels(vector<KERNEL> interp, int TotalPoints);
vector<KERNEL> readKernels(int total_nds);
vector<KERNEL> computeKernels(vector<POINT> points, vector<double> phi, int TotalPoints, int dim, KERNEL_PRMTRS kernel_prmtrs, double kernel_radius);
vector<double> computeGFD(double pnt1x, double pnt1y, double U1, vector<POINT> points, vector<double> U, vector<int> supp_pnt_nums, int totnbrpnts, double dm);

class MATERIALS
{
private:

public:
    vector<double> D;
	vector<double> Nu;
	vector<double> ro;
    int total;
    MATERIALS();
    MATERIALS(int);
    ~MATERIALS ();
};
vector<double> setVector(int Totalpoints, double val);
vector<double> setVector(vector<double> U, BOUNDARY_CONDITION BC);

//----------------------------------------------------------------
// INTERFACE EVOLUTION FUNCTIONS
//----------------------------------------------------------------
POINTSET generateInterfacePoints(vector<double> phi, vector<POINT> points, vector<KERNEL> kernels, int TotalPoints, double AvgPointSpacing);
vector<double> computeCurvature(vector<double> phi, int Totalpoints, vector<KERNEL> kernels, int dim);
vector<double> computeVolumeOfFluid(vector<KERNEL> kernels, vector<double> phi_now, int TotalPoints, vector<double> Vx, vector<double> Vy, vector<double> Vz, vector<double> Vn, int iterations, double dt);
vector<double> computeAllenCahn(double W, double M, vector<KERNEL> kernels, vector<double> phi, double dim, vector<POINT> points, int TotalPoints,  vector<double> Vx, vector<double> Vy, vector<double> Vz, vector<double> Vn, int iterations, double dt);
vector<double> computeAllenCahn2(double W, double M, vector<KERNEL> kernels, vector<double> phi_now, double dim, vector<POINT> points, int TotalPoints,  vector<double> Vx, vector<double> Vy, vector<double> Vz, vector<double> Vn, int iterations, double dt);
vector<double> computeCahnHilliard(double W, double M, vector<KERNEL> kernels, vector<double> phi, double dim, vector<POINT> points, int TotalPoints,  vector<double> Vx, vector<double> Vy, vector<double> Vz, vector<double> Vn, int iterations, double dt);
vector<double> computeCahnHilliard2(double W, double M, vector<KERNEL> kernels, vector<double> phi_now, double dim, vector<POINT> points, int TotalPoints,  vector<double> Vx, vector<double> Vy, vector<double> Vz, vector<double> Vn, int iterations, double dt);
vector<double> computeAllenCahnPhaseTransformation2D(int dim, vector<double> phi1, vector<double> T1, vector<double> D, INTERFACE_PRMTRS inter_param, double dt, vector<POINT> particles, int TotalParticles, vector<KERNEL> kernels);
vector<double> computeSharpenInterface(bool vof, double beta, double epsilon, vector<double> phi, int TotalPoints, vector<KERNEL> kernels, double dx,  double dt, int nt);

//----------------------------------------------------------------
// PREDEFINED INTERFACE ADVECTION VELOCITY FUNCTIONS
//----------------------------------------------------------------
vector<double> assignConstantExtendedVelocity(int TotalPoints, double b);
vector<double> assign4BranchExtendedVelocity(vector<POINT> points, int TotalPoints, double b);
vector<double> assign6BranchExtendedVelocity(vector<POINT> points, int TotalPoints, double b);
tuple<vector<double>,vector<double>,vector<double> > assignSingleVortexVelocity(vector<POINT> points, int TotalPoints);
tuple<vector<double>,vector<double>,vector<double> > assignSingleVortexVelocity(vector<POINT> points, int TotalPoints, double u, double v, double w);
tuple<vector<double>,vector<double>,vector<double> > assignSingleVortexVelocity2(vector<POINT> points, int TotalPoints, double u, double v, double w);
tuple<vector<double>,vector<double>,vector<double> > assignRotatingAdvectiveVelocity(vector<POINT> points, int TotalPoints, double u, double v, double w);
tuple<vector<double>,vector<double>,vector<double> > assignConstantAdvectiveVelocity(int TotalPoints, double u, double v, double w);
tuple<vector<double>,vector<double>,vector<double> > assignExtremeDeformationVelocity(vector<POINT> points, int TotalPoints);
tuple<vector<double>,vector<double>,vector<double> > assignAdvectionFlowVelocity(vector<POINT> points, int TotalPoints);
tuple<vector<double>,vector<double>,vector<double> > assign3DSingleVortexFlow(vector<POINT> points, int TotalPoints);

//----------------------------------------------------------------
// TRANSPORT EQUATION CORE SOLVERS
//----------------------------------------------------------------
tuple<vector<double>, vector<double> > computeTransport(int dim, vector<double> T, vector<double> Vx, vector<double> Vy, vector<double> Vz, vector<double> Q, vector<double> D, double dt, int t, vector<double> RHS, vector<POINT> points, vector<double> vof, int TotalPoints, vector<KERNEL> kernels, BOUNDARY_CONDITION BC_U);
vector<double> computeTransportGaussSeidel(int dim, vector<double> T, vector<double> Vx, vector<double> Vy, vector<double> Vz, vector<double> Q, vector<double> D, double dt, vector<POINT> points, vector<double> vof, int TotalPoints, vector<KERNEL> kernels, BOUNDARY_CONDITION BC_U);
vector<double> computeTransportUpwind(int dim, vector<double> T, vector<double> Vx, vector<double> Vy, vector<double> Vz, vector<double> Q, vector<double> vof, vector<double> diffusivity, int TotalPhases, double AV_factor, double dt, vector<POINT> points, int TotalPoints, vector<KERNEL> kernels, BOUNDARY_CONDITION BC_U, double upwind_ratio);
vector<double> computeTransportIterative(int dim, vector<double> T, vector<double> T_iter, vector<double> Vx, vector<double> Vy, vector<double> Vz, vector<double> Q, vector<double> vof, vector<double> diffusivity, int TotalPhases, double AV_factor, double dt, vector<POINT> points, int TotalPoints, vector<KERNEL> kernels, BOUNDARY_CONDITION BC_U);
vector<double> computeTransportGFD(double radius, vector<double> T, vector<double> Vx, vector<double> Vy, vector<double> Q, vector<double> vof, vector<double> diffusivity, int TotalPhases, double dt, vector<POINT> points, int TotalPoints, vector<KERNEL> kernels, BOUNDARY_CONDITION BC_U);
tuple<vector<double>,vector<double>,vector<double>,vector<double>,vector<double>>  computePeaksFunction(vector<POINT> points, int TotalPoints);

//----------------------------------------------------------------
// NAVIER-STOKES CORE SOLVERS
//----------------------------------------------------------------
tuple<vector<double>, vector<double>,vector<double>, vector<double>,vector<double>, vector<double> > solveNavierStokes(vector<double> P, vector<double> Vx1, vector<double> Vy1, vector<double> Vz1, vector<double> Fx, vector<double> Fy, vector<double> Fz, vector<double> phi, vector<POINT> points, int TotalPoints, vector<KERNEL> kernels, BOUNDARY_CONDITION BCs_Vx, BOUNDARY_CONDITION BCs_Vy, BOUNDARY_CONDITION BCs_Vz, BOUNDARY_CONDITION BCs_P, vector<double> Nu, vector<double> ro, int TotalPhases, SOLVER_SETTINGS settings, double AvgPntSpacing, double segma);

tuple<vector<double>, vector<double>,vector<double>, vector<double>,vector<double>, vector<double> > solveNavierStokesUpwind(vector<double> P0, vector<double> Vx1, vector<double> Vy1, vector<double> Vz1, vector<double> Fx, vector<double> Fy, vector<double> Fz, vector<double> phi, vector<POINT> points, int TotalPoints, vector<KERNEL> kernels, BOUNDARY_CONDITION BCs_Vx, BOUNDARY_CONDITION BCs_Vy, BOUNDARY_CONDITION BCs_Vz, BOUNDARY_CONDITION BCs_P, vector<double> Nu, vector<double> ro0, int TotalPhases, SOLVER_SETTINGS settings, double AvgPointSpacing, double segma, int dim);

vector<vector<double> > computeParticleShiftingVelocity(double factor, double srf_mrkr, double h, double dt, vector<POINT> pnts, int TotalPoints, vector<KERNEL> kernels);
vector<vector<double> > computeParticleShiftingVelocity_v2(double factor, double srf_mrkr, double h, double dt, vector<POINT> pnts, int TotalPoints, vector<KERNEL> kernels);

//----------------------------------------------------------------
// MISC FUNCTIONS
//----------------------------------------------------------------
vector<double> solve(vector<vector<double> > K, vector<double> F, int TotalPoints);
vector<vector<double> > inv(vector<vector<double> > a, int size);

template <class type>
type maximum(vector<type> V, int z)
{
    type max = V[1];
        #pragma omp parallel for
        for (int i = 2; i <= z; i++)
        {
            type number1 = V[i];
            if (number1 >= max)
            {max = number1;}
        }
    return max;
}

template <typename type>
type minimum(vector<type> V, int z)
{
    type min = V[1];
        #pragma omp parallel for
        for (int i = 2; i <= z; i++)
        {
            type number1 = V[i];
            if (number1 <= min)
            {min = number1;}
        }
    return min;
}

tuple<vector<int>, int> UniqueNums(vector<int> V, int z);

double mod(double number, double divnum);
double maxnum(double a, double b);
double minnum(double a, double b);

double checkConvergence(int iternum, vector<double> U_temp, vector<double> U_temp_old, int TotalPoints);

void randomSeed();
double random0to1();

template <class type>
double sign(type num)
{
    double num_sign=0;

    if (num > 0)
    {num_sign = 1;}
    else if (num < 0)
    {num_sign = -1;}
    return num_sign;
}
template <class type>
double heaveside(type phi)
{
    double num_sign=0;
    if (phi > 0)
    {num_sign = 1;}
    else if (phi < 0)
    {num_sign = 0;}
    else if (phi == 0)
    {num_sign = 0;}
    return num_sign;
}

//----------------------------------------------------------------
// MAIN SOLVERS
//----------------------------------------------------------------
void solveReactionDiffusion(POINTSET pointset, MATERIALS materials, BOUNDARY_CONDITION BC, SOLVER_SETTINGS settings, double f, double k, vector<double> U_A, vector<double> U_B);
void solveTransport(POINTSET pointset, MATERIALS materials, BOUNDARY_CONDITION BC, SOLVER_SETTINGS settings, vector<double> U, vector<double> Vx, vector<double> Vy, vector<double> Vz, vector<double> Q);
void solveTransportSteadyState(POINTSET pointset, MATERIALS materials, BOUNDARY_CONDITION BC, SOLVER_SETTINGS settings, vector<double> U, vector<double> Vx, vector<double> Vy, vector<double> Vz, vector<double> Q);
void solveTranportGFD(POINTSET pointset, MATERIALS materials, BOUNDARY_CONDITION BC, SOLVER_SETTINGS settings, vector<double> U, vector<double> Vx, vector<double> Vy, vector<double> Q, vector<double> phi);
void solveBurger(POINTSET pointset, MATERIALS materials, BOUNDARY_CONDITION BC, SOLVER_SETTINGS settings, vector<double> U);

void solveInterfaceEvolution(POINTSET pointset, INTERFACE_PRMTRS inter_param, SOLVER_SETTINGS settings, vector<double> Vn_ext, vector<double> Vx, vector<double> Vy, vector<double> Vz);
POINTSET solveCurvatureDrivenInterface(POINTSET pointset, INTERFACE_PRMTRS inter_param, SOLVER_SETTINGS settings);

void solveNavierStokes1Phase(POINTSET pointset, MATERIALS materials, BOUNDARY_CONDITIONS BC, SOLVER_SETTINGS settings, vector<double> U, vector<double> Vx, vector<double> Vy, vector<double> Vz, vector<double> Q, vector<double> Fx,vector<double> Fy,vector<double> Fz);
void solveNavierStokes2Phase(POINTSET pointset, MATERIALS materials, INTERFACE_PRMTRS inter_param, BOUNDARY_CONDITIONS BC, SOLVER_SETTINGS settings, vector<double> U, vector<double> Vx, vector<double> Vy, vector<double> Vz, vector<double> Q, vector<double> Fx, vector<double> Fy, vector<double> Fz);
void solveNavierStokesLagrangian(POINTSET pointset, POINTSET bndry, POINTSET bg_grid, MATERIALS materials,  INTERFACE_PRMTRS inter_param, BOUNDARY_CONDITIONS bc, SOLVER_SETTINGS settings, vector<double> phi0);

void solveAllenCahn(POINTSET pointset, MATERIALS materials, INTERFACE_PRMTRS inter_param, BOUNDARY_CONDITIONS BC, SOLVER_SETTINGS settings, vector<double> U, vector<double> Vx, vector<double> Vy, vector<double> Vz);
void solveAllenCahnMultipleSeeds(POINTSET pointset, MATERIALS materials, INTERFACE_PRMTRS inter_param, BOUNDARY_CONDITIONS bc, SOLVER_SETTINGS settings, vector<vector<double> > phi0, vector<double> U, vector<double> Vx, vector<double> Vy, vector<double> Vz);
void solveCahnHilliard(POINTSET pointset, INTERFACE_PRMTRS inter_param, SOLVER_SETTINGS settings, vector<double> Vn, vector<double> Vx, vector<double> Vy, vector<double> Vz);

#endif
