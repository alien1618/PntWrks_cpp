#include "../pntwrks.h"

vector<POINT> seed_center;

//********************************************************************************
// FUNCTIONS FOR COLLECTING THE MESHFREE INTERFACE POINTS
//********************************************************************************
INTERFACE_PRMTRS::INTERFACE_PRMTRS()
{
    delta_ratio = 1.5;
    surface_tension.function = 2;
    surface_tension.theta = 0;
    surface_tension.d0 = 0.0001;
    surface_tension.A = 15;
    surface_tension.b = 4;
    surface_tension.curv_effect = false;
    surface_tension.kinetic_mobility = false; 
    PFM.mobility = 1;
    PFM.tau0 = 0.0003;
    PFM.alpha = 0.9;         //smaller it is the smaller the solidification rate
    PFM.gamma = 10.0;        //the higher it is the more non-constant the interface thickness. smaller it is the more the surface tension
    TotalSeeds = 1;
    sharpen_interface = false;
    interp_at_bndry = false;
    surface_tension.surface_tension_co = 0;
}
INTERFACE_PRMTRS::~INTERFACE_PRMTRS()
{
    //destructor
}

void INTERFACE_PRMTRS::seeds(int num)
{
    seed_theta.resize(num+1);
    seed_center_x.resize(num+1);
    seed_center_y.resize(num+1);
    seed_center_z.resize(num+1);
    TotalSeeds = num;
}

POINTSET generateInterfacePoints(vector<double> phi, vector<POINT> points, vector<KERNEL> kernels, int TotalPoints, double AvgPointSpacing)
{
    //Function works for a VOF function (0 to 1)
    int TotalInterfacePoints = 0;
    vector<POINT> ip(1);
    POINT spacer;
    for(int i = 1; i <= TotalPoints; i++)
    {
        double phi1 = phi[i];
        for (int j = 1; j <= kernels[i].TotalNbrs; j++)
        {
            int pntnum = kernels[i].nbrs[j];
            double phi2 = phi[pntnum];
            double product = phi1 * phi2;
            if (product <= 1e-6)
            {
                double x1 = points[i].x;
                double y1 = points[i].y;
                double z1 = points[i].z;
                double x2 = points[pntnum].x;
                double y2 = points[pntnum].y;
                double z2 = points[pntnum].z;
                TotalInterfacePoints++;
                ip.push_back(spacer);
                ip[TotalInterfacePoints].x = (phi1*x2-phi2*x1)/(phi1-phi2);
                ip[TotalInterfacePoints].y = (phi1*y2-phi2*y1)/(phi1-phi2);
                ip[TotalInterfacePoints].z = (phi1*z2-phi2*z1)/(phi1-phi2);
            }
        }
    }

    int TotalUniquePoints = 0;
    vector<POINT> unique_p(2);
    vector<int> index(TotalInterfacePoints+1);
    for (int i = 1; i <= TotalInterfacePoints; i++)
    {index[i] = 0;}

    for (int i = 1; i <= TotalInterfacePoints; i++)
    {
        int m = 0;
        for (int j = 1; j <= TotalInterfacePoints; j++)
        {
            double d = ip[i].distance(ip[j]);
            if (index[j] == 0 && d <= 0.5*AvgPointSpacing)
            {
                m++;
                index[j] = m;
            }
        }
    }
    for (int i = 1; i <= TotalInterfacePoints; i++)
    {
        if (index[i] == 1)
        {
            TotalUniquePoints++;
            unique_p.push_back(spacer);
            unique_p[TotalUniquePoints] = ip[i];
        }
    }
    //cout << "Total unique interface points generation = " << TotalUniquePoints << endl;
    POINTSET interf;
    interf.TotalPoints = TotalUniquePoints;
    interf.points = unique_p;
    return interf;
}
