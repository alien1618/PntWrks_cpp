#include "../pntwrks.h"

POINTSET::POINTSET(string input_txt)
{
    vector<double> read_txt;
    read_txt.push_back(1);
    double x1;
    vector<int> spacer(1);
	lines.resize(1);
    surfaces.resize(1);
    volumes.resize(1);

    //----------------------------------------------------------------------
    // READ AND STORE POINTS
    //----------------------------------------------------------------------
    cout << "Reading " << input_txt << endl;
    ifstream inFile1;
    inFile1.open(input_txt.c_str());
    if (!inFile1)
    {
        cout << endl << "ERROR: UNABLE TO OPEN POINTS FILE" << endl;
        cout << "Quitting Program" << endl << endl;
        exit(0);
    }
    else
    {
        int inputs = 0;
        while (inFile1 >> x1)
        {
            inputs = inputs + 1;
            read_txt.push_back(1);
            read_txt[inputs] = x1;
        }
        inFile1.close();

        TotalPoints = int(read_txt[1]);
        points.resize(TotalPoints+1);
        int i = 1;
        for (int k = 3; k <= (TotalPoints*4+2); k=k+4)
        {
            points[i].x = read_txt[k+1];
            points[i].y = read_txt[k+2];
            points[i].z = read_txt[k+3];
            i++;
        }
        int l = 0;
        int f = 0;
        int e = 0;
        int iter = 0;
        for (int k = ((TotalPoints*4)+3); k <= inputs; k = k+iter+2)
        {
            int entitycode = int(read_txt[k+1]);
switch (entitycode)
            {
                case 102:
                {
                    l++;
                    lines.push_back(spacer);
                    LineNds = 2;
                    iter = LineNds;
                    lines[l].resize(LineNds+1);
                    ElementOrder = 1;
                    for (int m = 1; m <= LineNds; m++)
                    {lines[l][m] = int(read_txt[k+1+m]);}
                }break;
                case 103:
                {
                    l++;
                    lines.push_back(spacer);
                    LineNds = 3;
                    iter = LineNds;
                    lines[l].resize(LineNds+1);
                    ElementOrder = 2;
                    for (int m = 1; m <= LineNds; m++)
                    {lines[l][m] = int(read_txt[k+1+m]);}
                }break;
                case 204:
                {
                    f++;
                    surfaces.push_back(spacer);
                    SurfaceNds = 4;
                    iter = SurfaceNds;
                    surfaces[f].resize(SurfaceNds+1);
                    ElementOrder = 1;
                    SurfaceShape = 2;

                    for (int m = 1; m <= SurfaceNds; m++)
                    {surfaces[f][m] = int(read_txt[k+1+m]);}
                }break;
                case 208:
                {
                    f++;
                    surfaces.push_back(spacer);
                    SurfaceNds = 8;
                    iter = SurfaceNds;
                    surfaces[f].resize(SurfaceNds+1);
                    ElementOrder = 2;
                    SurfaceShape = 2;

                    for (int m = 1; m <= SurfaceNds; m++)
                    {surfaces[f][m] = int(read_txt[k+1+m]);}
                }break;
                case 203:
                {
                    f++;
                    surfaces.push_back(spacer);
                    ElementOrder = 1;
                    SurfaceShape = 1;
                    SurfaceNds = 3;
                    iter = SurfaceNds;
                    surfaces[f].resize(SurfaceNds+1);

                    for (int m = 1; m <= SurfaceNds; m++)
                    {surfaces[f][m] = int(read_txt[k+1+m]);}
                }break;
                case 206:
                {
                    f++;
                    surfaces.push_back(spacer);
                    SurfaceNds = 6;
                    iter = SurfaceNds;
                    surfaces[f].resize(SurfaceNds+1);
                    ElementOrder = 2;
                    SurfaceShape = 1;

                    for (int m = 1; m <= SurfaceNds; m++)
                    {surfaces[f][m] = int(read_txt[k+1+m]);}
                }break;
                case 308:
                {
                    e++;
                    volumes.push_back(spacer);
                    VolumeNds = 8;
                    iter = VolumeNds;
                    volumes[e].resize(VolumeNds+1); 
                    ElementOrder = 1;
                    VolumeShape = 2;

                    for (int m = 1; m <= VolumeNds; m++)
                    {volumes[e][m] = int(read_txt[k+1+m]);}
                }break;
                case 320:
                {
                    e++;
                    volumes.push_back(spacer);
                    VolumeNds = 20;
                    iter = VolumeNds;
                    volumes[e].resize(VolumeNds+1); 
                    ElementOrder = 2;
                    VolumeShape = 2;

                    for (int m = 1; m <= VolumeNds; m++)
                    {volumes[e][m] = int(read_txt[k+1+m]);}
                }break;
                case 304:
                {
                    e++;
                    volumes.push_back(spacer);
                    VolumeNds = 4;
                    iter = VolumeNds;
                    volumes[e].resize(VolumeNds+1); 
                    ElementOrder = 1;
                    VolumeShape = 1;

                    for (int m = 1; m <= VolumeNds; m++)
                    {volumes[e][m] = int(read_txt[k+1+m]);}
                }break;
                case 310:
                {
                    e++;
                    volumes.push_back(spacer);
                    VolumeNds = 10;
                    iter = VolumeNds;
                    volumes[e].resize(VolumeNds+1); 
                    ElementOrder = 2;
                    VolumeShape = 1;

                    for (int m = 1; m <= VolumeNds; m++)
                    {volumes[e][m] = int(read_txt[k+1+m]);}
                }break;
            }
        }
        TotalLines = l;
        TotalSurfaces = f;
        TotalVolumes = e;
        cout << "Total Points = " << TotalPoints << endl;
        cout << "Total Lines = " << TotalLines << endl;
        cout << "Total Surfaces = " << TotalSurfaces << endl;
        cout << "Total Volumes = " << TotalVolumes << endl;
        cout << "Mesh read successfuly..." << endl;
    }
    phi.resize(TotalPoints+1);
    for(int i = 1; i <= TotalPoints; i++)
    {phi[i] = 1;}
    read_txt.clear();
    if (TotalVolumes == 0)
    {dim = 2;}
    else
    {dim = 3;}
    if (TotalSurfaces > 0 || TotalVolumes > 0)
    {computeAvgPntSpacing();}
}
