#include "../pntwrks.h"

POINTSET::POINTSET(POINT p, DOMAIN l, RESOLUTION divss, string elemtype)
{
    cout << "Generating pointset with grid..." << endl;
    if (divss.z == 1)
    {
        if (divss.y > 1)
        {
            dim = 2;
            generateGridPnts2D(p, l, divss);
            if (elemtype == "T3")
            {generateGridTri(divss);}
            else if (elemtype == "Q4")
            {generateGridQuad(divss);}
            else
            {cout << "ERROR: element type not valid. Must be T3 or Q4 for 2D" << endl; exit(0);}
        }
        else if(divss.y == 1)
        {
            dim = 1;
            generateGridPnts1D(p, l, divss);
            TotalLines = 0;
            TotalSurfaces = 0;
            TotalVolumes = 0;
        }
    }
    else if (divss.z > 1)
    {
        dim = 3;
        generateGridPnts3D(p,l,divss);
        if (elemtype == "T4")
        {generateGridTet(p, l, divss);}
        else if (elemtype == "H8")
        {generateGridHex(p, l, divss);}
        else
        {cout << "ERROR: element type not valid. Must be T4 or H8 for 3D" << endl; exit(0);}
    }
    computeAvgPntSpacing();
    cout << "Dim = " << dim << endl;
    cout << "Total Points = " << TotalPoints << endl;
    cout << "Total Lines = " << TotalLines << endl;
    cout << "Total Surfaces = " << TotalSurfaces << endl;
    cout << "Total Volumes = " << TotalVolumes << endl;
    cout << "Pointset generated successfuly..." << endl;
}

POINTSET::POINTSET(POINT p, DOMAIN l, RESOLUTION divss)
{
    cout << "Generating pointset only..." << endl;
    if (divss.z == 1)
    {
        if (divss.y > 1)
        {
            dim = 2;
            generateGridPnts2D(p, l, divss);
        }
        else if(divss.y == 1)
        {
            dim = 1;
            generateGridPnts1D(p, l, divss);
        }
    }
    else if (divss.z > 1)
    {
        dim = 3;
        generateGridPnts3D(p,l,divss);
    }
    TotalLines = 0;
    TotalSurfaces = 0;
    TotalVolumes = 0;
    computeAvgPntSpacing();
    cout << "Dim = " << dim << endl;
    cout << "Total Points = " << TotalPoints << endl;
    cout << "Total Lines = " << TotalLines << endl;
    cout << "Total Surfaces = " << TotalSurfaces << endl;
    cout << "Total Volumes = " << TotalVolumes << endl;
    cout << "Pointset generated successfuly..." << endl;
}

//------------------------------------------------------------------
//FUNCTIONS RELATED TO GENERATING GRID POINTS
//------------------------------------------------------------------
void POINTSET::generateGridPnts1D(POINT p, DOMAIN s, RESOLUTION divs)
{
    points.resize(divs.x+1);
    double spacing_x = s.x/double(divs.x-1);
    points[1].x = p.x;
    points[1].y = 0;
    points[1].z = 0;
    int l = 1;
    for (int i = 1; i <= divs.x; i++)
    {
        points[l+1].x = points[l].x + spacing_x;
        points[l+1].y = 0;
        points[l+1].z = 0;
        l = l+1;
    }
    TotalPoints = l-1;
    phi.resize(TotalPoints+1);
    for (int i = 1; i <= TotalPoints; i++)
    {phi[i] = 1;}
}
void POINTSET::generateGridPnts2D(POINT p1, DOMAIN s, RESOLUTION divs)
{
    points.resize((divs.x+1)*(divs.y+1));
    double spacing_x = s.x/double(divs.x-1);
    double spacing_y = s.y/double(divs.y-1);
    points[1].x = p1.x;
    points[1].y = p1.y;
    points[1].z = 0;
    int l = 1;
    for (int j = 1; j <= divs.y; j++)
    {
        for (int i = 1; i <= divs.x; i++)
        {
            points[l+1].x = points[l].x + spacing_x;
            points[l+1].y = points[l].y;
            points[l+1].z = 0;
            l = l+1;
        }
        points[l].x = points[1].x;
        points[l].y = points[l].y + spacing_y;
        points[l].z = 0;
    }
    TotalPoints = l-1;
    phi.resize(TotalPoints+1);
    for(int i = 1; i <= TotalPoints; i++)
    {phi[i] = 1;}
}
void POINTSET::generateGridPnts3D(POINT p1, DOMAIN s, RESOLUTION divs)
{
    points.resize((divs.x+1)*(divs.y+1)*(divs.z+1));
    double spacing_x = s.x/double(divs.x-1);
    double spacing_y = s.y/double(divs.y-1);
    double spacing_z = s.z/double(divs.z-1);
    points[1].x = p1.x;
    points[1].y = p1.y;
    points[1].z = p1.z;
    int l = 1;
    for (int k = 1; k <= divs.z; k++)
    {
        for (int j = 1; j <= divs.y; j++)
        {
            for (int i = 1; i <= divs.x; i++)
            {
                points[l+1].x = points[l].x + spacing_x;
                points[l+1].y = points[l].y;
                points[l+1].z = points[l].z;
                l++;
            }
            points[l].x = points[1].x;
            points[l].y = points[l].y + spacing_y;
            points[l].z = points[l].z;
        }
        points[l].x = p1.x;
        points[l].y = p1.y;
        points[l].z = points[l].z+spacing_z;
    }
    TotalPoints = l-1;
    phi.resize(TotalPoints+1);
    for(int i = 1; i <= TotalPoints; i++)
    {phi[i] = 1;}
    cout << "Total points constructed = " << TotalPoints << endl;
}
//------------------------------------------------------------------
//FUNCTIONS RELATED TO GENERATING GRID ELEMENTS
//------------------------------------------------------------------
void POINTSET::generateGridQuad(RESOLUTION divs)
{
    int elems_x = divs.x-1;
    int elems_y = divs.y-1;
    dim = 2;
    ElementOrder = 1;
    SurfaceShape = 2;
    SurfaceNds = 4;

    int elem_num = elems_x * elems_y;
    surfaces.resize(elem_num+1);
    int node_num = (elems_x+1)* (elems_y+1);
    int node_n[node_num+1];
    for (int i = 1; i <= node_num; i++)
    {node_n[i] = i;}

    int k = 1;
    int l = 1;
    for (int j = 1; j<= elems_y; j++)
    {
        for (int i = 1; i <= elems_x; i++)
        {
            surfaces[l].resize(SurfaceNds+1);
            surfaces[l][1] = node_n[k];
            surfaces[l][2] = node_n[k+1];
            surfaces[l][3] = node_n[k+elems_x+2];
            surfaces[l][4] = node_n[k+elems_x+1];
            k = k+1;
            l = l+1;
        }
        k = k+1;
    }
    TotalLines = 0;
    TotalSurfaces = l-1;
    TotalVolumes = 0;
}
void POINTSET::generateGridTri(RESOLUTION divs)
{
    int elems_x = divs.x-1;
    int elems_y = divs.y-1;
    dim = 2;
    ElementOrder = 1;
    SurfaceShape = 1;
    SurfaceNds = 3;

    int elem_num = 2*elems_x * elems_y;
    surfaces.resize(elem_num+1);
    int node_num = (elems_x+1)*(elems_y+1);
    int node_n[node_num+1];
    for (int i = 1; i <= node_num; i++)
    {node_n[i] = i;}

    int k = 1;
    int l = 1;
    for (int j = 1; j<= elems_y; j++)
    {
        for (int i = 1; i <= elems_x; i++)
        {
            if (((mod(i,2) == 0) & (mod(j,2) != 0)) || ((mod(i,2) != 0) & (mod(j,2) == 0)))
            {
				surfaces[l].resize(SurfaceNds+1);
                surfaces[l][1] = node_n[k];
                surfaces[l][2] = node_n[k+1];
                surfaces[l][3] = node_n[k+elems_x+1];
				surfaces[l+1].resize(SurfaceNds+1);
                surfaces[l+1][1] = node_n[k+1];
                surfaces[l+1][2] = node_n[k+elems_x+2];
                surfaces[l+1][3] = node_n[k+elems_x+1];
            }
            else if (((mod(i,2) != 0) & (mod(j,2) != 0)) || ((mod(i,2) == 0) & (mod(j,2)== 0)))
            {
				surfaces[l].resize(SurfaceNds+1);
                surfaces[l][1] = node_n[k];
                surfaces[l][2] = node_n[k+1];
                surfaces[l][3] = node_n[k+elems_x+2];

				surfaces[l+1].resize(SurfaceNds+1);
                surfaces[l+1][1] = node_n[k];
                surfaces[l+1][2] = node_n[k+elems_x+2];
                surfaces[l+1][3] = node_n[k+elems_x+1];
            }
            k = k+1;
            l = l+2;
        }
        k = k+1;
    }
    TotalLines = 0;
    TotalSurfaces = l-1;
    TotalVolumes = 0;
}
void POINTSET::generateGridTri2(RESOLUTION divs)
{
    int elems_x = divs.x-1;
    int elems_y = divs.y-1;
    dim = 2;
    ElementOrder = 1;
    SurfaceShape = 1;
    SurfaceNds = 3;

    int elem_num = 2*elems_x * elems_y;
    surfaces.resize(elem_num+1);
    int node_num = (elems_x+1)*(elems_y+1);
    int node_n[node_num+1];
    for (int i = 1; i <= node_num; i++)
    {node_n[i] = i;}

    int k = 1;
    int l = 1;
    for (int j = 1; j<= elems_y; j++)
    {
        for (int i = 1; i <= elems_x; i++)
        {
			surfaces[l].resize(SurfaceNds+1);
            surfaces[l][1] = node_n[k];
            surfaces[l][2] = node_n[k+1];
            surfaces[l][3] = node_n[k+elems_x+1];
			surfaces[l+1].resize(SurfaceNds+1);
            surfaces[l+1][1] = node_n[k+1];
            surfaces[l+1][2] = node_n[k+elems_x+2];
            surfaces[l+1][3] = node_n[k+elems_x+1];
            l = l+2;
            k = k+1;
        }
        k = k+1;
    }
    TotalLines = 0;
    TotalSurfaces = elem_num;
    TotalVolumes = 0;
}
void POINTSET::generateGridTet(POINT p, DOMAIN size, RESOLUTION divs)
{
    dim = 3;
    ElementOrder = 1;
    VolumeShape = 1;
    SurfaceShape = 1;
    VolumeNds = 4;
    SurfaceNds = 3;
    
    vector<vector<vector<int> > > TotalPoints(divs.x+1, vector<vector<int> >(divs.y+1, vector<int>(divs.z+1)));
    TotalLines = 0;
    TotalVolumes = 5*(divs.x-1)*(divs.y-1)*(divs.z-1);
    volumes.resize(TotalVolumes+1);
    int e = 1;

    int l = 1;
    for (int k = 1; k <= divs.z; k++)
    {
        for (int j = 1; j <= divs.y; j++)
        {
            for (int i = 1; i <= divs.x; i++)
            {
                TotalPoints[i][j][k] = l;
                l++;
            }
        }
    }

    for (int k = 1; k <= divs.z-1; k++)
    {
        for (int j = 1; j <= divs.y-1; j++)
        {
            for (int i = 1; i <= divs.x-1; i++)
            {
                volumes[e].resize(VolumeNds+1);
                volumes[e][1] = TotalPoints[i][j][k];
                volumes[e][2] = TotalPoints[i+1][j][k];
                volumes[e][3] = TotalPoints[i][j+1][k];
                volumes[e][4] = TotalPoints[i][j][k+1];

                volumes[e+1].resize(VolumeNds+1);
                volumes[e+1][1] = TotalPoints[i+1][j][k];
                volumes[e+1][2] = TotalPoints[i+1][j+1][k];
                volumes[e+1][3] = TotalPoints[i][j+1][k];
                volumes[e+1][4] = TotalPoints[i+1][j+1][k+1];

                volumes[e+2].resize(VolumeNds+1);
                volumes[e+2][1] = TotalPoints[i][j][k+1];
                volumes[e+2][2] = TotalPoints[i][j+1][k];
                volumes[e+2][3] = TotalPoints[i][j+1][k+1];
                volumes[e+2][4] = TotalPoints[i+1][j+1][k+1];

                volumes[e+3].resize(VolumeNds+1);
                volumes[e+3][1] = TotalPoints[i][j][k+1];
                volumes[e+3][2] = TotalPoints[i+1][j][k];
                volumes[e+3][3] = TotalPoints[i+1][j+1][k+1];
                volumes[e+3][4] = TotalPoints[i+1][j][k+1];

                volumes[e+4].resize(VolumeNds+1);
                volumes[e+4][1] = TotalPoints[i][j+1][k];
                volumes[e+4][2] = TotalPoints[i+1][j][k];
                volumes[e+4][3] = TotalPoints[i][j][k+1];
                volumes[e+4][4] = TotalPoints[i+1][j+1][k+1];

                e = e+5;
            }
        }
    }
    TotalSurfaces = TotalVolumes*4;
    surfaces.resize(TotalSurfaces+1);
    int s = 1;
    for (int i = 1; i <= TotalVolumes; i++)
    {
        surfaces[s].resize(SurfaceNds+1);
        surfaces[s][1] = volumes[i][1];
        surfaces[s][2] = volumes[i][2];
        surfaces[s][3] = volumes[i][3];

        surfaces[s+1].resize(SurfaceNds+1);
        surfaces[s+1][1] = volumes[i][1];
        surfaces[s+1][2] = volumes[i][2];
        surfaces[s+1][3] = volumes[i][4];

        surfaces[s+2].resize(SurfaceNds+1);
        surfaces[s+2][1] = volumes[i][2];
        surfaces[s+2][2] = volumes[i][3];
        surfaces[s+2][3] = volumes[i][4];

        surfaces[s+3].resize(SurfaceNds+1);
        surfaces[s+3][1] = volumes[i][1];
        surfaces[s+3][2] = volumes[i][3];
        surfaces[s+3][3] = volumes[i][4];
        s = s + 4;
    }
}
void POINTSET::generateGridHex(POINT p, DOMAIN size, RESOLUTION divs)
{
	dim = 3;
    ElementOrder = 1;
    VolumeShape = 2;
    SurfaceShape = 2;
    VolumeNds = 8;
    SurfaceNds = 4;
    
    vector<vector<vector<int> > > TotalPoints(divs.x+1, vector<vector<int> >(divs.y+1, vector<int>(divs.z+1)));
    TotalLines = 0;
    TotalVolumes = (divs.x-1)*(divs.y-1)*(divs.z-1);
    volumes.resize(TotalVolumes+1);
    int e = 1;

    int l = 1;
    for (int k = 1; k <= divs.z; k++)
    {
        for (int j = 1; j <= divs.y; j++)
        {
            for (int i = 1; i <= divs.x; i++)
            {
                TotalPoints[i][j][k] = l;
                l++;
            }
        }
    }

    for (int k = 1; k <= divs.z-1; k++)
    {
        for (int j = 1; j <= divs.y-1; j++)
        {
            for (int i = 1; i <= divs.x-1; i++)
            {
                volumes[e].resize(VolumeNds+1);
                volumes[e][1] = TotalPoints[i][j][k];
                volumes[e][2] = TotalPoints[i+1][j][k];
                volumes[e][3] = TotalPoints[i+1][j+1][k];
                volumes[e][4] = TotalPoints[i][j+1][k];
                volumes[e][5] = TotalPoints[i][j][k+1];
                volumes[e][6] = TotalPoints[i+1][j][k+1];
                volumes[e][7] = TotalPoints[i+1][j+1][k+1];
                volumes[e][8] = TotalPoints[i][j+1][k+1];
                e = e+1;
            }
        }
    }
    TotalSurfaces = TotalVolumes*6;
    surfaces.resize(TotalSurfaces+1);
    int s = 1;
    for (int i = 1; i <= TotalVolumes; i++)
    {
        surfaces[s].resize(SurfaceNds+1);
        surfaces[s][1] = volumes[i][1];
        surfaces[s][2] = volumes[i][2];
        surfaces[s][3] = volumes[i][3];
		surfaces[s][4] = volumes[i][4];

        surfaces[s+1].resize(SurfaceNds+1);
        surfaces[s+1][1] = volumes[i][5];
        surfaces[s+1][2] = volumes[i][6];
        surfaces[s+1][3] = volumes[i][7];
        surfaces[s+1][4] = volumes[i][8];

        surfaces[s+2].resize(SurfaceNds+1);
        surfaces[s+2][1] = volumes[i][2];
        surfaces[s+2][2] = volumes[i][6];
        surfaces[s+2][3] = volumes[i][7];
		surfaces[s+2][4] = volumes[i][3];

        surfaces[s+3].resize(SurfaceNds+1);
        surfaces[s+3][1] = volumes[i][1];
        surfaces[s+3][2] = volumes[i][5];
        surfaces[s+3][3] = volumes[i][8];
        surfaces[s+3][4] = volumes[i][4];

		surfaces[s+4].resize(SurfaceNds+1);
        surfaces[s+4][1] = volumes[i][1];
        surfaces[s+4][2] = volumes[i][2];
        surfaces[s+4][3] = volumes[i][6];
        surfaces[s+4][4] = volumes[i][5];

        surfaces[s+5].resize(SurfaceNds+1);
        surfaces[s+5][1] = volumes[i][4];
        surfaces[s+5][2] = volumes[i][3];
        surfaces[s+5][3] = volumes[i][7];
        surfaces[s+5][4] = volumes[i][8];
        
        s = s + 6;
    }
}
