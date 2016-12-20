#include "LibUtilities/LibUtilities.h"

typedef struct{
	double x, y; // coordinate
	int ind; // No. of point
}point2d;

point2d PRI_VERT2d;

// Cross product function
double pt_crossProduct2d(point2d p1, point2d p2, point2d p3) {
	return ((p2.x - p1.x)*(p3.y - p1.y) - (p2.y - p1.y)*(p3.x - p1.x));
}

// distance of points
double pt_distance2d(point2d p1, point2d p2){
	return sqrt((p1.x - p2.x)*(p1.x - p2.x) + (p1.y - p2.y)*(p1.y - p2.y));
}

// Compare function
int pt_cmp2d(const void *p1, const void *p2){
	point2d *p3, *p4;
	double m;
	p3 = (point2d *)p1;
	p4 = (point2d *)p2;
	m = pt_crossProduct2d(PRI_VERT2d, *p3, *p4);
	if (m < 0) return 1;
	else if (m == 0 && (pt_distance2d(PRI_VERT2d, *p3) < pt_distance2d(PRI_VERT2d, *p4)))
		return 1;
	else return -1;
}

/**
 * @brief re-sort the vertex list to be anti-clockwise
 *
 * @param[in] Nvert number of vertex in vertlist
 * @param[in] vx coordinate of all the vertex
 * @param[in] vy coordinate of all the vertex
 * @param[in,out] vertlist array of vertex index
 */
void mr_resortEToV2d(int Nvert, double *vx, double *vy, int *vertlist){
	int i;

	point2d vertex[Nvert];

	for(i=0;i<Nvert;i++){
		int t = vertlist[i];
		vertex[i].x = vx[t];
		vertex[i].y = vy[t];
		vertex[i].ind = vertlist[i];
	}
	PRI_VERT2d.x  = vertex[0].x;
	PRI_VERT2d.y  = vertex[0].y;
	PRI_VERT2d.ind = vertex[0].ind;

    qsort(&vertex[1], Nvert - 1, sizeof(point2d), pt_cmp2d);

	for(i=0;i<Nvert;i++){
		vertlist[i] = vertex[i].ind;
	}
	return;
}