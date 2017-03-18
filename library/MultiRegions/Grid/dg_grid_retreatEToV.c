#include "Utility/utility.h"
#include "dg_grid.h"

/**
 * @brief
 */
typedef struct point2d{
	double x, y; // coordinate
	int ind; // No. of point
}point2d;

static point2d PRI_VERT2d;

/** cross product of two lines */
static double pt_crossProduct2d(point2d p1, point2d p2, point2d p3) {
	return ((p2.x - p1.x)*(p3.y - p1.y) - (p2.y - p1.y)*(p3.x - p1.x));
}

/** distance of 2d points */
static double pt_distance2d(point2d p1, point2d p2){
	return sqrt((p1.x - p2.x)*(p1.x - p2.x) + (p1.y - p2.y)*(p1.y - p2.y));
}

/** sort 2d point function */
static int pt_cmp2d(const void *p1, const void *p2){
	double m;
    point2d *p3 = (point2d *)p1;
    point2d *p4 = (point2d *)p2;
	m = pt_crossProduct2d(PRI_VERT2d, *p3, *p4);
	if (m < 0) return 1;
	else if (m == 0 && (pt_distance2d(PRI_VERT2d, *p3) < pt_distance2d(PRI_VERT2d, *p4)))
		return 1;
	else return -1;
}

/**
 * @brief
 * sort the vertex list in EToV to be anti-clockwise
 * @param[in] Nvert number of vertex in vertlist;
 * @param[in] vx coordinate of all the vertex;
 * @param[in] vy coordinate of all the vertex;
 * @param[in,out] vertlist array of vertex index;
 */
static void dg_grid_resortEToV2d(int Nvert, double *vx, double *vy, int *vertlist){

    register int i;
	point2d vertex[Nvert];

	for(i=0;i<Nvert;i++){
		int t = vertlist[i];
		vertex[i].x = vx[t];
		vertex[i].y = vy[t];
		vertex[i].ind = t;
	}
	PRI_VERT2d.x = vertex[0].x;
	PRI_VERT2d.y = vertex[0].y;
	PRI_VERT2d.ind = vertex[0].ind;

    qsort(&vertex[1], (size_t)(Nvert-1), sizeof(point2d), pt_cmp2d);
	for(i=0;i<Nvert;i++){
		vertlist[i] = vertex[i].ind;
	}
	return;
}
/**
 * @brief
 * @param grid
 */
void dg_grid_retreatEToV2d(dg_grid *grid){
	const int K = grid->K;
	const int Nv = dg_cell_Nv(grid->cell);
	int k;
	for(k=0;k<K;k++){
		dg_grid_resortEToV2d(Nv, grid->vx, grid->vy, grid->EToV[k]);
	}
	return;
}
/**
 * @brief
 * @param grid
 */
void dg_grid_retreatEToV3d(dg_grid *grid){
	return;
}