#include "LibUtilities/LibUtilities.h"
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

typedef struct{
	double x, y; // coordinate
	int nv; // No. of point
}POINT;

POINT stdvert;

// Cross product function
double Multiply(POINT p1, POINT p2, POINT p3) {
	return ((p2.x - p1.x)*(p3.y - p1.y) - (p2.y - p1.y)*(p3.x - p1.x));
}

// Distance of points
double Distance(POINT p1, POINT p2){
	return sqrt((p1.x - p2.x)*(p1.x - p2.x) + (p1.y - p2.y)*(p1.y - p2.y));
}

// Compare function
int cmpPoint(const void *p1, const void *p2){
	POINT *p3, *p4;
	double m;
	p3 = (POINT *)p1;
	p4 = (POINT *)p2;
	m = Multiply(stdvert, *p3, *p4);
	if (m < 0) return 1;
	else if (m == 0 && (Distance(stdvert, *p3) < Distance(stdvert, *p4)))
		return 1;
	else return -1;
}


void VertexSort(int Nvert, double *x, double *y, int *EToV){
	int i;
	POINT vertex[Nvert];

	for(i=0;i<Nvert;i++){
		int t = EToV[i];
		vertex[i].x = x[t];
		vertex[i].y = y[t];
		vertex[i].nv = EToV[i];
	}
	stdvert.x  = vertex[0].x;
	stdvert.y  = vertex[0].y;
	stdvert.nv = vertex[0].nv;

	qsort(&vertex[1], Nvert-1, sizeof(POINT), cmpPoint);

	for(i=0;i<Nvert;i++){
		EToV[i] = vertex[i].nv;
	}
	return;
}