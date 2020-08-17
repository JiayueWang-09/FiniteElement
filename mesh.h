#ifndef _MESH_ 
#define _MESH_ 

void BilinearGeneration(double L,double W,int NumL,int NumW,double *xi,double *eta,int **Node,double **Coordinate,double ***J);

void B_Matrix(int NumW,int NumL,double ***A,double ***G,double ***B,double ***J,double *xi,double *eta);

void ElementalStiffnessMatrix(int NumL,int NumW,double **D,double ***B,double ***J,double t,double ***Ke);

void GlobalStiffnessMatrix(double ***Ke,int NumW,int NumL,double **Kg,int **Nodes);

void Reduced_Matrix(double **Kg,double *F,double *Disp,double u,int NumL,int NumW);

void NodeForce(double *Force,double **Kg,double *Disp,int NumL,int NumW);

#endif