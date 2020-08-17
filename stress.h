#ifndef _STRESS_ 
#define _STRESS_ 

void StressStrain(double ***strain,double ***stress,double ***B,double *Disp,int **Nodes,int NumW,int NumL,double **D);

void VonMisesStress(double **Von_Mises,double ***stress,int NumW,int NumL);

#endif