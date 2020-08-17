#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Mesh.h"
#include"stress.h"
#include <iostream>
using namespace std;

void main()
{
	//Geometry//
	double t=1,L=100,W=100;
	int NumL=10,NumW=10;

	//Lamina Property//
	double E1=164300,E2=8977,v12=0.32,G12=5020,angle=45;
	double v21,Q11,Q12,Q21,Q22,Q66;
	v21=E2*v12/E1;
	Q11=E1/(1-v12*v21);
	Q12=Q21=v12*E2/(1-v12*v21);
	Q22=E2/(1-v12*v21);
	Q66=G12;

	double **Q;
	Q=(double**)malloc(3*sizeof(double*));
	for(int i=0;i<3;i++)
	{
		*(Q+i)=(double*)malloc(3*sizeof(double));
		for(int j=0;j<3;j++)
		{
			Q[i][j]=0.0;
		}
	}
	Q[0][0]=Q11*pow(cos(angle*3.1415/180),4)+Q22*pow(sin(angle*3.1415/180),4)+2*(Q12+2*Q66)*pow(cos(angle*3.1415/180),2)*pow(sin(angle*3.1415/180),2);
	Q[1][0]=Q[0][1]=(Q11+Q22-4*Q66)*pow(cos(angle*3.1415/180),2)*pow(sin(angle*3.1415/180),2)+Q12*(pow(cos(angle*3.1415/180),4)+pow(sin(angle*3.1415/180),4));
	Q[1][1]=Q11*pow(sin(angle*3.1415/180),4)+Q22*pow(cos(angle*3.1415/180),4)+2*(Q12+2*Q66)*pow(cos(angle*3.1415/180),2)*pow(sin(angle*3.1415/180),2);
	Q[0][2]=Q[2][0]=(Q11-Q12-2*Q66)*pow(cos(angle*3.1415/180),3)*sin(angle*3.1415/180)-(Q22-Q12-2*Q66)*pow(sin(angle*3.1415/180),3)*cos(angle*3.1415/180);
	Q[1][2]=Q[2][1]=(Q11-Q12-2*Q66)*pow(sin(angle*3.1415/180),3)*cos(angle*3.1415/180)-(Q22-Q12-2*Q66)*pow(cos(angle*3.1415/180),3)*sin(angle*3.1415/180);
	Q[2][2]=(Q11+Q22-2*Q12-2*Q66)*pow(cos(angle*3.1415/180),2)*pow(sin(angle*3.1415/180),2)+Q66*(pow(cos(angle*3.1415/180),4)+pow(sin(angle*3.1415/180),4));
	
	cout<<"Q Matrix for Lamina with fiber orientation "<<angle<<"бу (MPa)"<<endl;
	for(int i=0;i<3;i++)
	{
		for(int j=0;j<3;j++)
		{
			cout<<Q[i][j]<<" ";
		}
		cout<<endl;
	}
	cout<<endl;

	//Mesh and Stiffness Matrix Generation//
	double ***J,***B,***A,***G,***Ke,**Coordinate;
	double **Kg;
	int **Nodes;

	Coordinate=(double**)malloc((NumL+1)*(NumW+1)*sizeof(double*));
	for(int i=0;i<(NumL+1)*(NumW+1);i++)
	{
		*(Coordinate+i)=(double*)malloc(2*sizeof(double));
		for(int j=0;j<2;j++)
		{
			Coordinate[i][j]=0.0;
		}
	}

	Nodes=(int**)malloc(2*NumL*NumW*sizeof(int*));
	for(int i=0;i<2*NumL*NumW;i++)
	{
		*(Nodes+i)=(int*)malloc(4*sizeof(int));
		for(int j=0;j<4;j++)
			{
				Nodes[i][j]=0;
			}
	}

	double *eta,*xi;
	eta=(double*)malloc(2*sizeof(double));
	xi=(double*)malloc(2*sizeof(double));
	for(int i=0;i<2;i++)
	{
		eta[i]=0.0;
		xi[i]=0.0;
	}

	eta[0]=-1/sqrt(3.0),eta[1]=1/sqrt(3.0);
	xi[0]=-1/sqrt(3.0),xi[1]=1/sqrt(3.0);

	J=(double***)malloc(4*sizeof(double**));
	for(int i=0;i<4;i++)
	{
		*(J+i)=(double**)malloc(2*sizeof(double*));
		for(int j=0;j<2;j++)
		{
			*(*(J+i)+j)=(double*)malloc(2*sizeof(double));
				for(int n=0;n<2;n++)
				{
					J[i][j][n]=0.0;
				}
		}
	}

	A=(double***)malloc(4*sizeof(double**));
	for(int i=0;i<4;i++)
	{
		*(A+i)=(double**)malloc(3*sizeof(double*));
		for(int j=0;j<3;j++)
		{
			*(*(A+i)+j)=(double*)malloc(6*sizeof(double));
				for(int n=0;n<6;n++)
				{
					A[i][j][n]=0.0;
				}
		}
	}

	G=(double***)malloc(4*sizeof(double**));
	for(int i=0;i<NumW*NumL;i++)
	{
		*(G+i)=(double**)malloc(4*sizeof(double*));
		for(int j=0;j<4;j++)
		{
			*(*(G+i)+j)=(double*)malloc(8*sizeof(double));
			for(int m=0;m<8;m++)
			{
					G[i][j][m]=0.0;
			}
		}
	}

	B=(double***)malloc(4*sizeof(double**));
	for(int i=0;i<4;i++)
	{
		*(B+i)=(double**)malloc(3*sizeof(double*));
		for(int j=0;j<3;j++)
		{
			*(*(B+i)+j)=(double*)malloc(8*sizeof(double));
				for(int n=0;n<8;n++)
				{
					B[i][j][n]=0.0;
				}
		}
	}

	Ke = (double***)malloc(NumW*NumL*sizeof(double**));
	for(int i=0;i<NumW*NumL;i++)
	{
		*(Ke+i) = (double**)malloc(8*sizeof(double*));
		for(int j=0;j<8;j++)
		{
			*(*(Ke+i)+j) = (double*)malloc(8*sizeof(double));
			for(int k=0;k<8;k++)
			{
				Ke[i][j][k] = 0.0;
			}
		}
	}

	Kg=(double**)malloc(2*(NumL+1)*(NumW+1)*sizeof(double*));
	for(int i=0;i<2*(NumL+1)*(NumW+1);i++)
	{
		*(Kg+i)=(double*)malloc(2*(NumL+1)*(NumW+1)*sizeof(double));
		for(int j=0;j<2*(NumL+1)*(NumW+1);j++)
		{
			Kg[i][j]=0.0;
		}
	}

	BilinearGeneration(L,W,NumL,NumW,xi,eta,Nodes,Coordinate,J);
	B_Matrix(NumW,NumL,A,G,B,J,xi,eta);
	ElementalStiffnessMatrix(NumL,NumW,Q,B,J,t,Ke);
	GlobalStiffnessMatrix(Ke,NumW,NumL,Kg,Nodes);

	//Solve Displacement and Force//
	double *Disp,*Force;
	double u=0.5;
	Force=(double*)malloc(2*(NumL+1)*(NumW+1)*sizeof(double));
	Disp=(double*)malloc(2*(NumL+1)*(NumW+1)*sizeof(double));
	for(int i=0;i<2*(NumL+1)*(NumW+1);i++)
	{
		Force[i]=Disp[i]=0.0;
	}
	
	for(int k=0;k<NumW+1;k++)
	{
		Disp[2*(NumL+1)*k]=-u;
		Disp[2*(NumL+1)*k+2*NumL]=u;
	}

	Reduced_Matrix(Kg,Force,Disp,u,NumL,NumW);
	NodeForce(Force,Kg,Disp,NumL,NumW);


	// Stress and Strain//
	double ***strain,***stress;
	double **Von_Mises;
	strain=(double***)malloc(NumL*NumW*sizeof(double**));
	stress=(double***)malloc(NumL*NumW*sizeof(double**));
	for(int i=0;i<NumL*NumW;i++)
	{
		*(strain+i)=(double**)malloc(4*sizeof(double*));
		*(stress+i)=(double**)malloc(4*sizeof(double*));
		for(int j=0;j<4;j++)
		{
			*(*(strain+i)+j)=(double*)malloc(3*sizeof(double));
			*(*(stress+i)+j)=(double*)malloc(3*sizeof(double));
				for(int n=0;n<3;n++)
				{
					strain[i][j][n]=0.0;	
					stress[i][j][n]=0.0;
		
				}
		}
	}
	
	Von_Mises=(double**)malloc(NumL*NumW*sizeof(double*));
	for(int i=0;i<NumL*NumW;i++)
	{
		*(Von_Mises+i)=(double*)malloc(4*sizeof(double));
		for(int j=0;j<4;j++)
		{
			Von_Mises[i][j]=0.0;
		}
	}
	StressStrain(strain,stress,B,Disp,Nodes,NumW,NumL,Q);
	VonMisesStress(Von_Mises,stress,NumW,NumL);
	system("pause");
}