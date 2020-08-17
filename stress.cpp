#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include<math.h>
#include "MatrixMultiplication.h"
using namespace std;



void StressStrain(double ***strain,double ***stress,double ***B,double *Disp,int **Nodes,int NumW,int NumL,double **D)
{
	double *temp1;
	temp1=(double*)malloc(8*sizeof(double)); //´æ´¢µ¥ÔªÎ»ÒÆ
	for(int i=0;i<8;i++)
	{
		temp1[i]=0.0;
	}

	for(int i=0;i<NumW*NumL;i++)
	{
		for(int k=0;k<4;k++)
		{
			temp1[2*k]=Disp[2*Nodes[i][k]];
			temp1[2*k+1]=Disp[2*Nodes[i][k]+1];
		}
		for(int k=0;k<4;k++)
		{
			for(int m=0;m<3;m++)
			{
				double sum=0.0;
				for(int n=0;n<8;n++)
				{
					sum+=B[k][m][n]*temp1[n];
				}
				strain[i][k][m]=sum;
			}
		}
	}	

	for(int i=0;i<NumW*NumL;i++)
	{		
		for(int k=0;k<4;k++)
		{
			for(int m=0;m<3;m++)
			{
				double sum1=0.0;
				for(int n=0;n<3;n++)
				{
					sum1+=D[m][n]*strain[i][k][n];
				}
				stress[i][k][m]=sum1;
			}
		}
	}	


	cout<<"Ele IntegrateNumber   ex    ey   gxy"<<endl;
	for(int i=0;i<NumW*NumL;i++)
	{
		for(int j=0;j<4;j++)
		{
			cout<<i<<"       "<<j<<"       ";
			cout<<strain[i][j][0]<<"  "<<strain[i][j][1]<<"  "<<strain[i][j][2]<<endl;
		}
		
	}
	cout<<endl;
	cout<<"Ele IntegrateNumber   sx    sy   sxy"<<endl;
	for(int i=0;i<NumW*NumL;i++)
	{
		for(int j=0;j<4;j++)
		{
			cout<<i<<"       "<<j<<"       ";
			cout<<stress[i][j][0]<<"  "<<stress[i][j][1]<<"  "<<stress[i][j][2]<<endl;
		}
		
	}

}

void VonMisesStress(double **Von_Mises,double ***stress,int NumW,int NumL)
{
	for(int i=0;i<NumW*NumL;i++)
	{
		for(int j=0;j<4;j++)
		{
			Von_Mises[i][j]=sqrt(((stress[i][j][0]-stress[i][j][1])*(stress[i][j][0]-stress[i][j][1])+stress[i][j][0]*stress[i][j][0]+stress[i][j][1]*stress[i][j][1]+6*stress[i][j][2]*stress[i][j][2])/2);
		}
	}


	cout<<"Ele IntegrateNumber   Von_Mises Stress"<<endl;
	for(int i=0;i<NumW*NumL;i++)
	{
		for(int j=0;j<4;j++)
		{
			cout<<i<<"       "<<j<<"       ";
			cout<<Von_Mises[i][j]<<endl;
		}
		
	}
}
