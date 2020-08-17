#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include<math.h>
#include "MatrixMultiplication.h"
using namespace std;


void BilinearGeneration(double L,double W,int NumL,int NumW,double *xi,double *eta,int **Node,double **Coordinate,double ***J)
{
	cout<<"Element  Node1  Node2  Node3  Node4"<<endl;
	for(int i=0;i<NumL;i++)
	{
		for(int j=0;j<NumW;j++)
		{
			Node[i*NumW+j][0]=i*(NumL+1)+j;
			Node[i*NumW+j][1]=i*(NumL+1)+j+1;
			Node[i*NumW+j][2]=(i+1)*(NumL+1)+j+1;
			Node[i*NumW+j][3]=(i+1)*(NumL+1)+j;

			cout<<i*NumW+j<<"    "<<Node[i*NumL+j][0]<<"   "<<Node[i*NumL+j][1]<<"   "<<Node[i*NumL+j][2]<<"   "<<Node[i*NumL+j][3]<<endl;
		}
	}

	cout<<"Node   CoordinateX   CoordinateY"<<endl;
	for(int i=0;i<NumL+1;i++)
	{
		for(int j=0;j<NumW+1;j++)
		{
			Coordinate[i*(NumL+1)+j][0]=j*W/NumW;
			Coordinate[i*(NumL+1)+j][1]=i*L/NumL;

			cout<<i*(NumW+1)+j<<"    "<<Coordinate[i*(NumL+1)+j][0]<<"    "<<Coordinate[i*(NumL+1)+j][1]<<endl;
		}
	}

	cout<<NumW*NumL<<" elements are seperated."<<endl;

	
	
			for(int m=0;m<2;m++)
			{
				for(int n=0;n<2;n++)
				{
					J[2*m+n][0][0]=((eta[n]-1)*Coordinate[Node[0][0]][0]+(1-eta[n])*Coordinate[Node[0][1]][0]+(1+eta[n])*Coordinate[Node[0][2]][0]-(1+eta[n])*Coordinate[Node[0][3]][0])/4;
					J[2*m+n][0][1]=((eta[n]-1)*Coordinate[Node[0][0]][1]+(1-eta[n])*Coordinate[Node[0][1]][1]+(1+eta[n])*Coordinate[Node[0][2]][1]-(1+eta[n])*Coordinate[Node[0][3]][1])/4;
					J[2*m+n][1][0]=((xi[m]-1)*Coordinate[Node[0][0]][0]-(1+xi[m])*Coordinate[Node[0][1]][0]+(1+xi[m])*Coordinate[Node[0][2]][0]+(1-xi[m])*Coordinate[Node[0][3]][0])/4;
					J[2*m+n][1][1]=((xi[m]-1)*Coordinate[Node[0][0]][1]-(1+xi[m])*Coordinate[Node[0][1]][1]+(1+xi[m])*Coordinate[Node[0][2]][1]+(1-xi[m])*Coordinate[Node[0][3]][1])/4;
				}
			}		

		for(int j=0;j<4;j++)
		{
			cout<<endl;
			cout<<"Jacobian Matrix for Integrate point"<<j<<endl;
			for(int m=0;m<2;m++)
			{
				for(int n=0;n<2;n++)
				{
					cout<<J[j][m][n]<<" ";
				}
				cout<<endl;
			}
		}	
}

void B_Matrix(int NumW,int NumL,double ***A,double ***G,double ***B,double ***J,double *xi,double *eta)
{
	double detJ;

	for(int m=0;m<2;m++)
	{
		for(int n=0;n<2;n++)
		{
			G[2*m+n][0][0]=G[2*m+n][1][1]=(eta[m]-1)/4;
			G[2*m+n][0][2]=G[2*m+n][1][3]=(1-eta[m])/4;
			G[2*m+n][0][4]=G[2*m+n][1][5]=(eta[m]+1)/4;
			G[2*m+n][0][6]=G[2*m+n][1][7]=-(eta[m]+1)/4;
			G[2*m+n][2][0]=G[2*m+n][3][1]=(xi[n]-1)/4;
			G[2*m+n][2][2]=G[2*m+n][3][3]=(-xi[n]-1)/4;
			G[2*m+n][2][4]=G[2*m+n][3][5]=(xi[n]+1)/4;
			G[2*m+n][2][6]=G[2*m+n][3][7]=(1-xi[n])/4;

		}
	}
	G[0][0][1]=G[0][0][3]=G[0][0][5]=G[0][0][7]=G[0][1][0]=G[0][1][2]=G[0][1][4]=G[0][1][6]=0;
	G[0][2][1]=G[0][2][3]=G[0][2][5]=G[0][2][7]=G[0][3][0]=G[0][3][2]=G[0][3][4]=G[0][3][6]=0;

		for(int j=0;j<4;j++)
		{
			cout<<endl;
			cout<<"G Matrix for Integrate point"<<j<<endl;
			for(int m=0;m<4;m++)
			{
				for(int n=0;n<8;n++)
				{
					cout<<G[j][m][n]<<" ";
				}
				cout<<endl;
			}
		}	


		for(int j=0;j<4;j++)
		{
			detJ=J[j][1][1]*J[j][0][0]-J[j][0][1]*J[j][1][0];

			A[j][0][0]=J[j][1][1]/detJ;
			A[j][0][2]=-J[j][0][1]/detJ;
			A[j][1][1]=-J[j][1][0]/detJ;
			A[j][1][3]=J[j][0][0]/detJ;
			A[j][2][0]=-J[j][1][0]/detJ;
			A[j][2][1]=J[j][1][1]/detJ;
			A[j][2][2]=J[j][0][0]/detJ;
			A[j][2][3]=-J[j][0][1]/detJ;

			MatrixMultiply(A[j],3,4,G[j],4,8,B[j]);
		}

	

		for(int j=0;j<4;j++)
		{
			cout<<endl;
			cout<<"B Matrix for Integrate point"<<j<<endl;
			for(int m=0;m<3;m++)
			{
				for(int n=0;n<8;n++)
				{
					cout<<B[j][m][n]<<" ";
				}
				cout<<endl;
			}
		}	
	}

void ElementalStiffnessMatrix(int NumL,int NumW,double **D,double ***B,double ***J,double t,double ***Ke)
{
	double ***Temp,***Temp1;

	Temp= (double***)malloc(4*sizeof(double**));
	for(int i=0;i<4;i++)
	{
		*(Temp+i) = (double**)malloc(3*sizeof(double*));
		for(int j=0;j<3;j++)
		{
			*(*(Temp+i)+j) = (double*)malloc(8*sizeof(double));
			for(int k=0;k<8;k++)
			{
				Temp[i][j][k] = 0.0;
			}
		}
	}
	Temp1= (double***)malloc(4*sizeof(double**));
	for(int i=0;i<4;i++)
	{
		*(Temp1+i) = (double**)malloc(8*sizeof(double*));
		for(int j=0;j<8;j++)
		{
			*(*(Temp1+i)+j) = (double*)malloc(8*sizeof(double));
			for(int k=0;k<8;k++)
			{
				Temp1[i][j][k] = 0.0;
			}
		}
	}

	double detJ,sum=0.0;
	for(int i=0;i<NumL*NumW;i++)
	{
		for(int j=0;j<4;j++)
		{
			
			MatrixMultiply(D,3,3,B[j],3,8,Temp[j]);
			TransposeMatrixMultiply(B[j],3,8,Temp[j],3,8,Temp1[j]);		
		}
	}


		/*for(int j=0;j<4;j++)
		{
			cout<<endl;
			cout<<"D*B Matrix for Integrate point"<<j<<endl;
			for(int m=0;m<3;m++)
			{
				for(int n=0;n<8;n++)
				{
					cout<<Temp[j][m][n]<<" ";
				}
				cout<<endl;
			}
		}	*/

	for(int i=0;i<NumL*NumW;i++)
	{
		for(int m=0;m<8;m++)
		{
			for(int n=0;n<8;n++)
			{
				for(int j=0;j<4;j++)
				{
					detJ=J[j][1][1]*J[j][0][0]-J[j][0][1]*J[j][1][0];
					sum+=t*detJ*Temp1[j][m][n];
					detJ=0.0;
				}
					Ke[i][m][n]=sum;
					sum=0.0;
			}
		}
			
	}

	cout<<"Elemental stiffness Matrix"<<endl;
	for(int i=0;i<NumL*NumW;i++)
	{
		cout<<"Element "<<i<<endl;
		for(int m=0;m<8;m++)
		{
			for(int n=0;n<8;n++)
			{
				cout<<Ke[i][m][n]<<" ";
			}
			cout<<endl;
		}
			
	}
}

void GlobalStiffnessMatrix(double ***Ke,int NumW,int NumL,double **Kg,int **Nodes)
{
	for(int k=0;k<NumL*NumW;k++)
	{
		for(int i=0;i<4;i++)
		{
			for(int j=0;j<4;j++)
			{
				for(int m=0;m<2;m++)
				{
					for(int n=0;n<2;n++)
					{
						Kg[2*Nodes[k][i]+m][2*Nodes[k][j]+n]+=Ke[k][2*i+m][2*j+n];
					}
				}
				
			}
		}
	}
	cout<<"The Global stiffness matrix is"<<endl;
	for(int i=0;i<2*(NumL+1)*(NumW+1);i++)
	{
		for(int j=0;j<2*(NumL+1)*(NumW+1);j++)
		{
			cout<<Kg[i][j]<<" ";
		}
		cout<<endl;
	}
}

void Reduced_Matrix(double **Kg,double *F,double *Disp,double u,int NumL,int NumW)
{
	double **Temp;
	Temp=(double**)malloc(2*NumL*(NumW+1)*sizeof(double*));
	for(int i=0;i<2*NumL*(NumW+1);i++)
	{
		*(Temp+i)=(double*)malloc(2*(NumL+1)*(NumW+1)*sizeof(double));
		for(int j=0;j<2*(NumL+1)*(NumW+1);j++)
		{
			Temp[i][j]=0.0;
		}
	}

	double **Kg_Reduced;
	double *F_Reduced,*F_Temp,*Disp_Reduced,*Disp_Temp;
	F_Temp=(double*)malloc(2*NumL*(NumW+1)*sizeof(double));
	Disp_Temp=(double*)malloc(2*NumL*(NumW+1)*sizeof(double));
	for(int i=0;i<2*NumL*(NumW+1);i++)
	{
		F_Temp[i]=0.0;
		Disp_Temp[i]=0.0;
	}
	F_Reduced=(double*)malloc((2*NumL*(NumW+1)-1)*sizeof(double));
	Disp_Reduced=(double*)malloc((2*NumL*(NumW+1)-1)*sizeof(double));
	Kg_Reduced=(double**)malloc((2*NumL*(NumW+1)-1)*sizeof(double*));
	for(int i=0;i<(2*NumL*(NumW+1)-1);i++)
	{
		F_Reduced[i]=0.0;
		Disp_Reduced[i]=0.0;
		*(Kg_Reduced+i)=(double*)malloc((2*NumL*(NumW+1)-1)*sizeof(double));
		for(int j=0;j<(2*NumL*(NumW+1)-1);j++)
		{
			Kg_Reduced[i][j]=0.0;
		}
	}
	
	for(int k=0;k<NumW+1;k++)
	{
		for(int i=0;i<2*NumL;i++)
		{
			if(i<2*NumL-1)
			{
				for(int j=0;j<2*(NumL+1)*(NumW+1);j++)	
				{
					Temp[k*2*NumL+i][j]=Kg[k*2*NumL+2*k+1+i][j];//Remove Rows with unknown forces
				}
					Disp_Reduced[k*2*NumL+i]=Disp[k*2*NumL+2*k+1+i];
			}
			else
			{
					for(int j=0;j<2*(NumL+1)*(NumW+1);j++)	
				{
					Temp[k*2*NumL+i][j]=Kg[k*2*NumL+2*k+2+i][j];
				}
			}
		}
	}
	
	double sum=0.0;
	for(int i=0;i<2*NumL*(NumW+1);i++)
	{
		for(int k=0;k<NumW+1;k++)
		{
			sum+=-u*(Temp[i][2*(k+1)*NumL+2*k]-Temp[i][2*k*(NumL+1)]);
		}
		F_Temp[i]=sum;
		sum=0.0;
	}

	for(int i=0;i<(2*NumL*(NumW+1)-1);i++)
	{

				F_Reduced[i]=F_Temp[i+1];
	}

	/*cout<<"The Reduced Force is"<<endl;
	for(int i=0;i<(2*NumL*(NumW+1)-1);i++)
	{
			cout<<F_Reduced[i];
		cout<<endl;
	}*/

	for(int k=0;k<NumW+1;k++)
	{
		for(int j=0;j<2*NumL;j++)
		{
			if(j<2*NumL-1)
			{
				for(int i=0;i<2*NumL*(NumW+1);i++)	
				{
					Temp[i][k*2*NumL+j]=Temp[i][k*2*NumL+2*k+1+j];//Remove Columns with unknown forces
				}
			}
			else
			{
					for(int i=0;i<2*NumL*(NumW+1);i++)	
				{
					Temp[i][k*2*NumL+j]=Temp[i][k*2*NumL+2*k+2+j];
				}
			}
		}
	}
	

	/*cout<<"Temp stiffness matrix is"<<endl;
	for(int i=0;i<2*NumL*(NumW+1);i++)
	{
		for(int j=0;j<2*NumL*(NumW+1);j++)
		{
			cout<<Temp[i][j]<<" ";
		}
		cout<<endl;
	}
	cout<<endl;*/

	for(int i=0;i<(2*NumL*(NumW+1)-1);i++)
	{
		for(int j=0;j<(2*NumL*(NumW+1)-1);j++)
		{
				Kg_Reduced[i][j]=Temp[i+1][j+1];
		}
	}
	/*cout<<"The Reduced stiffness matrix is"<<endl;
	for(int i=0;i<(2*NumL*(NumW+1)-1);i++)
	{
		for(int j=0;j<(2*NumL*(NumW+1)-1);j++)
		{
			cout<<Kg_Reduced[i][j]<<" ";
		}
		cout<<endl;
	}*/


	Gauss(Kg_Reduced,Disp_Reduced,F_Reduced,(2*NumL*(NumW+1)-1));
	
	for(int i=0;i<2*NumL*(NumW+1)-1;i++)
	{
		Disp_Temp[i+1]=Disp_Reduced[i];
	}

	for(int i=0;i<NumW+1;i++)
	{
		for(int j=0;j<2*NumL;j++)
		{
			if(j<2*NumL-1)
			{
				Disp[i*2*NumL+2*i+1+j]=Disp_Temp[i*2*NumL+j];
			}
			else
			{
				Disp[i*2*NumL+2*i+2+j]=Disp_Temp[i*2*NumL+j];
			}
		}
	}


	cout<<"Node    Ux      Uy"<<endl;
	for(int i=0;i<(NumL+1)*(NumW+1);i++)
	{
		cout<<i<<"  "<<Disp[2*i]<<"  "<<Disp[2*i+1]<<endl;
	}
}

void NodeForce(double *Force,double **Kg,double *Disp,int NumL,int NumW)
{
	for(int i=0;i<2*(NumW+1)*(NumL+1);i++)
		{
			double sum=0.0;
			for(int j=0;j<2*(NumW+1)*(NumL+1);j++)
			{
				sum+=Kg[i][j]*Disp[j];
			}
			Force[i]=sum;
		}

	cout<<"Node          Fx            Fy"<<endl;
	for(int i=0;i<(NumW+1)*(NumL+1);i++)
	{
		cout<<i<<" "<<Force[2*i]<<" "<<Force[2*i+1]<<endl;
	}
}