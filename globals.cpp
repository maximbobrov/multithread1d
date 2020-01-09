
#include "globals.h"





double RHS[N_X][N_Y],RHS_p[N_X][N_Y],E_x[N_X][N_Y],E_y[N_X][N_Y];


//multigr--------------------------------------------

int num_thread = 1;
double Energy[N_X][N_Y];
double BoundaryLayer[N_X];

double div_[N_X][N_Y],div_J[N_X][N_Y],rho_[N_X][N_Y],phi_[N_X][N_Y],Ex[N_X][N_Y],Ey[N_X][N_Y], p_in[N_X][N_Y],rho_in[N_X][N_Y];
//double flow_1[N_X][N_Y];
//double flow_2[N_X][N_Y];

double Px_[N_X][N_Y]; //polaization
double Py_[N_X][N_Y];

double n_1[N_X][N_Y],n_1_prev[N_X][N_Y];
double n_2[N_X][N_Y];
double Jx_[N_X][N_Y],Jy_[N_X][N_Y];
double Ux_[N_X][N_Y];
double Uy_[N_X][N_Y];
double out_Ux[N_X][N_Y];
double out_Uy[N_X][N_Y];
double arr_X[N_X][N_Y];
double arr_Y[N_X][N_Y];
double p_in_x[N_X][N_Y];
double p_in_y[N_X][N_Y];
double Re=30.0;

 float* body_x, *body_y, *body_z;
double OMEGA=0.1;

double dx=2.5e-6/N_X;
double dy=0.5e-6/N_Y;
double dz=0.1e-6;
double U=1.0;

double dt=1e-13;
double D=0.01;
double b= 0.8;
double nu= 17.9*10e-6;
double alpha= 100;
double A=1500;
double rho=1000;
double eps = 1.0;

double m_e=9.1*10e-31;
double m_ip=2.32*10e-23;
double mu_1=4.0*10e-7; // mobility
double mu_2=2.0*10e-9;
double q_e= 1.6*10e-19;
double eps_0=8.85*10e-12;

int itn=0;
int clear_w=1.0;
//int t=0;


float rx=0;
float ry=0;
int mx0,my0;
int rotate=0;
float rx0=0;
float ry0=0;
double d_temp;
double mouse_x,mouse_y;

double r_2g=0.0;//cg res
