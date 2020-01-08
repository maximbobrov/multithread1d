

//just a test comment
#ifndef GLOBALS_H
#define GLOBALS_H

 //for multigrid
 #define N_X 81
 #define N_Y 81



 #define U_TAU 1.0


 #define W_WIDTH 800
 #define W_HEIGHT 800

 #define RES2_MIN 0.000001

#define VIEW_DIV0 7
#define VIEW_DIV  8
#define VIEW_PHI  9
#define VIEW_RHO  1
#define VIEW_RHO_2  13
#define VIEW_JX 2
#define VIEW_JY 3
#define VIEW_UX 4
#define VIEW_UY 5
#define VIEW_fieldUX 10
#define VIEW_fieldUY 11
#define VIEW_P 12


typedef struct
{
    double a,bm,bp,cm,cp;

    int w_bc_type,e_bc_type,n_bc_type,s_bc_type; //0--fixed value, 1--fixed gradient,2 --cyclic; 3 --init
    int w_bc_val,e_bc_val,n_bc_val,s_bc_val;

}INPUT_PARAM;

typedef struct
{
    int order;

    double C[8]; //C0*x + C1*x^2 + C2*x^3..

}poly;


void get_div(double f_x[N_X][N_Y],double f_y[N_X][N_Y],double out[N_X][N_Y]);

extern double RHS[N_X][N_Y],RHS_p[N_X][N_Y],E_x[N_X][N_Y],E_y[N_X][N_Y];


//multigr--------------------------------------------



extern double div_[N_X][N_Y],div_J[N_X][N_Y],rho_[N_X][N_Y],phi_[N_X][N_Y],Ex[N_X][N_Y],Ey[N_X][N_Y], p_in[N_X][N_Y],rho_in[N_X][N_Y];

//extern double flow_1[N_X][N_Y];
//extern double flow_2[N_X][N_Y];
extern double n_1[N_X][N_Y],n_1_prev[N_X][N_Y];
extern double n_2[N_X][N_Y];

extern double Jx_[N_X][N_Y],Jy_[N_X][N_Y];
extern double Ux_[N_X][N_Y];
extern double Uy_[N_X][N_Y];
extern double out_Ux[N_X][N_Y];
extern double out_Uy[N_X][N_Y];
extern double field_Ux[N_X][N_Y];
extern double field_Uy[N_X][N_Y];
extern double arr_X[N_X][N_Y];
extern double arr_Y[N_X][N_Y];
extern double p_in_x[N_X][N_Y];
extern double p_in_y[N_X][N_Y];

extern  double dx;
extern  double dy;
extern double U;


extern int num_thread;
extern double Energy[N_X][N_Y];
extern double BoundaryLayer[N_X];
extern  double dt;
extern double D;
extern double b;
extern double nu;
extern double rho;
extern double alpha;
extern double A;
//extern int t;
extern double m_e;
extern double m_ip;
extern double mu_1; // mobility
extern double mu_2;
extern double q_e;
extern double eps_0;

extern double OMEGA;

extern  int itn;
extern  int clear_w;


extern  float rx;
extern  float ry;
extern  int mx0,my0;
extern  int rotate;
extern  float rx0;
extern  float ry0;
extern double d_temp;
extern double mouse_x,mouse_y;

extern double r_2g;//cg res

extern float* body_x, *body_y, *body_z;
#endif
