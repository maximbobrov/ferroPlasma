

//just a test comment
#ifndef GLOBALS_H
#define GLOBALS_H

 //for multigrid
 //#define N_X 1025
 //#define N_Y 513
//#define N_Y_DIEL 255

#define M_PI 3.1415926535

#define N_X 257
#define N_Y 129
#define N_Y_DIEL 64

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

#define w_x0 -0.25e-6
#define w_x1 0.25e-6


#define w_y0 -0.5e-7
#define w_y1 0.5e-7

#define w_z0 -1.0e-7
#define w_z1 1.0e-7
//near_wall_layer width

#define DY_WALL ((w_y1-w_y0)/50)
#define Y_WALL (w_y0+DY_WALL)


 #define U_TAU 1.0

#define qe 1.6e-19
#define Me 9.11e-31
#define eps0 8.85e-12
#define pi4 12.5663706144


 #define W_WIDTH 600
 #define W_HEIGHT 300

 #define RES2_MIN 0.000001

#define VIEW_PHI  0
#define VIEW_E  1
#define VIEW_P 2
#define VIEW_PX 3
#define VIEW_DIV 4



const float inv4PI             = 0.25/M_PI;  // Laplace kernel coefficient

template<typename T>
class vec3 {
public:
  T x;
  T y;
  T z;
};

template<typename T>
class vec4 {
public:
  T x;
  T y;
  T z;
  T w;
};

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



double calc_poly(poly &p, double r, double x); //возвращает значение функции f(x) = x^2-2
double calc_d_poly(poly &p, double x); //возвращает значение производной
double calc_d2_poly(poly &p, double x); // значение второй производной
double solve_poly(poly &p,double _x0, double rhs,int itn);

void rand_init();


float my_rand(int i);
float getVms_from_Ev(float eps_in_ev);

void get_div(double f_x[N_X][N_Y],double f_y[N_X][N_Y],double out[N_X][N_Y]);



double delta_f(double E, double phi);

extern double RHS[N_X][N_Y],RHS_p[N_X][N_Y],E_x[N_X][N_Y],E_y[N_X][N_Y];


//multigr--------------------------------------------

extern double frequency;

extern double tt; //curr time;
extern bool move_particles;
extern double q[N_X],gau[N_Y];


extern double BoundaryLayer[N_X],WallEnergy[N_X];
extern double BoundaryLayerGauss[N_X];
extern int gaussL;

extern double Pins_top[N_X];  //nera electrode pinning
extern double Pins_bottom[N_Y];

extern double div_[N_X][N_Y],div_J[N_X][N_Y],rho_[N_X][N_Y],phi_[N_X][N_Y],Ex[N_X][N_Y],Ey[N_X][N_Y], Ey0[N_X][N_Y],p_in[N_X][N_Y],rho_in[N_X][N_Y];
extern int shift;
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
extern double sasign;

extern double Px_[N_X][N_Y],Px0_[N_X][N_Y]; //polaization
extern double Py_[N_X][N_Y],Py0_[N_X][N_Y];

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
#endif
