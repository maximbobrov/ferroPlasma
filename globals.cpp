
#include "globals.h"
#include<xmmintrin.h>



double RHS[N_X][N_Y],RHS_p[N_X][N_Y],E_x[N_X][N_Y],E_y[N_X][N_Y];

int num_thread=1;
uint32_t xr[32],yr[32],zr[32],wr[32];

bool g_emitElectrons=true;
double E_global=1e7;
double dtKoef = 1.;
double g_t = 0;
double g_phi = 5.;
double g_save_time = 0;
double g_save_time2 = 0;
double g_dphidt;
double get_time(void) {
    struct timeval tv;
    struct timezone tz;
    gettimeofday(&tv, &tz);
    return ((double)(tv.tv_sec+tv.tv_usec*1.0e-6));
}

void rand_init()
{
    for (int i=0;i<num_thread;i++)
    {
        xr[i] = 123456789+i;
        yr[i] = 362436069+2*i;
        zr[i] = 521288629+3*i;
        wr[i] = 88675123+5*i;
    }
}

float my_rand(int i) {
    uint32_t t;
    t = xr[i] ^ (xr[i] << 11);
    xr[i] = yr[i]; yr[i] = zr[i]; zr[i] = wr[i];
    wr[i] = wr[i] ^ (wr[i] >> 19) ^ (t ^ (t >> 8));
    return wr[i]*0.5/0x7FFFFFFF;
}


float getVms_from_Ev(float eps_in_ev)
{
    return sqrt(eps_in_ev)*5.9e5;
}

double delta_f(double E, double phi)
{
    double Emax=2.0;
    double E0=0.2;
    double E1=E-E0;

    if (E1<0) return -1; //just bounce back the electron

    double ks_d=1.0; //roughness in delta
    double ks_w=1.0; //roughness in w
    double deltam=3.0;
    double phi2pi=phi*phi/(2.0*M_PI);

    double w=(E1)/(Emax*(1.0+ks_w*phi2pi)-E0);


    double k= (w < 1) ? 0.56 : 0.25;


    double W= (w <= 3.6) ? pow(w*exp(1-w),k) : 1.125/pow(w,0.35);

    double delta=deltam*(1.0+ks_d*phi2pi)*W;

    return delta;


}


double calc_poly(poly &p, double r, double x) //возвращает значение функции f(x) = x^2-2
{
    double sum=0.0;
    double xn=x;
    for (int i=0;i<p.order;++i)
    {
        sum+=xn*p.C[i];
        xn*=x;
    }

    return sum-r;
}

double calc_d_poly(poly &p, double x) //возвращает значение производной
{
    double sum=p.C[0];
    double xn=x;
    for (int i=0;i<p.order-1;++i)
    {
        sum+=xn*p.C[i+1]*(i+2);
        xn*=x;
    }
    return sum;
}

double calc_d2_poly(poly &p, double x) // значение второй производной
{
    double sum=2.0*p.C[1];
    double xn=x;
    for (int i=0;i<p.order-2;++i)
    {
        sum+=xn*p.C[i+1]*(i+2)*(i+3);
        xn*=x;
    }
    return sum;
}

double solve_poly(poly &p,double _x0, double rhs,int itn)
{
    double x0,x,xn,eps;// вычисляемые приближения для корня
    double a, b,deltax;// границы отрезка и необходимая точность

    // deltax=deltax0;
    x=_x0;
    /* a=x0-deltax;
    b=x0+deltax;

    int nn=0;
    double sucs=f(a)*f(b);
    while ((sucs>0)&&(nn<10)) // если знаки функции на краях отрезка одинаковы
    {
        deltax=deltax*2;
        a=x0-deltax;
        b=x0+deltax;
        nn++;
        sucs=f(a)*f(b);
    }

    if (sucs>0)
    {
        printf("DIVERGENCEE in POLY \n",)
        return -1e30
    }
*/
    for (int i=0; i<itn; ++i)
    {
        x = x-calc_poly(p,rhs,x)/calc_d_poly(p,x); // считаем первое приближение
    }

    // printf("x = %lf f(x)=%e \n",x*1000.0,calc_poly(p,rhs,x));

    return x;

}

//multigr--------------------------------------------
double frequency = 2e7;

double tt=0.0; //curr time

double q[N_X],gau[N_Y]; //wall charge

double BoundaryLayer[N_X],WallEnergy[N_X];

double BoundaryLayerGauss[N_X];
int gaussL = 5;
double scale = 1;
bool move_particles=false;
int shift = 0;
double div_[N_X][N_Y],div_J[N_X][N_Y],rho_[N_X][N_Y],phi_[N_X][N_Y],Ex[N_X][N_Y],Ey[N_X][N_Y],Ey0[N_X][N_Y], p_in[N_X][N_Y],rho_in[N_X][N_Y];
//double flow_1[N_X][N_Y];
//double flow_2[N_X][N_Y];

double Px_[N_X][N_Y],Px0_[N_X][N_Y]; //polaization
double Py_[N_X][N_Y],Py0_[N_X][N_Y];

double Pins_top[N_X];  //nera electrode pinning
double Pins_bottom[N_Y];


double sasign = 1.0;
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

double OMEGA=0.1;

double dx=(w_x1-w_x0)/(N_X-1);//1.0e-4/N_X;
double dy=(w_y1-w_y0)/(N_Y-1);//0.3e-4/N_Y;
double U=1.0;

double dt=1.e-16;
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
