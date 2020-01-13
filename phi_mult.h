

#ifndef PHI_MULT_H
#define PHI_MULT_H

#include "globals.h"




double ddx(double var_[N_X][N_Y],int i, int j);
double ddy(double var_[N_X][N_Y],int i, int j);

double iter_phi_gs_bound(int in_);

double field_U(double U, double out_Ux[N_X][N_Y],double itn);
//double field_phi(double A, double phi_[N_X][N_Y],double itn);
//double jacobi(double phi[N_X][N_Y], double rhs[N_X][N_Y], int itn);

double jacobi(INPUT_PARAM par, double field[N_X][N_Y], double rhs[N_X][N_Y], int itn);
double multigrid(INPUT_PARAM par, double field[N_X][N_Y], double rhs[N_X][N_Y], int itn);
double multigrid_N(INPUT_PARAM par, double field[N_X][N_Y], double rhs[N_X][N_Y], int itn, int N);


double solve_current(double rho_in[N_X][N_Y],double Jx_in[N_X][N_Y],double Jy_in[N_X][N_Y],double out_Ux[N_X][N_Y],double out_Uy[N_X][N_Y], int itn);
double solve_current_1(double rho_in[N_X][N_Y],double Jx_in[N_X][N_Y],double Jy_in[N_X][N_Y],double out_Ux[N_X][N_Y],double out_Uy[N_X][N_Y], int itn);
double solve_current_new(double rho_in[N_X][N_Y],double Jx_in[N_X][N_Y],double Jy_in[N_X][N_Y],double out_Ux[N_X][N_Y],double out_Uy[N_X][N_Y], int itn);
double solve_current_old(double rho_in[N_X][N_Y],double Jx_in[N_X][N_Y],double Jy_in[N_X][N_Y],double out_Ux[N_X][N_Y],double out_Uy[N_X][N_Y], int itn);
double solve_current_phi(double rho_in[N_X][N_Y],double Jx_in[N_X][N_Y],double Jy_in[N_X][N_Y],double out_Ux[N_X][N_Y],double out_Uy[N_X][N_Y], int itn);


void get_div(double f_x[N_X][N_Y],double f_y[N_X][N_Y],double out[N_X][N_Y]);
double NavierStoks_solver(double p_in[N_X][N_Y], double out_Ux[N_X][N_Y],double out_Uy[N_X][N_Y], int itn);

double solve_poly(poly p,double _x0, double rhs,int itn);
int jacobi_polynomial(INPUT_PARAM par, poly pol,double field[N_X][N_Y],double rhs[N_X][N_Y], int itn);


double get_nearwall_potential(float x, float y);

#endif
