#include "pzsolver.h"


pzSolver::pzSolver()
{
    this->m_p_num=300;
    m_p=new pElem[m_p_num];
    m_dt=1e-11;

    double _dx,_dz;
    _dz=w_z1-w_z0;
    _dx=(w_z1-w_z0)/(m_p_num-1);
    m_dx=_dx;

    for (int i=0;i<m_p_num;i++) //first electrode
    {
        m_p[i].dl=100.0e-9; //100 nm width;
        double alpha=i*1.0/(m_p_num-1);
        m_p[i].r.x = (w_x0)*(1.0-alpha)+w_x1*alpha;
        m_p[i].r.y = (w_y0+m_p[i].dl+5e-9);
        m_p[i].p = 0.0;
        m_p[i].E=0.0;
        m_p[i].ds=_dx*_dz;
    }


    double kap=1.38e-10*0.15;
    m_par.a=(1.0/(m_dt))+(kap*2.0/(m_dx*m_dx));
    m_par.bp=-kap/(dx*dx);
    m_par.bm=-kap/(dx*dx);

    double alp,bet,gam,T,T0,rh;
    alp=3.324e5;
    bet=6.381e8;
    gam=7.89e9;
    T=300;
    T0=381;
    rh=0.0;

    double x0=-0.3;
    m_poly.order=5;
    m_poly.C[0]=2*alp*(T-T0) *0.5;//(T-T0); //x
    m_poly.C[1]=0.0;         //xx
    m_poly.C[2]=-4.0*bet *0.5;   //xxx
    m_poly.C[3]=0.0;        //x^4
    m_poly.C[4]=6.0*gam *0.5; //o.5 from crank-nikolson

    //-(fiy*0.0033- 2.0*alp*81*P1 - 4*bet*P1^3 +6*gam*P1^5)
    double emin=1e23;
    double emax=-1e23;
    for (int i=1; i<N_X-1; i++)
    {
        for (int j=1; j<N_Y-1; j++)
        {
            double lapl0 = kap *( - Py0_[i][j] * (2.0 / (dx * dx) + 2.0 / (dy * dy)) + 1.0 / (dx * dx) * (Py0_[i+1][j] + Py0_[i-1][j]) +  1.0 / (dy * dy) * (Py0_[i][j+1] + Py0_[i][j-1]));
            double poly0 = p.C[0] * Py0_[i][j] +
                    p.C[2] * Py0_[i][j] * Py0_[i][j] * Py0_[i][j] +
                    p.C[4] * Py0_[i][j] * Py0_[i][j] * Py0_[i][j] * Py0_[i][j] * Py0_[i][j];
            RHS_p[i][j]= /*1.4e7*(1.0-i*1.0/N_X)*/ -0.05*(Ey[i][j]+Ey0[i][j]) + lapl0  -poly0 +Py0_[i][j]/(dt_poly);

        }
    }

    //  printf("emin=%e emax=%e \n",emin,emax);
    time1 = get_time();
     //   printf("beforePolyn = %f\n", time1 - time);
    jacobi_polynomial( par_ferr, p,Py_,RHS_p, 4);

}


void pzSolver::solve_poly(int itn)
{

}

void pzSolver::get_q()
{

}
