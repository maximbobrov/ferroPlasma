#include "pzsolver.h"
#include "globals.h"


pzSolver::pzSolver()
{
    this->m_p_num=380;
    m_p=new pElem[m_p_num];
    m_dt=5e-11;//1e-11;

    double _dx,_dz;
    _dz=w_z1-w_z0;
    _dx=(w_z1-w_z0)/(m_p_num-1);
    m_dx=_dx;

    for (int i=0;i<m_p_num;i++) //first electrode
    {
        m_p[i].dl=50.0e-9; //100 nm width;
        double alpha=i*1.0/(m_p_num-1);
        m_p[i].r.x = (w_x0)*(1.0-alpha)+w_x1*alpha;
        m_p[i].r.y = (w_y0+(m_p[i].dl-1e-8)*0.5+5e-9);
        m_p[i].p = 0.0;//-0.26;//+rand()*0.043/RAND_MAX;
        m_p[i].p_prev = m_p[i].p;

        m_p[i].E=1.0e8;
        m_p[i].E_prev=m_p[i].E;
        m_p[i].ds=_dx*_dz;
    }


   /* m_p[0].p = 0.26;
    m_p[0].p_prev = m_p[0].p;*/

   /* for (int i=60;i<70;i++) //first electrode
    {
    m_p[i].p = 0.26;// 0.26;
        m_p[i].p_prev = m_p[i].p;
    }*/


    kappa=1.38e-10*0.15;//1.38e-10*0.15;
    m_par.a=(1.0/(m_dt))+(kappa*2.0/(m_dx*m_dx));
    m_par.bp=-kappa/(m_dx*m_dx);
    m_par.bm=-kappa/(m_dx*m_dx);

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
    //    jacobi_polynomial( par_ferr, p,Py_,RHS_p, 4);
}


void pzSolver::solvePz(int itn)
{



   //euler:

  /*  double alp,bet,gam,T,T0,rh;
    alp=3.324e5;
    bet=6.381e8;
    gam=7.89e9;
    T=300;
    T0=381;
    rh=0.0;

    double x0=-0.3;
    m_poly.order=5;
    m_poly.C[0]=2*alp*(T-T0);//(T-T0); //x
    m_poly.C[1]=0.0;         //xx
    m_poly.C[2]=-4.0*bet;   //xxx
    m_poly.C[3]=0.0;        //x^4
    m_poly.C[4]=6.0*gam; //o.5 from crank-nikolson*/


    m_par.a=(1.0/(m_dt))+(kappa*2.0/(m_dx*m_dx));
    m_par.bp=-kappa/(m_dx*m_dx);
    m_par.bm=-kappa/(m_dx*m_dx);
   ////////////




   getRHS();
   poly poly_new;

    double a,b_p,b_m,c_p,c_m;

    a=m_par.a;//((2.0)/(dx*dx)+2.0/(dy*dy));
    b_p=m_par.bp;//-1.0/(dx*dx);
    b_m=m_par.bm;//-1.0/(dx*dx);
    c_p=m_par.cp;//-1.0/(dy*dy);
    c_m=m_par.cm;//-1.0/(dy*dy);

    double rhs_;

    //field_x*a+....+pol*fieldx^..=rhs

    poly_new=m_poly;
    poly_new.C[0]+=a;
    for(int n=0;n<itn;n++)
    {
        for (int i=1; i<m_p_num-1; i++)
        {
            double f_xm=m_p[i-1].p;
            double f_xp=m_p[i+1].p;

            // if (i==1) f_xm=field[N_X-2][j];
            // if (i==N_X-2) f_xp=field[1][j];

            rhs_=m_p[i].RHS-(b_p*f_xp+b_m*f_xm);
            m_p[i].p=m_p[i].p*0.7+0.3*solve_poly(poly_new,m_p[i].p,rhs_,3);
        }
    }

    //step();

}

void pzSolver::getRHS()
{

    for (int i=1; i<m_p_num-1; i++)
    {
        double p_prev=m_p[i].p_prev;
        double pp_m,pp_p;
        pp_m=m_p[i-1].p_prev;
        pp_p=m_p[i+1].p_prev;

        double lapl0 = kappa *(- p_prev * (2.0 / (m_dx * m_dx)) + 1.0 / (m_dx * m_dx) * (pp_m + pp_p));

        double poly0 = m_poly.C[0] * p_prev +
                m_poly.C[2] * p_prev * p_prev * p_prev +
                m_poly.C[4] * p_prev * p_prev * p_prev * p_prev * p_prev;
        m_p[i].RHS = 0.05 * (m_p[i].E+m_p[i].E_prev) + lapl0 - poly0 + p_prev/(m_dt); //crank-nikolson
       // m_p[i].RHS = 0.05 * (m_p[i].E) + p_prev/(m_dt); //euler
    }
}

void pzSolver::get_q() //all charges are in elementary
{

    for (int i=0; i<m_p_num; i++)
    {
        m_p[i].q=-(m_p[i].p)*m_p[i].ds/qe;
    }


 /*   for (int i=60;i<70;i++) //first electrode
    {
        m_p[i].q=-(m_p[i].p)*m_p[i].ds/qe;
    }*/

}

void pzSolver::step()
{
    for (int i=1; i<m_p_num-1; i++)
    {
        m_p[i].p_prev=m_p[i].p;
        m_p[i].E_prev=m_p[i].E;
    }

}

vec3<double> pzSolver::getEdepol(double x, double y)
{
     vec3 <double> sum;
    sum.x=0.0; sum.y=0.0; sum.z=0.0;
    for (int i=0;i<m_p_num;i++)
    {
        double r2;
        double q;
        double dx,dy;
        double delta=1e-9;

        dx = m_p[i].r.x - x;
        dy = m_p[i].r.y + m_p[i].dl*0.5 - y;
        r2=(dx*dx+dy*dy);
        q=qe/eps0*m_p[i].q;

        sum.x+=q*dx*40000/(r2+delta*delta);
        sum.y+=q*dy*40000/(r2+delta*delta);

        dx = m_p[i].r.x - x;
        dy = m_p[i].r.y - m_p[i].dl*0.5 - y;
        r2=(dx*dx+dy*dy);

        //sum.x-=q*500.0*dx/(r2+delta);
        //sum.y-=q*500.0*dy/(r2+delta);


        sum.x-=q*dx*40000/(r2+delta*delta);
        sum.y-=q*dy*40000/(r2+delta*delta);
    }

    /*

    for (int i=0;i<m_elec_num;i++)
    {
        double r2;
        double q;
        double dx,dy;
        double delta=1e-9;

        dx = (m_electrodes[i].r0.x + m_electrodes[i].r1.x)*0.5 - x;
        dy = (m_electrodes[i].r0.y + m_electrodes[i].r1.y)*0.5 - y;
        r2=(dx*dx+dy*dy);
        q=(m_electrodes[i].rho1+m_electrodes[i].rho2)*0.5;

        sum.x+=q*dx/(r2+delta*delta)/log(delta);
        sum.y+=q*dy/(r2+delta*delta)/log(delta);
    }
    return sum;
     */

    return sum;
}

double pzSolver::getPhidepol(double x, double y)
{

    double sum=0.0;
    for (int i=0;i<m_p_num;i++)
    {
    //    int i=1;
        double r2;
        double q;
        double dx,dy;
        double delta=1e-9;

        dx = m_p[i].r.x - x;
        dy = m_p[i].r.y + m_p[i].dl*0.5 - y;
        r2=sqrt(dx*dx+dy*dy);
        q=qe/eps0*m_p[i].q;

      sum+=q*40000*log(r2+delta);
       // sum+=q/(r2+delta);


      dx = m_p[i].r.x - x;
      dy = m_p[i].r.y - m_p[i].dl*0.5 - y;
      r2=sqrt(dx*dx+dy*dy);

              sum-=q*40000*log(r2+delta);
             // sum-=q/(r2+delta);

    }
    return sum;
}
