#include "pzsolver.h"
#include "globals.h"
#include "pzsolver.h"


pzSolver::pzSolver()
{
    this->m_p_num=300;
    m_p=new pElem[m_p_num];
    m_rCentre=new vec2[2 * m_p_num];
    m_dt=1e-7;//15.0*45e-14;//1e-11;

    init();
}

void pzSolver::init()
{
    double _dx,_dz;
    _dz=w_z1-w_z0;
    _dx=1 * (w_x1/*-25e-6-*/- w_x0)/(m_p_num-1);
    m_dx=_dx;

    for (int i=0;i<m_p_num;i++) //first electrode
    {
        m_p[i].dl=dl_pz;// - 1e-6; //100 nm width;
        //double alpha=i*1.0/(m_p_num-1);
        m_p[i].r_top.x = w_x0 + m_dx * (i-g_i_wall_edge);// -18e-6;//(w_x0/*-25e-6*/)*(1.0-alpha)+w_x1*alpha;
        m_p[i].r_top.y = 0.0;//(w_y0+(m_p[i].dl-1e-8)*0.5+5e-9);

        //m_p[i].r_top.x = ;m_p[i].r.x;
        //m_p[i].r_top.y = m_p[i].r.y+m_p[i].dl*0.5;

        /*if (m_p[i].r.x<w_x0+50e-6 && m_p[i].r.x>w_x0+25e-6)
          m_p[i].p = 0.26;//-0.1*(rand()*1.0/RAND_MAX-0.5);//0.0;//-0.005;//-0.26;//+rand()*0.043/RAND_MAX;
          else*/
        m_p[i].p  =-0.26;
        if (i<g_i_wall_edge)
            m_p[i].p = 0.26;//0.0;//-0.26;

        if (i==g_i_wall_edge)
            m_p[i].p = 0.0;//0.0;//-0.26;


        m_p[i].p_prev = m_p[i].p;

        m_p[i].E=1.0e8;

        m_p[i].E_prev=m_p[i].E;
        m_p[i].ds=_dx*_dz;

        m_p[i].q_ext=0.0;
    }
#ifndef USE_MIRROR
    //m_p[0].p = 0.26;//0.005;//0.26;
    //m_p[1].p = 0.26;//0.005;//0.26;
    //m_p[2].p = 0.26;//0.005;//0.26;
#endif
    get_q();
    for (int i=0;i<m_p_num;i++) //first electrode
    {
        m_p[i].q_ext=-m_p[i].q;
        // printf("i=%d q=%e q_ext=%e \n",i,m_p[i].q,m_p[i].q_ext);
        m_p[i].q_0=m_p[i].q_ext;
    }
    //m_dx = 1e-9;
    kappa=1.1 * 1.38e-11;//1000000.1 * 1.38e-10*0.15;//1.38e-10*0.15;
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
    m_poly.C[1]=0.0;              //xx
    m_poly.C[2]=-4.0*bet *0.5;    //xxx
    m_poly.C[3]=0.0;              //x^4
    m_poly.C[4]=6.0*gam *0.5;     //o.5 from crank-nikolson
    get_q();
}


void  pzSolver::setWallPos(double a)
{
    for (int i=0;i<m_p_num;i++) //first electrode
    {
        if ((i*1.0)/(m_p_num-1)<a)
            m_p[i].p = 0.005;//0.26;//+rand()*0.043/RAND_MAX;
        else
            m_p[i].p = -0.005;//-0.26;
    }
    get_q();
}

void pzSolver::conduct(double sigma, double dt, int itn) //surface charge condudtion
{
    //true conduction:
    /*double a=-2e-11;//0.0001;//sigma*dt/m_dx;
    for (int i=0;i<m_p_num;i++)
    {
     m_p[i].q_tmp=m_p[i].q_ext;
    }
    double alpha=0.5;
    //for (int nn=0;nn<itn;nn++)
    {

        m_p[0].q_ext=m_p[0].q_tmp + a*(fmax(m_p[1].q_tmp+m_p[1].q,0)*m_p[1].Ex_s);


        for (int i=1;i<m_p_num-1;i++) //first electrode
        {
         //m_p[i].q_ext=m_p[i].q_ext*(1.0-alpha) + alpha*(m_p[i].q_tmp+a*m_p[i+1].Ex_s*m_p[i+1].q_ext)/(1.0+a*m_p[i].Ex_s);

            m_p[i].q_ext=m_p[i].q_tmp + a*( fmax(m_p[i+1].q_tmp+m_p[i+1].q,0)*m_p[i+1].Ex_s - fmax(m_p[i-1].q_tmp+m_p[i-1].q,0)*m_p[i-1].Ex_s);
        }
            m_p[m_p_num-1].q_ext=m_p[m_p_num-2].q_ext;//=m_p[m_p_num-1].q_tmp + a*( - m_p[m_p_num-2].q_tmp*m_p[m_p_num-2].Ex_s);


    }

    for (int i=0;i<m_p_num;i++)
    {
     m_p[i].q_tmp=0.0;
    }
    */

    double total_to_emit=0;
    int i=0;
    while (i<m_p_num)
    {
        double del_q=m_p[i].q+m_p[i].q_ext;
        if (del_q>=0)
        {
            m_p[i].q_ext-=del_q;
            total_to_emit+=del_q;

            i++;
        }else
        {
            if (total_to_emit+del_q<0)
            {
                m_p[i].q_ext+=total_to_emit;
                total_to_emit=0.0;
                i++;
                //break;
            }else
            {
                m_p[i].q_ext-=del_q;
                total_to_emit+=del_q;
                i++;
            }
        }
    }
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

        for (int i=0; i<m_p_num-1; i++)
        {
            double f_xm;
            if (i>0)
                f_xm=m_p[i-1].p;
            else
                f_xm=m_p[0].p;



            double f_xp=m_p[i+1].p;

            // if (i==1) f_xm=field[N_X-2][j];
            // if (i==N_X-2) f_xp=field[1][j];

            rhs_=m_p[i].RHS-(b_p*f_xp+b_m*f_xm);
            m_p[i].p=m_p[i].p*0.5+0.5*solve_poly(poly_new,m_p[i].p,rhs_,3);
        }
    }

    //step();

}

static const double a21 = 1.0 / 5.0,        a31 = 3.0 / 40.0,        a32 = 9.0 / 40.0
        ,                   a41 = 44.0 / 45.0,      a42 = -56.0 / 15.0,      a43 = 32.0 / 9.0
        ,                   a51 = 19372.0 / 6561.0, a52 = -25360.0 / 2187.0, a53 = 64448.0 / 6561.0
        ,                   a54 = -212.0 / 729.0,   a61 = 9017.0 / 3168.0,   a62 = -355.0 / 33.0
        ,                   a63 = 46732.0 / 5247.0, a64 = 49.0 / 176.0,      a65 = -5103.0 / 18656.0
        ,                   a71 = 35.0 / 384.0,     a72 = 0.0,               a73 = 500.0 / 1113.0
        ,                   a74 = 125.0 / 192.0,    a75 = -2187.0 / 6784.0,  a76 = 11.0 / 84.0;

void pzSolver::solvePzDOPRI(int itn)
{
    double  K1, K2, K3, K4, K5, K6, K7;
    double alp,bet,gam,T,T0,rh;
    alp=3.324e5;
    bet=6.381e8;
    gam=7.89e9;
    T=300;
    T0=381;
    rh=0.0;
    poly poly_new;
    poly_new.order=5;

    poly_new.C[1]=0.0;              //xx
    poly_new.C[2]=4.0*bet;    //xxx
    poly_new.C[3]=0.0;              //x^4
    poly_new.C[4]=-6.0*gam;     //
    double qepspi = (qe/(eps0*pi2))/(w_z1 - w_z0);
    //    int i=1;



    for (int i=0; i<m_p_num; i++)
    {
        double q= qepspi * m_p[i].ds/qe;

        double E_self=-q*(-log(m_p[i].dl+delta_phi)+log(delta_phi))/m_p[i].dl;
        double currStep = m_dt;
        double val = m_p[i].p;

        poly_new.C[0]=-2*alp*(T-T0) + E_self;//(T-T0); //x
        double rhs_;
        rhs_=-m_p[i].E+(m_p[i].p)*E_self+((m_p[i-1].p - 2 * m_p[i].p +m_p[i+1].p) * kappa/(m_dx*m_dx));
        K1 = calc_poly(poly_new,rhs_,val);
        K2 = val + currStep * (a21 * K1);
        K2 = calc_poly(poly_new,rhs_,K2);
        K3 = val + currStep * (a31 * K1 + a32 * K2);
        K3 = calc_poly(poly_new,rhs_,K3);
        K4 = val + currStep * (a41 * K1 + a42 * K2 + a43 * K3);
        K4 = calc_poly(poly_new,rhs_,K4);
        K5 = val + currStep * (a51 * K1 + a52 * K2 + a53 * K3 + a54 * K4);
        K5 = calc_poly(poly_new,rhs_,K5);
        K6 = val + currStep * (a61 * K1 + a62 * K2 + a63 * K3 + a64 * K4 + a65 * K5);
        K6 = calc_poly(poly_new,rhs_,K6);
        K7 = val + currStep * (a71 * K1 + a73 * K3 + a74 * K4 + a75 * K5 + a76 * K6);
        m_p[i].p = K7;
    }
}

void pzSolver::solvePz_steady(int itn) //1d steady state version
{
    //euler:

    // kappa=0.1 * 1.38e-10*0.15;//1.38e-10*0.15;
    m_par.a=(kappa*2.0/(m_dx*m_dx));
    m_par.bp=-kappa/(m_dx*m_dx);
    m_par.bm=-kappa/(m_dx*m_dx);

    double alp,bet,gam,T,T0,rh;
    alp=3.324e5;
    bet=6.381e8;
    gam=7.89e9;
    T=300;
    T0=381;
    rh=0.0;
    //E_self=(m_p[i].p) * log(delta) *m_p[i].ds/(eps0*pi2)/(w_z1 - w_z0)/m_p[i].dl;
    // double delta=1e-6;
    double mult=  log(delta_phi) *m_p[0].ds/(eps0*pi2)/(w_z1 - w_z0)/m_p[0].dl;
    m_poly.order=5;
    m_poly.C[0]=2*alp*(T-T0) -mult;//+1/eps_0 ;//(T-T0); //x
    m_poly.C[1]=0.0;              //xx
    m_poly.C[2]=-4.0*bet;    //xxx
    m_poly.C[3]=0.0;              //x^4
    m_poly.C[4]=6.0*gam;     //o.5 from crank-nikolson


    ////////////

    //getRHS();
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

    //poly_new.C[0]+=a;
    for(int n=0;n<itn;n++)
    {
        m_p[0].p=m_p[1].p;
        for (int i=1; i<m_p_num-1; i++)
        {
            double f_xm=m_p[i-1].p;
            double f_xp=m_p[i+1].p;

            // if (i==1) f_xm=field[N_X-2][j];
            // if (i==N_X-2) f_xp=field[1][j];



            rhs_=m_p[i].E;//-(b_p*f_xp+b_m*f_xm);
            m_p[i].p=m_p[i].p*0.999+0.001*solve_poly(poly_new,m_p[i].p,rhs_,3);
        }
        m_p[m_p_num-1].p=m_p[m_p_num-2].p;
    }

    /* for (int i=1; i<m_p_num-1; i+=30)
    {
    printf("i=%d p=%f E=%e \n",i,m_p[i].p, m_p[i].E);
    }*/
    //step();

}


void pzSolver::getRHS()
{
    double invDx2 = 1.0/ (m_dx * m_dx);
    for (int i=0; i<m_p_num-1; i++)
    {
        double p_prev=m_p[i].p_prev;
        double pp_m,pp_p;

        if (i>0) pp_m=m_p[i-1].p_prev;
        else pp_m=m_p[0].p_prev;

        pp_p=m_p[i+1].p_prev;

        double lapl0 = kappa *(- p_prev * (2.0 * invDx2) + 1.0 * invDx2 * (pp_m + pp_p));

        double pPrev2 = p_prev * p_prev;
        double pPrev3 = pPrev2 * p_prev;
        double poly0 = m_poly.C[0] * p_prev +
                m_poly.C[2] * pPrev3 +
                m_poly.C[4] * pPrev3 * pPrev2;
        m_p[i].RHS = (m_p[i].E+m_p[i].E_prev) + lapl0 - poly0 + p_prev/(m_dt); //crank-nikolson
        // m_p[i].RHS = 0.05 * (m_p[i].E) + p_prev/(m_dt); //euler
    }
}

void pzSolver::get_q() //all charges are in elementary
{

    for (int i=0; i<m_p_num; i++)
    {
        m_p[i].q=(m_p[i].p)*m_p[i].ds/qe;

        if(g_q_enable){
            if (g_i_wall>0)
            {
                if (i<g_i_wall)
                    m_p[i].q_ext=-m_p[i].q;

            }else
            {
                m_p[i].q_ext=-m_p[i].q;
            }
        }
    }

    for (int i=0;i<m_p_num;i+=MP_DIV)
    {
        m_p[i].r_top.charge=0.0;
        for (int j=0;j<MP_DIV;j++)
            m_p[i].r_top.charge += m_p[i + j].q + m_p[i + j].q_ext;

    }





    /*   for (int i=60;i<70;i++) //first electrode
    {
        m_p[i].q=-(m_p[i].p)*m_p[i].ds/qe;
    }*/

}

void pzSolver::step()
{
    for (int i=0; i<m_p_num; i++)
    {
        m_p[i].p_prev=m_p[i].p;
        m_p[i].E_prev=m_p[i].E;
    }

}



vec2 pzSolver::getEdepol(double x, double y, bool _is2d)
{
    vec2 sum(0.0,0.0,0.0);
    double d2=delta_phi*delta_phi;
    double qepspi = (qe/(eps0*pi2))/(w_z1 - w_z0);
    double lx=m_dx; //this should be changed
    double for_q_E=-0.5*qepspi/lx;

    for (int i=0;i<m_p_num;i+=MP_DIV)
    {
        if (!_is2d)
        {

            double r2;
            double q;
            double dx,dy;


            dx = m_p[i].r_top.x - x;
            dy = m_p[i].r_top.y - y;
            r2=(dx*dx+dy*dy);
            q=- qepspi *m_p[i].r_top.charge; // (m_p[i].q+m_p[i].q_ext);

            double c=q/((r2+d2));

            sum.x+=c*dx;
            sum.y+=c*dy;

            /*  dx = m_p[i].r.x - x;
        dy = m_p[i].r.y - m_p[i].dl*0.5 - y;
        r2=(dx*dx+dy*dy);
        q=-qe/(eps0*pi2) * (-m_p[i].q);

        c=q/((r2+delta*delta)*(w_z1 - w_z0));
        sum.x+=dx*c;
        sum.y+=dy*c; //zero charge at the bottom*/
        }else
        {
            double dx0 = m_p[i].r_top.x - lx*0.5 - x;
            double dx1 = m_p[i].r_top.x + lx*0.5 - x;
            dy = m_p[i].r_top.y - y;
            double dydy= dy*dy+1e-20;
            double ln_arg=(dx1*dx1+dydy)/(dx0*dx0+dydy);
            sum.x+=for_q_E*m_p[i].r_top.charge*log(ln_arg);
            sum.y+=2.0*for_q_E*m_p[i].r_top.charge*(atan(dx1/dy) - atan(dx0/dy));//dy*(atan(dx1/dy) - atan(dx0/dy))/(fabs(dy)+1e-10);
        }
    }
    return sum;
}


vec2 pzSolver::get_E_multiplier(double x,double y, int i, bool _is2d)
{
    vec2 sum(0.0,0.0,0.0);
    double d2=delta_phi*delta_phi;
    double qepspi = (qe/(eps0*pi2))/(w_z1 - w_z0);
    double lx=m_dx; //this should be changed
    double for_q_E=-0.5*qepspi/lx;


        if (!_is2d)
        {

            double r2;
            double q;
            double dx,dy;


            dx = m_p[i].r_top.x - x;
            dy = m_p[i].r_top.y - y;
            r2=(dx*dx+dy*dy);
            q=- qepspi /*m_p[i].r_top.charge*/; // (m_p[i].q+m_p[i].q_ext);

            double c=q/((r2+d2));

            sum.x+=c*dx;
            sum.y+=c*dy;


        }else
        {
            double dx0 = m_p[i].r_top.x - lx*0.5 - x;
            double dx1 = m_p[i].r_top.x + lx*0.5 - x;
            dy = m_p[i].r_top.y - y;
            double dydy= dy*dy+1e-20;
            double ln_arg=(dx1*dx1+dydy)/(dx0*dx0+dydy);
            sum.x+=for_q_E*log(ln_arg);
            sum.y+=2.0*for_q_E*(atan(dx1/dy) - atan(dx0/dy));//dy*(atan(dx1/dy) - atan(dx0/dy))/(fabs(dy)+1e-10);
        }

    return sum;

}

double pzSolver::getPhidepol(double x, double y)
{
    double sum=0.0;
    double r;
    double q;
    double dx,dy;
    //double delta=1e-6;

    double qepspi = (qe/(eps0*pi2))/(w_z1 - w_z0);
    for (int i=0;i<m_p_num;i+=MP_DIV)
    {
        //    int i=1;
        dx = m_p[i].r_top.x - x;
        dy = m_p[i].r_top.y - y;
        r=sqrt(dx*dx+dy*dy);
        q= qepspi * m_p[i].r_top.charge;//(m_p[i].q+m_p[i].q_ext);

        sum-=q*log(r+delta_phi);

        /*  dx = m_p[i].r.x - x;
        dy = m_p[i].r.y - m_p[i].dl*0.5 - y;
        r=sqrt(dx*dx+dy*dy);
        q=qe/(eps0*pi2) * (-m_p[i].q);

        sum+=q*log(r+delta)/(w_z1 - w_z0);*/
    }
    return sum;
}

double pzSolver::getPhiDiff(double x, double y, int i)
{
    double phi1d0 = getPhi1D(0,0,
                             m_p[i].r_top.x, m_p[i].r_top.y,
                             m_p[i].q);
    double phi1d = getPhi1D(x, y,
                            m_p[i].r_top.x, m_p[i].r_top.y,
                            m_p[i].q);
    double rig=(m_p[i].p+0.26)/0.52;

    double phi2d = getPhi2D(x, y,
                            m_p[i].r_top.x - m_dx/2,
                            m_p[i].r_top.x - m_dx/2 +  m_dx * rig,
                            m_p[i].r_top.y,
                            m_p[i].q/*m_p[i].r_top.charge*/);

    return phi2d - phi1d0;
}

vec2 pzSolver::getEDiff(double x, double y, int i, double delta)
{
    double phi1d = getPhi1D(x, y,
                            m_p[i].r_top.x, m_p[i].r_top.y,
                            m_p[i].r_top.charge);
    double phi1dx = getPhi1D(x + delta, y,
                             m_p[i].r_top.x, m_p[i].r_top.y,
                             m_p[i].r_top.charge);
    double phi1dy = getPhi1D(x, y + delta,
                             m_p[i].r_top.x, m_p[i].r_top.y,
                             m_p[i].r_top.charge);
    double rig=(m_p[i].p+0.26)/0.52;
    double phi2d = getPhi2D(x, y,
                            m_p[i].r_top.x - m_dx/2,
                            m_p[i].r_top.x - m_dx/2 +  m_dx * rig,
                            m_p[i].r_top.y,
                            m_p[i].r_top.charge);
    double phi2dx = getPhi2D(x + delta, y,
                             m_p[i].r_top.x - m_dx/2,
                             m_p[i].r_top.x - m_dx/2 +  m_dx * rig,
                             m_p[i].r_top.y,
                             m_p[i].r_top.charge);
    double phi2dy = getPhi2D(x, y + delta,
                             m_p[i].r_top.x - m_dx/2,
                             m_p[i].r_top.x - m_dx/2 +  m_dx * rig,
                             m_p[i].r_top.y,
                             m_p[i].r_top.charge);
    vec2 E(0,0,0);
    E.x = -(phi2dx - phi2d - phi1dx + phi1d)/delta;
    E.y = -(phi2dy - phi2d - phi1dy + phi1d)/delta;
    return E;
}


double pzSolver::getPhi2D(double x, double y, double x0,double x1,double y0,double charge) //assuming a horizontal patch of charge with constant density
{
    double sum=0.0;
    double r;
    double q;
    double dx0,dx1,dy,lx;

    lx=x1-x0;
    dx0 = x0 - x;
    dx1 = x1 - x;
    dy = y0 - y;
    double dydy= dy*dy;

    double qepspi = (qe/(eps0*pi2))/(w_z1 - w_z0);
    q= qepspi * charge;
    double part1=0.0;
    if (fabs(dy)>1e-10)
        part1 = dy*(atan(dx1/dy) - atan(dx0/dy));

    double part2 = 0.5*(dx1*log(dx1*dx1+dydy+1e-20) - dx0*log(dx0*dx0+dydy+1e-20));
    sum=q - (q/lx)*(part1 + part2);
    return sum;
}

double pzSolver::getPhi1D(double x, double y, double x0,double y0,double charge) //assuming a horizontal patch of charge with constant density
{
    double sum=0.0;
    double r;
    double q;
    double dx,dy,lx;

    double qepspi = (qe/(eps0*pi2))/(w_z1 - w_z0);
    q= qepspi * charge;


    dx = x0 - x;
    dy = y0 - y;
    sum-=q*log(sqrt(dx*dx+dy*dy));

    return sum;
}


double pzSolver::getE_self(int j)
{
    double sum=0.0;
    double r;
    double q;
    double dx,dy;
    //double delta=1e-6;

    double qepspi = (qe/(eps0*pi2))/(w_z1 - w_z0);
    double x=m_p[j].r_top.x;
    double y=m_p[j].r_top.y;

    for (int i=0;i<m_p_num;i+=MP_DIV)
    {
        //    int i=1;
        dx = m_p[i].r_top.x - x;
        dy = m_p[i].r_top.y - y;
        r=sqrt(dx*dx+dy*dy);
        q= qepspi * m_p[i].r_top.charge;//(m_p[i].q+m_p[i].q_ext);

        sum-=q*log(r+delta_phi);

        //  dx = m_p[i].r.x - x;
        dy = m_p[i].r_top.y - m_p[i].dl - y;
        r=sqrt(dx*dx+dy*dy);
        q= -qepspi * m_p[i].r_top.charge;

        sum-=q*log(r+delta_phi);
    }

    y=m_p[j].r_top.y - m_p[j].dl;

    double sum_b=0.0;
    for (int i=0;i<m_p_num;i+=MP_DIV)
    {
        //    int i=1;
        dx = m_p[i].r_top.x - x;
        dy = m_p[i].r_top.y - y;
        r=sqrt(dx*dx+dy*dy);
        q= qepspi * m_p[i].r_top.charge;//(m_p[i].q+m_p[i].q_ext);

        sum_b-=q*log(r+delta_phi);

        //  dx = m_p[i].r.x - x;
        dy = m_p[i].r_top.y - m_p[i].dl - y;
        r=sqrt(dx*dx+dy*dy);
        q= -qepspi * m_p[i].r_top.charge;

        sum_b-=q*log(r+delta_phi);
    }


    return (sum_b-sum)/m_p[j].dl;
}
