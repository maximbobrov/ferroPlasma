#include "pzsolver.h"
#include "globals.h"


pzSolver::pzSolver()
{
    this->m_p_num=380;
    m_p=new pElem[m_p_num];
    m_rCentre=new vec2[2 * m_p_num];
    m_dt=5.0*45e-14;//1e-11;

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
        m_p[i].dl=dl_pz; //100 nm width;
        //double alpha=i*1.0/(m_p_num-1);
        m_p[i].r.x = w_x0 + m_dx * i;//(w_x0/*-25e-6*/)*(1.0-alpha)+w_x1*alpha;
        m_p[i].r.y = w_y0+25e-6;//(w_y0+(m_p[i].dl-1e-8)*0.5+5e-9);

          /*if (m_p[i].r.x<w_x0+50e-6 && m_p[i].r.x>w_x0+25e-6)
          m_p[i].p = 0.26;//-0.1*(rand()*1.0/RAND_MAX-0.5);//0.0;//-0.005;//-0.26;//+rand()*0.043/RAND_MAX;
          else*/
        m_p[i].p = -0.26;
        #ifdef USE_MIRROR
        m_p[i].p=0.0;//-0.26;
        #endif

        m_p[i].p_prev = m_p[i].p;

        m_p[i].E=1.0e8;

        m_p[i].E_prev=m_p[i].E;
        m_p[i].ds=_dx*_dz;

        m_p[i].q_ext=0.0;
    }
    #ifndef USE_MIRROR
    m_p[0].p = 0.26;//0.005;//0.26;
    //m_p[1].p = 0.26;//0.005;//0.26;
    //m_p[2].p = 0.26;//0.005;//0.26;
    #endif
    get_q();
    for (int i=0;i<m_p_num;i++) //first electrode
    {
        m_p[i].q_ext=-m_p[i].q;
        printf("i=%d q=%e q_ext=%e \n",i,m_p[i].q,m_p[i].q_ext);
    }
    //m_p[100].q_ext+=100;
    //m_p[120].q_ext-=100;

    /* m_p[0].p = 0.26;
    m_p[0].p_prev = m_p[0].p;*/

    /* for (int i=60;i<70;i++) //first electrode
    {
    m_p[i].p = 0.26;// 0.26;
        m_p[i].p_prev = m_p[i].p;
    }*/

    //m_dx = 1e-9;
    kappa=0.1 * 1.38e-10*0.15;//1.38e-10*0.15;
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

    //-(fiy*0.0033- 2.0*alp*81*P1 - 4*bet*P1^3 +6*gam*P1^5)
    //    jacobi_polynomial( par_ferr, p,Py_,RHS_p, 4);
    //updateGridProp();
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
        m_p[i].RHS = (m_p[i].E+m_p[i].E_prev) + lapl0 - poly0 + p_prev/(m_dt); //crank-nikolson
        // m_p[i].RHS = 0.05 * (m_p[i].E) + p_prev/(m_dt); //euler
    }
}

void pzSolver::get_q() //all charges are in elementary
{

    for (int i=0; i<m_p_num; i++)
    {
        m_p[i].q=(m_p[i].p)*m_p[i].ds/qe;

        /*if (m_p[i].r.x<w_x0+50e-6 && m_p[i].r.x>w_x0+25e-6)
        {
            m_p[i].q=0.0;
            m_p[i].q_ext=0.0;
        }*/

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
    updateCharge();
}

void pzSolver::updateCharge()
{
    for( int i=0; i< m_p_num; i++ ){
        m_rCentre[i].charge = -m_p[i].q;
        m_rCentre[m_p_num + i].charge = m_p[i].q + m_p[i].q_ext;
    }
}

void pzSolver::updateGridProp()
{
    m_gridProp.NX = 30;
    m_gridProp.NY = 2;
    m_gridProp.startx = w_x0 - (w_x1-w_x0)*0.01;
    m_gridProp.starty = w_y0 - (w_y1-w_y0)*0.01;
    m_gridProp.endx = w_x1 + (w_x1-w_x0)*0.01;
    m_gridProp.endy =  w_y0 + (w_y1-w_y0)*0.501;
    m_gridProp.dx = (m_gridProp.endx - m_gridProp.startx) / (m_gridProp.NX  - 1);
    m_gridProp.dy = (m_gridProp.endy - m_gridProp.starty) / (m_gridProp.NY  - 1);

    for (int i = 0; i < m_gridProp.NX; i++)
        for (int j = 0; j < m_gridProp.NY; j++) {
            m_gridProp.gridCenters[i][j].x = 0.0;
            m_gridProp.gridCenters[i][j].y = 0.0;
            m_gridProp.gridCenters[i][j].charge = 0.0;
            m_gridProp.gridNeighbors[i][j].clear();
        }
    for( int i=0; i< m_p_num; i++ ){
        vec2 pos1(m_p[i].r.x, m_p[i].r.y - m_p[i].dl*0.5, -m_p[i].q);
        m_rCentre[i] = pos1;
        int x1Idx = fmin(int((pos1.x - m_gridProp.startx) / m_gridProp.dx), m_gridProp.NX - 1);
        int y1Idx = fmin(int((pos1.y - m_gridProp.starty) / m_gridProp.dy), m_gridProp.NY - 1);
        m_gridProp.gridCenters[x1Idx][y1Idx].charge += pos1.charge;
        m_gridProp.gridCenters[x1Idx][y1Idx].x += pos1.charge * pos1.x;
        m_gridProp.gridCenters[x1Idx][y1Idx].y += pos1.charge * pos1.y;
        m_gridProp.gridNeighbors[x1Idx][y1Idx].push_back(i);

        vec2 pos2(m_p[i].r.x, m_p[i].r.y + m_p[i].dl*0.5, m_p[i].q + m_p[i].q_ext);
        m_rCentre[m_p_num + i] = pos2;
        int x2Idx = fmin(int((pos2.x - m_gridProp.startx) / m_gridProp.dx), m_gridProp.NX - 1);
        int y2Idx = fmin(int((pos2.y - m_gridProp.starty) / m_gridProp.dy), m_gridProp.NY - 1);
        m_gridProp.gridCenters[x2Idx][y2Idx].charge += pos2.charge;
        m_gridProp.gridCenters[x2Idx][y2Idx].x += pos2.charge * pos2.x;
        m_gridProp.gridCenters[x2Idx][y2Idx].y += pos2.charge * pos2.y;
        m_gridProp.gridNeighbors[x2Idx][y2Idx].push_back(i);
    }
    for (int i = 0; i < m_gridProp.NX; i++)
        for (int j = 0; j < m_gridProp.NY; j++) {
            if(m_gridProp.gridNeighbors[i][j].size() == 0)
                continue;
            m_gridProp.gridCenters[i][j].x /= m_gridProp.gridCenters[i][j].charge;
            m_gridProp.gridCenters[i][j].y /= m_gridProp.gridCenters[i][j].charge;
        }
}
/*
vec2 pzSolver::getEField(const vec2& iCenterPos, const vec2& iFarPos)
{
    vec2 E;
    vec2 dist;

    double r2;
    double q;
    double dx,dy;
    double delta=1e-9;

    dx = iFarPos.x - iCenterPos.x;
    dy = iFarPos.y - iCenterPos.y;
    r2=(dx*dx+dy*dy);
    q=qe/eps0*(iFarPos.charge);

    E.x=q*dx*40000/(r2+delta*delta);
    E.y=q*dy*40000/(r2+delta*delta);

    return  E;
}

vec2 pzSolver::getPhiField(const vec2& iCenterPos, const vec2& iFarPos)
{
    vec2 E;
    vec2 dist;
    float invDist2;
    double delta=1e-9;

    dist.x = iCenterPos.x - iFarPos.x;
    dist.y = iCenterPos.y - iFarPos.y;

    invDist2 = iFarPos.charge / (dist.x*dist.x+dist.y*dist.y+1e-14)/log(delta);

    E.x = iFarPos.charge*log((dist.x*dist.x+dist.y*dist.y+1e-14)+delta)/log(delta);
    return  E;
}
*/
vec2 pzSolver::getEdepol(double x, double y)
{
    vec2 EField;
    vec2 pos(x,y,0.0);
    //getFieldFast(pos, m_rCentre, getEField, EField);
    //return EField;
    vec2 sum;
    sum.x=0.0; sum.y=0.0;
    for (int i=0;i<m_p_num;i++)
    {
        double r2;
        double q;
        double dx,dy;
        double delta=1e-6;

        dx = m_p[i].r.x - x;
        dy = m_p[i].r.y + m_p[i].dl*0.5 - y;
        r2=(dx*dx+dy*dy);
        q=-qe/(eps0*pi2) * (m_p[i].q+m_p[i].q_ext);

        double c=q/((r2+delta*delta)*(w_z1 - w_z0));

        sum.x+=c*dx;
        sum.y+=c*dy;

        /*  dx = m_p[i].r.x - x;
        dy = m_p[i].r.y - m_p[i].dl*0.5 - y;
        r2=(dx*dx+dy*dy);
        q=-qe/(eps0*pi2) * (-m_p[i].q);

        c=q/((r2+delta*delta)*(w_z1 - w_z0));
        sum.x+=dx*c;
        sum.y+=dy*c; //zero charge at the bottom*/
    }
    return sum;
}

double pzSolver::getPhidepol(double x, double y)
{
    vec2 PhiField;
    vec2 pos(x,y,0.0);
    //getFieldFast(pos, m_rCentre, getPhiField, PhiField);
    //return PhiField.x;
    double sum=0.0;
    for (int i=0;i<m_p_num;i++)
    {
        //    int i=1;
        double r;
        double q;
        double dx,dy;
        double delta=1e-6;

        dx = m_p[i].r.x - x;
        dy = m_p[i].r.y + m_p[i].dl*0.5 - y;
        r=sqrt(dx*dx+dy*dy);
        q=qe/(eps0*pi2) * (m_p[i].q+m_p[i].q_ext);

        sum-=q*log(r+delta)/(w_z1 - w_z0);

        /*  dx = m_p[i].r.x - x;
        dy = m_p[i].r.y - m_p[i].dl*0.5 - y;
        r=sqrt(dx*dx+dy*dy);
        q=qe/(eps0*pi2) * (-m_p[i].q);

        sum+=q*log(r+delta)/(w_z1 - w_z0);*/


    }
    return sum;
}
