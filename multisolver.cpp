#include "multisolver.h"

//for fast fields
double g_x_min=1e10;
double g_x_max=-1e10;
double g_y_min=1e10;
double g_y_max=-1e10;
double g_dx,g_dy;
int  g_grid_elec[50][50][200];
vec2 g_grid_mult_elec[50][50][ELEC_MAX];
vec2 g_grid_mult_pz[50][50][PZ_MAX];
int g_grid_elec_num[50][50];
int  g_grid_pz[50][50][200];
int g_grid_pz_num[50][50];

vec2 g_grid_Ephi[51][51];
vec2 g_grid_Ephi_near[51][51][2][2]; //last two indices are the in-cell coordinates

int g_Nx=40;
int g_Ny=20;

double to_elec_from_elec[ELEC_MAX][ELEC_MAX];
double to_elec_from_pz[ELEC_MAX][PZ_MAX];
double to_pzUp_from_elec[PZ_MAX][ELEC_MAX];
double to_pzUp_from_pz[PZ_MAX][PZ_MAX];

double to_pzDown_from_elec[PZ_MAX][ELEC_MAX];
double to_pzDown_from_pz[PZ_MAX][ELEC_MAX];



int traj_table_calculated=0;

vec2 multiSolver::trajectories[traj_num_max][2002];
int multiSolver::trajectories_num[traj_num_max];
int multiSolver::traj_num=0;

void multiSolver::solvePzAdaptive(double dtElec)
{
    double gamma = 1e4;
    double dt_loc = 1e-11;
    m_pzSolver->m_dt =  dt_loc;
    for (int aa=0;aa<int(dtElec *  gamma/ dt_loc);aa++)
    {
        m_pzSolver->get_q();
        updateEforPz();

        /* m_pzSolver->solvePz(4);
        m_pzSolver->step();*/

        //for (int i=0;i<100;i++)
        //{ m_pzSolver->solvePz(4);
        //     m_pzSolver->step();}
        //printf("hereeee-______------- \n");
        // m_pzSolver->solvePzDOPRI(0);
        solvePzDOPRI();
        //m_pzSolver->solvePz_steady(10);
    }
    m_pzSolver->m_dt = dtElec * gamma -  dt_loc * int(dtElec * gamma/ dt_loc);
    m_pzSolver->get_q();
    updateEforPz();

    solvePzDOPRI();
    //m_pzSolver->solvePzDOPRI(0);
    //m_pzSolver->solvePz(4);
    //m_pzSolver->step();
}

static const double a21 = 1.0 / 5.0,        a31 = 3.0 / 40.0,        a32 = 9.0 / 40.0
        ,                   a41 = 44.0 / 45.0,      a42 = -56.0 / 15.0,      a43 = 32.0 / 9.0
        ,                   a51 = 19372.0 / 6561.0, a52 = -25360.0 / 2187.0, a53 = 64448.0 / 6561.0
        ,                   a54 = -212.0 / 729.0,   a61 = 9017.0 / 3168.0,   a62 = -355.0 / 33.0
        ,                   a63 = 46732.0 / 5247.0, a64 = 49.0 / 176.0,      a65 = -5103.0 / 18656.0
        ,                   a71 = 35.0 / 384.0,     a72 = 0.0,               a73 = 500.0 / 1113.0
        ,                   a74 = 125.0 / 192.0,    a75 = -2187.0 / 6784.0,  a76 = 11.0 / 84.0;

void multiSolver::solvePzDOPRI()
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

    double zl=1.0/m_pzSolver->m_p[0].dl;


    for (int i=0; i<m_pzSolver->m_p_num; i++)
    {
        double q= /*qepspi */ m_pzSolver->m_p[i].ds/qe;

        double phi_up=to_pzUp_from_pz[i][i];
        double phi_down=to_pzDown_from_pz[i][i];


        double E_self=-q*(phi_up - phi_down)*zl;//-q*(-log(m_p[i].dl+delta_phi)+log(delta_phi))/m_p[i].dl;
        double currStep = m_pzSolver->m_dt;
        double val = m_pzSolver->m_p[i].p;

        poly_new.C[0]=-2*alp*(T-T0) + E_self;//(T-T0); //x
        double rhs_;
        rhs_=-m_pzSolver->m_p[i].E+(m_pzSolver->m_p[i].p)*E_self+((m_pzSolver->m_p[i-1].p - 2 * m_pzSolver->m_p[i].p +m_pzSolver->m_p[i+1].p) * m_pzSolver->kappa/(m_pzSolver->m_dx*m_pzSolver->m_dx));
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
        m_pzSolver->m_p[i].p = K7;
    }

    if (g_use_wall)
    {
        double r = (m_pzSolver->m_p[g_i_wall_tmp].p+0.26)/0.52;
        while( r>1.0 && g_i_wall_tmp<m_pzSolver->m_p_num-1)
        {
            g_i_wall_tmp++;
            r = (m_pzSolver->m_p[g_i_wall_tmp].p+0.26)/0.52;

        }
        if (rand()*1.0/RAND_MAX<0.01)
        {
            printf("g_i_wall_tmp = %d  r = %f dt_elec=%e \n",g_i_wall_tmp,r, dt_elec);
        }
    }
}

void multiSolver::solve(int itn)
{
    for (int i=0;i<m_Esolver->m_elec_num;i++)
    {
        g_phi=g_phi*0.999+g_phi_max*0.001;
        if (m_Esolver->m_electrodes[i].INDX==0) //first electrode -gphi
        {
            m_Esolver->m_electrodes[i].phi_fix=-g_phi;
        }else
        {
            m_Esolver->m_electrodes[i].phi_fix=g_phi;
        }
    }

    for (int nn=0;nn<itn;nn++)
    {
        double phi_depol0=m_pzSolver->getPhidepol(w_x0,w_y0,is2D);


        for (int i=0;i<m_Esolver->m_elec_num;i++)
        {
            m_Esolver->m_electrodes[i].phi_fix_charges = getPhi_at_electrode(i)-phi_depol0;
        }
        m_Esolver->solve_ls_fast();
        solvePzAdaptive(dt_elec);
    }
    //printf("dt_elec = %e\n", dt_elec);
    g_t+=dt_elec;
    g_save_time+=dt_elec;
    g_save_time2+=dt_elec;
    if (g_emitElectrons)
    {
        electronEmissionEndMoveToElectrode(dt_elec);
        pzEmission(dt_elec);
        // pzEmission_monte_carlo(15e-15,100,0.0003);//0.0003);
        //   pzEmission_monte_carlo(15e-15,200,0.12);//0.0003);
    }else
    {
        //pzEmission_monte_carlo(15e-15,100,0.0003);
        //   pzEmission(dt_elec);
        //pzEmission_monte_carlo(15e-15,200,0.12);//0.0003);
        //  printf("AAAA\n");
    }


    //pzEmissionHoriz(dt_elec);

    // double t1 = get_time();
    //printf("tall= %e t=%e\n", t1-t0, t11-t10);
}


multiSolver::multiSolver()
{
    dt_elec=1e-20;
}

void multiSolver::updateEforPz()
{
    double zl=1.0/m_pzSolver->m_p[0].dl;
    for (int i=0;i<m_pzSolver->m_p_num;i++)
    {
        double phi_up=getPhi_at_pz_up(i);
        double phi_down=getPhi_at_pz_down(i);
        double E_next=-(phi_up - phi_down)*zl;
        m_pzSolver->m_p[i].E = E_next;// - E_self;
    }
}

double  multiSolver::calcJ(double Ein)
{
    double E=Ein/100;//from V/m to V/cm
    double t2 = 1.1;
    double B = 145.0;//145.0
    double phi = 4.0;
    double y = 3.79 * 1e-4 * sqrt(fabs(B * E)) / phi;
    double tetta = 0.95 - 1.03 * y * y;
    return 1e4*(1.54 * 1.0e-6 * B * B * E * E / (t2  * phi)) * exp ( - 6.83 * 1.0e7 * pow(phi, 1.5) * tetta / fabs( B * E)); //in A/m^2
}

void multiSolver::updateTrajTable(bool leapfrog,double dt0,double dl0)
{

    fast_Fields_recalculate();
    slower_Fields_recalculate();

    double l_max = -1e10;
    double full_flux=1; //in electrons/s
    double max_flux=1;

    for (int i=0;i<m_Esolver->m_elec_num-1;i++)
    {
        m_Esolver->m_electrodes[i].EdotN=-1;
        if (m_Esolver->m_electrodes[i].canEmit)
        {


            vec2 E=get_slower_E(m_Esolver->m_electrodes[i].r.x,m_Esolver->m_electrodes[i].r.y);

            double ex=E.x;
            double ey=E.y;

            double l = ex * (-m_Esolver->m_electrodes[i].nx) + ey * (-m_Esolver->m_electrodes[i].ny);

            if (l>0)
            {
                m_Esolver->m_electrodes[i].EdotN=l;
                if(l>l_max)
                    l_max = l;
                double ds = m_Esolver->m_electrodes[i].dl*(w_z1-w_z0);
                double flux = calcJ(l)*ds/(fabs(qe));

                m_Esolver->m_electrodes[i].eCurrent=flux;

                full_flux+=flux;
                if (max_flux<flux) max_flux=flux;
            }
        }
    }


    dt_elec=fmin(electrons_in_pack/full_flux,1e-8);
    // dt_elec=fmin(electrons_in_pack/max_flux,1e-8);



    // printf("full_flux= %e maxL = %e el/s dt=%e t=%e phi=%f \n",full_flux,l_max,dt_elec,g_t,m_Esolver->m_electrodes[0].phi_fix);

    if (leapfrog) //fixed dt
    {
        traj_num=0;
        for (int i=0;i<m_Esolver->m_elec_num;i++) {
            vec2 r = m_Esolver->m_electrodes[i].r;
            vec2 v(0.0,0.0,0.0);
            vec2 acc(0.0, 0.0,0.0);
            bool inArea = true;

            double Dt=dt0;

            vec2 E_=get_slower_E(r.x,r.y);

            endPosTable[i] = -1;

            if (m_Esolver->m_electrodes[i].EdotN>0)
            {
                r.charge=m_Esolver->m_electrodes[i].eCurrent/full_flux;
                trajectories_num[traj_num]=1;
                trajectories[traj_num][0]=r;
                if (traj_num<traj_num_max) traj_num++;



                for (int j=0;j<200;j++)
                {

                    r.x+=Dt*v.x + acc.x*(Dt*Dt*0.5);
                    r.y+=Dt*v.y + acc.y*(Dt*Dt*0.5);

                    E_=get_slower_E(r.x,r.y);

                    double magn=qe/Me;//1e-1;

                    v.x += (magn*(E_.x)+acc.x)*Dt*0.5;
                    v.y += (magn*(E_.y)+acc.y)*Dt*0.5;

                    acc.x=magn*(E_.x);
                    acc.y=magn*(E_.y);

                    if (r.y+Dt*v.y<m_pzSolver->m_p[0].r_top.y) break;
                    r.x+=Dt*v.x;
                    r.y+=Dt*v.y;


                    trajectories[traj_num-1][trajectories_num[traj_num-1]]=r;
                    trajectories_num[traj_num-1]++;

                    if(r.y > w_y1){
                        inArea = false;
                        break;
                    }
                }


                int p_n=0;
                double xx = r.x;//(m_pzSolver->m_p[0].r_top.y - (r.y - r.x * v.y / v.x)) * v.x / v.y;
                r.x=xx;
                r.y=m_pzSolver->m_p[0].r_top.y;//last point
               // trajectories[traj_num-1][trajectories_num[traj_num-1]]=r;
               // trajectories_num[traj_num-1]++;

                p_n=(int) ((xx - m_pzSolver->m_p[0].r_top.x + 0.5 * m_pzSolver->m_dx) / m_pzSolver->m_dx);
                if(inArea && p_n >=0 && p_n<m_pzSolver->m_p_num)
                    endPosTable[i] = p_n;
                else
                    endPosTable[i] = -1;
            }
        }
    }else //fixed dl
    {
        traj_num=0;
        for (int i=0;i<m_Esolver->m_elec_num;i++) {
            vec2 r = m_Esolver->m_electrodes[i].r;
            vec2 v(0.0,0.0,0.0);
            bool inArea = true;
            double Dl=dl0;//0.5e-7;
            double Dt=1e-6;

            vec2 E_=get_slower_E(r.x,r.y);

            endPosTable[i] = -1;

            if (m_Esolver->m_electrodes[i].EdotN>0)
            {
                r.charge=m_Esolver->m_electrodes[i].eCurrent/full_flux;
                trajectories_num[traj_num]=1;
                trajectories[traj_num][0]=r;

                if (traj_num<traj_num_max) traj_num++;


                for (int j=0;j<2000;j++)
                {
                    E_=get_slower_E(r.x,r.y);
                    double magn=qe/Me;//1e-1;
                    double a_=fmax(fabs(magn*(E_.x)),fabs(magn*(E_.y)));
                    double v_=fmax(fabs(v.x),fabs(v.y));
                    Dt=(sqrt(v_*v_+2*a_*Dl)-v_)/(a_+1e-20);

                    v.x += magn*(E_.x)*Dt;
                    v.y += magn*(E_.y)*Dt;
                    if (r.y+Dt*v.y<m_pzSolver->m_p[0].r_top.y) break;
                    r.x+=Dt*v.x;
                    r.y+=Dt*v.y;

                    trajectories[traj_num-1][trajectories_num[traj_num-1]]=r;
                    trajectories_num[traj_num-1]++;

                    if(r.y > w_y1){
                        inArea = false;
                        break;
                    }
                }

                int p_n=0;
                double xx = r.x;//(m_pzSolver->m_p[0].r_top.y - (r.y - r.x * v.y / v.x)) * v.x / v.y;
                r.x=xx;
                r.y=m_pzSolver->m_p[0].r_top.y;//last point
               // trajectories[traj_num-1][trajectories_num[traj_num-1]]=r;
               // trajectories_num[traj_num-1]++;

                p_n=(int) ((xx - m_pzSolver->m_p[0].r_top.x + 0.5 * m_pzSolver->m_dx) / m_pzSolver->m_dx);
                if(inArea && p_n >=0 && p_n<m_pzSolver->m_p_num)
                    endPosTable[i] = p_n;
                else
                    endPosTable[i] = -1;
            }
        }
    }
}

int multiSolver::getEndPos(int i)
{
    return endPosTable[i];
}

void multiSolver::electronEmissionEndMoveToElectrode(double d_t)
{
   for (int i=0;i<m_Esolver->m_elec_num-1;i+=1)
    {
        if (m_Esolver->m_electrodes[i].canEmit)
        {
            vec2 E=get_slower_E(m_Esolver->m_electrodes[i].r.x,m_Esolver->m_electrodes[i].r.y);

            double ex=E.x;
            double ey=E.y;

            double l = ex * (-m_Esolver->m_electrodes[i].nx) + ey * (-m_Esolver->m_electrodes[i].ny);
            if (l>0)
            {
                double ds = m_Esolver->m_electrodes[i].dl*(w_z1-w_z0);
                double el_to_add = calcJ(l)*d_t*ds/(fabs(qe));
                m_Esolver->m_electrodes[i].eToEmit=el_to_add;
                m_Esolver->m_electrodes[i].eCurrent=el_to_add/d_t;
            }else
            {
                m_Esolver->m_electrodes[i].eToEmit=0.0;
                m_Esolver->m_electrodes[i].eCurrent=0.0;
            }
        }
    }
    int i=0;
    int i_left=i;
    int i_right=i+1;

    if (getEndPos(i)>getEndPos(i+1))
    {
        i_left=i+1;
        i_right=i;
    }

    int pzNumLeft = getEndPos(i_left);
    int pzNumRight = getEndPos(i_right);

    if(pzNumLeft >=0 && pzNumRight>=0 ){
        double q_all=0.0;
        for (int j = pzNumLeft; j <= pzNumRight; j++) {
            m_pzSolver->m_p[j].q_tmp = 1.0-fabs(pzNumLeft-j)*1.0/(pzNumRight-pzNumLeft+0.00001);
            q_all+=m_pzSolver->m_p[j].q_tmp;
        }

        for (int j = pzNumLeft; j <= pzNumRight; j++) {
            m_pzSolver->m_p[j].q_ext += m_pzSolver->m_p[j].q_tmp*(m_Esolver->m_electrodes[i].eToEmit/q_all);
        }
    }
    for (int i=1;i<m_Esolver->m_elec_num-1;i++)
    {
        if (m_Esolver->m_electrodes[i].canEmit)
        {
            int i_left=i-1;
            int i_right=i+1;

            if (getEndPos(i-1)>getEndPos(i+1))
            {
                i_left=i+1;
                i_right=i-1;
            }
            int pzNumLeft = getEndPos(i_left);
            int pzNumRight = getEndPos(i_right);
            int pzNumCenter= getEndPos(i);

            if(pzNumLeft >=0 && pzNumRight>=0 && pzNumCenter>=0){
                double q_all=0.0;
                for (int j = pzNumLeft; j < pzNumCenter; j++) {
                    m_pzSolver->m_p[j].q_tmp = 1.0-fabs(pzNumCenter-j)*1.0/(pzNumCenter-pzNumLeft+0.00001);
                    q_all+=m_pzSolver->m_p[j].q_tmp;
                }
                for (int j = pzNumCenter; j <= pzNumRight; j++) {
                    m_pzSolver->m_p[j].q_tmp = 1.0-fabs(pzNumCenter-j)*1.0/(pzNumRight-pzNumCenter+0.00001);
                    q_all+=m_pzSolver->m_p[j].q_tmp;
                }

                for (int j = pzNumLeft; j <= pzNumRight; j++) {
                    m_pzSolver->m_p[j].q_ext += m_pzSolver->m_p[j].q_tmp*(m_Esolver->m_electrodes[i].eToEmit/q_all);
                }

            }
        }
    }
    m_pzSolver->get_q();

}

void multiSolver::init()
{
    m_Esolver->init();

    m_pzSolver->init();
}

void multiSolver::pzEmission(double dt)
{
    double yy=0;//m_pzSolver->m_p[0].r_top.y;


    for (int itn=0;itn<20;itn++)
    {
        double b = 1e-12;//1e-12;
        static vec2 Ef[1000];

        for (int i=0;i<m_pzSolver->m_p_num;i++)
        {
            vec2 r;
            r.x=m_pzSolver->m_p[i].r_top.x;
            r.y=yy;
            //Ef[i] = get_slow_E(r.x,r.y);
        }

        for (int i=1;i<m_pzSolver->m_p_num-1;i++)
        {
            if(m_pzSolver->m_p[i].q_ext /*- m_pzSolver->m_p[i].q_0*/ > 0){
                vec2 r;
                r.x=m_pzSolver->m_p[i].r_top.x;
                r.y=yy;

                vec2 E_l,E_r,E_;
                E_l.x = -(getPhi_at_pz_up(i) - getPhi_at_pz_up(i-1))/(m_pzSolver->m_p[i-1].r_top.x - m_pzSolver->m_p[i].r_top.x);//get_slow_E(r.x,r.y);//Ef[i];
                E_r.x = -(getPhi_at_pz_up(i+1) - getPhi_at_pz_up(i))/(m_pzSolver->m_p[i+1].r_top.x - m_pzSolver->m_p[i].r_top.x);//get_slow_E(r.x,r.y);//Ef[i];

                if (fabs(E_l.x)>fabs(E_r.x))  E_.x=E_l.x;
                else E_.x=E_r.x;
                //double qq = fabs(E_.x * b + (rand()*2.0/RAND_MAX-1.0)*0.01) * (m_pzSolver->m_p[i].q_ext - m_pzSolver->m_p[i].q_0) ;

                int if_sum=(m_pzSolver->m_p[i].q_ext + m_pzSolver->m_p[i].q) > 0;
                int if_qext=m_pzSolver->m_p[i].q_ext > 0;
                double qq = fabs(E_.x * b ) * (m_pzSolver->m_p[i].q_ext + m_pzSolver->m_p[i].q) * if_sum * if_qext ;

                if(qq > (m_pzSolver->m_p[i].q_ext + m_pzSolver->m_p[i].q)*if_sum * if_qext)
                    qq  = (m_pzSolver->m_p[i].q_ext + m_pzSolver->m_p[i].q)* if_sum * if_qext;

                //  double qq = fabs(E_.x * b ) * (m_pzSolver->m_p[i].q_ext ) * if_qext ;
                m_pzSolver->m_p[i].q_ext-=qq;
                if(E_.x>0 && i != 0)
                    m_pzSolver->m_p[i-1].q_ext+=qq;
                if(E_.x<0 && i != m_pzSolver->m_p_num-1)
                    m_pzSolver->m_p[i+1].q_ext+=qq;
            }
        }
    }
}

void multiSolver::prepare_caches(bool _2D_elems)
{
    double  x,y;

    for (int n=0;n<m_Esolver->m_elec_num;n++)
    {
        x=m_Esolver->m_electrodes[n].r.x;
        y=m_Esolver->m_electrodes[n].r.y;

        for (int i=0;i<m_Esolver->m_chargeNum;i++)
        {
            to_elec_from_elec[n][i]=m_Esolver->getPhi_multiplier(x,y,i);
        }

        for (int i=0;i<m_pzSolver->m_p_num;i++)
        {
            to_elec_from_pz[n][i]=m_pzSolver->getPhi_multiplier(x,y,i,_2D_elems);
        }
    }

    for (int n=0;n<m_pzSolver->m_p_num;n++)
    {
        x=m_pzSolver->m_p[n].r_top.x;
        y=m_pzSolver->m_p[n].r_top.y;

        for (int i=0;i<m_Esolver->m_chargeNum;i++)
        {
            to_pzUp_from_elec[n][i]=m_Esolver->getPhi_multiplier(x,y,i);
        }

        for (int i=0;i<m_pzSolver->m_p_num;i++)
        {
            to_pzUp_from_pz[n][i]=m_pzSolver->getPhi_multiplier(x,y,i,_2D_elems);
        }
    }

    for (int n=0;n<m_pzSolver->m_p_num;n++)
    {
        x=m_pzSolver->m_p[n].r_top.x;
        y=m_pzSolver->m_p[n].r_top.y - m_pzSolver->m_p[n].dl;

        for (int i=0;i<m_Esolver->m_chargeNum;i++)
        {
            to_pzDown_from_elec[n][i]=m_Esolver->getPhi_multiplier(x,y,i);
        }

        for (int i=0;i<m_pzSolver->m_p_num;i++)
        {
            to_pzDown_from_pz[n][i]=m_pzSolver->getPhi_multiplier(x,y,i,_2D_elems);;
        }
    }
}
double multiSolver::getPhi_at_electrode(int n)
{
    double sum=0.0;
    /*for (int i=0;i<m_Esolver->m_chargeNum;i++)
    {
        sum+=m_Esolver->m_charges[i].charge*to_elec_from_elec[n][i];
    }*/
    for (int i=0;i<m_pzSolver->m_p_num;i+=MP_DIV)
    {
        sum+=(m_pzSolver->m_p[i].r_top.charge)*to_elec_from_pz[n][i];
    }

    if (g_use_wall)
    {
    sum+=m_pzSolver->getPhiDiff(m_Esolver->m_charges[n].x,m_Esolver->m_charges[n].y,g_i_wall_tmp);
    }

    return sum;
}

double multiSolver::getPhi_at_pz_up(int n)
{
    double sum=0.0;
    for (int i=0;i<m_Esolver->m_chargeNum;i++)
    {
        sum+=m_Esolver->m_charges[i].charge*to_pzUp_from_elec[n][i];
    }
    for (int i=0;i<m_pzSolver->m_p_num;i+=MP_DIV)
    {
        sum+=(m_pzSolver->m_p[i].r_top.charge)*to_pzUp_from_pz[n][i];
    }
    return sum;
}

double multiSolver::getPhi_at_pz_down(int n)
{
    double sum=0.0;
    for (int i=0;i<m_Esolver->m_chargeNum;i++)
    {
        sum+=m_Esolver->m_charges[i].charge*to_pzDown_from_elec[n][i];
    }
    for (int i=0;i<m_pzSolver->m_p_num;i+=MP_DIV)
    {
        sum+=(m_pzSolver->m_p[i].r_top.charge)*to_pzDown_from_pz[n][i];
    }
    return sum;
}



void multiSolver::fast_Fields_prepare()
{
    //1. determine the domain boundaries:
    g_x_min=1e10;
    g_x_max=-1e10;
    g_y_min=1e10;
    g_y_max=-1e10;

    for (int i=0;i<m_Esolver->m_chargeNum;i++)
    {
        if (g_x_max<m_Esolver->m_charges[i].x)
            g_x_max=m_Esolver->m_charges[i].x;

        if (g_y_max<m_Esolver->m_charges[i].y)
            g_y_max=m_Esolver->m_charges[i].y;

        if (g_x_min>m_Esolver->m_charges[i].x)
            g_x_min=m_Esolver->m_charges[i].x;

        if (g_y_min>m_Esolver->m_charges[i].y)
            g_y_min=m_Esolver->m_charges[i].y;
    }

    for (int i=0;i<m_pzSolver->m_p_num;i++)
    {
        if (g_x_max<m_pzSolver->m_p[i].r_top.x)
            g_x_max=m_pzSolver->m_p[i].r_top.x;

        if (g_y_max<m_pzSolver->m_p[i].r_top.y)
            g_y_max=m_pzSolver->m_p[i].r_top.y;

        if (g_x_min>m_pzSolver->m_p[i].r_top.x)
            g_x_min=m_pzSolver->m_p[i].r_top.x;

        if (g_y_min>m_pzSolver->m_p[i].r_top.y)
            g_y_min=m_pzSolver->m_p[i].r_top.y;

    }

    double small_x=(g_x_max-g_x_max)/100.0;
    g_x_min-=small_x;
    g_y_min-=small_x;
    g_x_max+=small_x;
    g_y_max+=small_x;

    printf("xmin=%e xmax=%e ymin=%e ymax=%e \n",g_x_min,g_x_max,g_y_min,g_y_max);

    //2. fillilng the arrays
    g_dx=(g_x_max-g_x_min)/g_Nx;
    g_Ny=(int)((g_y_max-g_y_min)/g_dx);
    g_dy=(g_y_max-g_y_min)/g_Ny;
    printf("dx=%e dy=%e NX=%d Ny=%d \n",g_dx,g_dy,g_Nx,g_Ny);

    for (int i=0;i<g_Nx;i++)
    {
        for (int j=0;j<g_Ny;j++)
        {
            g_grid_elec_num[i][j]=0;
            g_grid_pz_num[i][j]=0;

        }
    }
    for (int i=0;i<m_Esolver->m_chargeNum;i++)
    {
        int ni,nj;
        ni=(int)((m_Esolver->m_charges[i].x-g_x_min)/g_dx);
        nj=(int)((m_Esolver->m_charges[i].y-g_y_min)/g_dy);

        g_grid_elec[ni][nj][g_grid_elec_num[ni][nj]] = i;//&(m_Esolver->m_charges[i]);
        g_grid_elec_num[ni][nj]++;
    }
    for (int i=0;i<m_pzSolver->m_p_num;i++)
    {
        int ni,nj;
        ni=(int)((m_pzSolver->m_p[i].r_top.x-g_x_min)/g_dx);
        nj=(int)((m_pzSolver->m_p[i].r_top.y-g_y_min)/g_dy);

        g_grid_pz[ni][nj][g_grid_pz_num[ni][nj]] = i;//&(m_pzSolver->m_p[i].r_top);
        g_grid_pz_num[ni][nj]++;
    }


    printf("elec num=%d ELEC_MAX=%d \n",m_Esolver->m_chargeNum,ELEC_MAX);

    printf("pz num=%d pz_MAX=%d \n",m_pzSolver->m_p_num,PZ_MAX);

    if ((m_Esolver->m_chargeNum>ELEC_MAX)||(m_pzSolver->m_p_num>PZ_MAX))
    {
        printf ("\n\n\n\n EEEEEEEERRRRRRRRRRRRRRRRRRRRRRRRRORRRRRRRRRRRRRRRRRRRRRRR!!!!!!!!!!!!! ELEC_MAX or PZ_MAX \n\n\n\n");
    }

    for (int i=0;i<=g_Nx;i++)
    {
        for (int j=0;j<=g_Ny;j++)
        {

            double x,y;
            x=g_x_min+i*g_dx;
            y=g_y_min+j*g_dy;

            for (int k=0;k<m_Esolver->m_chargeNum;k++)
            {
                vec2 mult=m_Esolver->get_E_multiplier(x,y,k);
                g_grid_mult_elec[i][j][k]=mult;
            }

            for (int k=0;k<m_pzSolver->m_p_num;k++)
            {
                vec2 mult=m_pzSolver->get_E_multiplier(x,y,k,is2D);
                g_grid_mult_pz[i][j][k]=mult;
            }
        }
    }

    /*
    int sum1=0;
    for (int i=0;i<g_Nx;i++)
        for (int j=0;j<g_Ny;j++)
        {
            if (g_grid_num[i][j]>0)
                printf("i=%d j=%d gnum=%d \n",i,j,g_grid_num[i][j]);
            sum1+=g_grid_num[i][j];
        }

    printf("sum1=%d sum2=%d \n",sum1,m_pzSolver->m_p_num+m_Esolver->m_chargeNum);*/

}
vec2 multiSolver::getE_grid_ij(int i,int j)
{
    vec2 sum(0.0,0.0,0.0);
    for (int k=0;k<m_Esolver->m_chargeNum;k++)
    {

        sum.x+=g_grid_mult_elec[i][j][k].x*m_Esolver->m_charges[k].charge;
        sum.y+=g_grid_mult_elec[i][j][k].y*m_Esolver->m_charges[k].charge;
    }

    for (int k=0;k<m_pzSolver->m_p_num;k++)
    {
        sum.x+=g_grid_mult_pz[i][j][k].x*m_pzSolver->m_p[k].r_top.charge;
        sum.y+=g_grid_mult_pz[i][j][k].y*m_pzSolver->m_p[k].r_top.charge;
    }
    return sum;
}
void multiSolver::fast_Fields_recalculate()
{
    for (int i=0;i<=g_Nx;i++)
    {
        for (int j=0;j<=g_Ny;j++)
        {
            //double x,y;
            //x=g_x_min+i*g_dx;
            //y=g_y_min+j*g_dy;
            vec2 E =getE_grid_ij(i,j);// m_Esolver->getE(x,y);
            //m_pzSolver->getEdepol(x,y,is2D); //можно ускорить через кэш!
            g_grid_Ephi[i][j].x = E.x;
            g_grid_Ephi[i][j].y = E.y;
        }
    }
}

vec2 multiSolver::get_fast_E(double x, double y)
{
    int ni,nj;
    double a=((x-g_x_min)/g_dx);
    double b=((y-g_y_min)/g_dy);
    ni=(int)(a);
    nj=(int)(b);

    a-=ni;
    b-=nj;
    double xa=1.0-a;
    double xb=1.0-b;
    vec2 res;
    res.x=(g_grid_Ephi[ni][nj].x*(xa) + g_grid_Ephi[ni+1][nj].x*(a)) * (xb) +
            (g_grid_Ephi[ni][nj+1].x*(xa) + g_grid_Ephi[ni+1][nj+1].x*(a)) * (b);

    res.y=(g_grid_Ephi[ni][nj].y*(xa) + g_grid_Ephi[ni+1][nj].y*(a)) * (xb) +
            (g_grid_Ephi[ni][nj+1].y*(xa) + g_grid_Ephi[ni+1][nj+1].y*(a)) * (b);


    if (g_use_wall)
    {
    vec2 E_wall=m_pzSolver->getEDiff(x,y,g_i_wall_tmp);
    res.x+=E_wall.x;
    res.y+=E_wall.y;
    }

    return res;
}


void multiSolver::near_Fields_recalculate_cell(int i_, int j_)
{
    int im,ip,jm,jp;

    im= (i_-1>=0) ? i_-1 : 0;
    jm= (j_-1>=0) ? j_-1 : 0;

    ip= (i_+1<g_Nx) ? i_+1 : g_Nx-1;
    jp= (j_+1<g_Ny) ? j_+1 : g_Ny-1;



    vec2 E00(0.0,0.0,0.0);
    vec2 E10(0.0,0.0,0.0);
    vec2 E01(0.0,0.0,0.0);
    vec2 E11(0.0,0.0,0.0);

    double qepspi = (qe/(eps0*pi2))/(w_z1 - w_z0);
    //double delta=1e-6;
    double d2=delta_phi*delta_phi;

    double x00,y00, x10,y01;
    x00=g_x_min+i_*g_dx;
    y00=g_y_min+j_*g_dy;

    x10=g_x_min+(i_+1)*g_dx;
    y01=g_y_min+(j_+1)*g_dy;


    for (int i=im;i<=ip;i++)
    {
        for (int j=jm;j<=jp;j++)
        {
            for (int n=0; n<g_grid_elec_num[i][j]; n++)
            {
                int ind=g_grid_elec[i][j][n];

                E00.x+=g_grid_mult_elec[i_][j_][ind].x*m_Esolver->m_charges[ind].charge;
                E00.y+=g_grid_mult_elec[i_][j_][ind].y*m_Esolver->m_charges[ind].charge;

                E10.x+=g_grid_mult_elec[i_+1][j_][ind].x*m_Esolver->m_charges[ind].charge;
                E10.y+=g_grid_mult_elec[i_+1][j_][ind].y*m_Esolver->m_charges[ind].charge;

                E11.x+=g_grid_mult_elec[i_+1][j_+1][ind].x*m_Esolver->m_charges[ind].charge;
                E11.y+=g_grid_mult_elec[i_+1][j_+1][ind].y*m_Esolver->m_charges[ind].charge;

                E01.x+=g_grid_mult_elec[i_][j_+1][ind].x*m_Esolver->m_charges[ind].charge;
                E01.y+=g_grid_mult_elec[i_][j_+1][ind].y*m_Esolver->m_charges[ind].charge;
            }

            for (int n=0; n<g_grid_pz_num[i][j]; n++)
            {
                int ind=g_grid_pz[i][j][n];

                E00.x+=g_grid_mult_pz[i_][j_][ind].x*m_pzSolver->m_p[ind].r_top.charge;
                E00.y+=g_grid_mult_pz[i_][j_][ind].y*m_pzSolver->m_p[ind].r_top.charge;

                E10.x+=g_grid_mult_pz[i_+1][j_][ind].x*m_pzSolver->m_p[ind].r_top.charge;
                E10.y+=g_grid_mult_pz[i_+1][j_][ind].y*m_pzSolver->m_p[ind].r_top.charge;

                E11.x+=g_grid_mult_pz[i_+1][j_+1][ind].x*m_pzSolver->m_p[ind].r_top.charge;
                E11.y+=g_grid_mult_pz[i_+1][j_+1][ind].y*m_pzSolver->m_p[ind].r_top.charge;

                E01.x+=g_grid_mult_pz[i_][j_+1][ind].x*m_pzSolver->m_p[ind].r_top.charge;
                E01.y+=g_grid_mult_pz[i_][j_+1][ind].y*m_pzSolver->m_p[ind].r_top.charge;
            }
        }
    }

    g_grid_Ephi_near[i_][j_][0][0]=E00;
    g_grid_Ephi_near[i_][j_][1][0]=E10;
    g_grid_Ephi_near[i_][j_][0][1]=E01;
    g_grid_Ephi_near[i_][j_][1][1]=E11;
}


void multiSolver::slower_Fields_recalculate()
{
    for (int i=0;i<g_Nx;i++)
    {
        for (int j=0;j<g_Ny;j++)
        {
            near_Fields_recalculate_cell(i,j);
        }
    }
}

vec2 multiSolver::get_slower_E(double x, double y)
{
    int ni,nj;
    double a=((x-g_x_min)/g_dx);
    double b=((y-g_y_min)/g_dy);
    ni=(int)(a);
    nj=(int)(b);

    a-=ni;
    b-=nj;
    double xa=1.0-a;
    double xb=1.0-b;
    vec2 res;
    res.x=( (g_grid_Ephi[ni][nj].x - g_grid_Ephi_near[ni][nj][0][0].x)*(xa) + (g_grid_Ephi[ni+1][nj].x - g_grid_Ephi_near[ni][nj][1][0].x)*(a) ) * (xb) +
            ( (g_grid_Ephi[ni][nj+1].x - g_grid_Ephi_near[ni][nj][0][1].x)*(xa) + (g_grid_Ephi[ni+1][nj+1].x - g_grid_Ephi_near[ni][nj][1][1].x)*(a)) * (b);

    res.y=( (g_grid_Ephi[ni][nj].y - g_grid_Ephi_near[ni][nj][0][0].y)*(xa) + (g_grid_Ephi[ni+1][nj].y - g_grid_Ephi_near[ni][nj][1][0].y)*(a) ) * (xb) +
            ( (g_grid_Ephi[ni][nj+1].y - g_grid_Ephi_near[ni][nj][0][1].y)*(xa) + (g_grid_Ephi[ni+1][nj+1].y - g_grid_Ephi_near[ni][nj][1][1].y)*(a)) * (b);

    int im,ip,jm,jp;

    im= (ni-1>=0) ? ni-1 : 0;
    jm= (nj-1>=0) ? nj-1 : 0;

    ip= (ni+1<g_Nx) ? ni+1 : g_Nx-1;
    jp= (nj+1<g_Ny) ? nj+1 : g_Ny-1;

    for (int i=im;i<=ip;i++)
    {
        for (int j=jm;j<=jp;j++)
        {
            for (int n=0; n<g_grid_elec_num[i][j]; n++)
            {
                int ind=g_grid_elec[i][j][n];
                vec2 mult=m_Esolver->get_E_multiplier(x,y,ind);
                res.x+=mult.x*m_Esolver->m_charges[ind].charge;
                res.y+=mult.y*m_Esolver->m_charges[ind].charge;
            }
        }
    }

    for (int i=im;i<=ip;i++)
    {
        for (int j=jm;j<=jp;j++)
        {
            for (int n=0; n<g_grid_pz_num[i][j]; n++)
            {
                int ind=g_grid_pz[i][j][n];
                vec2 mult=m_pzSolver->get_E_multiplier(x,y,ind,is2D);
                res.x+=mult.x*m_pzSolver->m_p[ind].r_top.charge;
                res.y+=mult.y*m_pzSolver->m_p[ind].r_top.charge;
            }
        }
    }
    if (g_use_wall)
    {
    vec2 E_wall=m_pzSolver->getEDiff(x,y,g_i_wall_tmp);
    res.x+=E_wall.x;
    res.y+=E_wall.y;
    }

    return res;
}


vec2 multiSolver::get_slow_E(double x, double y)
{

    /*    double qepspi = (qe/(eps0*pi2))/(w_z1 - w_z0);
    //double delta=1e-6;
    double d2=delta_phi*delta_phi;
    vec2 res(0.0,0.0,0.0);
    for (int i=0;i<g_Nx;i++)
    {
        for (int j=0;j<g_Ny;j++)
        {
            for (int n=0; n<g_grid_num[i][j]; n++)
            {
                double r2;
                double q;
                double dx,dy;
                q=- qepspi *g_grid[i][j][n]->charge;

                dx = g_grid[i][j][n]->x - x;
                dy = g_grid[i][j][n]->y - y;
                r2=(dx*dx + dy*dy);
                res.x+=dx*q/(r2+d2);
                res.y+=dy*q/(r2+d2);
            }
        }
    }

    return res;*/
    vec2 res;
    vec2 Ee = m_Esolver->getE(x,y);
    vec2 Epz = m_pzSolver->getEdepol(x,y,is2D);//is2D); //можно ускорить через кэш!
    res.x = Ee.x + Epz.x;
    res.y = Ee.y +  Epz.y;

    if (g_use_wall)
    {
    vec2 E_wall=m_pzSolver->getEDiff(x,y,g_i_wall_tmp);
    res.x+=E_wall.x;
    res.y+=E_wall.y;
    }
    return res;
}


double multiSolver::get_slow_phi(double x, double y, bool _2d)
{
    /*double qepspi = (qe/(eps0*pi2))/(w_z1 - w_z0);
    //double delta=1e-6;
    double d2=delta_phi*delta_phi;
    double  res=0.0;
    for (int i=0;i<g_Nx;i++)
    {
        for (int j=0;j<g_Ny;j++)
        {
            for (int n=0; n<g_grid_num[i][j]; n++)
            {
                double r;
                double q;
                double dx,dy;
                q=- qepspi *g_grid[i][j][n]->charge;

                dx = g_grid[i][j][n]->x - x;
                dy = g_grid[i][j][n]->y - y;
                r=sqrt(dx*dx + dy*dy);
                sum+=q*log(r+delta_phi);
            }
        }
    }
    return res;
    */

    double phi_depol0=m_pzSolver->getPhidepol(w_x0,w_y0,_2d);
    double phi=m_pzSolver->getPhidepol(x,y,_2d)-phi_depol0+ m_Esolver->getPhi(x,y);
            if (g_use_wall)
            {
            phi+=m_pzSolver->getPhiDiff(x,y,g_i_wall_tmp);
            }
    return phi ;
}



