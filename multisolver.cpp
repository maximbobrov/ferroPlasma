#include "multisolver.h"

multiSolver::multiSolver()
{
    dt_elec=1e-20;
}

void multiSolver::updateEforPz()
{

    //double phi_down=m_Esolver->m_electrodes[m_Esolver->m_elec_num-1].phi_fix; //phi at lower electrode
    // printf("phi_down= %f \n",phi_down);

    double phi_depol0=m_pzSolver->getPhidepol(w_x0,w_y0);


    for (int i=0;i<m_pzSolver->m_p_num;i++)
    {
        //by Efield direct calculation
        /* double x,ym, y0,yp;
        x=m_pzSolver->m_p[i].r.x;
        ym=m_pzSolver->m_p[i].r.y+m_pzSolver->m_p[i].dl*0.22;
        y0=m_pzSolver->m_p[i].r.y;
        yp=m_pzSolver->m_p[i].r.y+m_pzSolver->m_p[i].dl*0.49;

        vec2 Em = m_Esolver->getE(x,ym);
        vec2 Epm = m_pzSolver->getEdepol(x,ym);
        vec2 Eem = m_elecSolver->getEe(x,ym);

        vec2 E0 = m_Esolver->getE(x,y0);
        vec2 Ep0 = m_pzSolver->getEdepol(x,y0);
        vec2 Ee0 = m_elecSolver->getEe(x,y0);

        vec2 Ep = m_Esolver->getE(x,yp);
        vec2 Epp = m_pzSolver->getEdepol(x,yp);
        vec2 Eep = m_elecSolver->getEe(x,yp);

        m_pzSolver->m_p[i].E =(Em.y + Epm.y+E0.y + Ep0.y+ Ep.y + Epp.y)/3.0 + (Eem.y +Ee0.y+Eep.y)/3;;
        */


        double x,yp,ym;
        x=m_pzSolver->m_p[i].r.x;
        yp=m_pzSolver->m_p[i].r.y+m_pzSolver->m_p[i].dl*0.5;
        ym=m_pzSolver->m_p[i].r.y-m_pzSolver->m_p[i].dl*0.5;

        double phi_er = m_Esolver->getPhi(x,yp);
        double phi_pz = m_pzSolver->getPhidepol(x,yp);

        double phi_up = (phi_er + (phi_pz-phi_depol0));

        phi_er = m_Esolver->getPhi(x,ym);
        phi_pz = m_pzSolver->getPhidepol(x,ym);

        double phi_down = (phi_er + (phi_pz-phi_depol0));


        double delta=1e-6; //this should be monitored
        //  m_p[i].q=(m_p[i].p)*m_p[i].ds/qe;

        //q=qe/(eps0*pi2) * (m_p[i].q+m_p[i].q_ext);
        double q_self=(m_pzSolver->m_p[i].p)*m_pzSolver->m_p[i].ds/(eps0*pi2);
        double phi_self=-q_self*log(delta)/(w_z1 - w_z0);
        double E_self=-phi_self/m_pzSolver->m_p[i].dl; //extract



        m_pzSolver->m_p[i].E = -(phi_up - phi_down)/m_pzSolver->m_p[i].dl - E_self;

        //now Ex_s Ey_s at the surface:

        vec2 Em = m_Esolver->getE(x,yp);
        vec2 Epm = m_pzSolver->getEdepol(x,yp);
        m_pzSolver->m_p[i].Ex_s =Em.x+Epm.x;
        m_pzSolver->m_p[i].Ey_s =Em.y+Epm.y;


    }
}

void multiSolver::updateEforElec()
{
    for (int i=0;i<m_elecSolver->m_numParticles;i++)
    {
        double x,y;
        x=m_elecSolver->m_bodyPos[i].x;
        y=m_elecSolver->m_bodyPos[i].y;
        
        vec2 E=m_Esolver->getE(x,y);
        vec2 Ep=m_pzSolver->getEdepol(x,y);
        vec2 Ee = m_elecSolver->getEe(x,y);
        m_elecSolver->m_bodyE[i].x=-(E.x+Ep.x+Ee.x);
        m_elecSolver->m_bodyE[i].y=-(E.y+Ep.y+Ee.y);
    }
}

/*void multiSolver::electronEmission(double d_t)
{
    double eMean=0.0;
    double emass=0.0;
    for (int i=0;i<m_Esolver->m_elec_num;i++)
    {
        if (m_Esolver->m_electrodes[i].canEmit)
        {
            vec2 E=m_Esolver->getE(m_Esolver->m_electrodes[i].r.x+1e-9,m_Esolver->m_electrodes[i].r.y);
            vec2 Ep=m_pzSolver->getEdepol(m_Esolver->m_electrodes[i].r.x+1e-9,m_Esolver->m_electrodes[i].r.y);
            vec2 Ee=m_elecSolver->getEe(m_Esolver->m_electrodes[i].r.x+1e-9,m_Esolver->m_electrodes[i].r.y);
            double ex=E.x+Ep.x+Ee.x;
            double ey=E.y+Ep.y+Ee.y;
            
            double l = ex * (-m_Esolver->m_electrodes[i].nx) + ey * (-m_Esolver->m_electrodes[i].ny);
            if (l>0)
            {
            eMean+=l;
            emass+=1.0;
            m_elecSolver->create_electron(m_Esolver->m_electrodes[i].r,l,d_t,m_Esolver->m_electrodes[i].dl*(w_z1-w_z0));
            }
        }
    }
    // printf("curr_elec_num=%d mass=%f E=%e\n",m_elecSolver->m_numParticles,emass, eMean/emass);
}
*/

void multiSolver::electronEmission(double d_t)
{
    for (int i=0;i<m_Esolver->m_elec_num;i++)
    {
        if (m_Esolver->m_electrodes[i].canEmit)
        {
            vec2 E=m_Esolver->getE(m_Esolver->m_electrodes[i].r.x,m_Esolver->m_electrodes[i].r.y);
            vec2 Ep=m_pzSolver->getEdepol(m_Esolver->m_electrodes[i].r.x,m_Esolver->m_electrodes[i].r.y);
            vec2 Ee=m_elecSolver->getEe(m_Esolver->m_electrodes[i].r.x,m_Esolver->m_electrodes[i].r.y);
            double ex=E.x+Ep.x+Ee.x;
            double ey=E.y+Ep.y+Ee.y;
            //m_Esolver->m_electrodes[i].charge*=0.6;
            
            double l = ex * (-m_Esolver->m_electrodes[i].nx) + ey * (-m_Esolver->m_electrodes[i].ny);
            if (l>0)
            {
                double ds = m_Esolver->m_electrodes[i].dl*(w_z1-w_z0);
                
                // double el_to_add = l/(E0 + 0.01);//m_elecSolver->calcJ(l)*d_t*ds/(fabs(qe)/**num_in_pack*/);
                //  for (int i = 0; i<200;i++) {
                //    el_to_add = el_to_add * 0.99 + 0.01 * m_elecSolver->calcJ(fmax(l - E0 * (/*m_Esolver->m_electrodes[i].eToEmit +*/ el_to_add),0.0))*d_t*ds/(fabs(qe)/**num_in_pack*/);
                // }
                // double E0=m_elecSolver->getEmult_dipole(2.0e-6);
                
                //double el_to_add = emis_tab.get_f(l, 1e-14) * d_t/1e-14;//m_elecSolver->calcJ(l)*d_t*ds/(fabs(qe));
                double el_to_add = m_elecSolver->calcJ(l)*d_t*ds/(fabs(qe));
                m_Esolver->m_electrodes[i].eToEmit+=el_to_add;
                m_Esolver->m_electrodes[i].eCurrent=l;
                
                /*if (i==0)
                    printf("E=%le J=%le to emit:%le\n", l, m_elecSolver->calcJ(l), el_to_add);*/
                
                if(m_Esolver->m_electrodes[i].eToEmit > 500)
                {
                    // m_Esolver->m_electrodes[i].charge-=m_Esolver->m_electrodes[i].eToEmit;
                    vec2 elecPos(m_Esolver->m_electrodes[i].r.x + 5e-7 * m_Esolver->m_electrodes[i].nx, m_Esolver->m_electrodes[i].r.y + 5e-7 * m_Esolver->m_electrodes[i].ny,0);
                    vec2 vel(0.00001 * 1.5e6 * m_Esolver->m_electrodes[i].nx, 0.00001 *  1.5e6 * m_Esolver->m_electrodes[i].ny,0);
                    //vec2 vel(0.0000 * 1.5e6 * m_Esolver->m_electrodes[i].nx, 0.0000 *  1.5e6 * m_Esolver->m_electrodes[i].ny,0);
                    
                    m_elecSolver->create_electrons(elecPos,vel,int(m_Esolver->m_electrodes[i].eToEmit));
                    m_Esolver->m_electrodes[i].eToEmit-=int(m_Esolver->m_electrodes[i].eToEmit);
                }
            }
        }
    }
    // printf("curr_elec_num=%d mass=%f E=%e\n",m_elecSolver->m_numParticles,emass, eMean/emass);
}

void multiSolver::updateTrajTable()
{


    double full_flux=1; //in electrons/s
    double electrons_in_pack=800;
    for (int i=0;i<m_Esolver->m_elec_num-1;i++)
    {
        if (m_Esolver->m_electrodes[i].canEmit)
        {
            vec2 E=m_Esolver->getE(m_Esolver->m_electrodes[i].r.x,m_Esolver->m_electrodes[i].r.y);
            vec2 Ep=m_pzSolver->getEdepol(m_Esolver->m_electrodes[i].r.x,m_Esolver->m_electrodes[i].r.y);
            double ex=E.x+Ep.x;
            double ey=E.y+Ep.y;
            double l = ex * (-m_Esolver->m_electrodes[i].nx) + ey * (-m_Esolver->m_electrodes[i].ny);
            if (l>0)
            {
                double ds = m_Esolver->m_electrodes[i].dl*(w_z1-w_z0);
                double flux = m_elecSolver->calcJ(l)*ds/(fabs(qe));
                full_flux+=flux;
            }
        }
    }


            dt_elec=fmin(electrons_in_pack/full_flux,1e-8);



    printf("full_flux= %e el/s dt=%e t=%e phi=%f \n",full_flux,dt_elec,g_t,m_Esolver->m_electrodes[0].phi_fix);


    for (int i=0;i<m_Esolver->m_elec_num;i++) {
        vec2 r = m_Esolver->m_electrodes[i].r;
        vec2 v(0.0,0.0,0.0);

        double Dl=0.5e-6;
        double Dt=1e-6;
        for (int j=0;j<1000;j++)
        {

            //vec2 Ee = elec_solver->getEe(r.x,r.y);
            vec2 Ed = m_Esolver->getE(r.x,r.y);
            vec2 Epz = m_pzSolver->getEdepol(r.x,r.y);
            vec2 E_;
            E_.x = Ed.x + Epz.x;
            E_.y = Ed.y + Epz.y;


            double magn=qe/Me;//1e-1;
            double a_=fmax(fabs(magn*(E_.x)),fabs(magn*(E_.y)));
            double v_=fmax(fabs(v.x),fabs(v.y));
            Dt=(sqrt(v_*v_+2*a_*Dl)-v_)/(a_+1e-20);

            v.x += magn*(E_.x)*Dt;
            v.y += magn*(E_.y)*Dt;
            if (r.y+Dt*v.y<0) break;
            r.x+=Dt*v.x;
            r.y+=Dt*v.y;
        }

        int p_n=0;
        double xx = (m_pzSolver->m_p[0].r.y + m_pzSolver->m_p[0].dl * 0.5 - (r.y - r.x * v.y / v.x)) * v.x / v.y;
        p_n=(int) ((xx - m_pzSolver->m_p[0].r.x + 0.5 * m_pzSolver->m_dx) / m_pzSolver->m_dx);
        if(p_n >=0 && p_n<m_pzSolver->m_p_num)
            endPosTable[i] = p_n;
        else
            endPosTable[i] = -1;
    }
}

int multiSolver::getEndPos(int i)
{
    return endPosTable[i];
}

void multiSolver::electronEmissionEndMoveToElectrode(double d_t)
{

    //find emitted charge
    for (int i=0;i<m_Esolver->m_elec_num-1;i++)
    {
        if (m_Esolver->m_electrodes[i].canEmit)
        {
            vec2 E=m_Esolver->getE(m_Esolver->m_electrodes[i].r.x,m_Esolver->m_electrodes[i].r.y);
            vec2 Ep=m_pzSolver->getEdepol(m_Esolver->m_electrodes[i].r.x,m_Esolver->m_electrodes[i].r.y);
            double ex=E.x+Ep.x;
            double ey=E.y+Ep.y;
            double l = ex * (-m_Esolver->m_electrodes[i].nx) + ey * (-m_Esolver->m_electrodes[i].ny);
            if (l>0)
            {
                double ds = m_Esolver->m_electrodes[i].dl*(w_z1-w_z0);
                double el_to_add = m_elecSolver->calcJ(l)*d_t*ds/(fabs(qe));
                m_Esolver->m_electrodes[i].eToEmit=el_to_add;
                m_Esolver->m_electrodes[i].eCurrent=m_elecSolver->calcJ(l);
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
            double e_t=m_Esolver->m_electrodes[i].eToEmit;

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

}



void multiSolver::init()
{
    m_Esolver->init();
    m_elecSolver->init();
    m_pzSolver->init();
}


void multiSolver::checkPotential()
{
    double phiMin = 1e100;
    double phiMax = -1e100;
    double phi_depol0=m_pzSolver->getPhidepol(w_x0,w_y0);
    double elec_depol0=m_elecSolver->getPhiSlow(w_x0,w_y0);
    for (int i=0;i<85;i++)
    {
        double x = m_Esolver->m_electrodes[i].r.x;
        double y = m_Esolver->m_electrodes[i].r.y;
        double phi = m_Esolver->getPhi(x, y)
                +m_elecSolver->getPhiSlow(x, y)
                +m_pzSolver->getPhidepol(x,y) - phi_depol0 - elec_depol0;
        phiMax = fmax(phiMax, phi);
        phiMin = fmin(phiMin, phi);
        //printf("phi1 = %f  phi2 = %f phi3 = %f\n", m_Esolver->getPhi(x, y),m_elecSolver->getPhiSlow(x, y),m_pzSolver->getPhidepol(x,y));
    }
    //printf("phiMax1 = %f phiMin1 = %f \n", phiMax, phiMin);
    
    phiMin = 1e100;
    phiMax = -1e100;
    for (int i=85;i<m_Esolver->m_elec_num;i++)
    {
        double x = m_Esolver->m_electrodes[i].r.x;
        double y = m_Esolver->m_electrodes[i].r.y;
        double phi = m_Esolver->getPhi(x, y)
                +m_elecSolver->getPhiSlow(x, y)
                +m_pzSolver->getPhidepol(x,y)- phi_depol0 - elec_depol0;
        phiMax = fmax(phiMax, phi);
        phiMin = fmin(phiMin, phi);
    }
    //printf("phiMax2 = %f phiMin2 = %f \n", phiMax, phiMin);
}
int traj_table_calculated=0;
void multiSolver::solve(int itn)
{

/*
    for (int i=0;i<m_Esolver->m_elec_num;i++)
    {
        if (m_Esolver->m_electrodes[i].phi_fix>0)
        {
            if (g_t<g_t_max)
                m_Esolver->m_electrodes[i].phi_fix=(g_t/g_t_max)*g_phi_max;
            else
                m_Esolver->m_electrodes[i].phi_fix=g_phi_max;
        }else
        {
            if (g_t<g_t_max)
                m_Esolver->m_electrodes[i].phi_fix=-(g_t/g_t_max)*g_phi_max;
            else
                m_Esolver->m_electrodes[i].phi_fix=-g_phi_max;
        }
    }*/

    for (int nn=0;nn<itn;nn++)
    {
#ifndef USE_MIRROR
        m_pzSolver->get_q();
#endif
        double phi_depol0=m_pzSolver->getPhidepol(w_x0,w_y0);

        for (int i=0;i<m_Esolver->m_elec_num;i++)
        {
            double x,y;
            x=m_Esolver->m_electrodes[i].r.x;
            y=m_Esolver->m_electrodes[i].r.y;

            m_Esolver->m_electrodes[i].phi_fix_charges = (m_pzSolver->getPhidepol(x,y)-phi_depol0
                                                          /*+ m_elecSolver->getPhiSlow(x, y) - elec_depol0*/
                                                          );
        }
        // m_Esolver->solve_ls_fast_PhiE();//solve_ls_fast();
        m_Esolver->solve_ls_fast();//solve_ls_fast();
        updateEforPz();
#ifndef USE_MIRROR
        //for (int aa=0;aa<2;aa++)
        {
        //updateEforPz();
        //m_pzSolver->conduct(1,-1e-12,5);
        //  m_pzSolver->solvePz(5);
        m_pzSolver->solvePz_steady(10);
        }
#endif
    }
    //updateEforPz();
    //m_pzSolver->conduct(1,-1e-12,5);
    if (g_emitElectrons)
    {


    if (traj_table_calculated==0)
    {
        updateTrajTable();
        traj_table_calculated=1;
    }else
    {
        //if (rand()*1.0/RAND_MAX<0.5)
        updateTrajTable();
    }

    g_t+=dt_elec;
    g_save_time+=dt_elec;
    g_save_time2+=dt_elec;

    //for(int n=0; n<elecSteps;n++)

        electronEmissionEndMoveToElectrode(dt_elec);
    }
    /*if (g_emitElectrons)
        electronEmission(dt_elec);*/

#ifndef USE_MIRROR
    // m_pzSolver->step();
#endif
    /*for(int n=0; n<elecSteps;n++)
    {
        updateEforElec();
        m_elecSolver->step(dt_elec/elecSteps);
        electronExchange(dt_elec/elecSteps);
    }*/
#ifndef USE_MIRROR
       pzEmission(dt_elec);
#endif
}

void multiSolver::step()
{
    
}

void multiSolver::preparePz()
{
    
}

void multiSolver::getExtChargeField()
{
    
}

void multiSolver::electronExchange(double dt)
{
    for (int i=0;i<m_elecSolver->m_numParticles;i++)
    {
        double x,y;
        x=m_elecSolver->m_bodyPos[i].x;
        y=m_elecSolver->m_bodyPos[i].y;
        
        if (y<m_pzSolver->m_p[0].r.y+m_pzSolver->m_p[0].dl*0.5)
        {
            double vx = m_elecSolver->m_bodyVel[i].x;
            double vy = m_elecSolver->m_bodyVel[i].y;
            double xx = (m_pzSolver->m_p[0].r.y + m_pzSolver->m_p[0].dl * 0.5 - (y - x * vy / vx)) * vx / vy;
            //int p_n=(int) ((x-m_pzSolver->m_p[0].r.x)/m_pzSolver->m_dx);
            int p_n=(int) ((xx - m_pzSolver->m_p[0].r.x + 0.5 * m_pzSolver->m_dx) / m_pzSolver->m_dx);
            //int p_n=(int) ((x-m_pzSolver->m_p[0].r.x-m_pzSolver->m_dx)/m_pzSolver->m_dx);
            if ((p_n>=0)&&(p_n<m_pzSolver->m_p_num))
            {
                //if (m_pzSolver->m_p[p_n].q_ext+m_elecSolver->m_bodyPos[i].charge<1.7e+2)
                {
                    double coef = 1.0;
#ifdef USE_MIRROR
                    coef = (1.0 - ((eps_pz-1.0)/(eps_pz+1.0)));
#endif
                    m_pzSolver->m_p[p_n].q_ext+= coef*m_elecSolver->m_bodyPos[i].charge;
                    m_elecSolver->delete_particle(i);
                }
                /*else
                {
                    m_elecSolver->m_bodyPos[i].y+=2.0*(m_pzSolver->m_p[0].r.y+m_pzSolver->m_p[0].dl*0.5-m_elecSolver->m_bodyPos[i].y);
                    m_elecSolver->m_bodyVel[i].y=fabs(m_elecSolver->m_bodyVel[i].y);
                }*/
            }else
            {
                m_elecSolver->delete_particle(i);
            }
        }
    }
    /*for (int i=0;i<50;i++) {
        //m_pzSolver->m_p[i].q_ext+=1;
        printf("m_pzSolver->m_p[%d] = %f\n",i,m_pzSolver->m_p[i].q_ext);
    }*/
}

void multiSolver::pzEmission(double dt)
{
   /* int succ=0;
    for (int i=0;i<m_pzSolver->m_p_num;i++)
    {
        int q_extra=int(m_pzSolver->m_p[i].q_ext+m_pzSolver->m_p[i].q);
        // if (g_emitElectrons==false)
        // printf("i=%d Qextra=%d q=%f qext=%f \n",i,q_extra,m_pzSolver->m_p[i].q,m_pzSolver->m_p[i].q_ext);
        if ((q_extra>0)&&(m_pzSolver->m_p[i].q_ext>0.0))
        {
            if ((rand()*1.0/RAND_MAX)>0.9)
            {
                if ((q_extra/2)>0)
                {
                    m_pzSolver->m_p[i].q_ext-=q_extra/2;
                    m_elecSolver->create_pz_electron(m_pzSolver->m_p[i].r.x,m_pzSolver->m_p[i].r.y+m_pzSolver->m_p[i].dl*0.5+0.75e-6,q_extra/2);
                    succ++;
                }
            }
        }
    }*/

    double yy=m_pzSolver->m_p[0].r.y + m_pzSolver->m_p[0].dl * 0.5;
    for (int i=0;i<m_pzSolver->m_p_num;i++)
    {
        int q_extra=int(m_pzSolver->m_p[i].q_ext+m_pzSolver->m_p[i].q);

        if ((q_extra>1)&&(m_pzSolver->m_p[i].q_ext-m_pzSolver->m_p[i].q_0>1.0))
        {
            vec2 r;
            r.x=m_pzSolver->m_p[i].r.x;
            r.y=0.0;
            vec2 Ed = m_Esolver->getE(r.x,r.y);
            vec2 Epz = m_pzSolver->getEdepol(r.x,r.y);
            vec2 E_;
            E_.x = Ed.x + Epz.x;
            E_.y = Ed.y + Epz.y;
            if (E_.y<0)
            {


            vec2 v(0.0,0.0,0.0);
            double Dl=0.5e-6;//
            double Dt=1e-6;
            for (int j=0;j<100;j++)
            {
                Ed = m_Esolver->getE(r.x,r.y);
                Epz = m_pzSolver->getEdepol(r.x,r.y);
                E_;
                E_.x = Ed.x + Epz.x;
                E_.y = Ed.y + Epz.y;
                double magn=qe/Me;//1e-1;
                double a_=fmax(fabs(magn*(E_.x)),fabs(magn*(E_.y)));
                double v_=fmax(fabs(v.x),fabs(v.y));
                Dt=(sqrt(v_*v_+2*a_*Dl)-v_)/(a_+1e-20);

                v.x += magn*(E_.x)*Dt;
                v.y += magn*(E_.y)*Dt;
                r.x+=Dt*v.x;
                r.y+=Dt*v.y;
                if (r.y<yy) break;
            }
                    int p_n=0;

                    double xx = (yy - (r.y - r.x * v.y / v.x)) * v.x / v.y;
                    p_n=(int) ((xx - m_pzSolver->m_p[0].r.x + 0.5 * m_pzSolver->m_dx) / m_pzSolver->m_dx);

                    double qq= q_extra*(rand()*0.5/RAND_MAX);

                    m_pzSolver->m_p[i].q_ext-=qq;

                    if (p_n >= 0 && p_n < m_pzSolver->m_p_num)   m_pzSolver->m_p[p_n].q_ext+=qq;
                    if (p_n < 0) m_pzSolver->m_p[0].q_ext+=qq;
            }

        }
    }

}
