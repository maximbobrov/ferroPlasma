#include "multisolver.h"

multiSolver::multiSolver()
{
}

void multiSolver::updateEforPz()
{

    //double phi_down=m_Esolver->m_electrodes[m_Esolver->m_elec_num-1].phi_fix; //phi at lower electrode
    // printf("phi_down= %f \n",phi_down);

    double phi_depol0=m_pzSolver->getPhidepol(w_x0,w_y0);
    double phi_e0=m_elecSolver->getPhiSlow(w_x0,w_y0);

    for (int i=0;i<m_pzSolver->m_p_num;i++)
    {
        //by Efield direct calculation
        double x,ym, y0,yp;
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


        /*
        //by PHI calculation for some reason its a bit underestmating the field values
        double x,yp;
        x=m_pzSolver->m_p[i].r.x;
        yp=m_pzSolver->m_p[i].r.y+m_pzSolver->m_p[i].dl*0.5;

        double phi_er = m_Esolver->getPhi(x,yp);
        double phi_pz = m_pzSolver->getPhidepol(x,yp);
        double phi_beta = m_elecSolver->getPhiSlow(x,yp);
        double phi_up = (phi_er + (phi_pz-phi_depol0) + (phi_beta-phi_e0));
        m_pzSolver->m_p[i].E = -(phi_up - phi_down)/m_pzSolver->m_p[i].dl;
*/
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
    double E0=m_elecSolver->getEmult_dipole(2.0e-6);
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

void multiSolver::solve(int itn)
{
    dt_elec = 1.5e-13 * dtKoef;
    static int elecSteps=10;
    
    g_t+=dt_elec;
    g_save_time+=dt_elec;
    g_save_time2+=dt_elec;

    //double tstart = get_time();
    //double t1,t2;
    /* for (int i=0;i<m_Esolver->m_elec_num && fabs(m_Esolver->m_electrodes[i].phi_fix) < g_phi;i++)
    {
        m_Esolver->m_electrodes[i].phi_fix +=2.0*((m_Esolver->m_electrodes[i].phi_fix>0)-0.5)*0.005*(dt_elec/0.5e-13);
        //  printf("i=%d phi=%e \n",m_Esolver->m_electrodes[i].phi_fix);
    }
*/
    //printf("phi=%e t=%f dt=%e\n",m_Esolver->m_electrodes[0].phi_fix, (g_t * 1e6),dt_elec);

    for (int nn=0;nn<itn;nn++)
    {
#ifndef USE_MIRROR
        m_pzSolver->get_q();
#endif


        double phi_depol0=m_pzSolver->getPhidepol(w_x0,w_y0);
        double elec_depol0=m_elecSolver->getPhiSlow(w_x0,w_y0);

        //double phi_fromCharges0 = m_Esolver->getPhiFromCharges(w_x0,w_y0);

        for (int i=0;i<m_Esolver->m_elec_num;i++)
        {
            double x,y;
            x=m_Esolver->m_electrodes[i].r.x;
            y=m_Esolver->m_electrodes[i].r.y;
            
            m_Esolver->m_electrodes[i].phi_fix_charges = (m_pzSolver->getPhidepol(x,y)-phi_depol0
                                                          + m_elecSolver->getPhiSlow(x, y) - elec_depol0
                                                          );
        }


        // m_Esolver->solve_ls_fast_PhiE();//solve_ls_fast();
        m_Esolver->solve_ls_fast();//solve_ls_fast();
#ifndef USE_MIRROR
        updateEforPz();
        m_pzSolver->solvePz(5);
#endif

    }


    if (g_emitElectrons)
        electronEmission(dt_elec);

#ifndef USE_MIRROR
    m_pzSolver->step();
#endif
    for(int n=0; n<elecSteps;n++)
    {
        updateEforElec();
        m_elecSolver->step(dt_elec/elecSteps);
        electronExchange(dt_elec/elecSteps);
    }
#ifndef USE_MIRROR
    pzEmission(dt_elec);
#endif

    /* //debug below
      
    for (int nn=0;nn<itn;nn++)
    {
        m_Esolver->solvePhi(20);
        updateEforPz();
        
        m_pzSolver->solvePz(5);
        
    }
    double dt_elec=15e-11;
    
    m_pzSolver->step();*/
    
    //update E for visualization purposes:
    /*double phi_depol0=m_pzSolver->getPhidepol(w_x0,w_y0);
    double elec_depol0=m_elecSolver->getPhiSlow(w_x0,w_y0);
  //t1 = get_time();

    for (int i=0;i<m_Esolver->m_elec_num;i++)
    {
        double x,y;
        x=m_Esolver->m_electrodes[i].r.x;
        y=m_Esolver->m_electrodes[i].r.y;
        m_Esolver->m_electrodes[i].phi_fix_charges=  (m_pzSolver->getPhidepol(x,y)-phi_depol0
                                                      + m_elecSolver->getPhiSlow(x, y) - elec_depol0
                                                      );
    }
    //t2 = get_time();
    m_Esolver->solve_ls_fast();*/

    //double tend = get_time();
    //printf("fullT = %e t2 = %e\n",(tend-tstart), (t2-t1));
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
    int succ=0;
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
    }
    /* if (succ>0)
        printf("succ electrons created %d \n",succ);*/
}
