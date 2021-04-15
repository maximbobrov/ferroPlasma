#include "multisolver.h"

double to_elec_from_elec[900][900];
double to_elec_from_pz[900][900];
double to_pzUp_from_elec[900][900];
double to_pzUp_from_pz[900][900];

double to_pzDown_from_elec[900][900];
double to_pzDown_from_pz[900][900];



multiSolver::multiSolver()
{
    dt_elec=1e-20;
}

void multiSolver::updateEforPz()
{

    //double phi_down=m_Esolver->m_electrodes[m_Esolver->m_elec_num-1].phi_fix; //phi at lower electrode
    // printf("phi_down= %f \n",phi_down);

    double phi_depol0=m_pzSolver->getPhidepol(w_x0,w_y0);

    double zl=1.0/m_pzSolver->m_p[0].dl;
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


        /* double x,yp,ym;
        x=m_pzSolver->m_p[i].r.x;
        yp=m_pzSolver->m_p[i].r.y+m_pzSolver->m_p[i].dl*0.5;
        ym=m_pzSolver->m_p[i].r.y-m_pzSolver->m_p[i].dl*0.5;

        double phi_er = m_Esolver->getPhi(x,yp);
        double phi_pz = m_pzSolver->getPhidepol(x,yp);

        double phi_up = (phi_er + (phi_pz-phi_depol0));

        phi_er = m_Esolver->getPhi(x,ym);
        phi_pz = m_pzSolver->getPhidepol(x,ym);

        double phi_down = (phi_er + (phi_pz-phi_depol0));
        */

        //double delta=1e-6; //this should be monitored
        //  m_p[i].q=(m_p[i].p)*m_p[i].ds/qe;

        //q=qe/(eps0*pi2) * (m_p[i].q+m_p[i].q_ext);
        // double q_self=(m_pzSolver->m_p[i].p)*m_pzSolver->m_p[i].ds/(eps0*pi2);
        //double phi_self=-q_self*log(delta)/(w_z1 - w_z0);
        //double E_self=-phi_self/m_pzSolver->m_p[i].dl; //extract

        double phi_up=getPhi_at_pz_up(i);
        double phi_down=getPhi_at_pz_down(i);

        m_pzSolver->m_p[i].E = -(phi_up - phi_down)*zl;// - E_self;

        //now Ex_s Ey_s at the surface:

        /*vec2 Em = m_Esolver->getE(x,yp);
        vec2 Epm = m_pzSolver->getEdepol(x,yp);
        m_pzSolver->m_p[i].Ex_s =Em.x+Epm.x;
        m_pzSolver->m_p[i].Ey_s =Em.y+Epm.y;*/


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
    double electrons_in_pack=200;
    for (int i=0;i<m_Esolver->m_elec_num-1;i++)
    {
        if (m_Esolver->m_electrodes[i].canEmit)
        {
            vec2 E=m_Esolver->getE(m_Esolver->m_electrodes[i].r.x,m_Esolver->m_electrodes[i].r.y);
            vec2 Ep=m_pzSolver->getEdepol(m_Esolver->m_electrodes[i].r.x,m_Esolver->m_electrodes[i].r.y);
            double ex=E.x+Ep.x;
            double ey=E.y+Ep.y;

            /*vec2 E=get_fast_E();

            double ex=E.x;
            double ey=E.y;*/

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
        bool inArea = true;
        double Dl=0.5e-7;
        double Dt=1e-6;
        for (int j=0;j<1000;j++)
        {
            /*vec2 Ed = m_Esolver->getE(r.x,r.y);
            vec2 Epz = m_pzSolver->getEdepol(r.x,r.y);*/
            vec2 E_=get_slower_E(r.x,r.y);
            /*E_.x = Ed.x + Epz.x;
            E_.y = Ed.y + Epz.y;*/


            double magn=qe/Me;//1e-1;
            double a_=fmax(fabs(magn*(E_.x)),fabs(magn*(E_.y)));
            double v_=fmax(fabs(v.x),fabs(v.y));
            Dt=(sqrt(v_*v_+2*a_*Dl)-v_)/(a_+1e-20);

            v.x += magn*(E_.x)*Dt;
            v.y += magn*(E_.y)*Dt;
            if (r.y+Dt*v.y<0) break;
            r.x+=Dt*v.x;
            r.y+=Dt*v.y;
            if(r.y > w_y1){
                inArea = false;
                break;
            }
        }

        int p_n=0;
        double xx = (m_pzSolver->m_p[0].r_top.y - (r.y - r.x * v.y / v.x)) * v.x / v.y;

        p_n=(int) ((xx - m_pzSolver->m_p[0].r_top.x + 0.5 * m_pzSolver->m_dx) / m_pzSolver->m_dx);
        if(inArea && p_n >=0 && p_n<m_pzSolver->m_p_num)
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
    for (int i=0;i<m_Esolver->m_elec_num-1;i+=1)
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
    m_pzSolver->get_q();

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
        double phi_depol0=m_pzSolver->getPhidepol(w_x0,w_y0);

        for (int i=0;i<m_Esolver->m_elec_num;i++)
        {
            /*            double x,y;
            x=m_Esolver->m_electrodes[i].r.x;
            y=m_Esolver->m_electrodes[i].r.y;
            m_Esolver->m_electrodes[i].phi_fix_charges = (m_pzSolver->getPhidepol(x,y)-phi_depol0);*/

            m_Esolver->m_electrodes[i].phi_fix_charges = getPhi_at_electrode(i)-phi_depol0;

        }
        m_Esolver->solve_ls_fast();
        for (int aa=0;aa<3;aa++)
        {
            m_pzSolver->get_q();
            updateEforPz();

            m_pzSolver->solvePz(4);
            m_pzSolver->step();
            // m_pzSolver->solvePz_steady(10);
        }
    }
    //double t11 = get_time();

    g_t+=dt_elec;
    g_save_time+=dt_elec;
    g_save_time2+=dt_elec;
    if (g_emitElectrons)
        electronEmissionEndMoveToElectrode(dt_elec);

    pzEmission(dt_elec);
    //pzEmissionHoriz(dt_elec);

    // double t1 = get_time();
    //printf("tall= %e t=%e\n", t1-t0, t11-t10);
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
        
        if (y<m_pzSolver->m_p[0].r_top.y)
        {
            double vx = m_elecSolver->m_bodyVel[i].x;
            double vy = m_elecSolver->m_bodyVel[i].y;
            double xx = (m_pzSolver->m_p[0].r_top.y - (y - x * vy / vx)) * vx / vy;
            //int p_n=(int) ((x-m_pzSolver->m_p[0].r.x)/m_pzSolver->m_dx);
            int p_n=(int) ((xx - m_pzSolver->m_p[0].r_top.x + 0.5 * m_pzSolver->m_dx) / m_pzSolver->m_dx);
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
    double yy=m_pzSolver->m_p[0].r_top.y;
    for (int i=0;i<m_pzSolver->m_p_num;i++)
    {
        if (m_pzSolver->m_p[i].q_ext-m_pzSolver->m_p[i].q_0>1.0)
        {
            vec2 r;
            r.x=m_pzSolver->m_p[i].r_top.x;
            r.y=0.0;
            vec2 Ed = m_Esolver->getE(r.x,r.y);
            vec2 Epz = m_pzSolver->getEdepol(r.x,r.y);
            vec2 E_=get_slower_E(r.x,r.y);
            double E0y = E_.y;
            if (E_.y<0)
            {
                vec2 v(0.0,0.0,0.0);
                double Dl=0.5e-7;//
                double Dt=1e-6;
                for (int j=0;j<200;j++)
                {
                    E_ = get_slower_E(r.x,r.y);
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
                p_n=(int) ((xx - m_pzSolver->m_p[0].r_top.x + 0.5 * m_pzSolver->m_dx) / m_pzSolver->m_dx);

                double qepspi = (qe/(eps0*pi2))/(w_z1 - w_z0);
                double delta=1e-6;
                double d2=delta*delta;
                double r2;
                double dy;
                dy = m_pzSolver->m_p[i].r_top.y;
                r2=(dy*dy);
                double qq =min(fabs(fabs(E0y/dy*(r2+d2)/qepspi)), m_pzSolver->m_p[i].q_ext - m_pzSolver->m_p[i].q_0);
                qq = 0.005 * (m_pzSolver->m_p[i].q_ext  - m_pzSolver->m_p[i].q_0);

                if(m_pzSolver->m_p[i].q_ext - qq - m_pzSolver->m_p[i].q_0>0){
                    m_pzSolver->m_p[i].q_ext -= qq;
                    if (p_n >= 0 && p_n < m_pzSolver->m_p_num)   m_pzSolver->m_p[p_n].q_ext += qq;
                    //if (p_n < 0) m_pzSolver->m_p[0].q_ext += qq;
                }else
                {
                    qq =  (m_pzSolver->m_p[i].q_ext - m_pzSolver->m_p[i].q_0);
                    m_pzSolver->m_p[i].q_ext -= qq;
                    if (p_n >= 0 && p_n < m_pzSolver->m_p_num)   m_pzSolver->m_p[p_n].q_ext += qq;
                    //if (p_n < 0) m_pzSolver->m_p[0].q_ext += qq;

                }
            }
        }
    }

/*
    for (int i=0;i<m_pzSolver->m_p_num;i++) //drift and diffusion
    {
        m_pzSolver->m_p[i].q_tmp=m_pzSolver->m_p[i].q_ext;
    }
    for (int itn=0;itn<2;itn++)
    {
        //m_pzSolver->m_p[0].q_ext =m_pzSolver->m_p[1].q_ext;

        vec2 r;
        r.x=m_pzSolver->m_p[1].r_top.x;
        r.y=yy;
        vec2 Ed = m_Esolver->getE(r.x,r.y);
        vec2 Epz = m_pzSolver->getEdepol(r.x,r.y);
        vec2 E_p;
        E_p.x = Ed.x + Epz.x;
        E_p.y = Ed.y + Epz.y;

        double a = 0.001;
        double b = 1e-11;

        m_pzSolver->m_p[0].q_ext=(m_pzSolver->m_p[0].q_tmp + a*(m_pzSolver->m_p[1].q_ext)  + b*(E_p.x*(m_pzSolver->m_p[1].q_ext-m_pzSolver->m_p[1].q_0)))/(1+2*a);

        for (int i=0;i<m_pzSolver->m_p_num;i++)
        {
            vec2 r;
            r.x=m_pzSolver->m_p[i+1].r_top.x;
            r.y=yy;
            vec2 Ed = m_Esolver->getE(r.x,r.y);
            vec2 Epz = m_pzSolver->getEdepol(r.x,r.y);
            vec2 E_p;
            E_p.x = Ed.x + Epz.x;
            E_p.y = Ed.y + Epz.y;

            r.x=m_pzSolver->m_p[i-1].r_top.x;
            r.y=yy;
            Ed = m_Esolver->getE(r.x,r.y);
            Epz = m_pzSolver->getEdepol(r.x,r.y);
            vec2 E_m;
            E_m.x = Ed.x + Epz.x;
            E_m.y = Ed.y + Epz.y;


            m_pzSolver->m_p[i].q_ext=(m_pzSolver->m_p[i].q_tmp+ a*(m_pzSolver->m_p[i+1].q_ext+m_pzSolver->m_p[i-1].q_ext)  + b*(E_p.x*(m_pzSolver->m_p[i+1].q_ext - m_pzSolver->m_p[i+1].q_0) - E_m.x*(m_pzSolver->m_p[i-1].q_ext-m_pzSolver->m_p[i-1].q_0)))/(1+2*a);

        }
        m_pzSolver->m_p[m_pzSolver->m_p_num-1].q_ext = m_pzSolver->m_p[m_pzSolver->m_p_num-2].q_ext;
    }*/
}

void multiSolver::pzEmissionHoriz(double dt)
{

    double yy=m_pzSolver->m_p[0].r_top.y;
    for (int i=0;i<m_pzSolver->m_p_num;i++)
    {
        if ((m_pzSolver->m_p[i].q_ext-m_pzSolver->m_p[i].q_0>1.0))
        {
            vec2 r;
            r.x=m_pzSolver->m_p[i].r_top.x;
            r.y=0.0;
            vec2 Ed = m_Esolver->getE(r.x,r.y);
            vec2 Epz = m_pzSolver->getEdepol(r.x,r.y);
            vec2 E_;
            E_.x = Ed.x + Epz.x;
            E_.y = Ed.y + Epz.y;
            double Emin_x = 10000;
            if (fabs(E_.x)>Emin_x)
            {
                int p_n=E_.x > 0 ? i-1 : i+1;

                double qepspi = (qe/(eps0*pi2))/(w_z1 - w_z0);
                double delta=1e-6;
                double d2=delta*delta;
                double r2;
                double dy,dx;
                dx = m_pzSolver->m_p[1].r_top.x - m_pzSolver->m_p[0].r_top.x;
                dy = m_pzSolver->m_p[i].r_top.y;
                r2=(dx*dx + dy*dy);
                double qq = min(fabs(0.01 * fabs(E_.x/dx*(r2+d2)/qepspi)), m_pzSolver->m_p[i].q_ext - m_pzSolver->m_p[i].q_0);

                printf("qq = %f\n", qq);
                //double qq= q_extra*(rand()*0.5/RAND_MAX);
                if(m_pzSolver->m_p[i].q_ext - m_pzSolver->m_p[i].q_0>0){

                    m_pzSolver->m_p[i].q_ext-=qq;
                    if (p_n >= 0 && p_n < m_pzSolver->m_p_num)   m_pzSolver->m_p[p_n].q_ext+=qq;
                    if (p_n < 0) m_pzSolver->m_p[0].q_ext+=qq;
                }
            }

        }
    }

}

void multiSolver::prepare_caches()
{
    double  x,y;
    static double qepspi = qe/(eps0*pi2);
    double r;
    double dx,dy;
    double delta=1e-6;

    for (int n=0;n<m_Esolver->m_elec_num;n++)
    {
        x=m_Esolver->m_electrodes[n].r.x;
        y=m_Esolver->m_electrodes[n].r.y;

        for (int i=0;i<m_Esolver->m_chargeNum;i++)
        {
            dx = m_Esolver->m_charges[i].x - x;
            dy =  m_Esolver->m_charges[i].y - y;
            r=sqrt(dx*dx+dy*dy);

            to_elec_from_elec[n][i]=-qepspi*log(r+delta)/(w_z1 - w_z0);
        }

        for (int i=0;i<m_pzSolver->m_p_num;i++)
        {
            dx = m_pzSolver->m_p[i].r_top.x - x;
            dy = m_pzSolver->m_p[i].r_top.y - y;
            r=sqrt(dx*dx+dy*dy);
            to_elec_from_pz[n][i]=-qepspi*log(r+delta)/(w_z1 - w_z0);
        }
    }

    for (int n=0;n<m_pzSolver->m_p_num;n++)
    {
        x=m_pzSolver->m_p[n].r_top.x;
        y=m_pzSolver->m_p[n].r_top.y;

        for (int i=0;i<m_Esolver->m_chargeNum;i++)
        {
            dx = m_Esolver->m_charges[i].x - x;
            dy =  m_Esolver->m_charges[i].y - y;
            r=sqrt(dx*dx+dy*dy);

            to_pzUp_from_elec[n][i]=-qepspi*log(r+delta)/(w_z1 - w_z0);
        }

        for (int i=0;i<m_pzSolver->m_p_num;i++)
        {
            dx = m_pzSolver->m_p[i].r_top.x - x;
            dy = m_pzSolver->m_p[i].r_top.y - y;
            r=sqrt(dx*dx+dy*dy);
            to_pzUp_from_pz[n][i]=-qepspi*log(r+delta)/(w_z1 - w_z0);
        }
    }

    for (int n=0;n<m_pzSolver->m_p_num;n++)
    {
        x=m_pzSolver->m_p[n].r_top.x;
        y=m_pzSolver->m_p[n].r_top.y - m_pzSolver->m_p[n].dl;

        for (int i=0;i<m_Esolver->m_chargeNum;i++)
        {
            dx = m_Esolver->m_charges[i].x - x;
            dy =  m_Esolver->m_charges[i].y - y;
            r=sqrt(dx*dx+dy*dy);

            to_pzDown_from_elec[n][i]=-qepspi*log(r+delta)/(w_z1 - w_z0);
        }

        for (int i=0;i<m_pzSolver->m_p_num;i++)
        {
            dx = m_pzSolver->m_p[i].r_top.x - x;
            dy = m_pzSolver->m_p[i].r_top.y - y;
            r=sqrt(dx*dx+dy*dy);
            to_pzDown_from_pz[n][i]=-qepspi*log(r+delta)/(w_z1 - w_z0);
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
    for (int i=0;i<m_pzSolver->m_p_num;i++)
    {
        sum+=(m_pzSolver->m_p[i].q+m_pzSolver->m_p[i].q_ext)*to_elec_from_pz[n][i];
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
    for (int i=0;i<m_pzSolver->m_p_num;i++)
    {
        sum+=(m_pzSolver->m_p[i].q+m_pzSolver->m_p[i].q_ext)*to_pzUp_from_pz[n][i];
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
    for (int i=0;i<m_pzSolver->m_p_num;i++)
    {
        sum+=(m_pzSolver->m_p[i].q+m_pzSolver->m_p[i].q_ext)*to_pzDown_from_pz[n][i];
    }
    return sum;
}


//for fast fields
double g_x_min=1e10;
double g_x_max=-1e10;
double g_y_min=1e10;
double g_y_max=-1e10;
double g_dx,g_dy;
vec2 * g_grid[50][50][200];
int g_grid_num[50][50];
vec2 g_grid_Ephi[51][51];
vec2 g_grid_Ephi_near[51][51][2][2]; //last two indices are the in-cell coordinates

int g_Nx=40;
int g_Ny=20;

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
            g_grid_num[i][j]=0;
        }
    }
    for (int i=0;i<m_Esolver->m_chargeNum;i++)
    {
        int ni,nj;
        ni=(int)((m_Esolver->m_charges[i].x-g_x_min)/g_dx);
        nj=(int)((m_Esolver->m_charges[i].y-g_y_min)/g_dy);

        g_grid[ni][nj][g_grid_num[ni][nj]] = &(m_Esolver->m_charges[i]);
        g_grid_num[ni][nj]++;
    }
    for (int i=0;i<m_pzSolver->m_p_num;i++)
    {
        int ni,nj;
        ni=(int)((m_pzSolver->m_p[i].r_top.x-g_x_min)/g_dx);
        nj=(int)((m_pzSolver->m_p[i].r_top.y-g_y_min)/g_dy);

        g_grid[ni][nj][g_grid_num[ni][nj]] = &(m_pzSolver->m_p[i].r_top);
        g_grid_num[ni][nj]++;
    }

    int sum1=0;
    for (int i=0;i<g_Nx;i++)
        for (int j=0;j<g_Ny;j++)
        {
            if (g_grid_num[i][j]>0)
                printf("i=%d j=%d gnum=%d \n",i,j,g_grid_num[i][j]);
            sum1+=g_grid_num[i][j];
        }

    printf("sum1=%d sum2=%d \n",sum1,m_pzSolver->m_p_num+m_Esolver->m_chargeNum);

}

void multiSolver::fast_Fields_recalculate()
{
    for (int i=0;i<=g_Nx;i++)
    {
        for (int j=0;j<=g_Ny;j++)
        {
            double x,y;
            x=g_x_min+i*g_dx;
            y=g_y_min+j*g_dy;
            vec2 Ee = m_Esolver->getE(x,y);
            vec2 Epz = m_pzSolver->getEdepol(x,y); //можно ускорить через кэш!
            g_grid_Ephi[i][j].x = Ee.x +  Epz.x;
            g_grid_Ephi[i][j].y = Ee.y +  Epz.y;
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
    double delta=1e-6;
    double d2=delta*delta;

    double x00,y00, x10,y01;
    x00=g_x_min+i_*g_dx;
    y00=g_y_min+j_*g_dy;

    x10=g_x_min+(i_+1)*g_dx;
    y01=g_y_min+(j_+1)*g_dy;


    for (int i=im;i<=ip;i++)
    {
        for (int j=jm;j<=jp;j++)
        {
            for (int n=0; n<g_grid_num[i][j]; n++)
            {
                double r2;
                double q;
                double dx,dy;
                q=- qepspi *g_grid[i][j][n]->charge;

                dx = g_grid[i][j][n]->x - x00;
                dy = g_grid[i][j][n]->y - y00;
                r2=(dx*dx + dy*dy);
                E00.x+=dx*q/(r2+d2);
                E00.y+=dy*q/(r2+d2);

                dx = g_grid[i][j][n]->x - x10;
                dy = g_grid[i][j][n]->y - y00;
                r2=(dx*dx + dy*dy);
                E10.x+=dx*q/(r2+d2);
                E10.y+=dy*q/(r2+d2);

                dx = g_grid[i][j][n]->x - x00;
                dy = g_grid[i][j][n]->y - y01;
                r2=(dx*dx + dy*dy);
                E01.x+=dx*q/(r2+d2);
                E01.y+=dy*q/(r2+d2);

                dx = g_grid[i][j][n]->x - x10;
                dy = g_grid[i][j][n]->y - y01;
                r2=(dx*dx + dy*dy);
                E11.x+=dx*q/(r2+d2);
                E11.y+=dy*q/(r2+d2);

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

    double qepspi = (qe/(eps0*pi2))/(w_z1 - w_z0);
    double delta=1e-6;
    double d2=delta*delta;

    for (int i=im;i<=ip;i++)
    {
        for (int j=jm;j<=jp;j++)
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

    return res;
}


vec2 multiSolver::get_slow_E(double x, double y)
{

    double qepspi = (qe/(eps0*pi2))/(w_z1 - w_z0);
    double delta=1e-6;
    double d2=delta*delta;
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

    return res;
}

