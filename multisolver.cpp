#include "multisolver.h"

multiSolver::multiSolver()
{
}

void multiSolver::updateEforPz()
{
    for (int i=0;i<m_pzSolver->m_p_num;i++)
    {
        double x,y;
        x=m_pzSolver->m_p[i].r.x;
        y=m_pzSolver->m_p[i].r.y;//+m_pzSolver->m_p[i].dl*0.4;

        vec2 E = m_Esolver->getE(x,y);
        vec2 Ep = m_pzSolver->getEdepol(x,y);
        vec2 Ee = m_elecSolver->getEe(x,y);
        m_pzSolver->m_p[i].E = E.y+/*E_global+*/Ep.y;
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
        m_elecSolver->m_bodyE[i].x=-(E.x+Ep.x);///10.0;
        m_elecSolver->m_bodyE[i].y=-(E.y+Ep.y);///10.0;
    }
}

void multiSolver::electronEmission(double d_t)
{
    double eMean=0.0;
    double emass=0.0;
    for (int i=0;i<m_Esolver->m_elec_num;i++)
    {
        if (m_Esolver->m_electrodes[i].canEmit)
        {
            vec2 E=m_Esolver->getE(m_Esolver->m_electrodes[i].r.x+1e-9,m_Esolver->m_electrodes[i].r.y);
            vec2 Ep=m_pzSolver->getEdepol(m_Esolver->m_electrodes[i].r.x+1e-9,m_Esolver->m_electrodes[i].r.y);
            double ex=E.x+Ep.x;
            double ey=E.y+Ep.y;
            double l=sqrt(ex*ex + ey*ey);
            eMean+=l;
            emass+=1.0;
            m_elecSolver->create_electron(m_Esolver->m_electrodes[i].r,l,d_t,m_Esolver->m_electrodes[i].dl*(w_z1-w_z0));
        }
    }
   // printf("curr_elec_num=%d mass=%f E=%e\n",m_elecSolver->m_numParticles,emass, eMean/emass);
}

void multiSolver::solve(int itn)
{


     double dt_elec=1.5e-16 * dtKoef;

    for (int i=0;i<m_Esolver->m_elec_num;i++)
    {
         m_Esolver->m_electrodes[i].phi_fix+=2.0*((m_Esolver->m_electrodes[i].phi_fix>0)-0.5)*0.0001*(dt_elec/1.5e-16);
      //  printf("i=%d phi=%e \n",m_Esolver->m_electrodes[i].phi_fix);
    }

printf("phi=%e \n",m_Esolver->m_electrodes[0].phi_fix);
    for (int nn=0;nn<itn;nn++)
    {
        m_pzSolver->get_q();

        double phi_depol0=m_pzSolver->getPhidepol(w_x0,w_y0);

        for (int i=0;i<m_Esolver->m_elec_num;i++)
        {
            double x,y;
            x=m_Esolver->m_electrodes[i].r.x;
            y=m_Esolver->m_electrodes[i].r.y;

            m_Esolver->m_electrodes[i].phi_fix_charges=(m_pzSolver->getPhidepol(x,y));//-phi_depol0);
        }
        m_Esolver->solve_ls_fast();
        updateEforPz();

       m_pzSolver->solvePz(5);

    }

    if (g_emitElectrons)
    electronEmission(dt_elec);

    printf("nele=%d \n",m_elecSolver->m_numParticles);
    updateEforElec();
    m_pzSolver->step();
    m_elecSolver->step(dt_elec);
    electronExchange(dt_elec);
    pzEmission(dt_elec);

    /* //debug below

    for (int nn=0;nn<itn;nn++)
    {
        m_Esolver->solvePhi(20);
        updateEforPz();

        m_pzSolver->solvePz(5);

    }
    double dt_elec=15e-11;

    m_pzSolver->step();*/

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
            int p_n=(int) ((x-m_pzSolver->m_p[0].r.x)/m_pzSolver->m_dx);
            if ((p_n>=0)&&(p_n<m_pzSolver->m_p_num))
            {
                //if (m_pzSolver->m_p[p_n].q_ext+m_elecSolver->m_bodyPos[i].charge<1.7e+2)
                {
                    m_pzSolver->m_p[p_n].q_ext+=m_elecSolver->m_bodyPos[i].charge;
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
                 if ((q_extra/4)>0)
                 {
                     m_pzSolver->m_p[i].q_ext-=q_extra/4;
                     m_elecSolver->create_pz_electron(m_pzSolver->m_p[i].r.x,m_pzSolver->m_p[i].r.y+m_pzSolver->m_p[i].dl*0.5+1.5e-9,q_extra/4);
                     succ++;
                 }
             }
         }
     }
     if (succ>0)
     printf("succ electrons created %d \n",succ);
}
