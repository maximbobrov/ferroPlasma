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
        y=m_pzSolver->m_p[i].r.y+m_pzSolver->m_p[i].dl*0.4;

        vec3<double> E=m_Esolver->getE(x,y);
        vec3<double> Ep=m_pzSolver->getEdepol(x,y);
        m_pzSolver->m_p[i].E=E.y+Ep.y;


    }
}

void multiSolver::updateEforElec()
{
    for (int i=0;i<m_elecSolver->m_numParticles;i++)
    {
        double x,y;
        x=m_elecSolver->m_bodyPos[i].x;
        y=m_elecSolver->m_bodyPos[i].y;

        vec3<double> E=m_Esolver->getE(x,y);
        vec3<double> Ep=m_pzSolver->getEdepol(x,y);
        m_elecSolver->m_bodyE[i].x=E.x+Ep.x;
        m_elecSolver->m_bodyE[i].y=E.y+Ep.y;


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
        vec3<double> E=m_Esolver->getE(m_Esolver->m_electrodes[i].r1.x+1e-9,m_Esolver->m_electrodes[i].r1.y);
        double l=sqrt(E.x*E.x +E.y*E.y);
        eMean+=l;
        emass+=1.0;
       m_elecSolver->create_electron(m_Esolver->m_electrodes[i].r1,l,d_t,m_Esolver->m_electrodes[i].dl*m_Esolver->m_dz);
        }

    }
    printf("curr_elec_num=%d mass=%f E=%e\n",m_elecSolver->m_numParticles,emass, eMean/emass);
}

void multiSolver::solve(int itn)
{
    for (int nn=0;nn<itn;nn++)
    {
        m_pzSolver->get_q();
        for (int i=0;i<m_Esolver->m_elec_num;i++)
        {
            double x,y;
            x=(m_Esolver->m_electrodes[i].r0.x+m_Esolver->m_electrodes[i].r1.x)*0.5;
            y=(m_Esolver->m_electrodes[i].r0.y+m_Esolver->m_electrodes[i].r1.y)*0.5;

            m_Esolver->m_electrodes[i].phi_fix_charges=m_pzSolver->getPhidepol(x,y);
        }
        m_Esolver->solvePhi(20);
        updateEforPz();

        m_pzSolver->solvePz(5);

    }
    double dt_elec=15e-11;
    electronEmission(dt_elec*10);
    updateEforElec();
    m_pzSolver->step();
    m_elecSolver->step(dt_elec);
    electronExchange(dt_elec);
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
                if (m_pzSolver->m_p[p_n].q_ext+m_elecSolver->m_bodyPos[i].w<1.7e+2)
                {
                    m_pzSolver->m_p[p_n].q_ext+=m_elecSolver->m_bodyPos[i].w;
                    m_elecSolver->delete_particle(i);
                }else
                {
                    m_elecSolver->m_bodyPos[i].y+=2.0*(m_pzSolver->m_p[0].r.y+m_pzSolver->m_p[0].dl*0.5-m_elecSolver->m_bodyPos[i].y);
                    m_elecSolver->m_bodyVel[i].y=fabs(m_elecSolver->m_bodyVel[i].y);
                }

            }else
            {
                m_elecSolver->delete_particle(i);
            }


        }



    }

}
