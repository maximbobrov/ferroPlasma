#include "multisolver.h"

multiSolver::multiSolver()
{
}

void multiSolver::updateEforPz()
{
    for (int i=0;i<m_pzSolver->m_p_num;i++)
    {
        vec3<double> E=m_Esolver->getE(m_pzSolver->m_p[i].r.x,m_pzSolver->m_p[i].r.y);
        vec3<double> Ep=m_pzSolver->getEdepol(m_pzSolver->m_p[i].r.x,m_pzSolver->m_p[i].r.y);
        m_pzSolver->m_p[i].E=E.y+Ep.y;
    }
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
        m_Esolver->solvePhi(50);
        updateEforPz();
        m_pzSolver->solvePz(5);
    }
    m_pzSolver->step();
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
