#ifndef MULTISOLVER_H
#define MULTISOLVER_H


#include "efieldlagrangian.h"
#include "pzsolver.h"
#include "electronlagrangian.h"

class multiSolver
{   
public:
    double dt_elec;
    eFieldLagrangian* m_Esolver;
    pzSolver* m_pzSolver;
    electronLagrangian* m_elecSolver;

    int endPosTable[1000];



    multiSolver();
    void checkPotential();
    void updateEforPz();
    void updateEforElec();
    void electronEmission(double dt);
    void electronEmissionEndMoveToElectrode(double dt);
    void updateTrajTable();
    int  getEndPos(int i);
    void solve(int itn);
    void init();
    void step();
    void preparePz();
    void getExtChargeField();
    void electronExchange(double dt);
    void pzEmission(double dt);
};

#endif // MULTISOLVER_H
