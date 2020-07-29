#ifndef MULTISOLVER_H
#define MULTISOLVER_H


#include "efieldlagrangian.h"
#include "pzsolver.h"
#include "electronlagrangian.h"

class multiSolver
{   
public:

    eFieldLagrangian* m_Esolver;
    pzSolver* m_pzSolver;
    electronLagrangian* m_elecSolver;



    multiSolver();
    void updateEforPz();
    void updateEforElec();
    void electronEmission(double dt);
    void solve(int itn);
    void step();
    void preparePz();
    void getExtChargeField();
    void electronExchange(double dt);
    void pzEmission(double dt);
};

#endif // MULTISOLVER_H
