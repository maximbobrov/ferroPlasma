#ifndef MULTISOLVER_H
#define MULTISOLVER_H

#include "sse_sum.h"
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
    void solve(int itn);
    void step();
    void preparePz();
    void getExtChargeField();

};

#endif // MULTISOLVER_H
