#ifndef PZSOLVER_H
#define PZSOLVER_H

//#include "globals.h"
#include "sse_sum.h"
class pzSolver
{
public:
    struct pElem   //linear electrode elem
    {
        double p,p_prev, ds, dl,E,E_prev,RHS;//C/m2, m, V
        vec3<double> r;
        double q1,q2;
    };
public:
    poly m_poly;

    INPUT_PARAM m_par;

    double m_dx;
    double m_dt;
    double kappa;
    pElem *m_p;
    int m_p_num;
    void solvePz(int itn);
    void getRHS();
    void get_q();


    pzSolver();

};

#endif // PZSOLVER_H
