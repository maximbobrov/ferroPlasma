#ifndef PZSOLVER_H
#define PZSOLVER_H

#include "globals.h"
class pzSolver
{
public:
    struct pElem   //linear electrode elem
    {
        double p, ds, dl,E;//C/m2, m, V
        vec3<double> r;
        double q1,q2;
    };
public:
    poly m_poly;

    INPUT_PARAM m_par;

    double m_dx;
    double m_dt;
    pElem *m_p;
    int m_p_num;
    void solve_poly(int itn);
    void get_q();

    pzSolver();

};

#endif // PZSOLVER_H
