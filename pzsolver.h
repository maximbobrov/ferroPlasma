#ifndef PZSOLVER_H
#define PZSOLVER_H

#include "globals.h"
#include "lagrangiansolver.h"

class pzSolver : lagrangianSolver
{
public:
    struct pElem   //linear electrode elem
    {
        double p,p_prev, ds, dl,E,E_prev,RHS;//C/m2, m, V
        vec2 r;
        double q, q_ext; //q is dipolar charge; q_ext is the attached charge
    };

public:
    poly m_poly;
    vec2* m_rCentre;
    INPUT_PARAM m_par;

    double m_dx;
    double m_dt;
    double kappa;
    pElem *m_p;
    int m_p_num;
    pzSolver();
    void solvePz(int itn);
    void getRHS();
    void get_q();
    void step();
    vec2 getEdepol(double x, double y);
    double getPhidepol(double x, double y);
    void updateCharge();
    void updateGridProp();
    void setWallPos(double a);

   // static vec2 getEField(const vec2& iCenterPos, const vec2& iFarPos);
   // static vec2 getPhiField(const vec2& iCenterPos, const vec2& iFarPos);
};

#endif // PZSOLVER_H
