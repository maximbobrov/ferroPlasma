#ifndef PZSOLVER_H
#define PZSOLVER_H

#include "globals.h"
#include "lagrangiansolver.h"

class pzSolver : lagrangianSolver
{
public:
    struct pElem   //linear electrode elem
    {
        double p,p_prev, ds, dl,E,E_prev,Ex_s,Ey_s,RHS;//C/m2, m, V
      //  vec2 r;
        vec2 r_top;
        double q, q_ext,q_tmp,q_0; //q is dipolar charge; q_ext is the attached charge
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
    void init();
    void solvePz(int itn);
    void solvePz_steady(int itn);
    void getRHS();
    void get_q();
    void step();
    vec2 getEdepol(double x, double y);
    double getPhidepol(double x, double y);
    void updateCharge();

    void setWallPos(double a);

    void conduct(double sigma, double dt, int itn);

    double getE_self(int j);
   // static vec2 getEField(const vec2& iCenterPos, const vec2& iFarPos);
   // static vec2 getPhiField(const vec2& iCenterPos, const vec2& iFarPos);
};

#endif // PZSOLVER_H
