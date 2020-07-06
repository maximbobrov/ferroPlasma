#ifndef EFIELDLAGRANGIAN_H
#define EFIELDLAGRANGIAN_H
#include "globals.h"
#include "lagrangiansolver.h"

class eFieldLagrangian : lagrangianSolver
{

public:
    struct eElem   //linear electrode elem
    {
        double rho1, rho2,dl,phi_ext,phi_fix,phi_fix_charges;//C/m2, m, V
        vec2 r0,r1;

        bool canEmit;
    };
    vec2* m_rCentre;
    eElem * m_electrodes;
    int m_elec_num;
    double m_dz;

public:

    eFieldLagrangian();
    void updatePhi(); //get external potential on electrodes
    void solvePhi(int itn); //update potential
    vec2 getE(double x,double y);
    double getPhi(double x,double y);
    void setElectrodeAngle(double deg);
    void updateGridProp();
    static vec2 getEField(const vec2& iPos1, const vec2& iPos2);
    static vec2 getPhiField(const vec2& iPos1, const vec2& iPos2);
};
#endif // EFIELDLAGRANGIAN_H
