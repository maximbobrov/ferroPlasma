#ifndef EFIELDLAGRANGIAN_H
#define EFIELDLAGRANGIAN_H
#include "globals.h"
#include "lagrangiansolver.h"

class eFieldLagrangian : lagrangianSolver
{

public:
    struct eElem   //linear electrode elem
    {
        double phi_ext,phi_fix,phi_fix_charges;//C/m2, m, V
        vec2 r;
        bool canEmit;
        double eToEmit;
        double dl;
        double nx, ny;
    };
    vec2* m_rCentre;
    eElem * m_electrodes;
    vec2 m_charges[1000];
    int m_chargeNum;
    int m_elec_num;

public:

    eFieldLagrangian();
    void init();
    void updatePhi(); //get external potential on electrodes
    void solvePhi(int itn); //update potential
    vec2 getE(double x,double y);
    double getPhi(double x,double y);
    void setElectrodeAngle(double deg);
    void updateGridProp();
    void addQuad(vec2 p[4], double dl, double phi, int emit);
    static vec2 getEField(const vec2& iPos1, const vec2& iPos2);
    static vec2 getPhiField(const vec2& iPos1, const vec2& iPos2);

    void initW();
    double getW(double s_x, double s_y, double t_x, double t_y);
    void solve_ls();
    void getInv();
    void solve_ls_fast();
};
#endif // EFIELDLAGRANGIAN_H
