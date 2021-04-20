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

    int points_distributed; //points are distributed over the grid

    multiSolver();

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
    void pzEmissionHoriz(double dt);
    void pzEmission_monte_carlo(double dt,int itn,double rat);
    double getPhi_at_electrode(int n);
    double getPhi_at_pz_down(int n);
    double getPhi_at_pz_up(int n);
    void prepare_caches();

    void fast_Fields_prepare(); //get fields on eul grid
    void fast_Fields_recalculate();
    vec2 get_fast_E(double x,double y);

    void near_Fields_recalculate_cell(int i_, int j_);
    void slower_Fields_recalculate();
    vec2 get_slower_E(double x, double y);
      vec2 get_slow_E(double x, double y);

};

#endif // MULTISOLVER_H
