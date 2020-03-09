#ifndef ELECTRONLAGRANGIAN_H
#define ELECTRONLAGRANGIAN_H
#include "globals.h"

class electronLagrangian
{
public:

    int m_threadIdx;
    int m_maxParticles;
    int m_numParticles;
    vec3<double>* m_bodyAccel;
    vec4<double>* m_bodyPos;
    vec3<double>* m_bodyVel;
    vec3<double>* m_bodyE;

    electronLagrangian();
    int create_electron(vec3<double> &pos, double Emag, double Dt, double ds); //creates electrons near the electrode form field emission
    void delete_particle(int particlesIdx);
    void wall_collision(int particlesIdx);
    double calcJ(double Ein);
    void step(double dt);
    void create_random_particles();
    void getEFromElectrons(vec3<double> &bodyAccel_, double x, double y, double z,  int n);

};

#endif // ELECTRONLAGRANGIAN_H
