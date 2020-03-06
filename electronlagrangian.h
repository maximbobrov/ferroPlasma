#ifndef ELECTRONLAGRANGIAN_H
#define ELECTRONLAGRANGIAN_H
#include "globals.h"

class electronLagrangian
{
public:

    int m_threadIdx;
    int m_maxParticles;
    int m_numParticles;
    vec3<float>* m_bodyAccel;
    vec4<float>* m_bodyPos;
    vec3<float>* m_bodyVel;

    electronLagrangian();
    void delete_particle(int particlesIdx);
    void wall_collision(int particlesIdx);
    void step(double dt);
    void create_random_particles();
    void getEFromElectrons(vec3<float> &bodyAccel_, double x, double y, double z,  int n);

};

#endif // ELECTRONLAGRANGIAN_H
