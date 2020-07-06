#ifndef LAGRANGIANSOLVER_H
#define LAGRANGIANSOLVER_H
#include "globals.h"

class lagrangianSolver
{
public:
    void getFieldFast(const vec2 iPos, vec2* iBodyPos, vec2(*getField)(const vec2&, const vec2&), vec2& oField);
    gridProp m_gridProp;
};

#endif // LAGRANGIANSOLVER_H
