#include "lagrangiansolver.h"


void lagrangianSolver::getFieldFast(const vec2 iPos, vec2* iBodyPos, vec2(*getField)(const vec2&, const vec2&), vec2& oField)
{
    oField.x = 0.0;
    oField.y = 0.0;
    int ixIdx = (iPos.x - m_gridProp.startx)/m_gridProp.dx;
    int iyIdx = (iPos.y - m_gridProp.starty)/m_gridProp.dy;

    vec2 field;
    //int num = 0;
    for (int i = 0; i < m_gridProp.NX; i++)
        for (int j = 0; j < m_gridProp.NY; j++) {
            if(i == ixIdx && j == iyIdx)
            {
                for (int k = 0; k < m_gridProp.gridNeighbors[i][j].size(); k++)
                {
                    vec2 pos(iBodyPos[m_gridProp.gridNeighbors[i][j].at(k)].x,iBodyPos[m_gridProp.gridNeighbors[i][j].at(k)].y, iBodyPos[m_gridProp.gridNeighbors[i][j].at(k)].charge);
                    field = getField(iPos, pos);
                    oField.x += field.x;
                    oField.y += field.y;
                }
                continue;
            }
            field = getField(iPos, m_gridProp.gridCenters[i][j]);
            oField.x += field.x;
            oField.y += field.y;
        }
}
