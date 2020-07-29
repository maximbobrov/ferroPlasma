
#include <stdio.h>
#include <stdlib.h>

#include  <GL/gl.h>
#include  <GL/glu.h>
#include  <GL/glut.h>/* glut.h includes gl.h and glu.h*/

//#include <my_include/gl.h>
//#include <my_include/glu.h>
//#include <my_include/glut.h>
#include  <math.h>
#include <time.h>
#include "globals.h"
#include <iostream>
#include <vector>
#include "efieldlagrangian.h"
#include "pzsolver.h"
#include "multisolver.h"
#include "electronlagrangian.h"




//#include <sys/time.h>



void save_fields();
void load_fields();




void display(void);
void sweep_init();
void init();

//void fmm_step(double dt);
//void sweep();
//void sweep_old();
int view=VIEW_PHI;

int redr=0;

double ck=1.0;

double cv=0.001;

double conv[100];
double conv1[100];

bool clearc=true;



std::vector <double>  avPy_;
std::vector <double>  avEy_;
std::vector <double>  area_;
std::vector <double>  T_;


bool view_px=true;

int i_tick=0;

eFieldLagrangian* lagr_solver;
pzSolver* pz_solver;
electronLagrangian* elec_solver;
multiSolver* multi_solver;


double wall_pos=0.0;
void updateEulFields()
{
    lagr_solver->updatePhi();
    lagr_solver->solvePhi(100);

     pz_solver->get_q();
    double phi_depol0=pz_solver->getPhidepol(w_x0,w_y0);


    multi_solver->updateEforPz();
    for (int i=0;i<N_X;i++)
    {
        for (int j=0;j<N_Y;j++)
        {
            double x,y;
            x=1.3*(w_x0+dx*(i));
            y=1.3*(w_y0+dy*j);

            Py_[i][j]=pz_solver->getPhidepol(x,y);//-phi_depol0;

            // printf("P=%e \n",Py_[i][j]);
            phi_[i][j]=lagr_solver->getPhi(x,y);//+Py_[i][j];



               Ey[i][j]=lagr_solver->getE(x,y).y;


            Ex[i][j]=pz_solver->getEdepol(x,y).y;//+Ey[i][j];



        }
    }
}


void display(void)
{
    double phi_max=1e-20;
    double e_max=1e-20;
    double p_max=1e-20;
    double px_max=1e-20;
    double div_max=1e-20;

    double Ec,Epc,phic,phipc,dphic,dphicp;
    for (int i=1;i<N_X-1;i++)
    {
        for (int j=1;j<N_Y-1;j++)
        {
            phi_max=fmax(phi_max,fabs(phi_[i][j]));
            p_max=fmax(p_max,fabs(Py_[i][j]));
            px_max=fmax(px_max,fabs(Px_[i][j]));
            e_max=fmax(e_max,fabs(Ey[i][j]));
            div_max=fmax(div_max,fabs(Ex[i][j]));


        }
    }

    double xc,yc,yp,ym;
    xc=(w_x0+w_x1)*0.5;
    yc=w_y0+(w_y1-w_y0)*0.25;
    yp=yc+0.5e-9;
    ym=yc-0.5e-9;

    Ec=lagr_solver->getE(xc,yc).y;
    phic=lagr_solver->getPhi(xc,yc);
    dphic=-(lagr_solver->getPhi(xc,yp)-lagr_solver->getPhi(xc,ym))/(yp-ym);


    Epc=pz_solver->getEdepol(xc,yc).y;
    phipc=pz_solver->getPhidepol(xc,yc);
    dphicp=-(pz_solver->getPhidepol(xc,yp)-pz_solver->getPhidepol(xc,ym))/(yp-ym);

    //printf("phi_max=%e Phi_p_max=%e Ey_max=%e Ep_max=%e\n",phi_max,p_max,e_max,div_max);
  //  printf("phi_c=%e Phi_p_c=%e Ey_c=%e fphic=%e Ep_c=%e ffp=%e \n",phic,phipc,Ec,dphic,Epc,dphicp);


    if (redr==1)
    {
        for (int i=0;i<10;i++)
            multi_solver->solve(1);
        double pzmax=0.0;
        double E_in=0.0;
        for (int i=0;i<pz_solver->m_p_num;i++)
        {
            if (fabs(pz_solver->m_p[i].p)>fabs(pzmax)) pzmax=(pz_solver->m_p[i].p);
            E_in+=pz_solver->m_p[i].E;
        }
        E_in/=pz_solver->m_p_num;

       // printf("pz_max=%e  E=%e E_in=%e \n",pzmax,E_global,E_in);
        updateEulFields();
        //  sweep();
    }



    int i,j;//,k,l;


    double l_2;//,tx,ty,tx0,ty0,vx,vy,v0x,v0y;
    /* clear window */


    if (clearc)
        glClear(GL_COLOR_BUFFER_BIT);


    glLoadIdentity();

    glRotatef(ry,1.0,0,0);
    glRotatef(rx,0.0,1.0,0);

    glColor3f(1,1,1);

    for (i=0;i<N_X-1;i++)
    {
        glBegin(GL_TRIANGLE_STRIP);
        for (j=0;j<N_Y;j++)
        {
            if (view==VIEW_PHI)
                l_2= ck*(phi_[i][j]);
            if (view==VIEW_P)
                l_2=ck*(Py_[i][j])/p_max;
            if (view==VIEW_PX)
                l_2=ck*(Px_[i][j])/px_max;
            if (view==VIEW_E)
                l_2=ck*(Ey[i][j])/e_max;
            if (view==VIEW_DIV)
                l_2=ck*(Ex[i][j])/div_max;

            glColor3f(l_2,l_2,-l_2);
            glVertex2f(1.3*(w_x0+dx*(i)),1.3*(w_y0+dy*j));

            if (view==VIEW_PHI)
                l_2=ck*(phi_[i+1][j]);
            if (view==VIEW_P)
                l_2=ck*(Py_[i+1][j])/p_max;
            if (view==VIEW_PX)
                l_2=ck*(Px_[i+1][j])/px_max;
            if (view==VIEW_E)
                l_2=ck*(Ey[i+1][j])/e_max;
            if (view==VIEW_DIV)
                l_2=ck*(Ex[i+1][j])/div_max;

            glColor3f(l_2,l_2,-l_2);
            glVertex2f(1.3*(w_x0+dx*(i+1)),1.3*(w_y0+dy*j));
        }
        glEnd();
    }


    glPointSize(3);
    glBegin(GL_POINTS);

    double rhomax=0.0;
    for (int i=0;i<lagr_solver->m_elec_num;i++)
    {
        if (fabs(lagr_solver->m_electrodes[i].phi_ext)>rhomax)
            rhomax=lagr_solver->m_electrodes[i].phi_ext;
    }


    //printf("phimax=%e \n",rhomax);
    for (int i=0;i<lagr_solver->m_elec_num;i++)
    {
        double c=ck*lagr_solver->m_electrodes[i].phi_ext/rhomax;
        //printf("i=%d phi=%e  phi_fix=%e\n",i,lagr_solver->m_electrodes[i].phi_ext,lagr_solver->m_electrodes[i].phi_fix);
        //glColor3f(c,c,-c);
        if (lagr_solver->m_electrodes[i].canEmit)
            glColor3f(0.25,0.25,1);
        else
            glColor3f(1,1,1);

        double x=(lagr_solver->m_electrodes[i].r.x);
        double y=(lagr_solver->m_electrodes[i].r.y);
        glVertex3f(x,y,0.0);
    }

    for (int i=0;i<lagr_solver->m_chargeNum;i++)
    {
        glColor3f(0,0.5,0);
        double x=(lagr_solver->m_charges[i].x);
        double y=(lagr_solver->m_charges[i].y);
        glVertex3f(x,y,0.0);
    }
    glEnd();

    glPointSize(5.5);
    glBegin(GL_POINTS);

    for( i=0; i<elec_solver->m_numParticles; i++ )
    {
        glColor3f(1.0,0.0,1.0);
        glVertex2f(elec_solver->m_bodyPos[i].x,elec_solver->m_bodyPos[i].y);
    }
    glEnd();

    double vel_scale=1900000.0;//sqrt(vel_scale/numParticles+0.0001);
    double leng_sacle=0.01*(w_x1-w_x0);
    glBegin(GL_LINES);

    for( i=0; i<elec_solver->m_numParticles; i++ )
    {
        // float pot=get_nearwall_potential(bodyPos[i].x,bodyPos[i].y);
        glColor3f(1,1,1);
        // glColor3f(ck*pot,ck*pot,-ck*pot);
        glVertex2f(elec_solver->m_bodyPos[i].x,elec_solver->m_bodyPos[i].y);
        glColor3f(0.0,0.0,0.0);
        glVertex2f(elec_solver->m_bodyPos[i].x+leng_sacle*elec_solver->m_bodyVel[i].x/vel_scale,
                   elec_solver->m_bodyPos[i].y+leng_sacle*elec_solver->m_bodyVel[i].y/vel_scale);
    }
    glEnd();


    /*   glColor3f(0.5,0.5,0.5);

    glBegin(GL_LINE_LOOP);

    glVertex3f(w_x0,w_y0,w_z0);
    glVertex3f(w_x1,w_y0,w_z0);
    glVertex3f(w_x1,w_y1,w_z0);
    glVertex3f(w_x0,w_y1,w_z0);
    glEnd();

    glColor3f(1,1,1);

    glBegin(GL_LINE_LOOP);

    glVertex3f(w_x0,w_y0,w_z1);
    glVertex3f(w_x1,w_y0,w_z1);
    glVertex3f(w_x1,w_y1,w_z1);
    glVertex3f(w_x0,w_y1,w_z1);
    glEnd();
*/

    if (view_px)
    {
        glBegin(GL_TRIANGLE_STRIP);
        for (int i=0;i<pz_solver->m_p_num;i++)
        {
            glColor3f(pz_solver->m_p[i].p/0.26,-pz_solver->m_p[i].p/0.26,0);
            glVertex2f(pz_solver->m_p[i].r.x/*-pz_solver->m_dx*0.5*/,pz_solver->m_p[i].r.y -pz_solver->m_p[i].dl*0.5);
            glVertex2f(pz_solver->m_p[i].r.x/*-pz_solver->m_dx*0.5*/,pz_solver->m_p[i].r.y +pz_solver->m_p[i].dl*0.5);
        }
        glEnd();
    }


    glBegin(GL_POINTS);
    glColor3f(1,0.5,0);
    glVertex2f(xc, yc);


    glEnd();

    glPointSize(3.0);

    glLineWidth(1.0);


    glDisable(GL_DEPTH_TEST);
    glLineWidth(1.5);
    glBegin(GL_LINE_STRIP);

    for( i=0; i < pz_solver->m_p_num; i++ )
    {
        glColor3f(0.5,0.5,1.0);
        glVertex2f(pz_solver->m_p[i].r.x, 1e-3 * scale * (pz_solver->m_p[i].q_ext+pz_solver->m_p[i].q) * (w_y1 - w_y0));
    }
    glEnd();

    glBegin(GL_LINE_STRIP);

    for( i=0; i < pz_solver->m_p_num; i++ )
    {
        glColor3f(0.1,0.5,0);
        glVertex2f(pz_solver->m_p[i].r.x, 1e-3 * scale * (pz_solver->m_p[i].q_ext) * (w_y1 - w_y0));
    }
    glEnd();

    /*glLineWidth(3.5);
    glBegin(GL_LINE_STRIP);

    for( i=0; i < pz_solver->m_p_num; i++ )
    {
        glColor3f(1.0,0.0,1.0);
        glVertex2f(pz_solver->m_p[i].r.x, 1e-10 * scale * pz_solver->m_p[i].E * (w_y1 - w_y0));
    }
    glEnd();

    glLineWidth(3.5);
    glBegin(GL_LINE_STRIP);
    glColor3f(1.0,1.0,1.0);
    glVertex2f(pz_solver->m_p[0].r.x, 0);
    glVertex2f(pz_solver->m_p[pz_solver->m_p_num - 1].r.x, 0);
    glEnd();*/

    glutSwapBuffers();
    if (redr==1) glutPostRedisplay();

}

void m_m(int x,int y) //mouse move
{

    if (rotate==1)
    {
        rx=rx0+0.5*(x-mx0);
        ry=ry0+0.5*(y-my0);
    }
    glutPostRedisplay();
}



void m_d(int button, int state,int x, int y)  //mouse down
{

    if (state==GLUT_UP)
    {
        rotate=0;
        rx0=rx;
        ry0=ry;
    }
    if (state==GLUT_DOWN)
    {
        rotate=1;
        mx0=x;
        my0=y;

    }

    mouse_x=(1.0*x)/W_WIDTH;

    mouse_y=(W_HEIGHT-(1.0*y))/W_HEIGHT;

    glutPostRedisplay();
}


double angle=0.0;
void kb(unsigned char key, int x, int y)
{
    int i,j,k,nn,n;
    double m,sum;
    double max_err=0.0;
    if (key=='.')
    {
        scale*=1.1;

        printf("scale=%e \n", scale);
    }

    if (key==',')
    {
        scale/=1.1;

        printf("scale=%e \n", scale);
    }
    if (key==']')
    {
        ck*=1.1;
        printf("ck=%e \n", ck);
    }

    if (key=='[')
    {
        ck/=1.1;
        printf("ck=%e \n", ck);
    }

    if (key=='1')
    {
        view=VIEW_PHI;
        printf("viewing PHI \n");
    }


    if (key=='2')
    {
        view=VIEW_P;
        printf("viewing phi_depol \n");
    }

    if (key=='3')
    {
        view=VIEW_E;
        printf("viewing E\n");
    }

    if (key=='4')
    {
        view=VIEW_DIV;

        printf("viewing E_depol\n");
    }


    if (key=='m')
    {
        move_particles=!move_particles;
    }

    if (key=='5')
    {
        // clearc=!clearc;
        view_px=!view_px;
    }


    if (key=='8')
    {
    //    wall_pos+=0.01;
    //    pz_solver->setWallPos(wall_pos);
     //  E_global=-E_global;
      //  printf("E_global=%e \n",E_global);

         g_emitElectrons=!g_emitElectrons;
    }

    if (key=='0')
    {
        //    wall_pos+=0.01;
        //    pz_solver->setWallPos(wall_pos);
        dtKoef*=1.1;
        printf("dtKoef=%e \n",dtKoef);
    }

    if (key=='9')
    {
        //  wall_pos-=0.01;
        // pz_solver->setWallPos(wall_pos);
        dtKoef/=1.1;
        printf("dtKoef=%e \n",dtKoef);

    }

    if (key=='=')
    {

        double phi_depol0=pz_solver->getPhidepol(w_x0,w_y0);
        printf("here i am ppp=%e \n",phi_depol0);
        for (int i=0;i<lagr_solver->m_elec_num;i++)
        {
            double x,y;
            x=lagr_solver->m_electrodes[i].r.x;
            y=lagr_solver->m_electrodes[i].r.y;

            lagr_solver->m_electrodes[i].phi_fix_charges=(pz_solver->getPhidepol(x,y)-phi_depol0);
        }
        printf("start ls \n");
        lagr_solver->solve_ls_fast();//solvePhi(20);
        printf("end ls \n");
        updateEulFields();
    }



    if (key=='d')
    {
        //  pz_solver->solvePz(100);
        multi_solver->solve(2);
        double pzmax=0.0;
        for (int i=0;i<pz_solver->m_p_num;i++)
        {
            if (fabs(pz_solver->m_p[i].p)>fabs(pzmax)) pzmax=(pz_solver->m_p[i].p);
        }
        printf(" pzmax =%e \n",pzmax);
    }

    if (key==' ')
    {
        redr=!redr;
    }

    glutPostRedisplay();
}




void init()
{
    glClearColor (0.0, 0.0, 0.0, 0.0);
    glColor3f(1.0, 1.0, 1.0);
    glMatrixMode (GL_PROJECTION);
    glLoadIdentity ();
    glOrtho(w_x0-(w_x1-w_x0)*0.3, w_x1*1.3, w_y0*1.3,w_y1*1.3, -10.0, 10.0);
    glMatrixMode (GL_MODELVIEW);

    lagr_solver= new eFieldLagrangian();
    lagr_solver->updatePhi();
    pz_solver= new pzSolver();

    elec_solver = new electronLagrangian();

    multi_solver = new multiSolver();
    multi_solver->m_Esolver = lagr_solver;
    multi_solver->m_pzSolver = pz_solver;
    multi_solver->m_elecSolver = elec_solver;

    //lagr_solver->solvePhi(10);
}



int main(int argc, char** argv)
{
    //srand(time(NULL));
    glutInit(&argc,argv);
    glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB);
    glutInitWindowSize(W_HEIGHT*(w_x1-w_x0)/(w_y1-w_y0),W_HEIGHT);
    glutInitWindowPosition(0,0);
    glutCreateWindow("simple");
    glutDisplayFunc(display);
    glutMotionFunc(m_m);
    glutMouseFunc(m_d);
    glutKeyboardFunc(kb);
    init();
    glutMainLoop();
}
