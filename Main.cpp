
#include <stdio.h>
#include <stdlib.h>

//#include  <GL/gl.h>
//#include  <GL/glu.h>
//#include  <GL/glut.h>/* glut.h includes gl.h and glu.h*/

#include <my_include/gl.h>
#include <my_include/glu.h>
#include <my_include/glut.h>
#include  <math.h>
#include <time.h>
#include "globals.h"


#include "phi_mult.h"


//#include <sys/time.h>



void save_fields();
void load_fields();




void display(void);
void sweep_init();
void init();

void fmm_step(double dt);
void sweep();
//void sweep_old();
int view=VIEW_PHI;

int redr=0;

double ck=2.0;

double cv=0.001;

double conv[100];
double conv1[100];

bool clearc=true;
/*#define MAIN
#include "fmm.h"
#undef MAIN
*/
#include "sse_sum.h"
const int maxParticles=8192;
int numParticles=0;//4096;

void filter_conv(int n,int nf, FILE* file,bool write)
{

    for (int i=0;i<100;i++)
    {
        conv[i]=0.0;
        conv1[i]=0.0;
    }



    for (int i=-nf;i<=nf;i++)
    {
        conv[50+i]=1;
        conv1[i]=0.0;
    }


    for (int nn=0;nn<n;nn++)
    {
        for (int i=1;i<99;i++)
        {
            conv1[i]=0.0;
            for (int j=-nf;j<=nf;j++)
            {
                //conv[50+i]=1;
                conv1[i]+=conv[j+i];0.0;
            }

            //  conv1[i]=(conv[i]+conv[i+1]+conv[i-1])/3.0;
        }

        for (int i=1;i<99;i++)
        {
            conv[i]=conv1[i];
        }
    }

    double nrm=conv[50];
    for (int i=0;i<100;i++)
    {
        conv[i]/=nrm;
    }

    {
        int ii=50;

        double di;
        while (conv[ii]>0.5)
            ii++;
        di=(ii-1)*(fabs(0.5-conv[ii]))/fabs(conv[ii]-conv[ii-1]) + ii*(fabs(0.5-conv[ii-1]))/fabs(conv[ii]-conv[ii-1]) ;

        double cmm=0.0;
        double mm=0.0;
        for (int i=0;i<100;i++)
        {
            cmm+=conv[i]*(i-50)*(i-50);
            mm+=conv[i];
        }
        //  fprintf(file,"n=%d disp=%f disp2=%f \n",n,di-50,sqrt(cmm/mm));
        if (write)
            fprintf(file,"%d %f %f \n",n,di-50,sqrt(cmm/mm));
        else
            printf("n=%d disp=%f disp2=%f \n",n,di-50,sqrt(cmm/mm));
    }



    //   glOrtho(-0.6, 0.6, -0.6,0.6, -10.0, 10.0);
    /*
    glColor3f(1,0,0);
    glBegin(GL_LINE_STRIP);

        for (int i=0;i<100;i++)
        {
    glVertex2f((i-50)/100.0*1.15,0.2+conv[i]*0.35);

        }
        glEnd();
        glPointSize(2);
        glColor3f(0,1,0);
        glBegin(GL_POINTS);

            for (int i=0;i<100;i++)
            {
        glVertex2f((i-50)/100.0*1.15,0.2+conv[i]*0.35);

            }
            glEnd();
            */



}



int i_tick=0;

void display(void)
{
    double phi_max=1e-20;
    double e_max=1e-20;
    double p_max=1e-20;
    double px_max=1e-20;
    double div_max=1e-20;

    for (int i=1;i<N_X-1;i++)
    {
        for (int j=1;j<N_Y-1;j++)
        {
            phi_max=fmax(phi_max,fabs(phi_[i][j]));
            p_max=fmax(p_max,fabs(Py_[i][j]));
            px_max=fmax(px_max,fabs(Px_[i][j]));
            e_max=fmax(e_max,fabs(Ey[i][j]));
            div_max=fmax(div_max,fabs(div_[i][j]));
        }
    }


    double bmax=0.0;
    double Emax=0.0001;


    for(int i=0; i<N_X; i++ )
    {
        if (fabs(q[i])>bmax)
        {
            bmax=fabs(q[i]);//BoundaryLayerGauss[i]);
        }

    }

    if (redr==1)
    {
        sweep();
        sweep();
        //fmm_step(dt*4000.0);
        //fmm_step(dt*4000.0);
        sweep();
        sweep();
        sweep();
        sweep();
        sweep();


    }



    i_tick++;
    if (i_tick>10){
        double emin=1e23;
        double emax=-1e23;
        for (int i=1; i<N_X-1; i++)
        {
            for (int j=1; j<N_Y-1; j++)
            {
                if (emin>Ey[i][j]) emin=Ey[i][j];
                if (emax<Ey[i][j]) emax=Ey[i][j];
            }
        }

        printf("emin=%e emax=%e \n",emin,emax);



        emin=1e23;
        emax=-1e23;
        for (int i=1; i<N_X-1; i++)
        {
            for (int j=1; j<N_Y-1; j++)
            {

                if (emin>Py_[i][j]) emin=Py_[i][j];
                if (emax<Py_[i][j]) emax=Py_[i][j];
            }
        }




        printf("pmin=%e pmax=%e \n",emin,emax);



        printf("bmax=%e divmax=%e \n",bmax,div_max);





        i_tick=0;
    }

    int i,j;//,k,l;


    double l_2;//,tx,ty,tx0,ty0,vx,vy,v0x,v0y;
    /* clear window */


    if (clearc)
        glClear(GL_COLOR_BUFFER_BIT);


    glLoadIdentity();

    // glTranslatef(0,0.4,0);
    // glRotatef(-90,0,0,1);
    glRotatef(ry,1.0,0,0);
    glRotatef(rx,0.0,1.0,0);

    glColor3f(1,1,1);

    //filter_conv(0);

    //filter_conv(2);
    //////////
    /*glColor3f(1,0,0);
glBegin(GL_LINE_STRIP);

    for (int i=0;i<100;i++)
    {
glVertex2f((i-50)/100.0*1.15,0.2+conv[i]*0.35);

    }
    glEnd();
    glPointSize(2);
    glColor3f(0,1,0);
    glBegin(GL_POINTS);

        for (int i=0;i<100;i++)
        {
    glVertex2f((i-50)/100.0*1.15,0.2+conv[i]*0.35);

        }
        glEnd();
*/
    /////////////



    //printf("Pxmax=%e \n",px_max);
    // glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);


    for (i=0;i<N_X-1;i++)
    {
        glBegin(GL_TRIANGLE_STRIP);


        for (j=0;j<N_Y;j++)
        {
            if (view==VIEW_PHI)
                l_2=ck*(phi_[i][j])/phi_max;
            if (view==VIEW_P)
                l_2=ck*(Py_[i][j])/p_max;
            if (view==VIEW_PX)
                l_2=ck*(Px_[i][j])/px_max;
            if (view==VIEW_E)
                l_2=ck*(Ey[i][j])/e_max;
            if (view==VIEW_DIV)
                l_2=ck*(div_[i][j])/div_max;

            glColor3f(l_2,l_2,-l_2);
            glVertex2f(dx*(i),dy*(j-N_Y/2));

            if (view==VIEW_PHI)
                l_2=ck*(phi_[i+1][j])/phi_max;
            if (view==VIEW_P)
                l_2=ck*(Py_[i+1][j])/p_max;
            if (view==VIEW_PX)
                l_2=ck*(Px_[i+1][j])/px_max;
            if (view==VIEW_E)
                l_2=ck*(Ey[i+1][j])/e_max;
            if (view==VIEW_DIV)
                l_2=ck*(div_[i+1][j])/div_max;

            glColor3f(l_2,l_2,-l_2);
            glVertex2f(dx*(i+1),dy*(j-N_Y/2));
        }


        glEnd();

    }


    //glEnable(GL_LINE_SMOOTH);

    /*    for (int ii=0;ii<N_X;ii++)
    {
        i=ii;
        glBegin(GL_LINES);
        for (int jj=0;jj<8;jj++)
        {

            j=jj*8+1;

            glColor3f(1,1,1);
            glVertex2f(dx*(i-N_X/2),dy*(j-N_Y/2));



            glColor3f(0.5,0.5,0.5);
            if (view_v==VIEW_JX)  glVertex2f(dx*(i-N_X/2)+cv*Jx_[i][j], dy*(j-N_Y/2)+cv*Jy_[i][j]);
            if (view_v==VIEW_UX)  glVertex2f(dx*(i-N_X/2)+cv*Ux_[i][j], dy*(j-N_Y/2)+cv*Uy_[i][j]);
        }
        glEnd();
    }

  */  glEnable(GL_BLEND);

    glEnable(GL_POINT_SMOOTH);
    glPointSize(1.5);



    double vel_scale=0.0;
    for( i=0; i<numParticles; i++ ) {
        vel_scale+=bodyVel[i].x*bodyVel[i].x+bodyVel[i].y*bodyVel[i].y+bodyVel[i].z*bodyVel[i].z;

    }

    vel_scale=sqrt(vel_scale/numParticles+0.0001);
    double leng_sacle=0.01*(w_x1-w_x0);
    glBegin(GL_LINES);

    for( i=0; i<numParticles; i++ )
    {
        // float pot=get_nearwall_potential(bodyPos[i].x,bodyPos[i].y);
        glColor3f(1.0*fabs(bodyVel[i].x)/vel_scale,4.0*fabs(bodyVel[i].y)/vel_scale,3.0*fabs(bodyVel[i].z)/vel_scale);
        // glColor3f(ck*pot,ck*pot,-ck*pot);
        glVertex3f(bodyPos[i].x,bodyPos[i].y,bodyPos[i].z);
        glColor3f(0.0,0.0,0.0);
        glVertex3f(bodyPos[i].x+leng_sacle*bodyVel[i].x/vel_scale,bodyPos[i].y+leng_sacle*bodyVel[i].y/vel_scale,bodyPos[i].z+leng_sacle*bodyVel[i].z/vel_scale);
    }
    glEnd();




    for( i=0; i<N_X; i++ )
    {

        /*        glColor3f((BoundaryLayerGauss[i])/bmax,-(BoundaryLayerGauss[i])/bmax,0);
        j = -2;
        glBegin(GL_QUADS);
        glVertex3f(w_x0+i*dx ,w_y0+j*dy ,0);
        glColor3f((BoundaryLayerGauss[i+1])/bmax,-(BoundaryLayerGauss[i+1])/bmax,0);
        glVertex3f(w_x0+(i+1)*dx ,w_y0+j*dy ,0);
        glVertex3f(w_x0+(i+1)*dx ,w_y0+(j-3)*dy ,0);
        glColor3f((BoundaryLayerGauss[i])/bmax,-(BoundaryLayerGauss[i])/bmax,0);
        glVertex3f(w_x0+i*dx ,w_y0+(j-3)*dy ,0);
        glEnd();
  */

        glColor3f((q[i])/bmax,-(q[i])/bmax,0);
        j = 0;
        glBegin(GL_QUADS);
        glVertex3f(w_x0+i*dx ,w_y0+j*dy ,0);
        glColor3f((q[i+1])/bmax,-(q[i+1])/bmax,0);
        glVertex3f(w_x0+(i+1)*dx ,w_y0+j*dy ,0);
        glVertex3f(w_x0+(i+1)*dx ,w_y0+(j-6)*dy ,0);
        glColor3f((q[i])/bmax,-(q[i])/bmax,0);
        glVertex3f(w_x0+i*dx ,w_y0+(j-6)*dy ,0);
        glEnd();

    }


    /*  for( i=0; i<N_X; i++ )
    {
        //glColor3f(ck*(BoundaryLayer[i]) * 1.0 / (maxB-minB),-ck*(BoundaryLayer[i]) * 1.0 / (maxB-minB),0);
        glColor3f((WallEnergy[i])/Emax,-(WallEnergy[i])/Emax,0);
        j = 0;
        glBegin(GL_QUADS);
        glVertex3f(w_x0+i*dx ,w_y0+j*dy ,0);
        glColor3f((WallEnergy[i+1])/Emax,-(WallEnergy[i+1])/Emax,0);
        glVertex3f(w_x0+(i+1)*dx ,w_y0+j*dy ,0);
        glVertex3f(w_x0+(i+1)*dx ,w_y0+(j-2)*dy ,0);
        glColor3f((WallEnergy[i])/Emax,-(WallEnergy[i])/Emax,0);
        glVertex3f(w_x0+i*dx ,w_y0+(j-2)*dy ,0);
        glEnd();
    }*/


    glColor3f(0.5,0.5,0.5);

    glBegin(GL_LINE_LOOP);

    /*glVertex3f(dx*(-N_X/2),dy*(-N_Y/2),0);
    glVertex3f(dx*(N_X-1-N_X/2),dy*(-N_Y/2),0);
    glVertex3f(dx*(N_X-1-N_X/2),dy*(N_Y-1-N_Y/2),0);
    glVertex3f(dx*(-N_X/2),dy*(N_Y-1-N_Y/2),0);*/

    glVertex3f(w_x0,w_y0,w_z0);
    glVertex3f(w_x1,w_y0,w_z0);
    glVertex3f(w_x1,w_y1,w_z0);
    glVertex3f(w_x0,w_y1,w_z0);
    glEnd();

    glColor3f(1,1,1);

    glBegin(GL_LINE_LOOP);

    /*glVertex3f(dx*(-N_X/2),dy*(-N_Y/2),0);
    glVertex3f(dx*(N_X-1-N_X/2),dy*(-N_Y/2),0);
    glVertex3f(dx*(N_X-1-N_X/2),dy*(N_Y-1-N_Y/2),0);
    glVertex3f(dx*(-N_X/2),dy*(N_Y-1-N_Y/2),0);*/

    glVertex3f(w_x0,w_y0,w_z1);
    glVertex3f(w_x1,w_y0,w_z1);
    glVertex3f(w_x1,w_y1,w_z1);
    glVertex3f(w_x0,w_y1,w_z1);
    glEnd();



    glPointSize(3.0);

    glLineWidth(1.0);

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

int iglob=0;



//FmmSystem tre;


vec3<float> getE(float x, float y)
{

    double a0=fmax(fmin(N_X-1,(N_X-1)*(x-w_x0)/(w_x1-w_x0)),0);
    int i0=(int)a0;
    a0-=i0;


    double b0=fmax(fmin(N_Y-1,(N_Y-1)*(y-w_y0)/(w_y1-w_y0)),0);
    int j0=(int)b0;
    b0-=j0;

    vec3<float> ret;

    ret.x=(1.0-b0)*((1.0-a0)*Ex[i0][j0]+(a0)*Ex[i0+1][j0]) + (b0)*((1.0-a0)*Ex[i0][j0+1]+(a0)*Ex[i0+1][j0+1]);
    ret.y=(1.0-b0)*((1.0-a0)*Ey[i0][j0]+(a0)*Ey[i0+1][j0]) + (b0)*((1.0-a0)*Ey[i0][j0+1]+(a0)*Ey[i0+1][j0+1]);
    ret.z;
    return ret;


}



vec4<float> init_pos()
{
    vec4<float> res;

    res.x = w_x0+0.001*(w_x1-w_x0);
    res.y = 0.0+0.0251*(w_y1)+my_rand(0)*my_rand(0)*(w_y1)*0.75;
    res.z = rand()/(float) RAND_MAX*(w_z1-w_z0)+w_z0;
    res.w = qe; //single electrons

    return res;

}

vec3<float> init_vel()
{
    vec3<float> vel;

    double angle = 0.98*acos(1-2*my_rand(0)) - M_PI / 2;
    double en = getVms_from_Ev(0.2);//0.1e2;
    vel.x = en *  cos(angle);
    vel.y = en * sin(angle);


    return vel;
}

/*void create_random_particles(int threadIdx, vec4<float> *bodyPos_, vec3<float> *bodyVel_,vec3<float> *bodyAccel_)
{
    int curNum = numParticles;
    int numToAdd = std::min(int(my_rand(threadIdx) * 5.0), maxParticles - numParticles-1);
    numParticles += numToAdd;
    for (int i = curNum; i < curNum + numToAdd; ++i)
    {
        vec4<float> pos=init_pos();//threadIdx);
        vec3<float> vel=init_vel();//(threadIdx);
        bodyPos_[i].x=pos.x;
        bodyPos_[i].y=pos.y;
        bodyPos_[i].z=pos.z;
        bodyPos_[i].w=pos.w;
        bodyVel_[i]=vel;
        bodyAccel_[i].x=0.0;
        bodyAccel_[i].y=0.0;
        bodyAccel_[i].z=0.0;
    }
}*/

double calcJ(double E)
{
    double t = 1.1;
    double B = 10.0;
    double phi = 4.0;
    double y = 3.79 * 1e-4 * sqrt(fabs(B * E)) / phi;
    double tetta = 0.95 - 1.03 * y * y;
    return (1.54 * 1e-6 * B * B * E * E / (t *t  * phi)) * exp( - 6.83 * 1e7 * pow(phi, 1.5) * tetta / fabs( B * E));
}


void create_random_particles(int threadIdx, vec4<float> *bodyPos_, vec3<float> *bodyVel_,vec3<float> *bodyAccel_)
{
    /*int curNum = numParticles[threadIdx];
    int numToAdd = std::min(int(my_rand(threadIdx) * 2.0), maxParticles - numParticles[threadIdx]-1);
    numParticles[threadIdx] += numToAdd;
    for (int i = curNum; i < curNum + numToAdd; ++i)
    {
        vec4 pos=init_pos(threadIdx);
        vec3 vel=init_vel(threadIdx);
        bodyPos_[i].x=pos.x;
        bodyPos_[i].y=pos.y;
        bodyPos_[i].z=pos.z;
        bodyPos_[i].w=pos.w;
        bodyVel_[i]=vel;
        bodyAccel_[i].x=0.0;
        bodyAccel_[i].y=0.0;
        bodyAccel_[i].z=0.0;
    }*/
    for (int i = (N_Y-1)/2; i < (N_Y-1); i+=1)
    {


        double dz = w_z1 - w_z0;
        double x_ = w_x0 + 0.001 * (w_x1-w_x0);//- N_X * dx / 2;
        double y_ = /*0.0251 * (w_y1) +*/w_y0 + i * (w_y1 - w_y0) / (N_Y-1);//i * dy - N_Y * dy / 2;
        double z_ = my_rand(0) * dz - dz / 2;
        //vec3 E = getE(x_, y_);
        vec3<float> E;
        getEFromElectrons(E, x_, y_, z_, numParticles);
        //printf("Ex=%e Ey = %e\n", E.x, E.y);
        E.x = Ex[0][i];
        E.y = Ey[0][i];
        E.x = fmax(0, -E.x);
        //printf("Ex=%e Ey = %e\n", E.x, E.y);
        double J =calcJ(0.01 * E.x);
        int chargeNum = 1e0;
        printf("JJJJ = %f jjj = %e E = %e \n", (dt * J * dy * dz * 1e4  * 6.24151 * 1e18 / chargeNum), J , 0.01 *E.x);
        int curNum = numParticles;
        int numToAdd =  std::min(int(dt * J * dy * dz * 1e4 * 6.24151 * 1e18 / chargeNum), maxParticles - numParticles-2);
        numToAdd = numToAdd > 5 ? 5 : numToAdd;
        if(my_rand(threadIdx) < (dt * J * dy * dz * 1e4 * 6.24151 * 1e18 / chargeNum) - int(dt * J * dy * dz * 1e4  * 6.24151 * 1e18 / chargeNum))
            numToAdd++;

        numParticles += numToAdd;
        for (int n = curNum; n < curNum + numToAdd; ++n)
        {
            bodyPos_[n].x = x_;
            bodyPos_[n].y = y_;
            bodyPos_[n].z = z_;
            bodyPos_[n].w = 1.6e-19 * chargeNum;
            double angle =  0.98 * acos(1 - 2 * my_rand(threadIdx)) - M_PI / 2;
            //printf("Ey = %f\n", angle);
            if (y_ < w_y0 + 0.1 * (w_y1 - w_y0) && angle < 0)
                angle *= -1;
            double en = getVms_from_Ev(0.2);
            bodyVel_[n].x = en * cos(angle);
            bodyVel_[n].y = en * sin(angle);
            bodyVel_[n].z = 0;
            bodyAccel_[n].x = 0.0;
            bodyAccel_[n].y = 0.0;
            bodyAccel_[n].z = 0.0;
        }
    }
}

void delete_particle(int threadIdx, int particlesIdx, vec3<float> *bodyAccel_, vec4<float> *bodyPos_, vec3<float> *bodyVel_)
{
    bodyPos_[particlesIdx] = bodyPos_[numParticles-1];
    bodyVel_[particlesIdx] = bodyVel_[numParticles-1];
    bodyAccel_[particlesIdx] = bodyAccel_[numParticles-1];
    numParticles -= 1;
    //float r = my_rand(threadIdx);
    //randMax  = randMax > r ? randMax : r;
    // randMin  = randMin < r ? randMin : r;
}

void wall_collision(int threadIdx, int particlesIdx, vec3<float> *bodyAccel_, vec4<float> *bodyPos_, vec3<float> *bodyVel_)
{
    int xIdx = (bodyPos[particlesIdx].x - w_x0)/dx;
    if ((xIdx>=0)&&(xIdx<N_X))
    {
        double me=9.1e-31*1e8;
        double Ev_in_J=1.6e-19;
        double v2=((bodyVel[particlesIdx].x*bodyVel[particlesIdx].x)+(bodyVel[particlesIdx].y*bodyVel[particlesIdx].y)+(bodyVel[particlesIdx].z*bodyVel[particlesIdx].z));
        // double v2_y=(bodyVel[particlesIdx].y*bodyVel[particlesIdx].y); //for a potential

        double E = 0.5*(me/Ev_in_J)*v2;






        //simple deletion
        //double v2y=bodyVel[particlesIdx].y*bodyVel[particlesIdx].y;
        //double phi



        q[xIdx]+=bodyPos_[particlesIdx].w/2.75;

        double addit=0.5;
        for (int i=xIdx-1;i>=fmax(xIdx-3,0);i--)
        {
            q[i]+=addit*bodyPos_[particlesIdx].w/2.75;
            addit*=0.5;
        }

        addit=0.5;
        for (int i=xIdx+1;i<=fmin(xIdx+3,N_X-1);i++)
        {
            q[i]+=addit*bodyPos_[particlesIdx].w/2.75;
            addit*=0.5;
        }

        if(my_rand(0)>0.9)
        {
            int curNum = numParticles;
            int numToAdd = std::min(int(my_rand(threadIdx) * 3.0), maxParticles - numParticles-1);
            numParticles += numToAdd;
            for (int i = curNum; i < curNum + numToAdd; ++i)
            {
                bodyPos_[i].x= bodyPos_[particlesIdx].x+(my_rand(threadIdx)*dx*2-dx);
                bodyPos_[i].y= 0.0+0.0151*(w_y1);//bodyPos_[particlesIdx].y;
                bodyPos_[i].z= bodyPos_[particlesIdx].z+(my_rand(threadIdx)*dx*2-dx);
                bodyPos_[i].w= bodyPos_[particlesIdx].w;
                double angleNew = 0.98*acos(1-2*my_rand(threadIdx))+0.01*M_PI;
                double En = my_rand(0)*0.4e2 ;
                bodyVel_[i].x =En * cos(angleNew);//0.5*sqrt(fabs(E/(num))) * cos(angleNew);
                bodyVel_[i].y =En * sin(angleNew); //0.5*sqrt(fabs(E/(num))) * sin(angleNew);

                bodyVel_[i].z = 0.0;
            }
            delete_particle( threadIdx, particlesIdx, bodyAccel_, bodyPos_, bodyVel_);
        }
        else
            delete_particle( threadIdx, particlesIdx, bodyAccel_, bodyPos_, bodyVel_);

        /*
        //BoundaryLayer[xIdx] += 1.0;


        //double delta= delta_f(E, phi);

///till here

        double phi  = fabs(atan2(bodyVel[particlesIdx].x, bodyVel[particlesIdx].y));

        double delta= delta_f(E, phi);

        int num = (int)(fabs(delta));

        if(my_rand(threadIdx) < (delta - 1.0 * num))
            num++;


        if (E<1e-5*WallEnergy[xIdx])
        {
            delete_particle( threadIdx, particlesIdx, bodyAccel_, bodyPos_, bodyVel_);
            BoundaryLayer[xIdx] += 1.0;
        }else
        {

            //  printf("collison occured xIdx=%d E=%e n=%d \n",xIdx,E,num);

            if((maxParticles - numParticles> num)&&(delta>0))
            {
                int size = std::min(maxParticles - numParticles, num);
                //for (int i = 0; i < size; ++i) {
                int i=numParticles;
                int n=0;
                while((i<maxParticles)&&(n<num))
                {

                    bodyPos_[i].x= bodyPos_[particlesIdx].x+3.0*(my_rand(threadIdx)*dx*2-dx);
                    bodyPos_[i].y= 0+0.001*(w_y1-0);//bodyPos_[particlesIdx].y;
                    bodyPos_[i].z= bodyPos_[particlesIdx].z+3.0*(my_rand(threadIdx)*dx*2-dx);
                    bodyPos_[i].w= bodyPos_[particlesIdx].w;
                    double angleNew = 0.98*acos(1-2*my_rand(threadIdx))+0.01*M_PI;

                    double v_new=0.6*sqrt(my_rand(threadIdx)*v2/(num));
                    bodyVel_[i].x =v_new * cos(angleNew);//0.5*sqrt(fabs(E/(num))) * cos(angleNew);
                    bodyVel_[i].y =v_new * sin(angleNew); //0.5*sqrt(fabs(E/(num))) * sin(angleNew);

                    bodyVel_[i].z = 0.0;

                    //BoundaryLayer[xIdx] -= 1.0;

                    q[xIdx]-=1.0/2.75;

                    double addit=0.5;
                    for (int i=xIdx-1;i>=fmax(xIdx-3,0);i--)
                    {
                        q[i]-=addit/2.75;
                        addit*=0.5;
                    }

                     addit=0.5;
                    for (int i=xIdx+1;i<=fmin(xIdx+3,N_X-1);i++)
                    {
                        q[i]-=addit/2.75;
                        addit*=0.5;
                    }

                    i++;
                    n++;

                    numParticles=i;
                }

                delete_particle( threadIdx, particlesIdx, bodyAccel_, bodyPos_, bodyVel_);


                        q[xIdx]+=1.0/2.75;

                        double addit=0.5;
                        for (int i=xIdx-1;i>=fmax(xIdx-3,0);i--)
                        {
                            q[i]+=addit/2.75;
                            addit*=0.5;
                        }

                         addit=0.5;
                        for (int i=xIdx+1;i<=fmin(xIdx+3,N_X-1);i++)
                        {
                            q[i]+=addit/2.75;
                            addit*=0.5;
                        }

            }else
            {
                bodyPos_[particlesIdx].y -= dt*bodyPos_[particlesIdx].y;
                bodyVel_[particlesIdx].y=fabs(bodyVel_[particlesIdx].y);
            }
        }

        /* double Em = 650;// in ev

        double phi  = atan2(bodyVel[particlesIdx].y, bodyVel[particlesIdx].x);
        double Z = 1.284 * E / (Em * (1 + phi * phi/(M_PI*M_PI)));
        double sigma_m = 6.4;
        double sigma = 1.526 * (1 + phi * phi/(4*M_PI*M_PI)) * (1 - exp(-pow(Z, 1.725))) * sigma_m /(pow(Z, 1.725)+0.1);

        int num = (int)(fabs(sigma));

        if(my_rand(threadIdx) < (sigma - 1.0 * num))
            num++;

         printf("eeee=%e num=%d \n",E,num);


        if(maxParticles - numParticles> num)
        {
            int size = std::min(maxParticles - numParticles, num);
            //for (int i = 0; i < size; ++i) {
            int i=numParticles;
            int n=0;
            while((i<maxParticles)&&(n<num))
            {

                bodyPos_[i].x= bodyPos_[particlesIdx].x+3.0*(my_rand(threadIdx)*dx*2-dx);
                bodyPos_[i].y= -0.15+0.001;//bodyPos_[particlesIdx].y;
                bodyPos_[i].z= bodyPos_[particlesIdx].z+3.0*(my_rand(threadIdx)*dx*2-dx);
                bodyPos_[i].w= bodyPos_[particlesIdx].w;
                double angleNew = 0.98*acos(1-2*my_rand(threadIdx))+0.01*M_PI;

                bodyVel_[i].x =sqrt(fabs(v2/(num))) * cos(angleNew);//0.5*sqrt(fabs(E/(num))) * cos(angleNew);
                bodyVel_[i].y =sqrt(fabs(v2/(num))) * sin(angleNew); //0.5*sqrt(fabs(E/(num))) * sin(angleNew);

                bodyVel_[i].z = 0.0;

                BoundaryLayer[xIdx] -= 1.0;
                i++;
                n++;

                numParticles=i;
            }

            delete_particle( threadIdx, particlesIdx, bodyAccel_, bodyPos_, bodyVel_);
            BoundaryLayer[xIdx] += 1.0;

        }else
        {
            bodyPos_[particlesIdx].y -= dt*bodyPos_[particlesIdx].y;
            bodyVel_[particlesIdx].y=fabs(bodyVel_[particlesIdx].y);
        }*/



    }
    else
    {
        delete_particle( threadIdx, particlesIdx, bodyAccel_, bodyPos_, bodyVel_);
    }
}

void fmm_step(double dt)
{
    int i;
    double tic,toc,timeFMM;

    static float BoundaryLayerGauss2[N_X];

    /*    for (i=1;i<N_X-1;i++)
    {
        BoundaryLayerGauss[i]=0.5*(BoundaryLayer[i]+0.5*(BoundaryLayer[i+1]+BoundaryLayer[i-1]));

    }

    BoundaryLayerGauss[0]=BoundaryLayerGauss[1];
    BoundaryLayerGauss[N_X-1]=BoundaryLayerGauss[N_X-2];

    for (int j=0;j<1;j++)
    {
        for (i=1;i<N_X-1;i++)
        {
            BoundaryLayerGauss2[i]=0.5*(BoundaryLayerGauss[i]+0.5*(BoundaryLayerGauss[i+1]+BoundaryLayerGauss[i-1]));

        }

        BoundaryLayerGauss2[0]=BoundaryLayerGauss2[1];
        BoundaryLayerGauss2[N_X-1]=BoundaryLayerGauss2[N_X-2];


        for (i=1;i<N_X-1;i++)
        {
            BoundaryLayerGauss[i]=0.5*(BoundaryLayerGauss2[i]+0.5*(BoundaryLayerGauss2[i+1]+BoundaryLayerGauss2[i-1]));

        }

        BoundaryLayerGauss[0]=BoundaryLayerGauss[1];
        BoundaryLayerGauss[N_X-1]=BoundaryLayerGauss[N_X-2];
    }*/


    /*

    for (i=0;i<N_X;i++)
    {
        float x=w_x0+dx*i;
        WallEnergy[i]=get_nearwall_potential(x,Y_WALL);
    }*/
   /* for (int i=1; i<7; i++)
    {
        printf("E =%e J=%e \n",i*2e7, calcJ(i*2e7));
    }*/
    create_random_particles(0, bodyPos, bodyVel,bodyAccel);
    //create_random_particles(0, bodyPos, bodyVel,bodyAccel);

    tic = get_time();
    //tre.fmmMain(numParticles,1);//treeOrFMM);
    direct_seq(numParticles);
    toc = get_time();
    timeFMM = toc-tic;

    double lay_num=0;
    for (i=0;i<N_X;i++)
    {
        lay_num+=fabs(BoundaryLayer[i]);
    }

    //  printf("fmm    : %g  pnum=%d bnum=%f ck=%e\n",timeFMM, numParticles,lay_num,ck);



    for( i=0; i<numParticles; i++ ) {


        float magn=qe/Me;//1e-1;
        vec3<float> ev=getE(bodyPos[i].x,bodyPos[i].y);
        bodyVel[i].x -= magn*(ev.x +bodyAccel[i].x)*dt;
        bodyVel[i].y -= magn*(ev.y +bodyAccel[i].y)*dt;
        bodyVel[i].z -=magn*dt*bodyAccel[i].z;

        /*
    bodyVel[i].x *= 0.9995;
    bodyVel[i].y *= 0.9995;
    bodyVel[i].z *= 0.9995;*/


    }

    for( i=0; i<numParticles; i++ )
    {
        bodyPos[i].x += dt*bodyVel[i].x;
        bodyPos[i].y += dt*bodyVel[i].y;
        bodyPos[i].z += dt*bodyVel[i].z;

        if (bodyPos[i].y<0)
        {
            // bodyPos[i].y -= dt*bodyVel[i].y;
            // bodyVel[i].y=fabs(bodyVel[i].y);
            wall_collision( 0, i, bodyAccel, bodyPos, bodyVel);
        }

        if (bodyPos[i].y>w_y1 || bodyPos[i].x < w_x0 || bodyPos[i].x>w_x1)
        {
            delete_particle( 0, i, bodyAccel, bodyPos, bodyVel);
        }

        /////////////
        if (bodyPos[i].z<w_z0) //periodic
        {
            bodyPos[i].z +=w_z1-w_z0 ;//dt*bodyVel[i].z;
            // bodyVel[i].z=fabs(bodyVel[i].z);
        }

        if (bodyPos[i].z>w_z1)
        {
            bodyPos[i].z -=w_z1-w_z0 ;//dt*bodyVel[i].z;
            //bodyVel[i].z=-fabs(bodyVel[i].z);
        }
    }
}

void fmm_init()
{
    int i;
    rand_init();
    bodyAccel = new vec3<float>[maxParticles];
    bodyVel = new vec3<float>[maxParticles];
    bodyPos = new vec4<float>[maxParticles];

    for( i=0; i<maxParticles; i++ ) {
        bodyPos[i].x = /*0.5*(w_x1-w_x0)+*/rand()/(float) RAND_MAX*(w_x1-w_x0)*0.95+w_x0;
        bodyPos[i].y = rand()/(float) RAND_MAX*(w_y1-0)*0.95+0;
        bodyPos[i].z = rand()/(float) RAND_MAX*(w_z1-w_z0)+w_z0;
        bodyPos[i].w = 1.0;

        bodyVel[i].x = 300.0*(0.4+rand()/(float) RAND_MAX*0.1);
        bodyVel[i].y = 300.0*(-0.005+rand()/(float) RAND_MAX*0.01-0.005);
        bodyVel[i].z = 300.0*(rand()/(float) RAND_MAX*0.01-0.005);
    }

    for (i=0; i<N_X;i++)
    {
        BoundaryLayer[i]=0.000001;
    }

    /*for( i=0; i<maxParticles; i++ ) {
  bodyPos[i].x = rand()/(float) RAND_MAX*2*M_PI-M_PI;
  bodyPos[i].y = rand()/(float) RAND_MAX*2*M_PI-M_PI;
  bodyPos[i].z = rand()/(float) RAND_MAX*2*M_PI-M_PI;
  bodyPos[i].w = rand()/(float) RAND_MAX;
}*/

}



/*
void fmm_init()
{
    int i,iteration,numParticles;
    double tic,toc,timeDirect,timeFMM,L2norm,difference,normalizer;
    vec3<float> *bodyAcceld;
    FmmKernel kernel;
    FmmSystem tree;
    std::fstream fid("time2.dat",std::ios::out);

    bodyAccel = new vec3<float>[maxParticles];
    bodyAcceld = new vec3<float>[maxParticles];
    bodyPos = new vec4<float>[maxParticles];

    for( i=0; i<maxParticles; i++ ) {
      bodyPos[i].x = rand()/(float) RAND_MAX*2*M_PI-M_PI;
      bodyPos[i].y = rand()/(float) RAND_MAX*2*M_PI-M_PI;
      bodyPos[i].z = rand()/(float) RAND_MAX*2*M_PI-M_PI;
      bodyPos[i].w = rand()/(float) RAND_MAX;
    }

   // for( iteration=0; iteration<25; iteration++ ) {
      numParticles =maxParticles; // int(pow(10,(iteration+32)/8.0));
      printf("N = %d\n",numParticles);

      tic = get_time();
      tree.fmmMain(numParticles,1);
      toc = get_time();
      timeFMM = toc-tic;
      printf("fmm    : %g\n",timeFMM);
      for ( i=0; i<9; i++ ) fid << t[i] << " ";
      fid << std::endl;
      for( i=0; i<numParticles; i++ ) {
        bodyAcceld[i].x = bodyAccel[i].x;
        bodyAcceld[i].y = bodyAccel[i].y;
        bodyAcceld[i].z = bodyAccel[i].z;
      }

      tic = get_time();
      kernel.direct(numParticles);
      toc = get_time();
      timeDirect = toc-tic;
      printf("direct : %g\n",timeDirect);

      L2norm = 0;
      for( i=0; i<numParticles; i++ ) {
        difference = (bodyAccel[i].x-bodyAcceld[i].x)*
               (bodyAccel[i].x-bodyAcceld[i].x)+
               (bodyAccel[i].y-bodyAcceld[i].y)*
               (bodyAccel[i].y-bodyAcceld[i].y)+
               (bodyAccel[i].z-bodyAcceld[i].z)*
               (bodyAccel[i].z-bodyAcceld[i].z);
        normalizer = bodyAccel[i].x*bodyAccel[i].x+
               bodyAccel[i].y*bodyAccel[i].y+
               bodyAccel[i].z*bodyAccel[i].z;
        L2norm += difference/normalizer/numParticles;
      }
      L2norm = sqrt(L2norm);
      printf("error  : %g\n\n",L2norm);

    fid.close();
}*/


void kb(unsigned char key, int x, int y)
{
    int i,j,k,nn,n;
    double m,sum;
    double max_err=0.0;
    if (key=='.')
    {


        ck*=1.1;


    }


    if (key==',')
    {


        ck/=1.1;


    }
    if (key==']')
    {


        dt*=1.1;
        printf("dt=%e \n",dt);

    }


    if (key=='[')
    {


        dt/=1.1;
        printf("dt=%e \n",dt);

    }


    if (key=='1')
    {

        view=VIEW_PHI;
        printf("viewing PHI \n");

    }


    if (key=='2')
    {

        view=VIEW_E;
        printf("viewing E \n");

    }

    if (key=='f')
    {
        move_particles=!move_particles;
        sasign*=-1.0;
        for(int i=0;i<20;i++)
          sweep();
        move_particles=!move_particles;
    }

    if (key=='m')
    {

        move_particles=!move_particles;

    }

    if (key=='3')
    {

        view=VIEW_P;
        printf("viewing P \n");

    }




    if (key=='4')
    {

        view=VIEW_PX;
        printf("viewing PX \n");

    }


    if (key=='5')
    {

        clearc=!clearc;
        // printf("viewing DIV \n");

    }





    if (key=='s')
    {

        //   for(int i=0;i<10;i++)
        sweep();


        //   fmm_step(0.0000001);

    }
    if (key=='d')
    {

        solve_current_1(rho_,Jx_,Jy_,Ux_,Uy_,80);
    }

    if (key=='u')
    {

        //  for(int i=0;i<100;i++)

        field_U(1.0, Ux_,1000);

    }

    if (key=='n')
    {

        //  for(int i=0;i<100;i++)

        NavierStoks_solver(p_in, Ux_, Uy_,1000);

    }


    if (key==' ')
    {
        redr=!redr;
        // filter_conv(iglob,1,NULL,false);
        //iglob++;
    }


    glutPostRedisplay();
}




void sweep_init()
{
    int i,j,k,n;
    double App,R_2,b,m;


    for (int i=1; i<N_X-1; i++)
    {
        /*        if (i<50)
            q[i]=-0.013/dy;
        else q[i]=0.0;*/
        q[i]=0.0;

    }

    float gsum=0.0;
    for (int i=0; i<N_Y; i++)
    {
        gau[i]=0.0;
        if (fabs(i-(N_Y_DIEL-1))<10)
        {
            int pow2n=1<<abs(i-(N_Y_DIEL-1));

            gau[i]=1.0/pow2n;
            gsum+=gau[i];
        }
    }


    for (int i=0; i<N_Y; i++)
    {
        gau[i]/=gsum;
    }



    for (i=0;i<N_X;i++)
    {

        for (j=0;j<N_Y;j++)
        {

            double x= (i-N_X/2)*dx;
            double y= (j-N_Y/2)*dy;

            div_[i][j]=0;
            rho_[i][j]=0.0;
            n_1[i][j]=exp(-(x*x+y*y)/((N_Y*0.1*dy)*(N_Y*0.1*dy)));
            n_1_prev[i][j]=n_1[i][j];

            n_2[i][j]=10e10;//10e10;
            p_in[i][j]=0;
            //  out_Ux[i][j]=sqrt(dy*j)*U;
            phi_[i][j]=0.0;//sin(i*6*M_PI/(N_X-1));

            if ((i==0)||(i==N_X-1)||(j==0)||(j==N_Y-1))
            {
                n_1[i][j]=0.0;
                n_2[i][j]=0.0;
            }


        }


    }

    //   field_phi(5, phi_, 1000);

    for (i=0;i<N_X;i++)
    {
        for (j=shift;j<N_Y_DIEL;j++)
        {

            phi_[i][j]=0;
            Py_[i][j]=0.26 * sin(6.28*i*15.0/N_X /*+ 6.28*/);//rand()*0.52/RAND_MAX - 0.26;
            Py0_[i][j]=Py_[i][j];

            Px_[i][j]=0.0;
            Px0_[i][j]=0.0;

        }
    }
    //  field_U(0.0, Ux_,500);

    /*  for (i=0; i<N_X; i++)
    {
        out_Ux[i][0]=0.0;
        out_Ux[j][N_Y]=0.0;
    }
  */
}
//double ddx()
void  get_rhs_NAVIER();



void advect()
{
    for (int i=1; i<N_X-1; i++)
    {
        for (int j=1; j<N_Y-1; j++)
        {
            n_1_prev[i][j]=n_1[i][j];
        }
    }

}

double tt=0.0;
void sweep()
{
    double rho_min=1e30;
    double rho_max=-1e30;

    tt+=1;


    double vol= fabs(dx*dy*(w_z1-w_z0));//volume of a cell
    for (int i=1; i<N_X-1; i++)
    {
        for (int j=1; j<N_Y-1; j++)
        {


            // rho_[i][j]=-n_1[i][j];//q_e*(-n_1[i][j])/eps_0;;//q_e*(-n_1[i][j]+n_2[i][j])/eps_0;
            RHS[i][j]=(Py_[i][j+1]-Py_[i][j-1])/(2.0*dy)/eps_0
                    /*+(Px_[i+1][j]-Px_[i-1][j])/(2.0*dx)/eps_0 ;*/ - gau[j]*q[i]/(eps_0*vol);


            div_[i][j]=RHS[i][j];
            if (div_[i][j]>rho_max) rho_max=div_[i][j];
            if (div_[i][j]<rho_min) rho_min=div_[i][j];

            // Py_[i][j]=0.0;
        }

    }



    //   printf("rho_min = %e rho_max= %e q=%e  \n", rho_min,rho_max,-0.26/dy/eps_0);


    for (int i=0;i<N_X;i++)
    {
        //phi_[i][N_Y-1]=0.0;
        phi_[i][0]=0;//*(1.0-i*1.0/N_X);

        //  phi_[i][N_Y-1]=sasign *50;//* sin(tt/3000);
        //if ((i<20)&&(i>0)) phi_[i][N_Y-1]=-A/4000;
        //if ((i>=60)&&(i<N_X-1)) phi_[i][0]=A/4000;
    }


    for (int j=N_Y_DIEL;j<N_Y;j++)
    {

        phi_[0][j]=-sasign*15.5;

    }

    /* for (int nnn=0;nnn<13;nnn++)
    {
        for (int i=1;i<N_X-1;i++)
        {
            // phi_[i][j]=A*sin(i*6*M_PI/(N_X-1)+1800.0*alpha*t*dt);
            phi_[i][N_Y-1]=(phi_[i][N_Y-1]+phi_[i-1][N_Y-1]+phi_[i+1][N_Y-1])/3.0;//-A;//.0;//A*sin(i*6*M_PI/(N_X-1)+1800.0*alpha*t*dt);
            phi_[i][0]=(phi_[i-1][0]+phi_[i][0]+phi_[i+1][0])/3.0;//A*sin(i*6*M_PI/(N_X-1)+1800.0*alpha*t*dt);
        }
    }*/


    INPUT_PARAM par;
    par.a=((2.0)/(dx*dx)+2.0/(dy*dy));
    par.bp=-1.0/(dx*dx);
    par.bm=-1.0/(dx*dx);
    par.cp=-1.0/(dy*dy);
    par.cm=-1.0/(dy*dy);

    par.w_bc_type=3;
    par.e_bc_type=1;
    par.n_bc_type=1;
    par.s_bc_type=3;

    par.w_bc_val=0.0;
    par.e_bc_val=0.0;
    par.n_bc_val=0.0;
    par.s_bc_val=0.0;

    //jacobi(par,phi_,div_,20);//phi
    multigrid_N(par,phi_,RHS,3,6);
    multigrid_N(par,phi_,RHS,3,6);
    /*multigrid_N(par,phi_,RHS,3,6);

    multigrid_N(par,phi_,RHS,3,6);
    multigrid_N(par,phi_,RHS,3,6);
    multigrid_N(par,phi_,RHS,3,6);*/


    //       multigrid_N(par,phi_,RHS,3,6);
    /*multigrid_N(par,phi_,RHS,8,3);
    multigrid_N(par,phi_,RHS,8,3);
    multigrid_N(par,phi_,RHS,8,3);
    multigrid_N(par,phi_,RHS,8,3);*/

    //

    //getE();



    /*n_1[i][j]*(1.0/dt+D*2.0/(dx*dx) + D*2.0/(dy*dy))
            - D*((n_1[i+1][j]+n_1[i-1][j])/(dx*dx) +(n_1[i][j+1]+n_1[i][j-1])/(dy*dy))
            = n_1_prev[i][j]/dt
                   +mu*(n_1[i][j]*rho_[i][j]+E x[i][j]*(n_1[i+1][j]-n_1[i-1][j])/(2.0*dx) + E y[i][j]*(n_1[i][j+1]-n_1[i][j-1])/(2.0*dy));*/


    double mu=1.0;
    double D=1.0;

    INPUT_PARAM par2;
    par2.a=(1.0/dt+D*2.0/(dx*dx) + D*2.0/(dy*dy));
    par2.bp=-D/(dx*dx);
    par2.bm=-D/(dx*dx);
    par2.cp=-D/(dy*dy);
    par2.cm=-D/(dy*dy);

    par2.w_bc_type=2;
    par2.e_bc_type=2;
    par2.n_bc_type=0;
    par2.s_bc_type=0;

    par2.w_bc_val=0.0;
    par2.e_bc_val=0.0;
    par2.n_bc_val=0.0;
    par2.s_bc_val=0.0;


    /*
    for (int i=1; i<N_X-1; i++)
    {
        for (int j=1; j<N_Y-1; j++)
        {

            int ip=(E_x[i][j]<0);
            int im=(E_x[i][j]>=0);
            int jp=(E_y[i][j]<0);
            int jm=(E_y[i][j]>=0);

            RHS_p[i][j]=n_1_prev[i][j]/dt
                    - mu*(n_1[i][j]*rho_[i][j]+ E_x[i][j]*(n_1[i+ip][j]-n_1[i-im][j])/(dx) + E_y[i][j]*(n_1[i][j+jp]-n_1[i][j-jm])/(dy));
        }
    }

    jacobi(par2,n_1,RHS_p,100);


    advect();
*/

    for (int i=0; i<N_X; i++)
    {
        for (int j=0; j<N_Y; j++)
        {
            Ex[i][j]=-ddx(phi_,i,j);
            Ey[i][j]=-ddy(phi_,i,j);
        }

    }

    double xix=1.17*eps_0;
    for (int i=1; i<N_X-1; i++)
    {
        for (int j=1; j<N_Y_DIEL; j++)
        {
            //Px0_[i][j]=Px_[i][j];
            Px_[i][j]=xix*Ex[i][j];
        }
    }

    INPUT_PARAM par_ferr;
    double kap=1.38e-10;
    par_ferr.a=(1.0/(1e5*dt))+(kap*2.0/(dx*dx) + kap*2.0/(dy*dy));
    par_ferr.bp=-kap/(dx*dx);
    par_ferr.bm=-kap/(dx*dx);
    par_ferr.cp=-kap/(dy*dy);
    par_ferr.cm=-kap/(dy*dy);

    par_ferr.w_bc_type=1;
    par_ferr.e_bc_type=0;
    par_ferr.n_bc_type=1;
    par_ferr.s_bc_type=1;

    par_ferr.w_bc_val=0.26;
    par_ferr.e_bc_val=0.0;
    par_ferr.n_bc_val=0.0;
    par_ferr.s_bc_val=0.0;

    poly p;

    double alp,bet,gam,T,T0,rh;
    alp=3.324e5;
    bet=6.381e8;
    gam=7.89e9;
    T=300;
    T0=381;
    rh=0.0;

    double x0=-0.3;
    p.order=5;
    p.C[0]=2*alp*(T-T0);//(T-T0); //x
    p.C[1]=0.0;         //xx
    p.C[2]=-4.0*bet;   //xxx
    p.C[3]=0.0;        //x^4
    p.C[4]=6.0*gam;

    //-(fiy*0.0033- 2.0*alp*81*P1 - 4*bet*P1^3 +6*gam*P1^5)


    double emin=1e23;
    double emax=-1e23;
    for (int i=1; i<N_X-1; i++)
    {
        for (int j=1; j<N_Y-1; j++)
        {

            RHS_p[i][j]= /*1.4e7*(1.0-i*1.0/N_X)*/ -Ey[i][j]+Py0_[i][j]/(1e5*dt);

        }
    }

    //  printf("emin=%e emax=%e \n",emin,emax);

    jacobi_polynomial( par_ferr, p,Py_,RHS_p, 5);

    emin=1e23;
    emax=-1e23;
    for (int i=1; i<N_X-1; i++)
    {
        for (int j=1; j<N_Y-1; j++)
        {

            if (emin>Py_[i][j]) emin=Py_[i][j];
            if (emax<Py_[i][j]) emax=Py_[i][j];
        }
    }

    for (int i=1; i<N_X-1; i++)
    {
        for (int j=1; j<N_Y-1; j++)
        {
            Py0_[i][j]=Py_[i][j];
        }
    }


    ///px
    ///
    ///

    INPUT_PARAM par_ferr_x;

    par_ferr.a=(1.0/dt);
    par_ferr.bp=0.0;
    par_ferr.bm=0.0;
    par_ferr.cp=0.0;
    par_ferr.cm=0.0;

    par_ferr.w_bc_type=1;
    par_ferr.e_bc_type=1;
    par_ferr.n_bc_type=1;
    par_ferr.s_bc_type=1;

    par_ferr.w_bc_val=0.0;
    par_ferr.e_bc_val=0.0;
    par_ferr.n_bc_val=0.0;
    par_ferr.s_bc_val=0.0;

    poly p_x;




    p_x.order=1;
    p.C[0]=1/xix;//0.0;//(T-T0); //x
    //p.C[1]=0.0;         //xx

    //-(fiy*0.0033- 2.0*alp*81*P1 - 4*bet*P1^3 +6*gam*P1^5)

    for (int i=1; i<N_X-1; i++)
    {
        for (int j=1; j<N_Y-1; j++)
        {

            RHS_p[i][j]=  -Ex[i][j]+Px0_[i][j]/dt;

        }
    }

    /* for (int i=0; i<N_X; i++)
    {
        RHS_p[i][0]=(Px_[i][1]-Px_[i][0])/dy;
        RHS_p[i][N_Y_DIEL-1]=(Px_[i][N_Y_DIEL-1]-Px_[i][N_Y_DIEL-2])/dy;
    }*/

    //  printf("emin=%e emax=%e \n",emin,emax);

    //  jacobi_polynomial( par_ferr_x, p_x,Px_,RHS_p, 5);




    //  printf("pmin=%e pmax=%e \n",emin,emax);
    if (move_particles)
    {
        //for (int i=1;i<10;i++)
        fmm_step(dt);
    }

}


void init()
{
    glClearColor (0.0, 0.0, 0.0, 0.0);
    glColor3f(1.0, 1.0, 1.0);
    glMatrixMode (GL_PROJECTION);
    glLoadIdentity ();
    glOrtho(w_x0-(w_x1-w_x0)*0.01, w_x1*1.01, w_y0*1.05,w_y1*1.05, -10.0, 10.0);
    glMatrixMode (GL_MODELVIEW);
    rand_init();
    sweep_init();

    fmm_init();

    poly p;

    double alp,bet,gam,T,T0,rh;
    alp=3.324e5;
    bet=6.381e8;
    gam=7.89e9;
    T=300;
    T0=381;



    double x1=0;
    double x0=0;
    double d_t=1e-11;
    for (int i=0;i<2000;i++)
    {
        rh=1.4e7+x0/d_t;

        p.order=5;
        p.C[0]=2*alp*(T-T0)+1.0/d_t; //x
        p.C[1]=0.0;         //xx
        p.C[2]=-4.0*bet;   //xxx
        p.C[3]=0.0;        //x^4
        p.C[4]=6.0*gam;

        x1=solve_poly(p,x1, rh,10);
        x0=x1;
        if (i%100==0)
            printf("i=%d x=%e \n",i,x1);
    }

    /*
    solve_poly(p,x0, rh,2);
    solve_poly(p,x0, rh,3);
    solve_poly(p,x0, rh,4);
    solve_poly(p,x0, rh,5);
    solve_poly(p,x0, rh,6);
    solve_poly(p,x0, rh,7);

    solve_poly(p,x0, rh,8);
    solve_poly(p,x0, rh,9);
    solve_poly(p,x0, rh,10);
    solve_poly(p,x0, rh,11);
*/



}



int main(int argc, char** argv)
{
    srand(time(NULL));
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
