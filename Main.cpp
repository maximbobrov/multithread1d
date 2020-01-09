
#include <stdio.h>
#include <stdlib.h>


#include <my_include/gl.h>
#include <my_include/glu.h>
#include <my_include/glut.h>
//#include <GL/glut.h>
#include  <math.h>
#include <time.h>
#include "globals.h"
#include <thread>
#include<xmmintrin.h>
#include <vector>
#include <iostream>
#include <algorithm>

#include "phi_mult.h"


//#include <sys/time.h>

#include "sse_sum.h"

void save_fields();
void load_fields();




void display(void);
void sweep_init();
void init();

void fmm_step(int threadIdx, vec3 *bodyAccel_, vec4 *bodyPos_, vec3 *bodyVel_, double dt);
void sweep();
//void sweep_old();
int view=VIEW_PHI;
int view_v=VIEW_JX;
int redr=0;

double ck=2.0;

double cv=0.001;

double conv[100];
double conv1[100];
float randMax = -1e30;
float randMin = 1e30;
/*#define MAIN
#include "fmm.h"
#undef MAIN
*/

const int maxParticles=8192;
//int numParticles=0;//4096;

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




void display(void)
{

    double sc = 0.4e6;
    if (redr==1)
    {
        /*    std::vector <std::thread> th_vec;
        for (int i = 0; i < num_thread; ++i)
        {
            th_vec.push_back(std::thread(fmm_step,i , bodyAccel[i], bodyPos[i], bodyVel[i],dt));
        }
        for (int i = 0; i < num_thread; ++i)
        {
            th_vec.at (i).join ();
        }*/
        int i=0;
        fmm_step(i , bodyAccel[i], bodyPos[i], bodyVel[i],dt);
    }
    int i,j;//,k,l;


    double l_2;


    // if (clear_w==1)
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


    for (i=0;i<N_X-1;i++)
    {
        glBegin(GL_TRIANGLE_STRIP);
        view=VIEW_PHI;
        for (j=0;j<N_Y;j++)
        {


            if (view==VIEW_RHO)
                l_2=ck*(n_1[i][j]);
            if (view==VIEW_RHO_2)
                l_2=ck*(n_2[i][j]/1.0e10);

            if (view==VIEW_PHI)
                l_2=ck*(phi_[i][j]/A);
            if (view==VIEW_JX)
                l_2=ck*(div_J[i][j]);
            if (view==VIEW_JY)
                l_2=ck*(Jy_[i][j]);
            if (view==VIEW_P)
                l_2=ck*(p_in[i][j]);
            /*  if (view==VIEW_fieldUX)
                l_2=ck*(field_Ux[i][j]);
            if (view==VIEW_fieldUY)
                l_2=ck*(field_Uy[i][j]);
            */
            glColor3f(l_2,l_2,-l_2);
            glVertex2f(sc*dx*(i-N_X/2),sc*dy*(j-N_Y/2));



            if (view==VIEW_RHO)
                l_2=ck*(n_1[i+1][j]);
            if (view==VIEW_RHO_2)
                l_2=ck*(n_2[i+1][j]/1.0e10);
            if (view==VIEW_PHI)
                l_2=ck*(phi_[i+1][j]/A);
            if (view==VIEW_JX)
                l_2=ck*(div_J[i+1][j]);
            if (view==VIEW_JY)
                l_2=ck*(Jy_[i+1][j]);
            if (view==VIEW_P)
                l_2=ck*(p_in[i+1][j]);
            /*if (view==VIEW_fieldUX)
                l_2=ck*(field_Ux[i+1][j]);
            if (view==VIEW_fieldUY)
                l_2=ck*(field_Uy[i+1][j]);
*/
            glColor3f(l_2,l_2,-l_2);
            glVertex2f(sc*dx*(i+1-N_X/2),sc*dy*(j-N_Y/2));
        }


        glEnd();

    }
    /* glBegin(GL_LINES);
        for (j=0;j<N_Y;j++)
        {


            glColor3f(1,1,1);
            glVertex2f(dx*i,dy*j);



            glColor3f(0,0,0);
            glVertex2f(dx*(i)+cv*Jx_[i][j],dy*j+cv*Jy_[i][j]);
        }
        glEnd();
      */

    //glEnable(GL_LINE_SMOOTH);

    for (int ii=0;ii<N_X;ii++)
    {
        i=ii;
        glBegin(GL_LINES);
        for (int jj=0;jj<8;jj++)
        {

            j=jj*8+1;

            glColor3f(1,1,1);
            glVertex2f(sc*dx*(i-N_X/2),sc*dy*(j-N_Y/2));


            glColor3f(0.5,0.5,0.5);
            if (view_v==VIEW_JX)  glVertex2f(sc*dx*(i-N_X/2)+cv*Jx_[i][j], sc*dy*(j-N_Y/2)+cv*Jy_[i][j]);
            if (view_v==VIEW_UX)  glVertex2f(sc*dx*(i-N_X/2)+cv*Ux_[i][j], sc*dy*(j-N_Y/2)+cv*Uy_[i][j]);
        }
        glEnd();
    }

    // glEnable(GL_BLEND);

    // glEnable(GL_POINT_SMOOTH);
    glPointSize(3.0);

    /*for( i=0; i<numParticles; i++ ) {
        bodyAccelCommon[i].x = 0;
        bodyAccelCommon[i].y = 0;
        bodyAccelCommon[i].z = 0;

        bodyVelCommon[i].x = 0;
        bodyVelCommon[i].y= 0;
        bodyVelCommon[i].z= 0;

        bodyPosCommon[i].x = 0;
        bodyPosCommon[i].y = 0;
        bodyPosCommon[i].z = 0;
        bodyPosCommon[i].w = 0;
    }*/
    for( i=0; i<N_X; i++ )
    {
        for( j=0; j<N_Y; j++ )
        {
            Energy[i][j] = 0.0;
        }
    }
    for( int k=0; k<num_thread; k++ )
    {
        for( i=0; i<numParticles[k]; i++ )
        {
            int x = (bodyPos[k][i].x + dx*(N_X/2))/dx;
            int y = (bodyPos[k][i].y + dy*(N_X/2))/dy;
            Energy[x][y] += (bodyVel[k][i].x*bodyVel[k][i].x)+(bodyVel[k][i].y*bodyVel[k][i].y)+(bodyVel[k][i].z*bodyVel[k][i].z);
            bodyAccelCommon[i] += bodyAccel[k][i];
            bodyVelCommon[i] +=  bodyVel[k][i];
            bodyPosCommon[i] +=  bodyPos[k][i];
        }
    }
    /*for( i=0; i<numParticles; i++ )
    {
        bodyAccelCommon[i] /= float(num_thread);
        bodyVelCommon[i] /= float(num_thread);
        bodyPosCommon[i] /= float(num_thread);

    }*/
    double maxEn=-1e30;
    double minEn=+1e30;
    double maxB=-1e30;
    double minB=+1e30;
    for( i=0; i<N_X; i++ )
    {
        minB = BoundaryLayer[i] < minB ? BoundaryLayer[i] : minB;
        maxB = BoundaryLayer[i] > maxB ? BoundaryLayer[i] : maxB;
        for( j=0; j<N_Y; j++ )
        {
            minEn = Energy[i][j] < minEn ? Energy[i][j] : minEn;
            maxEn = Energy[i][j] > maxEn ? Energy[i][j] : maxEn;
        }
    }

    glPointSize(5);


    /*for( i=0; i<N_X; i++ )
    {
        for( j=0; j<N_Y; j++ )
        {
            glColor3f(ck*(Energy[i][j]-minEn) * 1.0 / (maxEn-minEn),0,0);
            glBegin(GL_QUADS);

            glVertex3f(i*dx + dx*(-N_X/2),j*dy + dy*(-N_Y/2),0);
            glVertex3f((i+1)*dx + dx*(-N_X/2),j*dy + dy*(-N_Y/2),0);
            glVertex3f((i+1)*dx + dx*(-N_X/2),(j+1)*dy + dy*(-N_Y/2),0);
            glVertex3f(i*dx + dx*(-N_X/2),(j+1)*dy + dy*(-N_Y/2),0);

            glEnd();

        }
    }*/

    for( i=0; i<N_X; i++ )
    {
        //glColor3f(ck*(BoundaryLayer[i]/*-minB*/) * 1.0 / (maxB-minB),-ck*(BoundaryLayer[i]/*-minB*/) * 1.0 / (maxB-minB),0);
        glColor3f(ck*(BoundaryLayer[i]),-ck*(BoundaryLayer[i]),0);
        j = 0;
        glBegin(GL_QUADS);
        glVertex3f(sc*(i*dx + dx*(-N_X/2)),sc*(j*dy + dy*(-N_Y/2)),0);
        glVertex3f(sc*((i+1)*dx + dx*(-N_X/2)),sc*(j*dy + dy*(-N_Y/2)),0);
        glVertex3f(sc*((i+1)*dx + dx*(-N_X/2)),sc*((j-3)*dy + dy*(-N_Y/2)),0);
        glVertex3f(sc*(i*dx + dx*(-N_X/2)),sc*((j-3)*dy + dy*(-N_Y/2)),0);

        glEnd();
    }

    glPointSize(1);

    glBegin(GL_POINTS);

    //for( int k=0; k<num_thread; k++ )
    {
        int k=0;
        for( i=0; i<numParticles[k]; i++ ) {

            glColor3f(fabs(bodyVel[k][i].x)/400.0,fabs(bodyVel[k][i].y)/400.0,fabs(bodyVel[k][i].z)/400.0);
            glVertex3f(sc*bodyPos[k][i].x,sc*bodyPos[k][i].y,sc*bodyPos[k][i].z);
            //glColor3f(fabs(bodyVel[0][i].x)/400.0,fabs(bodyVel[0][i].y)/400.0,fabs(bodyVel[0][i].z)/400.0);
            //glVertex3f(bodyPos[0][i].x,bodyPos[0][i].y,bodyPos[0][i].z);

        }
    }
    glEnd();



    glColor3f(1,1,1);

    glBegin(GL_LINE_LOOP);

    glVertex3f(sc*dx*(-N_X/2),sc*dy*(-N_Y/2),0);
    glVertex3f(sc*dx*(N_X-1-N_X/2),sc*dy*(-N_Y/2),0);
    glVertex3f(sc*dx*(N_X-1-N_X/2),sc*dy*(N_Y-1-N_Y/2),0);
    glVertex3f(sc*dx*(-N_X/2),sc*dy*(N_Y-1-N_Y/2),0);
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

double w_x0=-dx * N_X /2;
double w_x1=dx * N_X /2;


double w_y0=-dy * N_Y /2;
double w_y1=dy * N_Y /2;

double w_z0=-dz/2;
double w_z1=dz/2;

//FmmSystem tre;


uint32_t xr[32],yr[32],zr[32],wr[32];

void rand_init()
{
    for (int i=0;i<num_thread;i++)
    {
        xr[i] = 123456789+i;
        yr[i] = 362436069+2*i;
        zr[i] = 521288629+3*i;
        wr[i] = 88675123+5*i;
    }
}

float my_rand(int i) {
    uint32_t t;
    t = xr[i] ^ (xr[i] << 11);
    xr[i] = yr[i]; yr[i] = zr[i]; zr[i] = wr[i];
    wr[i] = wr[i] ^ (wr[i] >> 19) ^ (t ^ (t >> 8));
    return wr[i]*0.5/0x7FFFFFFF;
}

vec3 getE(float x, float y)
{

    double a0=fmax(fmin(N_X-1,(N_X-1)*(x-w_x0)/(w_x1-w_x0)),0);
    int i0=(int)a0;
    a0-=i0;


    double b0=fmax(fmin(N_Y-1,(N_Y-1)*(y-w_y0)/(w_y1-w_y0)),0);
    int j0=(int)b0;
    b0-=j0;

    vec3 ret;

    ret.x=(1.0-b0)*((1.0-a0)*Ex[i0][j0]+(a0)*Ex[i0+1][j0]) + (b0)*((1.0-a0)*Ex[i0][j0+1]+(a0)*Ex[i0+1][j0+1]);
    ret.y=(1.0-b0)*((1.0-a0)*Ey[i0][j0]+(a0)*Ey[i0+1][j0]) + (b0)*((1.0-a0)*Ey[i0][j0+1]+(a0)*Ey[i0+1][j0+1]);
    ret.z;
    return ret;
}


vec4 init_pos(int threadIdx)
{
    vec4 res;
    res.x = dx * N_X / 2;
    res.y = my_rand(threadIdx)*dy * N_Y-dy * N_Y / 2;
    res.z = my_rand(threadIdx)*dz - dz / 2;
    res.w = 150;//my_rand(threadIdx);




    return res;
}

vec3 init_vel(int threadIdx)
{
    vec3 vel;
    double alpha = M_PI/4;
    double angle = my_rand(threadIdx) * alpha - alpha/2;
    double velMagn = 6e6;

    vel.x = velMagn * cos(angle);
    vel.y = velMagn * sin(angle);
    vel.z = my_rand(threadIdx)*10.0-5.0;


    return vel;
}

void create_random_particles(int threadIdx, vec4 *bodyPos_, vec3 *bodyVel_,vec3 *bodyAccel_)
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
    for (int i=0;i<(N_Y-1) * 1.0 / 2;i++)
    {
        double t = 1.1;
        double B = 1.0;
        double x_ = - N_X * dx / 2;
        double y_ = i * dy - N_Y * dy / 4;
        double z_ = my_rand(threadIdx) * dz - dz / 2;
        //vec3 E = getE(x_, y_);
        vec3 E;
        getEFromElectrons(E, bodyPos_, x_, y_, z_, numParticles[threadIdx]);
        printf("Ex=%e Ey = %e\n", E.x, E.y);
        E.x += Ex[0][i];
        E.y += Ey[0][i];
        E.x = fmax(0.1, -E.x);
        //printf("Ex=%e Ey = %e\n", E.x, E.y);
        double phi = 4.0;
        double y = 3.79 * 1e-4 * sqrt(fabs(B * E.x)) / phi;
        double tetta = 0.95 - 1.03 * y * y;
        double J = (1.54 * 1e-6 * B * B * E.x * E.x / (t * t * phi)) * exp( - 6.83 * 1e7 * pow(phi, 1.5) * tetta / fabs( B * E.x));
        //printf("JJJJ=%f jjj = %e E = %e \n",(dt * J * dy * dz * 1e-4 * 6.24151 * 1e18 / 1e5), J , E.x);
        int curNum = numParticles[threadIdx];
        int numToAdd = std::min(int(dt * J * dy * dz * 1e-4 * 6.24151 * 1e18 / 1e3), maxParticles - numParticles[threadIdx]-2);
        if(my_rand(threadIdx) < (dt * J * dy * dz * 1e-4 * 6.24151 * 1e18 / 1e3) - int(dt * J * dy * dz * 1e-4 * 6.24151 * 1e18 / 1e3))
          numToAdd++;
        numParticles[threadIdx] += numToAdd;
        for (int n = curNum; n < curNum + numToAdd; ++n)
        {
            bodyPos_[n].x=x_;
            bodyPos_[n].y=y_;
            bodyPos_[n].z=z_;
            bodyPos_[n].w=1.6e-19 * 1e3;
            double angle = 0.98*acos(1-2*my_rand(threadIdx)) - M_PI / 2;
            double en = 1.1e5;
            bodyVel_[n].x = en * cos(angle);
            bodyVel_[n].y = en * sin(angle);
            bodyVel_[n].z = 0;
            bodyAccel_[n].x=0.0;
            bodyAccel_[n].y=0.0;
            bodyAccel_[n].z=0.0;
        }
    }
}

void delete_particle(int threadIdx, int particlesIdx, vec3 *bodyAccel_, vec4 *bodyPos_, vec3 *bodyVel_)
{
    bodyPos_[particlesIdx] = bodyPos_[numParticles[threadIdx]-1];
    bodyVel_[particlesIdx] = bodyVel_[numParticles[threadIdx]-1];
    bodyAccel_[particlesIdx] = bodyAccel_[numParticles[threadIdx]-1];
    numParticles[threadIdx] -= 1;
    //float r = my_rand(threadIdx);
    //randMax  = randMax > r ? randMax : r;
    // randMin  = randMin < r ? randMin : r;
}

void wall_colision(int threadIdx, int particlesIdx, vec3 *bodyAccel_, vec4 *bodyPos_, vec3 *bodyVel_)
{
    int xIdx = (bodyPos[threadIdx][particlesIdx].x + dx*(N_X/2))/dx;


    if ((xIdx>=0)&&(xIdx<N_X))
    {
        double me=9.1e-31;
        double Ev_in_J=1.6e-19;

        double v2=((bodyVel[threadIdx][particlesIdx].x*bodyVel[threadIdx][particlesIdx].x)+(bodyVel[threadIdx][particlesIdx].y*bodyVel[threadIdx][particlesIdx].y)+(bodyVel[threadIdx][particlesIdx].z*bodyVel[threadIdx][particlesIdx].z));
        double E = 0.5*(me/Ev_in_J)*v2;
        double Em = 650;// in ev

        double phi  = atan2(bodyVel[threadIdx][particlesIdx].y, bodyVel[threadIdx][particlesIdx].x);
        double Z = 1.284 * E / (Em * (1 + phi * phi/(M_PI*M_PI)));
        double sigma_m = 6.4;
        double sigma = 1.526 * (1 + phi * phi/(4*M_PI*M_PI)) * (1 - exp(-pow(Z, 1.725))) * sigma_m /(pow(Z, 1.725)+0.1);

        int num = (int)(fabs(sigma));

        if(my_rand(threadIdx) < (sigma - 1.0 * num))
            num++;

        //printf("eeee=%e num=%d \n",E,num);


        if(maxParticles - numParticles[threadIdx]> num)
        {
            int size = std::min(maxParticles - numParticles[threadIdx], num);
            //for (int i = 0; i < size; ++i) {
            int i=numParticles[threadIdx];
            int n=0;
            while((i<maxParticles)&&(n<num))
            {

                bodyPos_[i].x= bodyPos_[particlesIdx].x+3.0*(my_rand(threadIdx)*dx*2-dx);
                bodyPos_[i].y= bodyPos_[particlesIdx].y;
                bodyPos_[i].z= bodyPos_[particlesIdx].z+3.0*(my_rand(threadIdx)*dx*2-dx);
                bodyPos_[i].w= bodyPos_[particlesIdx].w;
                double angleNew = 0.98*acos(1-2*my_rand(threadIdx))+0.01*M_PI;

                bodyVel_[i].x =sqrt(fabs(v2/(num))) * cos(angleNew);//0.5*sqrt(fabs(E/(num))) * cos(angleNew);
                bodyVel_[i].y =sqrt(fabs(v2/(num))) * sin(angleNew); //0.5*sqrt(fabs(E/(num))) * sin(angleNew);

                bodyVel_[i].z = 0.0;

                BoundaryLayer[xIdx] -= 1.0;
                i++;
                n++;

                numParticles[threadIdx]=i;
            }

            delete_particle( threadIdx, particlesIdx, bodyAccel_, bodyPos_, bodyVel_);
            BoundaryLayer[xIdx] += 1.0;

        }else
        {
            bodyPos_[particlesIdx].y -= dt*bodyPos_[particlesIdx].y;
            bodyVel_[particlesIdx].y=fabs(bodyVel_[particlesIdx].y);
        }
        /*bodyPos_[particlesIdx].y -= dt*bodyPos_[particlesIdx].y;
        bodyVel_[particlesIdx].y=fabs(bodyVel_[particlesIdx].y);
        double angle =acos(1-2*my_rand(threadIdx));// * M_PI;//asin(bodyVel_[particlesIdx].y/bodyVel_[particlesIdx].x);
        bodyVel_[particlesIdx].x=sqrt(E-E_excit)*cos(angle);
        bodyVel_[particlesIdx].y=sqrt(E-E_excit)*sin(angle);*/
        //bodyVel_[numParticles[threadIdx]-1].z = my_rand(threadIdx)*10.0-5.0;


        /*numParticles[threadIdx] ++;
        bodyPos_[numParticles[threadIdx]-1].x= bodyPos_[particlesIdx].x;
        bodyPos_[numParticles[threadIdx]-1].y= bodyPos_[particlesIdx].y;
        bodyPos_[numParticles[threadIdx]-1].z= bodyPos_[particlesIdx].z;
        bodyPos_[numParticles[threadIdx]-1].w= bodyPos_[particlesIdx].w;
        double angleNew = my_rand(threadIdx) * M_PI;
        bodyVel_[numParticles[threadIdx]-1].x = sqrt(E) * cos(angleNew);
        bodyVel_[numParticles[threadIdx]-1].y = -sqrt(E) * sin(angleNew);
        bodyVel_[numParticles[threadIdx]-1].z = my_rand(threadIdx)*10.0-5.0;

        BoundaryLayer[xIdx] -= 1.0;*/
        /*
    else
    {
        if(my_rand(threadIdx)>0.9)
        {
            delete_particle( threadIdx, particlesIdx, bodyAccel_, bodyPos_, bodyVel_);
            BoundaryLayer[xIdx] += 1.0;
        }
        else
        {
        bodyPos_[particlesIdx].y -= dt*bodyPos_[particlesIdx].y;
        bodyVel_[particlesIdx].y=0.98*fabs(bodyVel_[particlesIdx].y);
        }
    }*/
    }
}

void fmm_step(int threadIdx, vec3 *bodyAccel_, vec4 *bodyPos_, vec3 *bodyVel_, double dt)
{
    int i;
    double tic,toc,timeFMM;


    create_random_particles(threadIdx, bodyPos_, bodyVel_,bodyAccel_);
    tic = get_time();
    direct(bodyAccel_, bodyPos_, bodyVel_, numParticles[threadIdx]);
    toc = get_time();
    timeFMM = toc-tic;
    printf("fmm    : %g\n",timeFMM);

    //printf("aa = %f  bb = %f\n",randMin, randMax);

    double q_over_m = 1.0;//1.6e-19/9.1e-31;
    for( i=0; i<numParticles[threadIdx]; i++ )
    {
        vec3 ev=getE(bodyPos_[i].x,bodyPos_[i].y);
        bodyVel_[i].x -= -dt*q_over_m*(ev.x+ bodyAccel_[i].x);
        bodyVel_[i].y -= -dt*q_over_m*(ev.y +bodyAccel_[i].y);
        bodyVel_[i].z -= dt*bodyAccel_[i].z;
        /*bodyVel_[i].x *= 0.9995;
        bodyVel_[i].y *= 0.9995;
        bodyVel_[i].z *= 0.9995;*/
    }

    for( i=0; i<numParticles[threadIdx]; i++ )
    {
        bodyPos_[i].x += dt*bodyVel_[i].x;
        bodyPos_[i].y += dt*bodyVel_[i].y;
        bodyPos_[i].z += dt*bodyVel_[i].z;

        if (bodyPos_[i].y<w_y0)
        {
            /*bodyPos_[i].y -= dt*bodyVel_[i].y;
            bodyVel_[i].y=fabs(bodyVel_[i].y);*/
            bodyPos_[i].y -= dt*bodyVel_[i].y;
            wall_colision( threadIdx, i, bodyAccel_, bodyPos_, bodyVel_);
        }

        if (bodyPos_[i].y>w_y1 || bodyPos_[i].x < w_x0 || bodyPos_[i].x>w_x1)
        {
            delete_particle( threadIdx, i, bodyAccel_, bodyPos_, bodyVel_);
        }

        if (bodyPos_[i].z<w_z0)
        {
            bodyPos_[i].z -= dt*bodyVel_[i].z;
            bodyVel_[i].z=fabs(bodyVel_[i].z);
        }

        if (bodyPos_[i].z>w_z1)
        {
            bodyPos_[i].z -= dt*bodyVel_[i].z;
            bodyVel_[i].z=-fabs(bodyVel_[i].z);
        }
    }
}

void fmm_init()
{
    rand_init();
    bodyAccelCommon = new vec3[maxParticles];
    bodyVelCommon = new vec3[maxParticles];
    bodyPosCommon = new vec4[maxParticles];

    numParticles = new int[num_thread];

    bodyAccel = new vec3*[num_thread];
    bodyVel = new vec3*[num_thread];
    bodyPos = new vec4*[num_thread];
    for (int k = 0; k < num_thread; ++k)
    {
        numParticles[k] = 0;
        bodyAccel[k] = new vec3[maxParticles];
        bodyVel [k] = new vec3[maxParticles];
        bodyPos[k] = new vec4[maxParticles];
    }

    for( int i=0; i<N_X; i++ )
    {
        BoundaryLayer[i] = 0.0;
    }
}

void kb(unsigned char key, int x, int y)
{
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
        cv*=1.1;
    }

    if (key=='[')
    {
        cv/=1.1;
    }

    if (key=='s')
    {
        sweep();
    }

    if (key==' ')
    {
        redr=!redr;
    }
    glutPostRedisplay();
}

void sweep_init()
{
    int i,j,k,n;
    double App,R_2,b,m;

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

            n_2[i][j]=10e10;
            p_in[i][j]=0;
            phi_[i][j]=sin(i*6*M_PI/(N_X-1));

            if ((i==0)||(i==N_X-1)||(j==0)||(j==N_Y-1))
            {
                n_1[i][j]=0.0;
                n_2[i][j]=0.0;
            }
        }
    }

    for (i=0;i<N_X;i++)
    {
        for (j=0;j<N_Y-1;j++)
        {

            phi_[i][j]=0;

        }
    }
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





void sweep()
{
    double rho_min=1e30;
    double rho_max=-1e30;

    for (int i=1; i<N_X-1; i++)
    {
        for (int j=1; j<N_Y-1; j++)
        {
            rho_[i][j]=-n_1[i][j];//q_e*(-n_1[i][j])/eps_0;;//q_e*(-n_1[i][j]+n_2[i][j])/eps_0;
            RHS[i][j]=-20000.5*n_1[i][j];

            if (rho_[i][j]>rho_max) rho_max=rho_[i][j];
            if (rho_[i][j]<rho_min) rho_min=rho_[i][j];
        }

    }

    printf("rho_min = %e rho_max= %e \n", rho_min,rho_max);


    for (int i=0;i<N_X;i++)
    {
        phi_[i][N_Y-1]=0.0;
        phi_[i][0]=0.0;
        if ((i<20)&&(i>0)) phi_[i][N_Y-1]=-A;
        if ((i>=60)&&(i<N_X-1)) phi_[i][0]=A;
    }

    for (int nnn=0;nnn<13;nnn++)
    {
        for (int i=1;i<N_X-1;i++)
        {
            // phi_[i][j]=A*sin(i*6*M_PI/(N_X-1)+1800.0*alpha*t*dt);
            phi_[i][N_Y-1]=(phi_[i][N_Y-1]+phi_[i-1][N_Y-1]+phi_[i+1][N_Y-1])/3.0;//-A;//.0;//A*sin(i*6*M_PI/(N_X-1)+1800.0*alpha*t*dt);
            phi_[i][0]=(phi_[i-1][0]+phi_[i][0]+phi_[i+1][0])/3.0;//A*sin(i*6*M_PI/(N_X-1)+1800.0*alpha*t*dt);
        }
    }


    INPUT_PARAM par;
    par.a=((2.0)/(dx*dx)+2.0/(dy*dy));
    par.bp=-1.0/(dx*dx);
    par.bm=-1.0/(dx*dx);
    par.cp=-1.0/(dy*dy);
    par.cm=-1.0/(dy*dy);

    par.w_bc_type=2;
    par.e_bc_type=2;
    par.n_bc_type=3;
    par.s_bc_type=3;

    par.w_bc_val=0.0;
    par.e_bc_val=0.0;
    par.n_bc_val=0.0;
    par.s_bc_val=0.0;

    //jacobi(par,phi_,div_,20);//phi
    multigrid_N(par,phi_,RHS,8,3);
    multigrid_N(par,phi_,RHS,8,3);
    multigrid_N(par,phi_,RHS,8,3);
    multigrid_N(par,phi_,RHS,8,3);
    multigrid_N(par,phi_,RHS,8,3);

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


    for (int i=0; i<N_X; i++)
    {
        for (int j=0; j<N_Y; j++)
        {
            Ex[i][j]=ddx(phi_,i,j);
            Ey[i][j]=ddy(phi_,i,j);
        }

    }

    //    (n_1[i][j]-n_1_prev[i][j])/dt = D*((n_1[i+1][j]+n_1[i-1][j]-2.0*n_1[i][j])/(dx*dx) +(n_1[i][j+1]+n_1[i][j-1]-2.0*n_1[i][j])/(dy*dy))
    //                 +mu*(n_1[i][j]*rho_[i][j]+E x[i][j]*(n_1[i+1][j]-n_1[i-1][j])/(2.0*dx) + E y[i][j]*(n_1[i][j+1]-n_1[i][j-1])/(2.0*dy));


    /*
    for (int i=1; i<N_X-1; i++)
    {
        for (int j=1; j<N_Y-1; j++)
        {

            RHS_p[i][j]=rho_[i][j]/eps_0;
        }

    }


    a=((2.0*D)/(dx*dx)+(2.0*D)/(dy*dy)+ b*(d2phidx+d2phidy));
    b_p=(-1.0*D)/(dx*dx)- b*dphidx/(2.0*dx);
    b_m=(-1.0*D)/(dx*dx)+ b*dphidx/(2.0*dx);
    c_p=(-1.0*D)/(dy*dy)- b*dphidy/(2.0*dy);
    c_m=(-1.0*D)/(dy*dy)+ b*dphidy/(2.0*dy);


    rho_in[i][j]=rho_in[i][j]*0.9+0.1*(-(b_p*rho_in[i+1][j]+b_m*rho_in[i-1][j]+c_p*rho_in[i][j+1]+c_m*rho_in[i][j-1]))/a;
*/


    /*  multigrid_N(par,phi_,RHS,8,3);
    multigrid_N(par,phi_,RHS,8,3);
    multigrid_N(par,phi_,RHS,8,3);*/
    //  get_rhs_NAVIER();


    //solve_current(rho_,Jx_,Jy_,Ux_,Uy_,1000);
    // solve_current_1(rho_,Jx_,Jy_,Ux_,Uy_,150);
    // solve_current_phi(rho_,Jx_, Jy_,Ux_,Uy_,100);

}


void init()
{
    glClearColor (0.0, 0.0, 0.0, 0.0);
    glColor3f(1.0, 1.0, 1.0);
    glMatrixMode (GL_PROJECTION);
    glLoadIdentity ();
    glOrtho(-0.6, 0.6, -0.6,0.6, -10.0, 10.0);
    glMatrixMode (GL_MODELVIEW);

    sweep_init();

    fmm_init();

}

double delta_f(double E, double phi)
{
    double Emax=1.0;
    double E0=0.0;
    double E1=E-E0;

    if (E1<0) return -1; //just bounce back the electron

    double ks_d=1.0; //roughness in delta
    double ks_w=1.0; //roughness in w
    double deltam=1.0;
    double phi2pi=phi*phi/(2.0*M_PI);

    double w=(E1)/(Emax*(1.0+ks_w*phi2pi)-E0);


    double k= (w < 1) ? 0.56 : 0.25;


    double W= (w <= 3.6) ? pow(w*exp(1-w),k) : 1.125/pow(w,0.35);

    double delta=deltam*(1.0+ks_d*phi2pi)*W;

    return delta;


}


double delta2_f(double E, double phi)
{
    double Emax=1.0;
    double E0=0.0;
    double E1=E-E0;

    if (E1<0) return -1; //just bounce back the electron

    double ks_d=1.0; //roughness in delta
    double ks_w=1.0; //roughness in w
    double deltam=1.0;
    double phi2pi=phi*phi/(2.0*M_PI);

    double w=(E1)/(Emax*(1.0+ks_w*phi2pi)-E0);


    double k= (w < 1) ? 0.56 : 0.25;

    double alpha=std::min(w-3,1.0);
    alpha=std::max(alpha,0.0);

    double W=  pow(w*exp(1-w),k)*(1.0-alpha)+alpha* 1.125/pow(w,0.35);

    double delta=deltam*(1.0+ks_d*phi2pi)*W;

    return delta;


}


double Sample_energies_sec() //полная энергия должна сохраняться...
{

}

int main(int argc, char** argv)
{

    FILE* f=fopen("delta.txt","w");
    for (int i=0;i<200;i++)
    {
        fprintf(f,"%f %f %f\n",i*9.0/200,delta_f(i*9.0/200,0.0),delta2_f(i*9.0/200,0.0));
    }
    fclose(f);

    srand(time(NULL));
    glutInit(&argc,argv);
    glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB);
    glutInitWindowSize(W_WIDTH,W_HEIGHT);
    glutInitWindowPosition(0,0);
    glutCreateWindow("simple");
    glutDisplayFunc(display);
    glutMotionFunc(m_m);
    glutMouseFunc(m_d);
    glutKeyboardFunc(kb);
    init();
    glutMainLoop();
}
