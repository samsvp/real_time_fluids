#include <stdlib.h>
#include <stdio.h>
#include <GL/glut.h>
#include <chrono>

#include "include/solver.hpp"


// window
static int win_id;
static int win_x = 512;
static int win_y = 512;

// mouse
static std::vector<int> mouse_down = {0, 0, 0};
static int last_mx = 0;
static int last_my = 0;

auto start = std::chrono::system_clock::now();

Fluid* fluid;


void render_density()
{
    int i, j;
	float x, y, h, d00, d01, d10, d11;

    float N = fluid->N; 
	h = 1.0f/N;

	glBegin(GL_QUADS);

    for(i=0 ; i<=N ; i++) {
        x =(i-0.5f)*h;
        for(j=0 ; j<=N ; j++) {
            y =(j-0.5f)*h;

            d00 = fluid->density[IX(i,j,N)];
            d01 = fluid->density[IX(i,j+1,N)];
            d10 = fluid->density[IX(i+1,j,N)];
            d11 = fluid->density[IX(i+1,j+1,N)];

            glColor3f(d00, 5 * d00, 10 * d00); glVertex2f(x, y);
            glColor3f(d10, 10 * d10, d10); glVertex2f(x+h, y);
            glColor3f(d11, d11, 25 * d11); glVertex2f(x+h, y+h);
            glColor3f(d01, 35 * d01, d01); glVertex2f(x, y+h);
        }
    }

	glEnd();
}

/*
  ----------------------------------------------------------------------
   GLUT callback routines
  ----------------------------------------------------------------------
*/

static void key_func(unsigned char key, int x, int y)
{
    switch(key)
	{

		case 'q':
		case 'Q':
        case 27: // Escape key
			glutDestroyWindow(win_id);
			break;

		case 'c':
		case 'C':
			fluid->clear();
			break;
	}
}


static void mouse_func(int button, int state, int x, int y)
{
	mouse_down[button] = state == GLUT_DOWN;

    if (state == GLUT_DOWN && button == GLUT_LEFT_BUTTON)
    {
        fluid->add_density(x, fluid->N-y, 100.0f);
    }
}


static void motion_func(int x, int y)
{
    if (mouse_down[0])
    {
        float delta_x = x - last_mx;
        float delta_y = y - last_my;
        fluid->add_velocity(x, fluid->N-y, 5 * delta_x, -5 * delta_y);

        fluid->add_density(x, fluid->N-y, 10.0f);
    }

    last_mx = x;
    last_my = y;
}


static void reshape_func(int width, int height)
{
	glutSetWindow(win_id);
	glutReshapeWindow(width, height);

	win_x = width;
	win_y = height;
}


static void pre_display()
{
	glViewport(0, 0, win_x, win_y);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(0.0, 1.0, 0.0, 1.0);
	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT);
}


static void post_display()
{
	glutSwapBuffers();
}


static void idle_func()
{
	glutSetWindow(win_id);
	glutPostRedisplay();
}

static void display_func()
{
	auto now = std::chrono::system_clock::now();
    std::chrono::duration<double> diff = now - start;
	double dt = diff.count();
	printf("%g FPS\n", 1/dt);
	start = now;

	pre_display();

    fluid->step(dt);
    render_density();

	post_display();
}


static void open_glut_window()
{
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);

	glutInitWindowPosition(0, 0);
	glutInitWindowSize(win_x, win_y);
	win_id = glutCreateWindow("Alias | wavefront");

	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT);
	glutSwapBuffers();
	glClear(GL_COLOR_BUFFER_BIT);
	glutSwapBuffers();

	pre_display();

	glutKeyboardFunc(key_func);
	glutMouseFunc(mouse_func);
	glutMotionFunc(motion_func);
	glutReshapeFunc(reshape_func);
	glutIdleFunc(idle_func);
	glutDisplayFunc(display_func);
}


int main(int argc, char ** argv)
{
    glutInit(&argc, argv);

    int N = win_x + 1;
    float dt = 1.0/60.;
    float diff = 0.0f;
    float visc = 0.00005f;
    fluid = new Fluid(N, dt, diff, visc);


    fprintf(stderr, "Using : N=%d dt=%g diff=%g visc=%g\n",
        N, dt, diff, visc);

    printf("\n\nHow to use this demo:\n\n");
	printf("\t Add densities with the right mouse button\n");
	printf("\t Add velocities with the left mouse button and dragging the mouse\n");
	printf("\t Toggle density/velocity display with the 'v' key\n");
	printf("\t Clear the simulation by pressing the 'c' key\n");
	printf("\t Quit by pressing the 'q' key\n");

	float dvel = 0;

	open_glut_window();

	glutMainLoop();

    delete fluid;

    return 0;
}