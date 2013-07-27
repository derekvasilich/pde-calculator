/*
 *  main.c
 *  Assignment Two
 *
 *  Created by Derek Williams on 10-01-30.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 *	This code numerically solves hyperbolic PDEs of the form: 
 *  
 *		Dt[u(x, t)] + a(t, x) * Dx[u(x, t)] = F(x, t)
 *
 *		where Dt[] and Dx[] are the differential operators for t and x
 *
 *  The following finite difference schemes are implemented:
 *		
 *		Forward-Time Forward-Space
 *		Forward-Time Backward-Space
 *		Forward-Time Central-Space
 *		Lax Fredrichs
 *		Leapfrog
 *		Equilibrium
 *      Lax Wendroff
 *
 *  The solutions are animated in a window. Data is saved to text files when 's' is pressed.
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <OpenGL/OpenGL.h>
#include <GLUT/GLUT.h>

#include "main.h"

int activeScheme = DEFAULT_SCHEME;
int lastFrameTime = 0;
int blowup = 0;
double dx = 0.1;
double dt = 0.0;
int Xm = 0;
int Tn = 0;
double lambda = 0.8;
float time = 0;
double **data = 0;
double *maxv = 0;
double maxmaxv = 0;
double *error = 0;
int maxVmN = 0;
int pbc = 1;				// periodic boundary conditions
int pause = 0;				// start/pause animation
int right = 1;				// direction of wave propagation
int zeroRightBoundary = 0;

char *schemeNames[NUM_SCHEMES] = {
	FORWARD_SPACE_NAME,
	BACKWARD_SPACE_NAME,
	CENTRAL_SPACE_NAME,
	LAX_FREDRICHS_NAME,
	LEAPFROG_NAME,
	EQUILIBRIUM_NAME,
	LAX_WENDROFF_NAME
};

char *aFuncNames[2] = {
	"",
	"[1 + 1/4(3-x)(1+x)] "
};

int activeA = 0;


double FZero(double t, double x)
{
	return 0;
}

double FOne(double t, double x)
{
	return sin(M_PI*(x-t));
}

double FAOne(double t, double x)
{
	return t*sin(M_PI*(x-t));
}

double (*FActual)(double t, double x) = &FZero;
double (*FFunc)(double t, double x) = &FZero;

double BCOne(double t, double x)
{
	return -(1+t)*sin(M_PI*t);
}

double (*boundCondFunc)(double t, double x) = &FZero;

/* functions to specify a(x, t) parameter */

double aConst(double t, double x)
{
	return ((right) ? 1 : -1);
}

double aSpace(double t, double x)
{
	return 1 + 0.25 * (3 - x) * (1 + x);
}

double (*aFunc)(double t, double x) = &aConst;

/* functions to specify initial conditions */

double initOne (double t, double x)
{
	if ( fabs(x) <= 0.5)
		return pow(cos(M_PI * x), 2);
	return 0;
}

double initTwo (double t, double x)
{
	double xabs = fabs(x);
	if ( xabs <= 1.0 )
		return 1.0 - xabs;
	return 0;
}

double initThree (double t, double x)
{
	return sin(2*M_PI*x);
}

double initFour (double t, double x)
{
	return sin(M_PI*x);
}

double (*initFunc)(double t, double x) = &initOne;

double solution (double t, double x)
{
	double a = (aFunc)(t, x);
	return (initFunc)(0, x - a*t) + FActual(x, t);
}

int initialize (void) 
{
	int i;
	
	dt = dx*lambda;

	Xm = ceil((BOUNDS_RIGHT - BOUNDS_LEFT) / dx + 1);
	Tn = (TIME_END - TIME_START) / dt;
		
	if (data)
		free(data);
	
	data = malloc(Tn * sizeof(double *));
	if (data == NULL)
		perror("Error allocating memory.\n");

	maxv = malloc(Tn * sizeof(double));
	if (maxv == NULL)
		perror("Error allocating memory.\n");

	error = malloc(Tn * sizeof(double));
	if (error == NULL)
		perror("Error allocating memory.\n");
	
	for (i = 0; i < Tn; i++) {
		maxv[i] = -1e100;
		error[i] = 0;
		
		data[i] = malloc(Xm * sizeof (double));
		if (data[i] == NULL)
			perror("Error allocating memory.\n");
	}
}

double applyScheme(int scheme, int n, int m)
{
	double result = 0.0;	
	double lambda2 = lambda*lambda;		
	double t = n * dt;
	double x = BOUNDS_LEFT + m * dx;

	double a = (aFunc)(t, x);
	
	int mm = m - 1;
	int mp = m + 1;

	if (pbc)
	{
		mp = mp % Xm;
		if (mm < 0) 
			mm = Xm - 1;
	}
	else if (m == (Xm-1) && zeroRightBoundary) {
		return 0.0;
	}
	else if (m == (Xm-1) && activeScheme != BACKWARD_SPACE) {
		mp = m;
	}

	switch (scheme) {
		case FORWARD_SPACE:
			result += data[n-1][m] * ( 1 + a * lambda );
			result += data[n-1][mp] * ( - a * lambda );
			break;

		case BACKWARD_SPACE:
			result += data[n-1][mm] * ( a * lambda );
			result += data[n-1][m] * ( 1 - a * lambda );
			break;

		case CENTRAL_SPACE:
			result += data[n-1][mm] * ( a * 0.5 * lambda );
			result += data[n-1][m];
			result += data[n-1][mp] * ( - a * 0.5 * lambda );
			break;

		case LEAPFROG:
			result += data[n-2][m];
			result += data[n-1][mm] * ( a * lambda );
			result += data[n-1][mp] * ( - a * lambda );
			break;

		case LAX_FREDRICHS:
			result += data[n-1][mm] * ( 0.5 * (1 + a * lambda) );
			result += data[n-1][mp] * ( 0.5 * (1 - a * lambda) );
			break;			
					
		case EQUILIBRIUM:
			result += ( data[n-1][mm] + dt*(FFunc)(t-dt, x-dx) ) * ( 0.5 * a * ( lambda + a*lambda2) );
			result += ( data[n-1][m]  + dt*(FFunc)(t-dt, x) ) * ( 1 - a*lambda2 );
			if (m != (Xm-1))
				result += ( data[n-1][mp] + dt*(FFunc)(t-dt, x+dx) )* ( 0.5 * a * (-lambda + a*lambda2) );
			else {
				result += 2 * ( data[n-1][m] + dt*(FFunc)(t-dt, x) )* ( 0.5 * a * (-lambda + a*lambda2) );
				result -= 2 * ( data[n-1][mm] + dt*(FFunc)(t-dt, x-dx) )* ( 0.5 * a * (-lambda + a*lambda2) );
			}
			break;
			
		case LAX_WENDROFF:
			result += data[n-1][m] * ( 1 - a * a * lambda2 );
			result += data[n-1][mm] * ( 0.5 * a * lambda * ( 1 + a * lambda) );
			result += data[n-1][mp] * ( 0.5 * a * lambda * ( -1 + a * lambda) );
			break;						

		default:
			break;
	}
	
	if (abs(result) > 5 && !blowup)
		blowup = n;

	if (result > maxv[n])
		maxv[n] = result;
	
	if (result > maxmaxv)
		maxmaxv = result;
	
	if (!pbc && m != (Xm-1)) // exclude duplicate from error calculation
		error[n] += pow(solution(t, x) - result, 2);
	
	return result;
}

void save (void)
{
	char t[1000];
	int m, n;
	
	FILE *o;
	
	sprintf(t, "V-%s.dat", schemeNames[activeScheme]);
	o = fopen(t, "w");
	
	for (n = 0; n < Tn; n++) {		
		for (m = 1; m < Xm; m++)	
			fprintf(o, "%.4f ", data[n][m]);
		fprintf(o, "\n");
	}
	fflush(o);
	fclose(o);
	
	sprintf(t, "Max-%s.dat", schemeNames[activeScheme]);
	o = fopen(t, "w");

	for (n = 0; n < Tn; n++) {		
		fprintf (o, "%.4f ", maxv[n]);
	}
	fprintf(o, "\n");

	fflush(o);
	fclose(o);
	
	sprintf(t, "Error-%d-%s.dat", (int)(1/dx), schemeNames[activeScheme]);
	o = fopen(t, "w");
	
	for (n = 0; n < Tn; n++) {		
		fprintf (o, "%.4f ", sqrt(dx*error[n]));
	}
	fprintf(o, "\n");
	
	fclose(o);
}

void finiteDifference (void)
{
	int m, n, l, i;
	double v0m;
		
	for (i = 0; i < Tn; i++) {
		maxv[i] = - 1e100;
	}
	maxmaxv = -1e100;

	/* apply initial condition */
	for (m = 0; m < Xm; m++)
	{
		v0m = (initFunc)(0, BOUNDS_LEFT + m*dx);

		data[0][m] = v0m;
		
		if (v0m > maxv[0])
			maxv[0] = v0m;
		if (v0m > maxmaxv)
			maxmaxv = v0m;
	}
	
	l = 1;
	if (activeScheme == LEAPFROG)
	{
		l = 2;
		if (pbc)
			data[1][0] = applyScheme(CENTRAL_SPACE, 1, 0);
		else
			data[1][0] = (boundCondFunc)(dt, BOUNDS_LEFT);				// boundary condition

		for (m = 1; m < Xm; m++)
			data[1][m] = applyScheme(CENTRAL_SPACE, 1, m);
	}

	for (n = l; n < Tn; n++) {
		if (pbc)			
			data[n][0] = applyScheme(activeScheme, n, 0);
		else
			data[n][0] = (boundCondFunc)(n*dt, BOUNDS_LEFT);			// boundary condition
			
		for (m = 1; m < Xm; m++)	
			data[n][m] = applyScheme(activeScheme, n, m);
	} 
	
	blowup = 0;
	time = 0;
}

void display(void)
{
	int x, t, i;
	int linespacing = 5;
	char c[100], *pc;
	int width, wwidth;
	wwidth = (BOUNDS_RIGHT - BOUNDS_LEFT);
	
	if (lastFrameTime == 0)
    {
        lastFrameTime = glutGet(GLUT_ELAPSED_TIME);
    }
    
	int now = glutGet(GLUT_ELAPSED_TIME);
    int elapsedMilliseconds = now - lastFrameTime;
	float elapsedTime = elapsedMilliseconds / 1000.0f / dt;
    lastFrameTime = now;
	
	if (!pause)
		time += elapsedTime;		
	
	t = time;

	if (t >= Tn) 
	{
		t = Tn - 1;
		time = TIME_END/dt;
	}
	
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	/* draw title */
	
	sprintf(c, "Dt[u(x, t)] + %sDx[u(x, t)] = 0", aFuncNames[activeA]);
	pc = c;
	glColor3f(0.0f, 0.0f, 0.0f);
	width = glutBitmapLength(GLUT_BITMAP_HELVETICA_18, (unsigned char *)pc);
	glRasterPos2f(BOUNDS_LEFT + 0.5*wwidth*(1 - (float)width/glutGet(GLUT_WINDOW_WIDTH)), 2.7);
	pc = c;
	while (*pc != '\0')
		glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, *pc++);
	
	sprintf(c, "%s (Time %.2f)", schemeNames[activeScheme], TIME_START + time*dt);
	pc = c;
	glColor3f(0.0f, 0.0f, 0.0f);
	width = glutBitmapLength(GLUT_BITMAP_HELVETICA_18, (unsigned char *)pc);
	glRasterPos2f(BOUNDS_LEFT + 0.5*wwidth*(1 - (float)width/glutGet(GLUT_WINDOW_WIDTH)), 2.5);
	pc = c;
	while (*pc != '\0')
		glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, *pc++);
	
	if (blowup)
	{
		sprintf(c, "Blow-up (Time %.2f)", TIME_START + blowup*dt);
		pc = c;
		glColor3f(1.0f, 0.0f, 0.0f);
		width = glutBitmapLength(GLUT_BITMAP_HELVETICA_18, (unsigned char *)pc);
		glRasterPos2f(BOUNDS_LEFT + 0.5*wwidth*(1 - (float)width/glutGet(GLUT_WINDOW_WIDTH)), 2.3);
		pc = c;
		while (*pc != '\0')
			glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, *pc++);
	}

	/* draw grid */

	glLineWidth (1);

	glBegin(GL_LINES);
	for (x = 0; x <= Xm*linespacing; x++) {
		if (x % linespacing == 0)
			glColor3f(0.5f, 0.5f, 0.5f);
		else
			glColor3f(0.9f, 0.9f, 0.9f);
		glVertex2f(BOUNDS_LEFT + x/(double)linespacing, -1);
		glVertex2f(BOUNDS_LEFT + x/(double)linespacing, 2);
	}
	glEnd();	

	glBegin(GL_LINES);
	for (x = 0; x <= 3*linespacing; x++) {
		if (x % linespacing == 0)
			glColor3f(0.5f, 0.5f, 0.5f);
		else
			glColor3f(0.9f, 0.9f, 0.9f);
		glVertex2f(BOUNDS_LEFT, -1.0 + x/(double)linespacing);
		glVertex2f(BOUNDS_RIGHT, -1.0 + x/(double)linespacing);
	}
	glEnd();	

	for (x = 0; x <= (BOUNDS_RIGHT - BOUNDS_LEFT); x++) {
		sprintf(c, "%d", BOUNDS_LEFT + x);
		pc = c;
		glColor3f(0.0f, 0.0f, 0.0f);
		glRasterPos2f(BOUNDS_LEFT + x - 0.01, -1.1);
		while (*pc != '\0')
			glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10, *pc++);		
		pc = c;
		glRasterPos2f(BOUNDS_LEFT + x - 0.01, 2.05);
		while (*pc != '\0')
			glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10, *pc++);		
		
	}
	
	/* draw error max(v_m(n)) */

	if (maxVmN)
	{
		glColor3f(0.0f, 0.0f, 1.0f);
		glBegin(GL_LINE_STRIP);
		for (i = 0; i < Tn; i++) {
			glVertex2f(BOUNDS_LEFT + i * (BOUNDS_RIGHT - BOUNDS_LEFT)/(double)Tn, maxv[i]*3/maxmaxv - 1);
		}
		glEnd();	
	}
	else {
	
		/* draw actual solution */

		glPushAttrib(GL_LINE_BIT);
		glLineStipple (3, 0xAAAA);
		glColor3f(0.0f, 0.0f, 0.0f);
		glBegin(GL_LINE_STRIP);
		for (x = 0; x < Xm; x++) {
			glVertex2f(BOUNDS_LEFT + x*dx, solution((t % Tn) * dt, BOUNDS_LEFT + x*dx));
		}
		glEnd();	
		glPopAttrib();
	
		/* draw solution */
	
		glColor3f(0.0f, 0.0f, 0.0f);
		glBegin(GL_LINE_STRIP);
		for (x = 0; x < Xm; x++) {
			glVertex2f(BOUNDS_LEFT + x*dx, data[t % Tn][x]);
		}
		glEnd();		
	}
	
	glutSwapBuffers();

}

void reshape(int width, int height)
{
    glViewport(0, 0, width, height);
    
	glEnable(GL_LINE_SMOOTH);
	glEnable(GL_LINE_STIPPLE);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(BOUNDS_LEFT - 0.1, BOUNDS_RIGHT + 0.1, -1.2, 3);
    glMatrixMode(GL_MODELVIEW);
}

void idle(void)
{
    glutPostRedisplay();
}

void Menu (int value) {
	int dirty = 0;
	
	switch (value) {
		case FORWARD_SPACE:
		case BACKWARD_SPACE:
		case CENTRAL_SPACE:
		case LAX_FREDRICHS:
		case LEAPFROG:
		case EQUILIBRIUM:
		case LAX_WENDROFF:
			activeScheme = value;
			break;

		case H_10:
		case H_20:
		case H_40:
		case H_80:
			dirty = (dx!=1/(double)value); 
			dx = 1/(double)value;
			break;
			
		case LAMBDA_8:
			dirty = (lambda != 0.8);
			lambda = 0.8;
			break;

		case LAMBDA_16:
			dirty = (lambda != 1.6);
			lambda = 1.6;
			break;
						
		case A_CONST:
			aFunc = &aConst;
			activeA = 0;
			break;		

		case A_SPACE:
			aFunc = &aSpace;
			activeA = 1;
			break;		

		case INIT_ONE:
			initFunc = &initOne;
			break;
			
		case INIT_TWO:
			initFunc = &initTwo;
			break;
			
		case INIT_THREE:
			initFunc = &initThree;
			break;		
			
		case INIT_FOUR:
			initFunc = &initFour;
			break;		

		case F_ZERO:
			FFunc = &FZero;
			FActual = &FZero;
			break;
			
		case F_ONE:
			FFunc = &FOne;
			FActual = &FAOne;
			break;
			
		case MAX_VM_N:
			maxVmN = !maxVmN;
			
		default:
			return;
	}

	if (dirty) {
		initialize();
	}
	
	finiteDifference();
}

void buildMenu (void)
{
	int i, schemeMenuId, hMenuId, lambdaMenuId, aMenuId, initMenuId, fMenuId, menuId;
	
	schemeMenuId = glutCreateMenu(Menu);
	
	for (i = 0; i < NUM_SCHEMES; i++)
		glutAddMenuEntry(schemeNames[i], i);

		
	hMenuId = glutCreateMenu(Menu);
	
	glutAddMenuEntry(H_10_NAME, H_10);
	glutAddMenuEntry(H_20_NAME, H_20);
	glutAddMenuEntry(H_40_NAME, H_40);
	glutAddMenuEntry(H_80_NAME, H_80);

	lambdaMenuId = glutCreateMenu(Menu);

	glutAddMenuEntry(LAMBDA_8_NAME, LAMBDA_8);
	glutAddMenuEntry(LAMBDA_16_NAME, LAMBDA_16);

	aMenuId = glutCreateMenu(Menu);
	
	glutAddMenuEntry(A_CONST_NAME, A_CONST);
	glutAddMenuEntry(A_SPACE_NAME, A_SPACE);

	initMenuId = glutCreateMenu(Menu);

	glutAddMenuEntry(INIT_ONE_NAME, INIT_ONE);
	glutAddMenuEntry(INIT_TWO_NAME, INIT_TWO);
	glutAddMenuEntry(INIT_THREE_NAME, INIT_THREE);
	glutAddMenuEntry(INIT_FOUR_NAME, INIT_FOUR);

	fMenuId = glutCreateMenu(Menu);
	
	glutAddMenuEntry(F_ZERO_NAME, F_ZERO);
	glutAddMenuEntry(F_ONE_NAME, F_ONE);
	
	menuId = glutCreateMenu(Menu);
	
	glutAddSubMenu("Scheme", schemeMenuId);
	glutAddSubMenu("Grid Spacing (h)", hMenuId);
	glutAddSubMenu("Lambda (k/h)", lambdaMenuId);
	glutAddSubMenu("a(t, x)", aMenuId);
	glutAddSubMenu("Initial Condition (u0(x))v", initMenuId);
	glutAddSubMenu("F(t, x)", fMenuId);	
	glutAddMenuEntry("max(v_m(n))", MAX_VM_N);
	
	glutAttachMenu(GLUT_LEFT_BUTTON);
	
}

void Mouse (int button, int state, int x, int y)
{
	
}

void Keyboard (unsigned char key, int x, int y)
{
	switch (key) {
		case ' ':
			pause = !pause;
			break;

		case 'd':
			right = !right;
			finiteDifference();
			break;
			
		case 's':
			save();
			break;
			
		case 'p':
			pbc = !pbc;
			finiteDifference();
			break;

		case 'r':
			zeroRightBoundary = !zeroRightBoundary;
			finiteDifference();
			break;			
			
		case 'x':
			boundCondFunc = &BCOne;
			finiteDifference();
			break;
			
		default:
			break;
	}
}

int main(int argc, char** argv)
{
	int t, x;

	initialize();
	finiteDifference();

    glutInit(&argc, argv);
    
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
    glutInitWindowSize(640, 640);
    
    glutCreateWindow("Numeric Hyperbolic PDE Solver");
    
	glClearColor(1.0f, 1.0f, 1.0f, 0.0f);
	
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutIdleFunc(idle);
	glutMouseFunc(Mouse);
	glutKeyboardFunc(Keyboard);
	
	buildMenu();
	
    glutMainLoop();
	
    return EXIT_SUCCESS;
}