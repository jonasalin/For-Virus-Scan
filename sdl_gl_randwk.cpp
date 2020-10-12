#include <windows.h>
#include <wingdi.h>
#include <SDL.h>
#include <GL\gl.h>
#include <GL\glu.h>

#include <iostream.h>
#include <fstream.h>
#include <stdlib.h>

#include "utils//math3d.cpp"
#define RADIUS 1.53
#define NUMBER_OF_STEPS 10000
#define NDEBUG
#include <assert.h>

#ifndef NDEBUG
#define TRACE(STRING) { cout << endl << STRING; }
#define TRACEV(STRING, value) { cout << endl << STRING << value; }
#else
#define TRACE(STRING) (void)0;
#define TRACEV(STRING, value) (void)0;
#endif

class Program {
public:
	Program(double ss, int ns); 
	~Program() { delete[] points, steps, max, pts; };
	void Init(void);
	void UpdateScreen(void);
	double GetSquaredDist(void) { return sum.lengthSq(); };
	double GetStepSize(void) { return stepsize; };
	Vect_3d& GetMaxVec(void) { return *max; };
	int GetNumSteps(void) { return numsteps; };
	Uint32 TimeCheck(void);
	void SetTime(Uint32 time) { next_time = time; };
	void IncrTime(Uint32 time) { next_time += time; };
	void handle_key_down( SDL_keysym* keysym );
	GLboolean should_rotate; 
private:
	Uint32 next_time;
	int numsteps;
	double stepsize;
	Vect_3d *points, *steps, sum, *max;	
	GLfloat *pts, rot_angle;
	Matrix Ry;
};

Program::Program(double ss, int ns)
: sum(0, 0, 0), Ry(3), rot_angle(0.0), stepsize(ss), numsteps(ns), should_rotate(GL_TRUE)
{
	pts = new GLfloat[numsteps*3];
	points = new Vect_3d[numsteps];
	steps = new Vect_3d[numsteps];
	max = &steps[0];
}

void Program::UpdateScreen(void)
{
	/* Clear the color and depth buffers. */
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	
	/* We don't want to modify the projection matrix. */
	glMatrixMode( GL_MODELVIEW );
	glLoadIdentity( );
		
	/* Move down the z-axis. */
	glTranslatef( 0.0, 0.0, -250.0 );

	/* Rotate. */
	
	glRotatef( rot_angle, 0.0, 1.0, 0.0 );
	

	if( should_rotate ) 
	{
		if( (rot_angle+=0.5) > 360.0f ) 
			rot_angle = 0.0f;
	}
	
	/* Send our coordinate data to the pipeline. */

	glBegin( GL_LINE_STRIP );
	glColor3f(1.0, 0.0, 0.0);
	for(int i = 0, j = 0; i < numsteps; i++, j += 3)
		glVertex3fv(&pts[j]);
	glEnd( );
	
	glBegin( GL_POINTS );
	glColor3f(0.0, 1.0, 0.0);
	glVertex3fv(&pts[0]);
	glColor3f(1.0, 1.0, 1.0);
	glVertex3fv(&pts[(numsteps-1)*3]);
	glEnd();

	glFlush();
	SDL_GL_SwapBuffers( );
}

void Program::Init(void)
{
	double r, angle;
	Matrix My(3), Mz(3);

//	using namespace rand_rng32;
//	seed(time(NULL));
	srand(time(NULL));
	
	for(unsigned long i = 0; i < numsteps; i++)
	{
		Vect_3d v(stepsize, 0, 0);
		/* Create rot matrices for rotating v random angles around y and z axes */
		rot_matr(My, y, rand() % 360);
		rot_matr(Mz, z, rand() % 360);
		/* Rotate v */
		v = My * v;
		v = Mz * v;
		/* Add to the sum rand-walk vector */
		sum += v;
		steps[i] = sum;	// Store each step
	
		if(steps[i].lengthSq() > max->lengthSq())
			max = &steps[i];
	}
	
	/* Calculate center of mass */
	Vect_3d cm(0, 0, 0);
	for (unsigned long i = 0; i < numsteps; i++)
		cm += steps[i];
	cm = cm / numsteps;
	
	/* Move origo to the center of mass and further translate the points back from viewpoint by tz */
	for(unsigned long i = 0; i < numsteps; i++)
		points[i] = steps[i] - cm;
	
	/* Store the coordinates in OpenGL-friendly array for later drawing on screen */
	for(unsigned long i=0, j=0; i < numsteps; i++, j+=3)
	{
		pts[j] = points[i](0);
		pts[j+1] = points[i](1);
		pts[j+2] = points[i](2);
	}
}

Uint32 Program::TimeCheck(void)
{
	Uint32 now = SDL_GetTicks();
	if(next_time <= now)
		return 0;
	else
		return next_time - now;
}


static void quit_sdl( int code )
{
	/*
	* Quit SDL so we can release the fullscreen
	* mode and restore the previous video settings,
	* etc.
	*/
	SDL_Quit( );
	
	/* Exit program. */
	exit( code );
}

void Program::handle_key_down( SDL_keysym* keysym )
{
	
	/*
	* We're only interested if 'Esc' has
	* been presssed.
	*
	* EXERCISE:
	* Handle the arrow keys and have that change the
	* viewing position/angle.
	*/
	switch( keysym->sym ) 
	{
		case SDLK_ESCAPE:
		quit_sdl( 0 );
		break;
		case SDLK_SPACE:
		should_rotate = !should_rotate;
		break;
		default:
		break;
	}	
}

static void process_events( Program &pgm )
{
	/* Our SDL event placeholder. */
	SDL_Event event;
	
	/* Grab all the events off the queue. */
	while( SDL_PollEvent( &event ) ) 
	{
		
		switch( event.type ) 
		{
			case SDL_KEYDOWN:
			/* Handle key presses. */
			pgm.handle_key_down( &event.key.keysym );
			break;
			case SDL_QUIT:
			/* Handle quit requests (like Ctrl-c). */
			quit_sdl( 0 );
			break;
		}
		
	}
	
}


static void setup_opengl( int width, int height )
{
	float ratio = (float) width / (float) height;
	
	/* Our shading model--Gouraud (smooth). */
	glShadeModel( GL_SMOOTH );
	
	/* Culling. */
//	glCullFace( GL_BACK );
//	glFrontFace( GL_CCW );
//	glEnable( GL_CULL_FACE );
	
	/* Set the clear color. */
	glClearColor( 0, 0, 0, 0 );
	
	/* Setup our viewport. */
	glViewport( 0, 0, width, height );
	
	/*
	* Change to the projection matrix and set
	* our viewing volume.
	*/
	glMatrixMode( GL_PROJECTION );
	glLoadIdentity( );
	/*
	* EXERCISE:
	* Replace this with a call to glFrustum.
	*/
	gluPerspective( 30.0, ratio, 1.0, 1024.0 );
	
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glEnable(GL_LINE_SMOOTH);
	glEnable(GL_POINT_SMOOTH);
	glLineWidth(1.0);
	glPointSize(5.0);
}

int main( int argc, char* argv[] )
{
	/* Information about the current video settings. */
	const SDL_VideoInfo* info = NULL;
	/* Dimensions of our window. */
	int width = 0;
	int height = 0;
	/* Color depth in bits of our window. */
	int bpp = 0;
	/* Flags we will pass into SDL_SetVideoMode. */
	int flags = 0;
/*
	ofstream outf("randw.txt");
	if(!outf)
	{
		cerr << "Could not create output file." << endl;
		exit(1);
	}
*/	
	/* First, initialize SDL's video subsystem. */
	if( SDL_Init( SDL_INIT_VIDEO ) < 0 ) 
	{
		/* Failed, exit. */
		cerr << "Video initialization failed: " << SDL_GetError( ) << endl;
		quit_sdl( 1 );
	}
	
	/* Let's get some video information. */
	info = SDL_GetVideoInfo( );
	
	if( !info ) 
	{
		/* This should probably never happen. */
		cerr << "Video query failed: " << SDL_GetError( ) << endl;
		quit_sdl( 1 );
	}
	
	width = 800;
	height = 600;
	bpp = info->vfmt->BitsPerPixel;
	
	SDL_GL_SetAttribute( SDL_GL_RED_SIZE, 5 );
	SDL_GL_SetAttribute( SDL_GL_GREEN_SIZE, 5 );
	SDL_GL_SetAttribute( SDL_GL_BLUE_SIZE, 5 );
	SDL_GL_SetAttribute( SDL_GL_DEPTH_SIZE, 16 );
	SDL_GL_SetAttribute( SDL_GL_DOUBLEBUFFER, 1 );
	
	flags = SDL_OPENGL; // | SDL_FULLSCREEN;
	
	if(SDL_SetVideoMode(width, height, bpp, flags) == 0)
	{
		/*
		* This could happen for a variety of reasons,
		* including DISPLAY not being set, the specified
		* resolution not being available, etc.
		*/
		cerr << "Video mode set failed: " << SDL_GetError( ) << endl;
		quit_sdl( 1 );
	}
	
	Program pgm(RADIUS, NUMBER_OF_STEPS); 
	
	pgm.Init();
/*
	outf << "Step size: " << pgm.GetStepSize() << endl;
	outf << "Number of steps: " << pgm.GetNumSteps() << endl;
	outf << "Final vector distance from origo: " << sqrt(pgm.GetSquaredDist()) << endl;
	outf << "Maximum distance from origo: " << pgm.GetMaxVec().length() 
		<< ", coordinates: (" << pgm.GetMaxVec()(0) << ", " << pgm.GetMaxVec()(1) << ", " << pgm.GetMaxVec()(2) << ")" << endl;
	outf.close();
*/	
	
	/*
	* At this point, we should have a properly setup
	* double-buffered window for use with OpenGL.
	*/
	setup_opengl( width, height );
	
	/*
	* Now we want to begin our normal app process--
	* an event loop with a lot of redrawing.
	*/
	
	SDL_Event event;

	Uint32 start = SDL_GetTicks();
	pgm.SetTime(start + 10);
	
	while( 1 ) 
	{
		/* Process incoming events. */
		process_events(pgm);
		/* Draw the screen. */
		pgm.UpdateScreen();
		SDL_Delay(pgm.TimeCheck());
		pgm.IncrTime(10);
	}
	
	return 0;
}
