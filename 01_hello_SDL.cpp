/*This source code copyrighted by Lazy Foo' Productions 2004-2024
and may not be redistributed without written permission.*/

//Using SDL and standard IO
#include <SDL.h>
#include <stdio.h>
#include <cstdlib>
#include <time.h>
#include <SDL_ttf.h>
#include <string>
#undef main
#define _GBORDER 1
#define _GFILLINTERIOR 0

//Screen dimension constants
const int SCREEN_WIDTH = 320;
const int SCREEN_HEIGHT = 240;

//The window we'll be rendering to
SDL_Window* window = NULL;

//The surface contained by the window
SDL_Surface* screenSurface = NULL;


SDL_Renderer* renderer = NULL;

SDL_Texture* texture, * text;
TTF_Font* font;


int color = 0;
int global_x, global_y;
int text_x, text_y;
int isMap = 0;

void _setcolor(int new_color);
void _moveto(int x1, int y1);
int _setpixel(int x, int y);
int _lineto(int x2, int y2);
int _rectangle(int order, int x1, int y1, int x2, int y2);


// D E F I N E S /////////////////////////////////////////////////////////////

// #define DEBUG 1

#define OVERBOARD          48 // the absolute closest a player can get to a wall

#define INTERSECTION_FOUND 1

// constants used to represent angles

#define ANGLE_0     0
#define ANGLE_1     5
#define ANGLE_2     10
#define ANGLE_4     20
#define ANGLE_5     25
#define ANGLE_6     30
#define ANGLE_15    80
#define ANGLE_30    160
#define ANGLE_45    240
#define ANGLE_60    320
#define ANGLE_90    480
#define ANGLE_135   720
#define ANGLE_180   960
#define ANGLE_225   1200
#define ANGLE_270   1440
#define ANGLE_315   1680
#define ANGLE_360   1920

#define WORLD_ROWS    16        // number of rows in the game world
#define WORLD_COLUMNS 16        // number of columns in the game world
#define CELL_X_SIZE   64        // size of a cell in the gamw world
#define CELL_Y_SIZE   64

// size of overall game world

#define WORLD_X_SIZE  (WORLD_COLUMNS * CELL_X_SIZE)
#define WORLD_Y_SIZE  (WORLD_ROWS    * CELL_Y_SIZE)

// G L O B A L S /////////////////////////////////////////////////////////////

//unsigned int * clock = (unsigned int *)0x0000046C; // pointer to internal
// 18.2 clicks/sec


// world map of nxn cells, each cell is 64x64 pixels

char world[WORLD_ROWS][WORLD_COLUMNS] = { 
	{1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1},
	{1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1},
	{1,0,0,1,1,1,1,1,0,0,0,0,0,0,0,1},
	{1,0,0,1,0,0,0,1,0,1,0,1,0,1,0,1},
	{1,0,0,1,0,0,0,1,0,0,0,0,0,0,0,1},
	{1,0,0,1,0,0,1,1,0,0,0,0,0,0,0,1},
	{1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1},
	{1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1},
	{1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1},
	{1,0,0,1,1,0,0,1,1,1,1,1,1,0,0,1},
	{1,0,0,1,0,0,0,0,0,0,0,0,1,0,0,1},
	{1,0,0,1,1,1,0,0,0,0,0,0,1,0,0,1},
	{1,0,0,1,0,0,0,0,0,0,0,0,1,0,0,1},
	{1,0,0,1,1,1,1,1,1,1,1,0,1,0,0,1},
	{1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1},
	{1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1}
};       // pointer to matrix of cells that make up
// world

float tan_table[1921];
float inv_tan_table[1921];

float y_step[1921];
float x_step[1921]; 

float cos_table[1921]; 

float inv_cos_table[1921]; 
float inv_sin_table[1921];
float cos_local[361];
float sin_local[361];

// F U N C T I O N S /////////////////////////////////////////////////////////

float fabs(float f)
{
	return f >= 0.0 ? f : -f;
}
//Cordic in 32 bit signed fixed point math
//Function is valid for arguments in range -pi/2 -- pi/2
//for values pi/2--pi: value = half_pi-(theta-half_pi) and similarly for values -pi---pi/2
//
// 1.0 = 1073741824
// 1/k = 0.6072529350088812561694
// pi = 3.1415926536897932384626

// https://www.dcs.gla.ac.uk/~jhw/cordic/
// https://www.dcs.gla.ac.uk/~jhw/cordic/cordic-32bit.h

#define M_PI 3.1415926536
#define cordic_1K 0x26DD3B6A
#define half_pi 0x6487ED51
#define MUL 1073741824.000000f
#define CORDIC_NTAB 32

int cordic_ctab[] = {
	0x3243F6A8,
	0x1DAC6705,
	0x0FADBAFC,
	0x07F56EA6,
	0x03FEAB76,
	0x01FFD55B,
	0x00FFFAAA,
	0x007FFF55,
	0x003FFFEA,
	0x001FFFFD,
	0x000FFFFF,
	0x0007FFFF,
	0x0003FFFF,
	0x0001FFFF,
	0x0000FFFF,
	0x00007FFF,
	0x00003FFF,
	0x00001FFF,
	0x00000FFF,
	0x000007FF,
	0x000003FF,
	0x000001FF,
	0x000000FF,
	0x0000007F,
	0x0000003F,
	0x0000001F,
	0x0000000F,
	0x00000008,
	0x00000004,
	0x00000002,
	0x00000001,
	0x00000000,
};

void cordic(int theta, int* s, int* c, int n)
{
	int k, d, tx, ty, tz;
	int x = cordic_1K, y = 0, z = theta;
	n = (n > CORDIC_NTAB) ? CORDIC_NTAB : n;
	for (k = 0; k < n; ++k)
	{
		d = z >> 31;
		// get sign. for other architectures, you might want to use the more portable version
		// d = z>=0 ? 0 : -1;
		tx = x - (((y >> k) ^ d) - d);
		ty = y + (((x >> k) ^ d) - d);
		tz = z - ((cordic_ctab[k] ^ d) - d);
		x = tx;
		y = ty;
		z = tz;
	}
	*c = x;
	*s = y;
}



void cordic_grad(float grad, int* sin, int* cos)
{
	float p = grad * 3.14159265f / 180;
	cordic((p * MUL), sin, cos, 32);
}

//rad_angle = (3.272e-4) + ang * 2 * 3.141592654 / ANGLE_360;

float sinc_rad(float r)
{
	int s, c;
	float rad = r;
	if (rad >= 2 * M_PI) {
		int p = rad / (2 * M_PI);
		rad = rad - (p * 2 * M_PI);
	}

	float c_rad = rad;

	if (rad > M_PI/2 && rad <= M_PI)
	{
		c_rad = M_PI - rad;
	}

	if (rad > M_PI && rad <= 3*M_PI/2)
	{
		c_rad = rad - M_PI;
	}
	if (rad > 3*M_PI/2 && rad <= 2*M_PI)
	{
		c_rad = 2*M_PI - rad;
	}

	cordic(c_rad*MUL, &s, &c,32);
	
	float sn = s / MUL;

	if (rad > M_PI && rad <= 2*M_PI)
	{
		sn = sn * (-1.0f);
	}
	return sn;
}


float cosc_rad(float r)
{
	int s, c;
	float rad = r;
	if (rad >= 2 * M_PI) {
		int p = rad / (2 * M_PI);
		rad = rad - (p * 2 * M_PI);
	}
	float c_rad = rad;

	if (rad > M_PI/2 && rad <= M_PI)
	{
		c_rad = M_PI - rad;
	}

	if (rad > M_PI && rad <= 3*M_PI/2)
	{
		c_rad = rad - M_PI;
	}
	if (rad > 3*M_PI/2 && rad <= 2*M_PI)
	{
		c_rad = 2*M_PI - rad;
	}

	cordic(c_rad*MUL, &s, &c,32);
	float cs = c / MUL;

	if ((rad > M_PI/2 && rad <= M_PI) || (rad > M_PI && rad <= 3 * M_PI / 2) )
	{
		cs = cs * (-1.0f);
	}
	return cs;
}




float sinc(float g)
{
	int s, c;
	float grad = g;
	float c_grad = grad;
	if (grad > 90 && grad <= 180)
	{
		c_grad = 180 - grad;
	}

	if (grad > 180 && grad <= 270)
	{
		c_grad = grad - 180;
	}
	if (grad > 270 && grad <= 380)
	{
		c_grad = 360 - grad;
	}

	cordic_grad(c_grad, &s, &c);
	float sn = s / MUL;

	if (grad > 180 && grad <= 360)
	{
		sn = sn * (-1.0f);
	}
	return sn;
}

float cosc(float g)
{
	int s, c;
	float grad = g;
	float c_grad = grad;

	if (grad > 90 && grad <= 180)
	{
		c_grad = 180 - grad;
	}

	if (grad > 180 && grad <= 270)
	{
		c_grad = grad - 180;
	}
	if (grad > 270 && grad <= 380)
	{
		c_grad = 360 - grad;
	}

	cordic_grad(c_grad, &s, &c);
	float cs = c / MUL;

	if (grad > 180 && grad <= 360)
	{
		cs = cs * (-1.0f);
	}
	return cs;
}

float tanc_rad(float r)
{
	int s, c;
	float rad = r ;
	if (rad > 2 * M_PI) {
		int p = rad / (2 * M_PI);
		rad = rad - (p * 2 * M_PI);
	}
	float c_rad = rad;
	if (rad > M_PI/2 && rad <= 2*M_PI)
	{
		c_rad = M_PI - rad;
	}

	if (rad > M_PI && rad <= 3*M_PI/2)
	{
		c_rad = rad - M_PI;
	}
	if (rad > 3*M_PI/2 )
	{
		c_rad = 2*M_PI - rad;
	}

	cordic(c_rad*MUL, &s, &c,32);
	float cs = c;
	float sn = s;

	if (rad > M_PI/2 && rad <= M_PI)
	{
		sn = sn * (-1.0f);
	}

	if (rad > 3*M_PI/2  && rad <= 2*M_PI)
	{
		sn = sn * (-1.0f);
	}


	return cs == 0 ? 0 : (float) sn / cs;
}

void Build_Tables(void)
{

	int ang;
	float rad_angle;



	// create tables, sit back for a sec!

	for (ang = ANGLE_0; ang <= ANGLE_360; ang++)
	{

		rad_angle = (3.272e-4) + ang * 2 * 3.141592654 / ANGLE_360;

		tan_table[ang] = tanc_rad(rad_angle);
		inv_tan_table[ang] = 1 / tan_table[ang];


		// tangent has the incorrect signs in all quadrants except 1, so
		// manually fix the signs of each quadrant since the tangent is
		// equivalent to the slope of a line and if the tangent is wrong
		// then the ray that is case will be wrong

		if (ang >= ANGLE_0 && ang < ANGLE_180)
		{
			y_step[ang] = fabs(tan_table[ang] * CELL_Y_SIZE);
		}
		else
			y_step[ang] = -fabs(tan_table[ang] * CELL_Y_SIZE);

		if (ang >= ANGLE_90 && ang < ANGLE_270)
		{
			x_step[ang] = -fabs(inv_tan_table[ang] * CELL_X_SIZE);
		}
		else
		{
			x_step[ang] = fabs(inv_tan_table[ang] * CELL_X_SIZE);
		}

		// create the sin and cosine tables to copute distances

		inv_cos_table[ang] = 1 / cosc_rad(rad_angle);

		inv_sin_table[ang] = 1 / sinc_rad(rad_angle);
		//printf("%d %g %g %g \n", ang, rad_angle, sin(rad_angle), sinc_rad(rad_angle));


	} // end for ang

// create view filter table.  There is a cosine wave modulated on top of
// the view distance as a side effect of casting from a fixed point.
// to cancell this effect out, we multiple by the inverse of the cosine
// and the result is the proper scale.  Without this we would see a
// fishbowl effect, which might be desired in some cases?

	for (ang = -ANGLE_30; ang <= ANGLE_30; ang++)
	{

		rad_angle = (3.272e-4) + ang * 2 * 3.141592654 / ANGLE_360;

		cos_table[ang + ANGLE_30] = 1 / cosc_rad(rad_angle);

	} // end for

} // end Build_Tables

/////////////////////////////////////////////////////////////////////////////
 
void Timer(int clicks)
{
	// this function uses the internal time keeper timer i.e. the one that goes
	// at 18.2 clicks/sec to to a time delay.  You can find a 322 bit value of
	// this timer at 0000:046Ch

	unsigned int now;

	// get current time

	now = clock();

	// wait till time has gone past current time plus the amount we eanted to
	// wait.  Note each click is approx. 55 milliseconds.

	while (abs((long)(clock() - now)) < clicks) {}

} // end Timer

/////////////////////////////////////////////////////////////////////////////

void sline(long x1, long y1, long x2, long y2, int color)
{

	// used a a diagnostic function to draw a scaled line

	x1 = x1 / 4;
	y1 = 256 - (y1 / 4);

	x2 = x2 / 4;
	y2 = 256 - (y2 / 4);

	_setcolor(color);
	_moveto((int)x1, (int)y1);
	_lineto((int)x2, (int)y2);

} // end sline

/////////////////////////////////////////////////////////////////////////////

void splot(long x, long y, int color)
{
	// used as a diagnostic function to draw a scaled point

	x = x / 4;
	y = 256 - (y / 4);

	_setcolor(color);

	_setpixel((int)x, (int)y);
	_setpixel((int)x + 1, (int)y);
	_setpixel((int)x, (int)y + 1);
	_setpixel((int)x + 1, (int)y + 1);

} // end splot

/////////////////////////////////////////////////////////////////////////////

void Draw_2D_Map(void)
{
	// draw 2-D map of world

	int row, column, block, t, done = 0;

	for (row = 0; row < WORLD_ROWS; row++)
	{
		for (column = 0; column < WORLD_COLUMNS; column++)
		{

			block = world[row][column];

			// test if there is a solid block there

			if (block == 0)
			{

				_setcolor(253);
				_rectangle(_GBORDER, column * CELL_X_SIZE / 4, row * CELL_Y_SIZE / 4,
					column * CELL_X_SIZE / 4 + CELL_X_SIZE / 4 - 1, row * CELL_Y_SIZE / 4 + CELL_Y_SIZE / 4 - 1);

			}
			else
			{

				_setcolor(152);
				_rectangle(_GFILLINTERIOR, column * CELL_X_SIZE / 4, row * CELL_Y_SIZE / 4,
					column * CELL_X_SIZE / 4 + CELL_X_SIZE / 4 - 1, row * CELL_Y_SIZE / 4 + CELL_Y_SIZE / 4 - 1);

			}

		} // end for column

	} // end for row

} // end Draw_2D_Map

/////////////////////////////////////////////////////////////////////////////

void Ray_Caster(long x, long y, long view_angle)
{
	// This function casts out 320 rays from the viewer and builds up the video
	// display based on the intersections with the walls. The 320 rays are
	// cast in such a way that they all fit into a 60 degree field of view
	// a ray is cast and then the distance to the first horizontal and vertical
	// edge that has a cell in it is recorded.  The intersection that has the
	// closer distance to the user is the one that is used to draw the bitmap.
	// the distance is used to compute the height of the "sliver" of texture
	// or line that will be drawn on the screen

	// note: this function uses floating point (slow), no optimizations (slower)
	// and finally it makes calls to Microsofts Graphics libraries (slowest!)
	// however, writing it in this manner makes it many orders of magnitude
	// easier to understand.

	int rcolor = 55;

	long xray = 0,        // tracks the progress of a ray looking for Y interesctions
		yray = 0,        // tracks the progress of a ray looking for X interesctions
		next_y_cell,   // used to figure out the quadrant of the ray
		next_x_cell,
		cell_x,        // the current cell that the ray is in
		cell_y,
		x_bound,       // the next vertical and horizontal intersection point
		y_bound,
		xb_save,       // storage to record intersections cell boundaries
		yb_save,
		x_delta,       // the amount needed to move to get to the next cell
		y_delta,       // position
		ray,           // the current ray being cast 0-320
		casting = 2,     // tracks the progress of the X and Y component of the ray
		x_hit_type,    // records the block that was intersected, used to figure
		y_hit_type,    // out which texture to use

		top,           // used to compute the top and bottom of the sliver that
		bottom;        // is drawn symetrically around the bisecting plane of the
	// screens vertical extents


	float xi,           // used to track the x and y intersections
		yi,
		xi_save,      // used to save exact x and y intersection points
		yi_save,
		dist_x,       // the distance of the x and y ray intersections from
		dist_y,       // the viewpoint
		scale;        // the final scale to draw the "sliver" in

	// S E C T I O N  1 /////////////////////////////////////////////////////////v

	// initialization

	// compute starting angle from player.  Field of view is 60 degrees, so
	// subtract half of that current view angle

	if ((view_angle -= ANGLE_30) < 0)
	{
		// wrap angle around
		view_angle = ANGLE_360 + view_angle;
	} // end if

	rcolor = 1 + rand() % 14;

	// loop through all 320 rays

	// section 2

	
	

	for (ray = 0; ray < SCREEN_WIDTH; ray++)
	{
		
		// S E C T I O N  2 /////////////////////////////////////////////////////////

			// compute first x intersection

			// need to know which half plane we are casting from relative to Y axis
		float rad = (3.272e-4) + view_angle * 2 * 3.141592654 / ANGLE_360;
		
		if (view_angle >= ANGLE_0 && view_angle < ANGLE_180)
		{
			
			// compute first horizontal line that could be intersected with ray
			// note: it will be above player

			y_bound = CELL_Y_SIZE + CELL_Y_SIZE * (y / CELL_Y_SIZE);

			// compute delta to get to next horizontal line

			y_delta = CELL_Y_SIZE;

			// based on first possible horizontal intersection line, compute X
			// intercept, so that casting can begin

			//xi = inv_tan_table[view_angle] * (y_bound - y) + x;
			

			//printf("%d %g %g %g \n", view_angle, rad, tanc_rad(rad), tan_table[view_angle]);
			//printf("%d %g %g %g \n\n", view_angle, rad, 1.0f / tanc_rad(rad), inv_tan_table[view_angle]);

			//xi = inv_tan_table[view_angle] * (y_bound - y) + x;
			xi = inv_tan_table[view_angle] * (y_bound - y) + x;


			// set cell delta

			next_y_cell = 0;

		} // end if upper half plane
		else
		{

			// compute first horizontal line that could be intersected with ray
			// note: it will be below player

			y_bound = CELL_Y_SIZE * (y / CELL_Y_SIZE);

			// compute delta to get to next horizontal line

			y_delta = -CELL_Y_SIZE;

			// based on first possible horizontal intersection line, compute X
			// intercept, so that casting can begin
			

			xi = inv_tan_table[view_angle] * (y_bound - y) + x;

			// set cell delta

			next_y_cell = -1;

		} // end else lower half plane


 // S E C T I O N  3 /////////////////////////////////////////////////////////

	 // compute first y intersection

	 // need to know which half plane we are casting from relative to X axis

		if (view_angle < ANGLE_90 || view_angle >= ANGLE_270)
		{

			// compute first vertical line that could be intersected with ray
			// note: it will be to the right of player

			x_bound = CELL_X_SIZE + CELL_X_SIZE * (x / CELL_X_SIZE);

			// compute delta to get to next vertical line

			x_delta = CELL_X_SIZE;

			// based on first possible vertical intersection line, compute Y
			// intercept, so that casting can begin

			yi = tan_table[view_angle]* (x_bound - x) + y;
			//printf("%d, ", x / CELL_X_SIZE);
			// set cell delta

			next_x_cell = 0;

		} // end if right half plane
		else
		{

			// compute first vertical line that could be intersected with ray
			// note: it will be to the left of player

			x_bound = CELL_X_SIZE * (x / CELL_X_SIZE);

			// compute delta to get to next vertical line

			x_delta = -CELL_X_SIZE;

		//	// based on first possible vertical intersection line, compute Y
		//	// intercept, so that casting can begin

			yi = tan_table[view_angle] * (x_bound - x) + y;

		//	// set cell delta

			next_x_cell = -1;

		} // end else right half plane

 // begin cast

		casting = 2;                // two rays to cast simultaneously
		xray = yray = 0;                // reset intersection flags


		// S E C T I O N  4 /////////////////////////////////////////////////////////

		while (casting)
		{

			// continue casting each ray in parallel

			if (xray != INTERSECTION_FOUND)
			{

				// test for asymtotic ray

				// if (view_angle==ANGLE_90 || view_angle==ANGLE_270)

				if (fabs(y_step[view_angle]) == 0)
				{
					xray = INTERSECTION_FOUND;
					casting--;
					dist_x = 1e+8;

				} // end if asymtotic ray

			 // compute current map position to inspect

				cell_x = ((x_bound + next_x_cell) / CELL_X_SIZE);
				cell_y = (long)(yi / CELL_Y_SIZE);

				// test if there is a block where the current x ray is intersecting

				//if ((x_hit_type = world[(WORLD_ROWS - 1) - cell_y][cell_x]) != 0)
				if (cell_x < 0) cell_x = 0;
				if (cell_x > WORLD_COLUMNS) cell_x = WORLD_COLUMNS - 1; 
				if (cell_y < 0) cell_y = 0; 
				if (cell_y > WORLD_ROWS) cell_y = WORLD_ROWS - 1;
				if ((x_hit_type = world[(WORLD_ROWS - 1) - cell_y][cell_x]) != 0)
				{
					// compute distance

					dist_x = (yi - y) * inv_sin_table[view_angle];
					yi_save = yi;
					xb_save = x_bound;

					// terminate X casting

					xray = INTERSECTION_FOUND;
					casting--;

				} // end if a hit
				else
				{
					// compute next Y intercept

					yi += y_step[view_angle];

				} // end else

			} // end if x ray has intersected

// S E C T I O N  5 /////////////////////////////////////////////////////////

			if (yray != INTERSECTION_FOUND)
			{

				// test for asymtotic ray

				// if (view_angle==ANGLE_0 || view_angle==ANGLE_180)

				if (fabs(x_step[view_angle]) == 0)
				{
					yray = INTERSECTION_FOUND;
					casting--;
					dist_y = 1e+8;

				} // end if asymtotic ray

			 // compute current map position to inspect

				cell_x = (long)(xi / CELL_X_SIZE);
				cell_y = ((y_bound + next_y_cell) / CELL_Y_SIZE);
								if (cell_x < 0) cell_x = 0;
				if (cell_x > WORLD_COLUMNS) cell_x = WORLD_COLUMNS - 1; 
				if (cell_y < 0) cell_y = 0; 
				if (cell_y > WORLD_ROWS) cell_y = WORLD_ROWS - 1;
				// test if there is a block where the current y ray is intersecting

				if ((y_hit_type = world[(WORLD_ROWS - 1) - cell_y][cell_x]) != 0)
				{
					// compute distance

					dist_y = (xi - x) * inv_cos_table[view_angle];
					xi_save = xi;
					yb_save = y_bound;

					// terminate Y casting

					yray = INTERSECTION_FOUND;
					casting--;

				} // end if a hit
				else
				{
					// compute next X intercept

					xi += x_step[view_angle];

				} // end else

			} // end if y ray has intersected

		 // move to next possible intersection points


			x_bound += x_delta;
			y_bound += y_delta;


			// _settextposition(38,40);
			// printf("x_bound = %ld, y_bound = %ld    ",x_bound,y_bound);

		} // end while not done


// S E C T I O N  6 /////////////////////////////////////////////////////////

	// at this point, we know that the ray has succesfully hit both a
	// vertical wall and a horizontal wall, so we need to see which one
	// was closer and then render it

	// note: latter we will replace the crude monochrome line with a sliver
	// of texture, but this is good enough for now

		if (dist_x < dist_y)
		{


			// there was a vertical wall closer than the horizontal

			// compute actual scale and multiply by view filter so that spherical
			// distortion is cancelled

			scale = cos_table[ray] * 15000 / (1e-10 + dist_x);

			// compute top and bottom and do a very crude clip

			if ((top = 100 - scale / 2) < 1)
				top = 1;

			if ((bottom = top + scale) > 199)
				bottom = 199;

			// draw wall sliver and place some dividers up

			if (((long)yi_save) % CELL_Y_SIZE <= 1)
				_setcolor(222);
			else
				_setcolor(221);

			
			if(isMap == 1) {
				sline(x, y, (long)xb_save, (long)yi_save, rcolor);
			}
			else {
				_moveto((int)(SCREEN_WIDTH - ray - 1), (int)top);
				_lineto((int)(SCREEN_WIDTH - ray - 1), (int)bottom);
			}
		}
		else // must of hit a horizontal wall first
		{


			// compute actual scale and multiply by view filter so that spherical
			// distortion is cancelled

			scale = cos_table[ray] * 15000 / (1e-10 + dist_y);

			// compute top and bottom and do a very crude clip

			if ((top = 100 - scale / 2) < 1)
				top = 1;

			if ((bottom = top + scale) > 199)
				bottom = 199;

			// draw wall sliver and place some dividers up

			if (((long)xi_save) % CELL_X_SIZE <= 1)
				_setcolor(215);
			else
				_setcolor(210);
			if (isMap == 1) {
				sline(x, y, (long)xi_save, (long)yb_save, rcolor);
			}
			else {
				_moveto((int)(SCREEN_WIDTH - ray - 1), (int)top);
				_lineto((int)(SCREEN_WIDTH - ray - 1), (int)bottom);
			}

		} // end else

 // S E C T I O N  7 /////////////////////////////////////////////////////////

	 // cast next ray

		if (++view_angle >= ANGLE_360)
		{
			// reset angle back to zero

			view_angle = 0;

		} // end if

	} // end for ray

} // end Ray_Caster



int _lineto(int x2, int y2) {
	int res = SDL_RenderDrawLine(renderer, global_x, global_y, x2, y2);
	if (res < 0) {
		SDL_ShowSimpleMessageBox(SDL_MESSAGEBOX_ERROR, "Couldn't create renderer!", SDL_GetError(), window);
	}

	_moveto(x2, y2);
	return 0;
}



void _setcolor(int new_color) {
	int sol = new_color % 3;

	int	b = (new_color >> 16) & 0xFF;
	int r = (new_color >> 8) & 0xFF;
	int g = new_color & 0xFF;

	if (sol == 1) {
		g = (new_color >> 16) & 0xFF;
		b = (new_color >> 8) & 0xFF;
		r = new_color & 0xFF;
	}
	if (sol == 2) {
		r = (new_color >> 16) & 0xFF;
		g = (new_color >> 8) & 0xFF;
		b = new_color & 0xFF;
	}

	SDL_SetRenderDrawColor(renderer, r, g, b, 255);
}

int _setpixel(int x, int y) {
	SDL_RenderDrawPoint(renderer, x, y);
	return 0;
}

void _moveto(int x1, int y1) {
	global_x = x1;
	global_y = y1;
}

int _rectangle(int border, int x1, int y1, int x2, int y2) {

	if (border) {
		SDL_Rect rect;
		rect.x = x1;
		rect.y = y1;
		rect.h = y2 - y1;
		rect.w = x2 - x1;
		return SDL_RenderDrawRect(renderer, &rect);
	}
	else
	{
		SDL_Rect rect;
		rect.x = x1;
		rect.y = y1;
		rect.h = y2 - y1;
		rect.w = x2 - x1;
		return SDL_RenderFillRect(renderer, &rect);
	}
}

int  main(void)
{

	int row, column, block, t, done = 0;

	long x, y, view_angle, x_cell, y_cell, x_sub_cell, y_sub_cell;

	float dx, dy;



	srand(13);
	printf("Build tables");
	Build_Tables();

	 //Initialize SDL
	if (SDL_Init(SDL_INIT_VIDEO) < 0)
	{
		printf("SDL could not initialize! SDL_Error: %s\n", SDL_GetError());
		return 0;
	}
	else
	{


		//Create window
		window = SDL_CreateWindow("SDL Tutorial", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, SCREEN_WIDTH, SCREEN_HEIGHT, SDL_WINDOW_SHOWN);
		if (window == NULL)
		{
			printf("Window could not be created! SDL_Error: %s\n", SDL_GetError());
			return 0;
		}
		else
		{

			SDL_CreateWindowAndRenderer(SCREEN_WIDTH, SCREEN_HEIGHT, SDL_WINDOW_RESIZABLE, &window, &renderer);

			if (!renderer) {
				printf("Failed to create window and renderer:   %s\n", SDL_GetError());
				return 0;
			}
			// Initialize SDL_ttf
			if (TTF_Init() < 0) {
				printf("Error intializing SDL_ttf: %s\n", SDL_GetError());
				return false;
			}
			// Load font
			font = TTF_OpenFont("C:\\_Development\\sdl\\Arial.ttf", 12);
			if (!font) {
				printf("TTF  %s\n", SDL_GetError());
				return false;
			}


			// draw top view of world

			if (isMap == 1) {
				Draw_2D_Map();
			}
			else
			{
				_setcolor(115);
				_rectangle(_GBORDER, 0, 0, 319, 199);
			}

			SDL_RenderPresent(renderer);

			// draw window around view port

		

			x = 8 * 64 + 25;
			y = 3 * 64 + 25;


			view_angle = ANGLE_60;

			// render initial view
			if (isMap == 0) {
				_setcolor(118);
				_rectangle(_GFILLINTERIOR, 0, 100, 319, 199);
			}
			
			Ray_Caster(x, y, view_angle);
			SDL_RenderPresent(renderer);

			bool KEYS[322];  // 322 is the number of SDLK_DOWN events

			for (int i = 0; i < 322; i++) { // init them all to false
				KEYS[i] = false;
			}



			// wait for user to press q to quit

			while (!done)
			{

				// has keyboard been hit?


				SDL_Event event; bool quit = false;
				while (quit == false)
				{
					while (SDL_PollEvent(&event))
					{
						switch (event.type) {
							// exit if the window is closed
						case SDL_QUIT:
							quit = true; // set game state to done,(do what you want here)
							break;
							// check for keypresses
						case SDL_KEYDOWN:
							KEYS[event.key.keysym.sym] = true;
							break;
						case SDL_KEYUP:
							KEYS[event.key.keysym.sym] = false;
							break;
						default:
							break;
						}


						if (event.type == SDL_KEYDOWN) {
							// reset deltas

							dx = dy = 0;

							// clear viewport
							if (isMap == 1) {
								Draw_2D_Map();
							}
							else {
								_setcolor(0);
								_rectangle(_GFILLINTERIOR, 0, 0, 320, 200);
								_setcolor(118);
								_rectangle(_GFILLINTERIOR, 0, 100, 320, 200);
							}
							// what is user doing

							if (KEYS[SDLK_d]) {

								if ((view_angle -= ANGLE_6) < ANGLE_0)
									view_angle = ANGLE_360;

							}

							if (KEYS[SDLK_a])
							{
								if ((view_angle += ANGLE_6) >= ANGLE_360)
									view_angle = ANGLE_0;

							}


							if (KEYS[SDLK_w])
							{

								// move player along view vector foward


								dx = cosc_rad(6.28 * view_angle / ANGLE_360) * 10;
								dy = sinc_rad(6.28 * view_angle / ANGLE_360) * 10;

							}

							if (KEYS[SDLK_s])
							{
								// move player along view vector backward

								dx = -cosc_rad(6.28 * view_angle / ANGLE_360) * 10;
								dy = -sinc_rad(6.28 * view_angle / ANGLE_360) * 10;

								// test if player is bumping into a wall

							}

							if (KEYS[SDLK_m]) {
								_setcolor(0);
								_rectangle(_GFILLINTERIOR, 0, 0, 320, 240);

								if (isMap == 0) {
									isMap = 1;
									Draw_2D_Map();
								}
								else
								{
									isMap = 0;
									_setcolor(118);
									_rectangle(_GFILLINTERIOR, 0, 100, 320, 200);
								}
							}


							if (KEYS[SDLK_q]) { done = 1; quit = true; }


							// move player

							x += dx;
							y += dy;

							// test if user has bumped into a wall i.e. test if there
							// is a cell within the direction of motion, if so back up !

							// compute cell position

							x_cell = x / CELL_X_SIZE;
							y_cell = y / CELL_Y_SIZE;

							// compute position relative to cell

							x_sub_cell = x % CELL_X_SIZE;
							y_sub_cell = y % CELL_Y_SIZE;


							// resolve motion into it's x and y components

							if (dx > 0)
							{
								// moving right

								if ((world[(WORLD_ROWS - 1) - y_cell][x_cell + 1] != 0) &&
									(x_sub_cell > (CELL_X_SIZE - OVERBOARD)))
								{
									// back player up amount he steped over the line

									x -= (x_sub_cell - (CELL_X_SIZE - OVERBOARD));

								} // end if need to back up

							}
							else
							{
								// moving left

								if ((world[(WORLD_ROWS - 1) - y_cell][x_cell - 1] != 0) &&
									(x_sub_cell < (OVERBOARD)))
								{
									// back player up amount he steped over the line

									x += (OVERBOARD - x_sub_cell);

								} // end if need to back up

							} // end else

							if (dy > 0)
							{
								// moving up

								if ((world[(WORLD_ROWS - 1) - (y_cell + 1)][x_cell] != 0) &&
									(y_sub_cell > (CELL_Y_SIZE - OVERBOARD)))
								{
									// back player up amount he steped over the line

									y -= (y_sub_cell - (CELL_Y_SIZE - OVERBOARD));

								} // end if need to back up
							}
							else
							{
								// moving down

								if ((world[(WORLD_ROWS - 1) - (y_cell - 1)][x_cell] != 0) &&
									(y_sub_cell < (OVERBOARD)))
								{
									// back player up amount he steped over the line

									y += (OVERBOARD - y_sub_cell);

								} // end if need to back up

							} // end else

						 // render the view

							Ray_Caster(x, y, view_angle);

							// display status

							SDL_RenderPresent(renderer);
						}

					}
				}




			} //end while

				//Destroy window
			SDL_DestroyWindow(window);

			//Quit SDL subsystems
			SDL_Quit();

			return 0;
		}
	}
} // end main

