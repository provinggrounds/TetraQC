/*
 *	quasi.c			Eric Weeks, 6-15-95
 *
 * formatted for tab stop = 5 (in vi: type  ":set ts=5 <return>")

Program:  quasi.c

link to this program:
    http://www.physics.emory.edu/~weeks/software/quasi.html
email: weeks@physics.emory.edu

This program is public domain, but please leave my name, email,
and web link attached.



Started:  2-8-94 by Eric Weeks, for problem set #3 of P392K
Designed to output lattice points for a quasi-crystal, following method
outlined in homework assignment (creating five sets of parallel lines
offset by 72 degrees, and assigning lattice points based on the empty
spaces between lines.)

V01: 2-08-94:	started from shell.c
	:	added in dot product routine to determine k values
V02: 6-15-95:  added in postscript code from fieldln.c - v13
V03: 6-21-95:  added in polygon fill; -M option
V05: 6-25-95:  tried to change for 7-fold symmetry
V06: 9-05-95:	tried to redo order of loops
     8-29-96:  added in -nz options; added maxmax variable
V07: 9-14-96:  add in -M 456 options
V08: 9-14-96:  add in -M 78 options
V09: 9-14-96:  add in -N option, change -M option
     2-04-98:  slight change to postscript code


 *
 * FLOW OF PROGRAM:
 * main : processes cmd line options, calls postscript stuff
 *	printheader  : prints CGLE postscript header
 *	quasi	   : does quasi-crystal
 *		plot    : subroutine to do postscript lines
 *   printtrailer : final postscript commands needed to draw page
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define PI 3.14159265358979323846264338328
#define EE 2.71828182845904523536
/* following two parameters are PostScript linewidths (in cm) */
#define LINEWIDTH 0.015
#define MAX 100				/* Increase this to make more tiles */

/* --- function declarations --- */
void quasi();
void plot();
double getdx();
void printheader();
void printtrailer();
void segment();

/* --- variables used only by 'plot' subroutine --- */
double window;
double oldx,oldy;
double color;
int fillon;
double gray;

/* --- variables used by whole program --- */
int midon[2];
double scale;
char *prgname;
float offsetx,offsety;
int oldflag,rotate;
double xcen,ycen;			/* center of picture */
int symmetry;
int zfill;
int maxmax;

/*  M   M   AAA   IIIII  N   N
 *  MM MM  A   A    I    NN  N
 *  M M M  AAAAA    I    N N N   ************************************
 *  M   M  A   A    I    N  NN
 *  M   M  A   A  IIIII  N   N
 */
main(argc,argv)
int argc;
char **argv;
{
	int c;
	double magnifya;
	extern int optind;
	extern char *optarg;

	/* --- set up defaults --- */
	window = 20.0;
	scale = 15.0;				/* 15 cm wide boxes				*/
	magnifya = 1.0;			/* no magnification				*/
	rotate = 0;				/* don't rotate picture			*/
	fillon = 0;				/* don't do polygon fill			*/
	midon[0] = 0;				/* don't connect midpoints		*/
	midon[1] = 0;				/* don't connect midpoints		*/
	symmetry = 5;				/* five-fold symmetry			*/
	zfill = 0;				/* don't use zfill option		*/
	maxmax = 30;
	while ((c = getopt(argc, argv, "hfs:m:FM:S:zn:N:")) != EOF)
		switch (c)  {
		case 'h':  fprintf(stderr,"Usage: %s [options]\n",argv[0]);
fprintf(stderr," -h   : this help message\n");
fprintf(stderr," -s # : size of output (width) in cm [%.1f]\n",scale);
fprintf(stderr," -f   : fill polygons\n");
fprintf(stderr," -z   : fill color is according to polygon type\n");
fprintf(stderr," -F   : flip picture 90 degrees\n");
fprintf(stderr," -m # : magnification factor \n");
fprintf(stderr," -M # : midpoint type for skinny diamonds\n");
fprintf(stderr," -N # : midpoint type for fat diamonds\n");
fprintf(stderr,"          1 = acute angle sides joined\n");
fprintf(stderr,"          2 = obtuse angle sides joined\n");
fprintf(stderr,"          3 = opposite sides joined to make cross\n");
fprintf(stderr,"          4 = all sides joined to make rectangle\n");
fprintf(stderr,"          5 = randomly choose 1 or 2\n");
fprintf(stderr,"          6 = randomly choose 1, 2, or 4\n");
fprintf(stderr," -S # : degrees of symmetry [%d]\n",symmetry);
fprintf(stderr," -n # : number of lines to use\n",maxmax);
				exit(1);
				break;
		case 's': scale = atof(optarg);
				break;
		case 'm': magnifya = atof(optarg);
				/* fprintf(stderr,"%f\n",magnifya);	HELP */
				break;
		case 'F': rotate = 1;
				break;
		case 'f': fillon = 1;
				break;
		case 'M': midon[0] = atoi(optarg);
				break;
		case 'N': midon[1] = atoi(optarg);
				break;
		case 'S': symmetry = atoi(optarg);
				break;
		case 'z': zfill = 1;
				fillon = 1;
				break;
		case 'n': maxmax = atoi(optarg);
				break;
		}

	xcen = 0.0; ycen = 0.0;
	window = window/magnifya;
	/* --- done with initialization stuff --- */

	printheader();

	offsetx = 2.0;					/* lower left corner of picture */
	offsety = 2.0;
	quasi();

	printtrailer();
	exit(0);
} /* END OF MAIN /* END OF MAIN /* END OF MAIN /* END OF MAIN /* END OF MAIN */






/*******************************************************************


 *******************************************************************/


void quasi()
{
	int t,r,index[50],i,j,m,n,flag;
	double vx[50],vy[50],mm[50],b[50][MAX],phi,x0,y0,x1,y1,dx;
	/* v's are vectors, m's are slopes, b's are intercepts */
	double midx1,midx2,midx3,midx4,midy1,midy2,midy3,midy4;
	double dx1,dx2,dy1,dy2,dist1,dist2;
	int themin,themax,nnn,mmm,flag2;
	double minmin,rad1,rad2,rad;
	int halfmax;
	int midsix = 0;
	int type,segtype;

	halfmax = maxmax/2;


	for (t=0;t<symmetry;t++)  {
		phi = (t*2.0)/(1.0*symmetry)*PI;
		vx[t] = cos(phi);
		vy[t] = sin(phi);
		mm[t] = vy[t]/vx[t];
		for (r=0;r<maxmax;r++)  {
			y1 = vy[t]*(t*0.1132) - vx[t]*(r-halfmax);  /* offset */
			x1 = vx[t]*(t*0.2137) + vy[t]*(r-halfmax);
			b[t][r] = y1 - mm[t]*x1;		/* intercept */
		}
	}

	/* t is 1st direction, r is 2nd.  look for intersection between pairs
	 * of lines in these two directions. (will be x0,y0) */

	color = 0.2;
	themax = (maxmax-1);
	themin = themax/2;
	for (minmin=0.0;minmin<=(double)(themax);minmin+=0.4)  {
		rad1 = minmin*minmin;
		rad2 = (minmin+0.4)*(minmin+0.4);
		for (n=1;n<themax;n++)  {
			for (m=1;m<themax;m++)  {
				rad = (double)((n-themin)*(n-themin)+(m-themin)*(m-themin));
				/* rad = (double)(n*n+m*m); */
				if ((rad>=rad1)&&(rad<rad2))  {

			for (t=0;t<(symmetry-1);t++)  {
				for (r=t+1;r<symmetry;r++)  {
					x0 = (b[t][n] - b[r][m])/(mm[r]-mm[t]);
					y0 = mm[t]*x0 + b[t][n];
					flag = 0;
					for (i=0;i<symmetry;i++)  {
						if ((i!=t) && (i!=r))  {
							dx = -x0*vy[i]+(y0-b[i][0])*vx[i];
							index[i] = -dx;
							if ((index[i]>(maxmax-3))||(index[i]<1))
								flag=1;
						}
					}
					if (flag==0)  {
						index[t] = n-1;    index[r] = m-1;
						x0 = 0.0;          y0 = 0.0;
						for (i=0;i<symmetry;i++)  {
							x0 += vx[i]*index[i];
							y0 += vy[i]*index[i];
						}
						if (midon[0]>0)  gray = 0.8;	/* faint lines */
						/* color of tile unless zfill==1 */
						color += .05;
						if (color>1.0)  color = 0.2;
if (zfill==1)  {
	color = 0.0;
	for (i=0;i<symmetry;i++)  {
		color += index[i];
	}
	while (color>((symmetry-1.0)/2.0))  
		color -= ((symmetry-1.0)/2.0);
	color = color/((symmetry-1.0)/2.0)*.8+.1;
	color += fabs(vx[t]*vx[r]+vy[t]*vy[r]);	/* dot product */
	if (color>1.0) color-=1.0;
}
						plot(x0,y0,0);
						x0 += vx[t];  y0 += vy[t];
						plot(x0,y0,1);
						x0 += vx[r];  y0 += vy[r];
						plot(x0,y0,1);
						x0 -= vx[t];  y0 -= vy[t];
						plot(x0,y0,1);
						x0 -= vx[r];  y0 -= vy[r];
						plot(x0,y0,2);
						if (midon[0]>0)  {
							midx1 = x0 + vx[t]*0.5;
							midy1 = y0 + vy[t]*0.5;
							midx2 = x0 + vx[t] + vx[r]*0.5;
							midy2 = y0 + vy[t] + vy[r]*0.5;
							midx3 = x0 + vx[r] + vx[t]*0.5;
							midy3 = y0 + vy[r] + vy[t]*0.5;
							midx4 = x0 + vx[r]*0.5;
							midy4 = y0 + vy[r]*0.5;
							dx1 = midx1-midx2;  dy1 = midy1-midy2;
							dist1 = dx1*dx1+dy1*dy1;
							dx2 = midx2-midx3;  dy2 = midy2-midy3;
							dist2 = dx2*dx2+dy2*dy2;
							gray = 0.0;		/* dark lines */
							if (dist1*dist2<.1)  
								type = 0;
							else
								type = 1;
							segtype = midon[type];
							if ((segtype==1)||(segtype==2))  {
								if (dist1>dist2)  segtype = 3-segtype;
							} else if (segtype==5)  {
								midsix = 1 - midsix;
								segtype = midsix + 1;
							} else if (segtype==6)  {
								midsix++;
								if (midsix>2)  midsix = 0;
								segtype = midsix + 1;
							}

							if (segtype==3)  {
								/* X's */
								segment(midx1,midy1,midx3,midy3);
								segment(midx2,midy2,midx4,midy4);
							} else if (segtype==1)  {
								segment(midx1,midy1,midx2,midy2);
								segment(midx3,midy3,midx4,midy4);
							} else if (segtype==2)  {
								segment(midx1,midy1,midx4,midy4);
								segment(midx2,midy2,midx3,midy3);
							} else if (segtype==4)  {
								/* boxes */
								segment(midx1,midy1,midx2,midy2);
								segment(midx3,midy3,midx4,midy4);
								segment(midx1,midy1,midx4,midy4);
								segment(midx2,midy2,midx3,midy3);
							}
						}
					}
					}
					}
				}
			}
		}
	}
}




void segment(double x1, double y1, double x2, double y2)
{
	plot (x1,y1,0);
	plot (x2,y2,2);
}







/*  PPPP   L       OOO   TTTTT
 *  P   P  L      O   O    T
 *  PPPP   L      O   O    T
 *  P      L      O   O    T
 *  P      LLLLL   OOO     T
 */

void plot(double x,double y,int plotflag)  {
	double dx,dy,swap;
	double cmx,cmy;		/* x,y in centimeters */
	/* flag variable:  0 = start line; 1 = lineto; 2 = endpoint */

	dx = getdx(x,xcen);
	dy = getdx(y,ycen);
	if (rotate==1)  { swap = dx; dx = dy; dy = swap; }

	if ((dx<1.3)&&(dy<1.0)&&(dx>0)&&(dy>0))  {		/* in window */
		cmx = dx * scale  +  offsetx;
		cmy = dy * scale  +  offsety;
		if (plotflag<1)  {
			printf("%.2f %.2f m\n",cmx,cmy);
		} else {
			if (oldflag==1)  printf("%.2f %.2f m\n",cmx,cmy);
			printf("%.2f %.2f l\n",cmx,cmy);
			if (plotflag==2)  {
				printf("cp\n");
				if (fillon==1)  {
					printf("g\n");
					printf("%.1f sg\n",color);
					printf("fill\n");
					printf("h\n");
				}
				if (midon>0)  printf("%.1f sg\n",gray);
				printf("a\n");
			}
		}
		oldflag = 0;
	} else {
		oldflag = 1;
	}

}


















/* Subroutine for plot */
double getdx(double x,double center)
{
	double dx;

	dx = (x - center) / window;
	dx = 0.5 * (dx + 1.0);
	/*  dx  : 0 = left/bottom, +1 = right/top */
	return dx;
}






























/* VARIOUS POSTSCRIPT STUFF */




void printheader()  {
/* ================================================================ */
/*  PostScript Header (taken from CGLE output)
 * ================================================================ */

printf("%%!PS-Adobe-1.0 \n");
printf("%%%%BoundingBox: -1 -1 766.354 567.929 \n");
printf("%%%%EndComments \n");
printf("%%%%EndProlog \n");
printf("gsave \n");
printf(" \n");
printf("/f {findfont exch scalefont setfont} bind def \n");
printf("/s {show} bind def \n");
printf("/ps {true charpath} bind def \n");
printf("/l {lineto} bind def \n");
printf("/m {newpath moveto} bind def \n");
printf("/sg {setgray} bind def\n");
printf("/a {stroke} bind def\n");
printf("/cp {closepath} bind def\n");
printf("/g {gsave} bind def\n");
printf("/h {grestore} bind def\n");
printf("matrix currentmatrix /originmat exch def \n");
printf("/umatrix {originmat matrix concatmatrix setmatrix} def \n");
printf(" \n");
 printf("%% Flipping coord system \n");
printf("[8.35928e-09 28.3465 -28.3465 8.35928e-09 609.449 28.6299] umatrix \n");
printf("[] 0 setdash \n");
printf("0 0 0 setrgbcolor \n");
printf("0 0 m \n");
printf("%f setlinewidth \n",LINEWIDTH);
/* END OF HEADER =============================================== */

}











void printtrailer()  {

printf("showpage grestore \n");
printf("%%%%Trailer\n");

}





