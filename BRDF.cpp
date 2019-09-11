/***author: Samia Kabir*****/

#include <GL/glut.h>
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <math.h>

using namespace std;

#define ImageW 600
#define ImageH 400
#define pi 3.1416

struct Coord3D { double x, y, z; };
void GetColor(Coord3D view,Coord3D normal, Coord3D light, int SphNum, float& R, float& G, float& B);
double Diffuse(double n);
double Specular(Coord3D view,Coord3D normal, Coord3D light, int SphNum,double n);
float framebuffer[ImageH][ImageW][3];
double M[6]={0.3,0.6,1.0,0.5,0.9,0.4}; //surface roughness of 6 spheres where 1st three=Plastic,2nd two=Fe and 3rd=diamond
int wavelength[6]={436,455,510,587,635,695}; //6 samples
double energy[6]={1.106,1.2655,1.3497,1.3708,1.3065,1.1538};  //spectral irradiance for sunlight
double eta[6][6]={                                      //refractive index for 3 materials,6 spheres
                  {1.502,1.502,1.502,2.336,2.336,2.444},
                  {1.500,1.500,1.500,2.606,2.606,2.438},
                  {1.495,1.495,1.495,2.807,2.807,2.428},
                  {1.490,1.490,1.490,2.931,2.931,2.417},
                  {1.488,1.488,1.488,2.897,2.897,2.410},
                  {1.486,1.486,1.486,2.872,2.872,2.406}
                 };

double match_func[6][3]={
                          {0.32271,0.04337,1.78750},
                          {0.28266,0.07238,1.66444},    //color matching function
                          {0.01557,0.52049,0.11768},
                          {1.07842,0.84150,0.00011},
                          {0.57554,0.24169,0.00000},
                          {0.01575,0.00602,0.00000}
                        };

double sRGB[3][3]= {                                    //conversion matrix from XYZ to sRGB
                    {3.2406,-1.5372,-0.4986},
                    {-0.9689,1.8758,0.0415},
                    {0.0557,-0.2040,1.0570}
                   };

double RGB[3][1];        //R,G,B values
double col_val[3][1];    //X,Y,Z values

Coord3D Light;
Coord3D SphCent[6];
double SphRad[6];

double F,F0,F90=1;     //Fresnel terms
double n;              //refractive index
float theta, alpha;   //incident angle and angle between Normal and H
float c,g;            //terms to compute F
double D,Geo;
double Rs,Rd;          //Rs=specular component Rd=diffuse component
Coord3D H;             //Bisector of light and view vector
double m;
double s=0.7,d=0.3;    //fraction of light that is specular and diffuse
double BR;             //bidirectional reflectance
int L=1;



// Normalizes the vector passed in
void normalize(double& x, double& y, double& z) {
	float temp = sqrt(x*x + y*y + z*z);
	if (temp > 0.0) {
		x /= temp;
		y /= temp;
		z /= temp;
	}
	else {
		x = 0.0;
		y = 0.0;
		z = 0.0;
	}
}

// Returns dot product of two vectors
float dot(double x1, double y1, double z1, double x2, double y2, double z2) {
	return (x1*x2 + y1*y2 + z1*z2);
}

// Returns angle between two vectors
float angle(double x1, double y1, double z1, double x2, double y2, double z2) {
	normalize(x1, y1, z1);
	normalize(x2, y2, z2);
	return  acos(dot(x1, y1, z1, x2, y2, z2));
}
//return the minimum of three numbers
double min(double x,double y,double z)
{
    if(x<y)
    {
        if(z<x)
          return z;
        else
          return x;
    }
    else
    {
        if(z<y)
          return z;
        else
          return y;

    }
}


// Get Color for a point on the surface
void GetColor(Coord3D view,   // Normalized Vector pointing FROM eye TO surface
	Coord3D normal, // Normalized Vector giving surface normal
	Coord3D light,  // Normalized Vector pointing FROM surface TO light
	int SphNum,     // Sphere Number (0-5)
	float& R,       // Return these values for surface color.
	float& G,
	float& B)
{

   int count;

   col_val[0][0]= 0;
   col_val[1][0]= 0;
   col_val[2][0]= 0;

  for(count=0;count<6;count++)        //iterate for all 6 samples
  {
    L=count;
    n=eta[L][SphNum];
    m=M[SphNum];


  /************calculating diffuse component***********************/
    F0= ((n-1)/(n+1))*((n-1)/(n+1)); // fresnel term at 0 incident
    Rd=(F0/pi);                       //bidirectional reflectance of material at normal
  //  cout<< Rd << endl;

  /************calculateing spceular component *********************/

  //calculate H

   H.x=(view.x+light.x)/2;
   H.y=(view.y+light.y)/2;
   H.z=(view.z+light.z)/2;

   normalize(H.x,H.y,H.z);

   //calculate alpha and theta

   alpha=angle(H.x,H.y,H.z,normal.x,normal.y,normal.z);

   theta=angle(H.x,H.y,H.z,view.x,view.y,view.z);

   //Calculate F at theta
   c=cos(theta);
   float temp= (n*n) +(c*c) -1;
   g=sqrt(temp);

   float t1= ((g-c)*(g-c))/((g+c)*(g+c));
   float t2= (c*(g+c))-1;
   float t3= (c*(g-c))+1;
   float t4= ((t2*t2)/(t3*t3))+1;
   F= t1*t4*0.5;

   //calculate D for rougness m

    float tmp1= (tan(alpha)/m);
    float param=tmp1*tmp1;
    double tmp2= exp(-1*param);
    double tmp3= pow(cos(alpha),4.0);
    D= tmp2/(m*m*tmp3);

    //calculate Geo
    double temp1=dot(normal.x,normal.y,normal.z,H.x,H.y,H.z);
    double temp2=dot(view.x,view.y,view.z,H.x,H.y,H.z);
    double term1=2*temp1*dot(normal.x,normal.y,normal.z,view.x,view.y,view.z)/temp2;
    double term2=2*temp1*dot(normal.x,normal.y,normal.z,light.x,light.y,light.z)/temp2;
    Geo=min(1.0,term1,term2);

    //calculating Rs
    Rs=(F/pi)*D*Geo/(dot(normal.x,normal.y,normal.z,view.x,view.y,view.z)*dot(normal.x,normal.y,normal.z,light.x,light.y,light.z));

  /*************calculating BRDF****************************/
    BR= (s*Rs)+(d*Rd);

   col_val[0][0]+= (0.01*energy[L]+(BR*energy[L]))*match_func[L][0];  //X //(0.001*energy) is a part for ambient component
   col_val[1][0]+= (0.01*energy[L]+(BR*energy[L]))*match_func[L][1];  //Y
   col_val[2][0]+= (0.01*energy[L]+(BR*energy[L]))*match_func[L][2];  //Z

  }

   /****matrix multiplication to convert XYZ to sRGB****/

   int i,j,k;
   for ( i = 0; i < 3; i++)
  {
     for ( j = 0; j < 1; j++)
     {
        RGB[i][j] = 0;
        for ( k = 0; k < 3; k++)
        {
          RGB[i][j] = RGB[i][j] + sRGB[i][k] * col_val[k][j];
        }
     }

  }
  /*****Gamma Correction******/
  int nl,nn;

  for(nl=0;nl<3;nl++)
  {
      for(nn=0;nn<1;nn++)
      {
          if(RGB[nl][nn]<=0.00313)
          {
              RGB[nl][nn]=RGB[nl][nn]*12.92;
          }

          else
          {
              RGB[nl][nn]=(pow(RGB[nl][nn],0.4167)*1.055)-0.055;
          }
      }
  }
  R=dot(light.x, light.y, light.z, normal.x, normal.y, normal.z)*RGB[0][0];
  G=dot(light.x, light.y, light.z, normal.x, normal.y, normal.z)*RGB[1][0];
  B=dot(light.x, light.y, light.z, normal.x, normal.y, normal.z)*RGB[2][0];

  return;
}

// Draws the scene
void drawit(void)
{
	glDrawPixels(ImageW, ImageH, GL_RGB, GL_FLOAT, framebuffer);
	glFlush();
}

// Sets pixel x,y to the color RGB
void setFramebuffer(int x, int y, float R, float G, float B)
{
	if (R <= 1.0)
		if (R >= 0.0)
			framebuffer[y][x][0] = R;
		else
			framebuffer[y][x][0] = 0.0;
	else
		framebuffer[y][x][0] = 1.0;
	if (G <= 1.0)
		if (G >= 0.0)
			framebuffer[y][x][1] = G;
		else
			framebuffer[y][x][1] = 0.0;
	else
		framebuffer[y][x][1] = 1.0;
	if (B <= 1.0)
		if (B >= 0.0)
			framebuffer[y][x][2] = B;
		else
			framebuffer[y][x][2] = 0.0;
	else
		framebuffer[y][x][2] = 1.0;
}


void display(void)
{
	int i, j, k;
	float R, G, B;
	Coord3D refpt;
	Coord3D view;
	Coord3D normal;
	Coord3D light;
	Coord3D intpt;
	double xstep = 12.0 / ImageW;
	double ystep = 8.0 / ImageH;
	double t;
	double a, b, c;
	int intsphere;

	refpt.x = -6.0 + xstep / 2.0;
	refpt.y = -4.0 + ystep / 2.0;
	refpt.z = -10.0;

	for (i = 0; i<ImageW; i++, refpt.x += xstep) {
		for (j = 0; j<ImageH; j++, refpt.y += ystep) {
			// Compute the view vector
			view.x = refpt.x; view.y = refpt.y; view.z = refpt.z;
			normalize(view.x, view.y, view.z);

			// Find intersection with sphere (if any) - only 1 sphere can intesect.
			intsphere = -1;
			for (k = 0; (k<6) && (intsphere == -1); k++) {
				a = 1.0;  // Since normalized;
				b = 2.0*view.x*(-SphCent[k].x) + 2.0*view.y*(-SphCent[k].y) + 2.0*view.z*(-SphCent[k].z);
				c = SphCent[k].x*SphCent[k].x + SphCent[k].y*SphCent[k].y + SphCent[k].z*SphCent[k].z -
					SphRad[k] * SphRad[k];
				if ((b*b - 4 * a*c) >= 0.0) {  // We have an intersection with that sphere
											   // Want nearest of two intersections
					t = (-b - sqrt(b*b - 4 * a*c)) / 2.0;
					intsphere = k;
				}
			}

			if (intsphere != -1) { // We had an intersection with a sphere
				intpt.x = t*view.x; intpt.y = t*view.y; intpt.z = t*view.z;
				normal.x = (intpt.x - SphCent[intsphere].x) / SphRad[intsphere];
				normal.y = (intpt.y - SphCent[intsphere].y) / SphRad[intsphere];
				normal.z = (intpt.z - SphCent[intsphere].z) / SphRad[intsphere];
				normalize(normal.x, normal.y, normal.z);

				light.x = Light.x - intpt.x;
				light.y = Light.y - intpt.y;
				light.z = Light.z - intpt.z;
				normalize(light.x, light.y, light.z);
				GetColor(view, normal, light, intsphere, R, G, B);

			}
			else {
				R = G = B = 0.0;
			}
			setFramebuffer(i, j, R, G, B);
		}
		refpt.y = -4.0 + ystep / 2.0;
	}

	drawit();
}

void init(void)
{
	int i, j;

	// Initialize framebuffer to clear
	for (i = 0; i<ImageH; i++) {
		for (j = 0; j<ImageW; j++) {
			framebuffer[i][j][0] = 0.0;
			framebuffer[i][j][1] = 0.0;
			framebuffer[i][j][2] = 0.0;
		}
	}

	// Create Sphere data
	SphCent[0].x = -3.0;
	SphCent[0].y = 1.5;
	SphCent[0].z = -10.0;
	SphCent[1].x = 0.0;
	SphCent[1].y = 1.5;
	SphCent[1].z = -10.0;
	SphCent[2].x = 3.0;
	SphCent[2].y = 1.5;
	SphCent[2].z = -10.0;
	SphCent[3].x = -3.0;
	SphCent[3].y = -1.5;
	SphCent[3].z = -10.0;
	SphCent[4].x = 0.0;
	SphCent[4].y = -1.5;
	SphCent[4].z = -10.0;
	SphCent[5].x = 3.0;
	SphCent[5].y = -1.5;
	SphCent[5].z = -10.0;
	for (i = 0; i<6; i++) SphRad[i] = 1.0;

	// Set Light Position
	Light.x = 0.0;
	Light.y = 6.0;
	Light.z = 0.0;

	// Eye is at origin, looking down -z axis, y axis is up,
	// Looks at 8x6 window centered around z = -10.
}

int main(int argc, char** argv)
{

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
	glutInitWindowSize(ImageW, ImageH);
	glutInitWindowPosition(100, 100);
	glutCreateWindow("<Samia Kabir> - 641 Assignment 2");
	init();
	glutDisplayFunc(display);
	glutMainLoop();
	return 0;
}
