// Mandelbrot set example
// Adam Sampson <a.sampson@abertay.ac.uk>

#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <complex>
#include <fstream>
#include <iostream>
#include <amp.h>
#include <amp_math.h>
#include "mandelbrot.h"
#include <memory.h>

// Import things we need from the standard library
using std::chrono::duration_cast;
using std::chrono::milliseconds;

using std::complex;
using std::cout;
using std::endl;
using std::ofstream;
using namespace concurrency;

// Define the alias "the_clock" for the clock type we're going to use.
typedef std::chrono::steady_clock the_clock;

// The size of the image to generate.
const int WIDTH = 640;
const int HEIGHT = 480;
static const int TileSize = 16;
// The number of times to iterate before we assume that a point isn't in the
// Mandelbrot set.
// (You may need to turn this up if you zoom further into the set.)
const int MAX_ITERATIONS = 500;

// The image data.
// Each pixel is represented as 0xRRGGBB.
uint32_t image[HEIGHT][WIDTH];
uint32_t filteredImage[HEIGHT][WIDTH];
//Define kernel
float kernel[7][7]{
	{1,2,3,4,3,2,1},
	{2,3,4,5,4,3,2},
	{ 3,4,5,6,5,4,3 },
	{4,5,6,7,6,5,4},
	{3,4,5,6,5,4,3},
	{ 2,3,4,5,4,3,2 },
	{ 1,2,3,4,3,2,1 }
};

// using our own Complex number structure and definitions as the Complex type is not available
struct Complex1 {
	float x;
	float y;
};

Complex1 c_add(Complex1 c1, Complex1 c2) restrict(cpu, amp) 
{
	Complex1 tmp;
	float a = c1.x;
	float b = c1.y;
	float c = c2.x;
	float d = c2.y;
	tmp.x = a + c;
	tmp.y = b + d;
	return tmp;
} // c_add

float c_abs(Complex1 c) restrict(cpu, amp)
{
	return concurrency::fast_math::sqrt(c.x*c.x + c.y*c.y);
} // c_abs
Complex1 c_mul(Complex1 c1, Complex1 c2) restrict(cpu, amp)
{
	Complex1 tmp;
	float a = c1.x;
	float b = c1.y;
	float c = c2.x;
	float d = c2.y;
	tmp.x = a*c - b*d;
	tmp.y = b*c + a*d;
	return tmp;
}


// Write the image to a TGA file with the given name.
// Format specification: http://www.gamers.org/dEngine/quake3/TGA.txt
void write_tga(const char *filename)
{
	ofstream outfile(filename, ofstream::binary);

	uint8_t header[18] = {
		0, // no image ID
		0, // no colour map
		2, // uncompresse d 24-bit image
		0, 0, 0, 0, 0, // empty colour map specification
		0, 0, // X origin
		0, 0, // Y origin
		WIDTH & 0xFF, (WIDTH >> 8) & 0xFF, // width
		HEIGHT & 0xFF, (HEIGHT >> 8) & 0xFF, // height
		24, // bits per pixel
		0, // image descriptor
	};
	outfile.write((const char *)header, 18);

	for (int y = 0; y < HEIGHT; ++y)
	{
		for (int x = 0; x < WIDTH; ++x)
		{
			uint8_t pixel[3] = {
				image[y][x] & 0xFF, // blue channel
				(image[y][x] >> 8) & 0xFF, // green channel
				(image[y][x] >> 16) & 0xFF, // red channel
			};
			outfile.write((const char *)pixel, 3);
		}
	}

	outfile.close();
	if (!outfile)
	{
		// An error has occurred at some point since we opened the file.
		cout << "Error writing to " << filename << endl;
		exit(1);
	}
}


// Render the Mandelbrot set into the image array.
// The parameters specify the region on the complex plane to plot.
void compute_mandelbrot(float left, float right, float top, float bottom)
{
	
	uint32_t *pImage = &(image[0][0]);

	array_view<uint32_t, 2> a(HEIGHT, WIDTH, pImage);
	a.discard_data();

	the_clock::time_point start = the_clock::now();

	parallel_for_each(a.extent, [=](concurrency::index<2> idx) restrict(amp) {
		// compute Mandelbrot here i.e. Mandelbrot kernel / shader
		//USE THREAD ID/INDEX TO MAP INTO THE COMPLEX PLANE
		int y = idx[0];
		int x = idx[1];

		Complex1 c;
		c.x = left + (x * (right - left) / WIDTH);
		c.y = top + (y * (bottom - top) / HEIGHT);

		// Work out the point in the complex plane that
		// corresponds to this pixel in the output image.

		// Start off z at (0, 0).
		Complex1 z;
		z.x = 0.0;
		z.y = 0.0;
				

		// Iterate z = z^2 + c until z moves more than 2 units
		// away from (0, 0), or we've iterated too many times.
		int iterations = 0;
		while (c_abs(z) < 2.0 && iterations < MAX_ITERATIONS)
		{
			z = c_mul(z, z);
			z = c_add(z, c);

			++iterations;
		}

		if (iterations == MAX_ITERATIONS)
		{
			// z didn't escape from the circle.
			// This point is in the Mandelbrot set.
			a[y][x] = 0x000000;
		}

		else
		{
			float colour = ((float)iterations / (float)MAX_ITERATIONS) * (float)0xFFFFFF;
			a[y][x] = colour;
		}
				
	});

	the_clock::time_point end = the_clock::now();
	auto time_taken = duration_cast<milliseconds>(end - start).count();
	//cout << "Time taken to compute mandlebrot " << time_taken << "ms" << endl;

}

void AMPWarmUpRoutine() {

	uint32_t *pImage = &(image[0][0]);

	array_view<uint32_t, 2> fi(HEIGHT, WIDTH, pImage);
	parallel_for_each(fi.extent.tile<TileSize, TileSize>(), [=](tiled_index<TileSize, TileSize> tidx) restrict(amp) {});

}

//Function that genreates filter to apply over image pixels
void generateFilter() {

	//This function based upon equation in "Computer Vision"(2001) Linda G. Shapiro & George C. Stockman, p.149 

	float stDev = 1.0;
	float r = 2.0 * stDev * stDev;
	float s = 2.0 * stDev * stDev;
	float sum = 0.0;

	the_clock::time_point start = the_clock::now();

	//Generate Kernel

	for (int x = -3; x <= 3; x++)
	{
		for (int y = -3; y <= 3; y++)
		{
			r = fast_math::sqrt(x*x + y*y);
			kernel[x + 3][y + 3] = (fast_math::exp(-(r*r) / s)) / (3.14159265358979323846 * s);
			sum += kernel[x + 3][y + 3];
		}
	}

	//Normalise Kernel
	for (int i = 0; i < 7; ++i)
	{
		for (int j = 0; j < 7; ++j)
		{
			kernel[i][j] /= sum;
		}
	}

	the_clock::time_point end = the_clock::now();
	auto time_taken = duration_cast<milliseconds>(end - start).count();

	//Output Kernel
	//for (int i = 0; i < 7; ++i)
	//{
	//	for (int j = 0; j < 7; ++j)
	//	{
	//		cout << kernel[i][j]<< "\t";
	//	}
	//	cout << endl;
	//}

	//cout << "Time taken to compute kernel " << time_taken << "ms" << endl;
	//getchar();

}

void applyFilter(){

	int kernelRadius = 7 / 2;
	uint8_t r, g, b;
	float kernelTotal = 0.0;
	float redTotal = 0.0, blueTotal = 0.0, greenTotal = 0.0;

	the_clock::time_point start = the_clock::now();

	for (int y = 0; y < HEIGHT; y++) { //loop through image height
		for (int x = 0; x < WIDTH; x++) {//loop through image width

			float redTotal = 0.0, blueTotal = 0.0, greenTotal = 0.0;
			float kernelTotal = 0.0;

			for (int v = 0; v < 7; v++) { //loop through for kernel height
				for (int u = 0; u < 7; u++) { //loop through for kernel width

					// Current position
					int cX = x + u - kernelRadius;
					int cY = y + v - kernelRadius;

					//Make sure we stay in boundries
					if (cX < 0 || cX > WIDTH - 1 || cY < 0 || cY > HEIGHT - 1)
					{
						continue;
					}

					//Get colour channels of current pixel
					r = (image[cY][cX] >> 16);
					g = (image[cY][cX] >> 8);
					b = (image[cY][cX]);

					//Calculate Totals
					redTotal += r *kernel[v][u];
					greenTotal += g *kernel[v][u];
					blueTotal += b *kernel[v][u];
					kernelTotal += kernel[v][u];

				}
			}

			//Calculate new pixel values
			r = (redTotal / kernelTotal);
			g = (greenTotal / kernelTotal);
			b = (blueTotal / kernelTotal);

			//Generate new weighted pixel
			image[y][x] = r << 16 | g << 8 | b;

		}
	}

	the_clock::time_point end = the_clock::now();
	auto time_taken = duration_cast<milliseconds>(end - start).count();
	cout << "Time taken to apply kernel(sequential) " << time_taken << "ms" << endl;

}

void applyFilterAMP() {

	uint32_t *pImage = &(image[0][0]);
	uint32_t *pFilteredImage = &(filteredImage[0][0]);
	float *pKernel = &(kernel[0][0]);

	array_view<float, 2>k(7, 7, pKernel);
	array_view<uint32_t, 2> m(HEIGHT, WIDTH, pImage);
	array_view<uint32_t, 2> fi(HEIGHT, WIDTH, pFilteredImage);

	fi.discard_data();

	the_clock::time_point start = the_clock::now();

	//parallel_for_each(fi.extent, [=](concurrency::index<2> idx) restrict(amp) {
	parallel_for_each(fi.extent.tile<TileSize, TileSize>(), [=](tiled_index<TileSize, TileSize> tidx) restrict(amp) {

		int kernelRadius = 7 / 2;
		//int y = idx[0];
		//int x = idx[1];
		int y = tidx.global[0];
		int x = tidx.global[1];

			float kernelTotal = 0.0;
			float redTotal = 0.0, blueTotal = 0.0, greenTotal = 0.0;

			uint32_t r, g, b, newPixel;

			for (int v = 0; v < 7; v++) { //loop through for kernel height
				for (int u = 0; u < 7; u++) { //loop through for kernel width

					// Current position
					int cX = tidx.global[1] + u - kernelRadius;
					int cY = tidx.global[0] + v - kernelRadius;

					//Make sure we stay in boundries
					if (cX < 0 || cX > WIDTH - 1 || cY < 0 || cY > HEIGHT - 1)
					{
						continue;
					}

					//Get colour channels of pixel
					r = (m[cY][cX] >> 16) & 0x000000FF;
					g = (m[cY][cX] >> 8) & 0x000000FF;
					b = m[cY][cX] & 0x000000FF;

					//Calculate Totals
					redTotal += r *k[v][u];
					greenTotal += g *k[v][u];
					blueTotal += b *k[v][u];
					kernelTotal += k[v][u];

				}
			}

			//Calculate new pixel values
			r = (redTotal / kernelTotal);
			g = (greenTotal / kernelTotal);
			b = (blueTotal / kernelTotal);

			//Generate new weighted pixel
			fi[y][x] = (((r & 0x000000FF) << 16) | ((g & 0x000000FF) << 8) | ((b & 0x000000FF)));

	});

	//Sync back to CPU
	fi.synchronize();

	the_clock::time_point end = the_clock::now();
	auto time_taken = duration_cast<milliseconds>(end - start).count();
	cout << "Time taken to apply kernel(parallel) " << time_taken << "ms" << endl;

	//Copy from filtered image, to mandelbrot
	errno_t err = memcpy_s(image, sizeof(uint32_t)*(HEIGHT * WIDTH), filteredImage, sizeof(uint32_t)*(HEIGHT * WIDTH));
	if (err)
	{
		printf("Error executing memcpy_s.\n");
	}

}

int main(int argc, char *argv[])
{

	cout << "Working on Image of Size " << WIDTH << " x " << HEIGHT << " using a tile size of " << TileSize << " x " << TileSize << endl;

	cout << "Please wait..." << endl;

	AMPWarmUpRoutine();//Get AMP ready

	generateFilter(); // Function call to create a filter

	// This shows the whole set.
	//compute_mandelbrot(-2.0, 1.0, 1.125, -1.125);

	// This zooms in on an interesting bit of detail.
	compute_mandelbrot(-0.751085, -0.734975, 0.118378, 0.134488);

	write_tga("before.tga"); //Output before filter

	applyFilter();//Apply filter

	write_tga("SequentialAfter.tga"); //Output after filter
	
	// This shows the whole set
	//compute_mandelbrot(-2.0, 1.0, 1.125, -1.125);

	// This zooms in on an interesting bit of detail.
	compute_mandelbrot(-0.751085, -0.734975, 0.118378, 0.134488);
	
	applyFilterAMP();//Apply filter with AMP

	write_tga("ParallelAfter.tga"); //Output after filter

	cout << "Finished..." << endl;

	getchar();
}