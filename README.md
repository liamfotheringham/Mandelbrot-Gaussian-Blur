# Mandelbrot-Gaussian-Blur
Application generates a Mandelbrot Set and applies a gaussian blur to the generated image whilst utilising parallel processing.

## The Application
* Applies a 7x7 Generated Gaussian Filter on a Mandelbrot Image
* Uses AMP GPGPU Processing to Apply Filter and Calculate New Weighted Pixel
* Task Decomposed Using Tiled GPGPU Processing
## Calculation Process
![alt text](https://github.com/liamfotheringham/Mandelbrot-Gaussian-Blur/blob/master/ReadMe%20images/CalculationProcess.png)
## Task Decomposition
![alt text](https://github.com/liamfotheringham/Mandelbrot-Gaussian-Blur/blob/master/ReadMe%20images/TeskDecomposition.png)
## Generated Outputs
![alt text](https://github.com/liamfotheringham/Mandelbrot-Gaussian-Blur/blob/master/ReadMe%20images/GeneratedOutputs.png)
## Critical Evaluation
Rethink Use of Tiles

Utilise tile_static:
  * Reduce Global Memory Reads
  * Would Require Barriers
  * Memory Overhead (Cache)
  * Wasteful Use of Threads

Separable Filters:
  * O(Kernel Width * Image Width * Image Height) + O(Kernel Height * Image Width * Image Height) < O(Kernel Width * Kernel Height * Image Width * Image Height)

Further Developments:
  * Sobel Filter (Edge Detection)
  * Load Images Into Application
