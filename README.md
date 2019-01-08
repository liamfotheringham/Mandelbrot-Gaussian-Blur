# Mandelbrot-Gaussian-Blur
Application generates a Mandelbrot Set and applies a gaussian blur to the generated image whilst utilising parallel processing.

Gaussian Blur will be applied to the image using sequential processing aswell as parallel processing, with both images being exported in .tga image format.

The time for each process is calculated and displayed to the console.

## The Application
* Applies a 7x7 Generated Gaussian Filter on a Mandelbrot Image
* Uses AMP GPGPU Processing to Apply Filter and Calculate New Weighted Pixel
* Task Decomposed Using Tiled GPGPU Processing
## Calculation Process
![alt text](https://github.com/liamfotheringham/Mandelbrot-Gaussian-Blur/blob/master/ReadMe%20images/CalculationProcess.png)
## Task Decomposition
![alt text](https://github.com/liamfotheringham/Mandelbrot-Gaussian-Blur/blob/master/ReadMe%20images/TeskDecomposition.png)
* Breaks Image into 'Tiles'
* Works Through the Extent, Calculating Each Weighted Pixel
## Generated Outputs
![alt text](https://github.com/liamfotheringham/Mandelbrot-Gaussian-Blur/blob/master/ReadMe%20images/GeneratedOutputs.png)
