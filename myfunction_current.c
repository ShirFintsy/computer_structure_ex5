#include <stdbool.h> 

typedef struct {
   unsigned char red;
   unsigned char green;
   unsigned char blue;
} pixel;

typedef struct {
    int red;
    int green;
    int blue;
    // int num;
} pixel_sum;


/* Compute min and max of two integers, respectively */
int min(int a, int b) { return (a < b ? a : b); }
int max(int a, int b) { return (a > b ? a : b); }

int calcIndex(int i, int j, int n) { ///////////////////:)
	//return ((i)*(n)+(j));
    int dup = i*n;
    int sum = dup+j;
    return sum;
}

/*
 * initialize_pixel_sum - Initializes all fields of sum to 0
 */
void initialize_pixel_sum(pixel_sum *sum) {
	sum->red = sum->green = sum->blue = 0;
	// sum->num = 0;
	return;
}

/*
 * assign_sum_to_pixel - Truncates pixel's new value to match the range [0,255]
 */
static void assign_sum_to_pixel(pixel *current_pixel, pixel_sum sum, int kernelScale) {

	// divide by kernel's weight
    sum.red /= kernelScale;
    sum.green /= kernelScale;
    sum.blue /= kernelScale;
//	sum.red = sum.red / kernelScale;
//	sum.green = sum.green / kernelScale;
//	sum.blue = sum.blue / kernelScale;

	// truncate each pixel's color values to match the range [0,255]
	current_pixel->red = (unsigned char) (min(max(sum.red, 0), 255));
	current_pixel->green = (unsigned char) (min(max(sum.green, 0), 255));
	current_pixel->blue = (unsigned char) (min(max(sum.blue, 0), 255));
	return;
}

/*
* sum_pixels_by_weight - Sums pixel values, scaled by given weight
*/
static void sum_pixels_by_weight(pixel_sum *sum, pixel p, int weight) {
	sum->red += ((int) p.red) * weight;
	sum->green += ((int) p.green) * weight;
	sum->blue += ((int) p.blue) * weight;
	// sum->num++;
	return;
}

/*
 *  Applies kernel for pixel at (i,j)
 */
static pixel applyKernel(int dim, int i, int j, pixel *src, int kernelSize, int kernel[kernelSize][kernelSize], int kernelScale, bool filter) {

	int ii, jj;
	//int currRow, currCol;
	pixel_sum sum;
	pixel current_pixel;
	int min_intensity = 766; // arbitrary value that is higher than maximum possible intensity, which is 255*3=765
	int max_intensity = -1; // arbitrary value that is lower than minimum possible intensity, which is 0
	int min_row, min_col, max_row, max_col;
	pixel loop_pixel;

	initialize_pixel_sum(&sum);

    // The calculate is out of the for loop- so it won't calculate every iteration:
    int mini = min(i+1, dim-1);
    int minij = min(j+1, dim-1);
    for(ii = max(i-1, 0); ii <= mini; ii++) {
        for(jj = max(j-1, 0); jj <= minij; jj++) {
            // Initialize the parameters so we could save an if conditional:
			int kRow, kCol;

			// compute row index in kernel
//            if (ii < i) {
//                kRow = 0;
//            } else if (ii > i) {
//				kRow = 2;
//            }

			if (ii < i) {
				kRow = 0;
			} else if (ii > i) {
				kRow = 2;
			} else {
				kRow = 1;
			}

			// compute column index in kernel
			if (jj < j) {
				kCol = 0;
			} else if (jj > j) {
				kCol = 2;
			}
            else {
				kCol = 1;
			}

			// apply kernel on pixel at [ii,jj]
			sum_pixels_by_weight(&sum, src[calcIndex(ii, jj, dim)], kernel[kRow][kCol]);
		}
	}

	if (filter) {
		// find min and max coordinates
        // The calculate is out of the for loop- so it won't calculate every iteration:
        int mini = min(i+1, dim-1), minij = min(j+1, dim-1);
        for(ii = max(i-1, 0); ii <= mini; ii++) {
            for(jj = max(j-1, 0); jj <= minij; jj++) {
				// check if smaller than min or higher than max and update
				loop_pixel = src[calcIndex(ii, jj, dim)];
                //int sumPixel = ((int) loop_pixel.red) + ((int) loop_pixel.green) + ((int) loop_pixel.blue); ////:)
                int sumPixel = loop_pixel.red + loop_pixel.green + loop_pixel.blue;
				if (sumPixel <= min_intensity) {
					min_intensity = sumPixel;
					min_row = ii;
					min_col = jj;
				}
				if (sumPixel > max_intensity) {
					max_intensity = sumPixel;
					max_row = ii;
					max_col = jj;
				}
			}
		}
		// filter out min and max
		sum_pixels_by_weight(&sum, src[calcIndex(min_row, min_col, dim)], -1);
		sum_pixels_by_weight(&sum, src[calcIndex(max_row, max_col, dim)], -1);
	}

	// assign kernel's result to pixel at [i,j]
	assign_sum_to_pixel(&current_pixel, sum, kernelScale);
	return current_pixel;
}

/*
* Apply the kernel over each pixel.
* Ignore pixels where the kernel exceeds bounds. These are pixels with row index smaller than kernelSize/2 and/or
* column index smaller than kernelSize/2
*/
void smooth(int dim, pixel *src, pixel *dst, int kernelSize, int kernel[kernelSize][kernelSize], int kernelScale, bool filter) {

	int divKernel = kernelSize / 2;
    int dis = dim - divKernel;
	for (int i = divKernel; i < dis; i++) {
		for (int j = divKernel; j < dis ; j++) {
			dst[calcIndex(i, j, dim)] = applyKernel(dim, i, j, src, kernelSize, kernel, kernelScale, filter);
		}
	}
}

void charsToPixels(Image *charsImg, pixel* pixels) {

	int row, col;
	for (row = 0 ; row < m ; row++) {
		for (col = 0 ; col < n ; col++) {
            int rowNCol = row*n + col;
            int dupRowCol = 3*rowNCol;
			pixels[rowNCol].red = image->data[dupRowCol];
			pixels[rowNCol].green = image->data[dupRowCol + 1];
			pixels[rowNCol].blue = image->data[dupRowCol + 2];
		}
	}
}

void pixelsToChars(pixel* pixels, Image *charsImg) {

	int row, col;
	for (row = 0 ; row < m ; row++) {
		for (col = 0 ; col < n ; col++) {
            int rowNCol = n*row + col;
            int dupRowCol = 3*rowNCol; //sum insted of duplicate in 3.
			image->data[dupRowCol] = pixels[rowNCol].red;
			image->data[dupRowCol + 1] = pixels[rowNCol].green;
			image->data[dupRowCol + 2] = pixels[rowNCol].blue;
		}
	}
}

void copyPixels(pixel* src, pixel* dst) { ///////////:)

	int row, col;
    int tempN = n/2;
    int rowN = 0; // no duplicate
    for (row = 0 ; row < m ; row++) {
        for (col = 0 ; col < tempN ; col++) {
            int rowNCol = rowN + col;
            dst[rowNCol].red = src[rowNCol].red;
            dst[rowNCol].green = src[rowNCol].green;
            dst[rowNCol].blue = src[rowNCol].blue;
            //loop unrolling:
            dst[rowNCol+tempN].red = src[rowNCol+tempN].red;
            dst[rowNCol+tempN].green = src[rowNCol+tempN].green;
            dst[rowNCol+tempN].blue = src[rowNCol+tempN].blue;
        }
        rowN +=n; // insted we have adding
    }
//	for (row = 0 ; row < m ; row++) {
//		for (col = 0 ; col < tempN ; col++) {
//            int rowNCol = row*n + col;
//			dst[rowNCol].red = src[rowNCol].red;
//			dst[rowNCol].green = src[rowNCol].green;
//			dst[rowNCol].blue = src[rowNCol].blue;
//            //loop unrolling:
//            dst[rowNCol+tempN].red = src[rowNCol+tempN].red;
//            dst[rowNCol+tempN].green = src[rowNCol+tempN].green;
//            dst[rowNCol+tempN].blue = src[rowNCol+tempN].blue;
//		}
//	}
}

void doConvolution(Image *image, int kernelSize, int kernel[kernelSize][kernelSize], int kernelScale, bool filter) {

	pixel* pixelsImg = malloc(m*n*sizeof(pixel));
	pixel* backupOrg = malloc(m*n*sizeof(pixel));

	charsToPixels(image, pixelsImg);
	copyPixels(pixelsImg, backupOrg);

	smooth(m, backupOrg, pixelsImg, kernelSize, kernel, kernelScale, filter);

	pixelsToChars(pixelsImg, image);

	free(pixelsImg);
	free(backupOrg);
}

void myfunction(Image *image, char* srcImgpName, char* blurRsltImgName, char* sharpRsltImgName, char* filteredBlurRsltImgName, char* filteredSharpRsltImgName, char flag) {

	/*
	* [1, 1, 1]
	* [1, 1, 1]
	* [1, 1, 1]
	*/
	int blurKernel[3][3] = {{1, 1, 1}, {1, 1, 1}, {1, 1, 1}};
//    int blurKernel[3][3];
//    for (int i = 0; i < 3; ++i) {
//        blurKernel[i][0] = 1;
//        blurKernel[i][1] = 1;
//        blurKernel[i][2] = 1;
//    }

	/*
	* [-1, -1, -1]
	* [-1, 9, -1]
	* [-1, -1, -1]
	*/
	int sharpKernel[3][3] = {{-1,-1,-1},{-1,9,-1},{-1,-1,-1}};

	if (flag == '1') {	
		// blur image
		doConvolution(image, 3, blurKernel, 9, false);

		// write result image to file
		writeBMP(image, srcImgpName, blurRsltImgName);	

		// sharpen the resulting image
		doConvolution(image, 3, sharpKernel, 1, false);
		
		// write result image to file
		writeBMP(image, srcImgpName, sharpRsltImgName);	
	} else {
		// apply extermum filtered kernel to blur image
		doConvolution(image, 3, blurKernel, 7, true);

		// write result image to file
		writeBMP(image, srcImgpName, filteredBlurRsltImgName);

		// sharpen the resulting image
		doConvolution(image, 3, sharpKernel, 1, false);

		// write result image to file
		writeBMP(image, srcImgpName, filteredSharpRsltImgName);	
	}
}

