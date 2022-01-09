#include <stdbool.h> 

#define IMG_SIZE 3*m*n
#define min(a,b)  (a < b ? a : b)
#define max(a,b)  (a > b ? a : b)
#define LOCATION(i,j,n) (i*n+j)
typedef struct {
   unsigned char red;
   unsigned char green;
   unsigned char blue;
} pixel;

typedef struct {
    int red;
    int green;
    int blue;
} pixel_sum;


/* Compute min and max of two integers, respectively */
//int min(int a, int b) { return (a < b ? a : b); }
//int max(int a, int b) { return (a > b ? a : b); }

int calcIndex(int i, int j, int n) {
	return ((i)*(n)+(j));
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
static void assign_sum_to_pixel(pixel *current_pixel, pixel_sum sum, int kernelScale) { //divied by 9 or 7 and copy sum to current pixel

	// divide by kernel's weight
	sum.red = sum.red / kernelScale;
	sum.green = sum.green / kernelScale;
	sum.blue = sum.blue / kernelScale;

	// truncate each pixel's color values to match the range [0,255]
	current_pixel->red = (unsigned char) (min(max(sum.red, 0), 255));
	current_pixel->green = (unsigned char) (min(max(sum.green, 0), 255));
	current_pixel->blue = (unsigned char) (min(max(sum.blue, 0), 255));
	return;
}

/*
* sum_pixels_by_weight - Sums pixel values, scaled by given weight
*/
static void sum_pixels_by_weight(pixel_sum *sum, pixel p, int weight) { // add to sum the pixel
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
    int currRow, currCol;
    pixel_sum sum;
    pixel current_pixel;
    int min_intensity = 766; // arbitrary value that is higher than maximum possible intensity, which is 255*3=765
    int max_intensity = -1; // arbitrary value that is lower than minimum possible intensity, which is 0
    int min_row, min_col, max_row, max_col;
    pixel loop_pixel;

    initialize_pixel_sum(&sum);

    for(ii = max(i-1, 0); ii <= min(i+1, dim-1); ii++) {
        for(jj = max(j-1, 0); jj <= min(j+1, dim-1); jj++) {

            int kRow, kCol;

            // compute row index in kernel
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
            } else {
                kCol = 1;
            }

            // apply kernel on pixel at [ii,jj]
            sum_pixels_by_weight(&sum, src[calcIndex(ii, jj, dim)], kernel[kRow][kCol]);
        }
    }

    if (filter) {
        // find min and max coordinates
        for(ii = max(i-1, 0); ii <= min(i+1, dim-1); ii++) {
            for(jj = max(j-1, 0); jj <= min(j+1, dim-1); jj++) {
                // check if smaller than min or higher than max and update
                loop_pixel = src[calcIndex(ii, jj, dim)];
                if ((((int) loop_pixel.red) + ((int) loop_pixel.green) + ((int) loop_pixel.blue)) <= min_intensity) {
                    min_intensity = (((int) loop_pixel.red) + ((int) loop_pixel.green) + ((int) loop_pixel.blue));
                    min_row = ii;
                    min_col = jj;
                }
                if ((((int) loop_pixel.red) + ((int) loop_pixel.green) + ((int) loop_pixel.blue)) > max_intensity) {
                    max_intensity = (((int) loop_pixel.red) + ((int) loop_pixel.green) + ((int) loop_pixel.blue));
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
void smoothBlurNoFilter(pixel *src, pixel *dst) {
    //int dimLoop =dim-1;
    int dimLoop =m-1;
    int i,j;
    for (i = 1 ; i < dimLoop; i++) {
        for (j =  1 ; j < dimLoop; j++) {
            //int location = calcIndex(i, j, m);
            int location = LOCATION(i,j,m);
            pixel_sum sum;
            sum.red = 0, sum.green = 0, sum.blue = 0;
            int k;
            // sum all pixels around current pixel:
            for (k = i-1; k <= i+1; ++k) { // for a row in src
                //int currLoc = k*m +j; //not sure about it - we can do with the func calcIndex
                int currLoc = LOCATION(k,j,m);
                sum.red += src[currLoc-1].red;
                sum.red += src[currLoc].red;
                sum.red += src[currLoc+1].red;
                sum.green += src[currLoc-1].green;
                sum.green += src[currLoc].green;
                sum.green += src[currLoc+1].green;
                sum.blue += src[currLoc-1].blue;
                sum.blue += src[currLoc].blue;
                sum.blue += src[currLoc+1].blue;
            }
            dst[location].red = sum.red /9;
            dst[location].green = sum.green /9;
            dst[location].blue = sum.blue /9;
//            dst[location].red = sum.red * 0.1111111111;
//            dst[location].green = sum.green * 0.1111111111;
//            dst[location].blue = sum.blue * 0.1111111111;
        }
    }
}

void smoothBlurWithFilter(pixel *src, pixel *dst) {
    int dimLoop =m-1; /// maybe in a different way
    int i,j;
    for (i = 1 ; i < dimLoop; i++) {
        for (j =  1 ; j < dimLoop; j++) {
            int location = LOCATION(i, j, m);
            //int location = calcIndex(i, j, m);
            pixel_sum sum;
            sum.red = 0, sum.green = 0, sum.blue = 0;
            int maxi = 776, mini = -1;
            int sumTotal, minRow, maxRow, minCol, maxCol,k;
            // sum all pixels around current pixel:
            for (k = i-1; k <= i+1; ++k) { // for a row in src
                //int currLoc = k*m +j; //not sure about it - we can do with the func calcIndex
                int currLoc = LOCATION(k,j,m);

                sum.red += src[currLoc-1].red;
                sum.green += src[currLoc-1].green;
                sum.blue += src[currLoc-1].blue;
                sumTotal = (int) (src[currLoc-1].red + src[currLoc-1].green + src[currLoc-1].blue);
                if (sumTotal <= maxi) {
                    maxi = sumTotal;
                    maxRow = k, maxCol = j-1;
                }
                if (sumTotal > mini) {
                    mini = sumTotal;
                    minRow = k, minCol = j-1;
                }

                sum.red += src[currLoc].red;
                sum.green += src[currLoc].green;
                sum.blue += src[currLoc].blue;
                sumTotal = (int) (src[currLoc].red + src[currLoc].green + src[currLoc].blue);
                if (sumTotal <= maxi) {
                    maxi = sumTotal;
                    maxRow = k, maxCol = j;
                }
                if (sumTotal > mini) {
                    mini = sumTotal;
                    minRow = k, minCol = j;
                }

                sum.red += src[currLoc+1].red;
                sum.green += src[currLoc+1].green;
                sum.blue += src[currLoc+1].blue;
                sumTotal = (int) (src[currLoc+1].red + src[currLoc+1].green + src[currLoc+1].blue);
                if (sumTotal <= maxi) {
                    maxi = sumTotal;
                    maxRow = k, maxCol = j+1;
                }
                if (sumTotal > mini) {
                    mini = sumTotal;
                    minRow = k, minCol = j+1;
                }
            }
            int minLoc = minRow*m +minCol;
            int maxLoc = maxRow*m +maxCol;
            sum.red += (-1) * src[minLoc].red;
            sum.green += (-1) * src[minLoc].green;
            sum.blue += (-1) * src[minLoc].blue;
            sum.red += (-1) * src[maxLoc].red;
            sum.green += (-1) * src[maxLoc].green;
            sum.blue += (-1) * src[maxLoc].blue;
            //divide by 7:
            dst[location].red = sum.red * 0.1428571429;
            dst[location].green = sum.green * 0.1428571429;
            dst[location].blue = sum.blue * 0.1428571429;
        }
    }
}

void smoothSharpNoFilter(pixel *src, pixel *dst) {
    int dimLoop =m-1; /// maybe in a different way
    int i,j,k;
    for (i = 1 ; i < dimLoop; i++) {
        for (j =  1 ; j < dimLoop; j++) {
            int location = LOCATION(i, j, m);
            pixel_sum sum;
            sum.red = 0, sum.green = 0, sum.blue = 0;
            // sum all pixels around current pixel:
            for (k = i-1; k <= i+1; k+=2) { // for a row in src-  2 rows- first and last (1,3)
                int currLoc = LOCATION(k,j,m);
                //int currLoc = k*m +j; //not sure about it - we can do with the func calcIndex
                sum.red += src[currLoc-1].red;
                sum.green += src[currLoc-1].green;
                sum.blue +=  src[currLoc-1].blue;

                sum.red +=  src[currLoc].red;
                sum.green +=  src[currLoc].green;
                sum.blue +=  src[currLoc].blue;

                sum.red +=  src[currLoc+1].red;
                sum.green +=  src[currLoc+1].green;
                sum.blue +=  src[currLoc+1].blue;
            }
            // instead of duplicate every iteration in k loop:
            sum.red *= (-1);
            sum.green *= (-1);
            sum.blue *= (-1);
            //add the middle row to sum:
            int currLoc = LOCATION(i,j,m);
            sum.red += (-1) * src[currLoc-1].red;
            sum.green += (-1) * src[currLoc-1].green;
            sum.blue += (-1) * src[currLoc-1].blue;

            sum.red += 9 * src[currLoc].red;
            sum.green += 9 * src[currLoc].green;
            sum.blue += 9 * src[currLoc].blue;

            sum.red += (-1) * src[currLoc+1].red;
            sum.green += (-1) * src[currLoc+1].green;
            sum.blue += (-1) * src[currLoc+1].blue;

            dst[location].red = (unsigned char) (min(max(sum.red, 0), 255));
            dst[location].green = (unsigned char) (min(max(sum.green, 0), 255));
            dst[location].blue = (unsigned char) (min(max(sum.blue, 0), 255));
        }
    }
}

void smooth(int dim, pixel *src, pixel *dst, int kernelSize, int kernel[kernelSize][kernelSize], int kernelScale, bool filter) {

    int i, j;
    for (i = kernelSize / 2 ; i < dim - kernelSize / 2; i++) {
        for (j =  kernelSize / 2 ; j < dim - kernelSize / 2 ; j++) {
            dst[calcIndex(i, j, dim)] = applyKernel(dim, i, j, src, kernelSize, kernel, kernelScale, filter);
        }
    }
}
void charsToPixels(Image *charsImg, pixel* pixels) {

	int row, col;
	for (row = 0 ; row < m ; row++) {
		for (col = 0 ; col < n ; col++) {

			pixels[row*n + col].red = image->data[3*row*n + 3*col];
			pixels[row*n + col].green = image->data[3*row*n + 3*col + 1];
			pixels[row*n + col].blue = image->data[3*row*n + 3*col + 2];
		}
	}
}

void pixelsToChars(pixel* pixels, Image *charsImg) {

	int row, col;
	for (row = 0 ; row < m ; row++) {
		for (col = 0 ; col < n ; col++) {

			image->data[3*row*n + 3*col] = pixels[row*n + col].red;
			image->data[3*row*n + 3*col + 1] = pixels[row*n + col].green;
			image->data[3*row*n + 3*col + 2] = pixels[row*n + col].blue;
		}
	}
}

void copyPixels(pixel* src, pixel* dst) {

	int row, col;
	for (row = 0 ; row < m ; row++) {
		for (col = 0 ; col < n ; col++) {

			dst[row*n + col].red = src[row*n + col].red;
			dst[row*n + col].green = src[row*n + col].green;
			dst[row*n + col].blue = src[row*n + col].blue;
		}
	}
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

void doConvolutionBlurNoFilter(Image *image) {

	pixel* pixelsImg = malloc(m*n*sizeof(pixel));
	pixel* backupOrg = malloc(m*n*sizeof(pixel));

    //coping image:
    memcpy(pixelsImg, image->data, IMG_SIZE);
    memcpy(backupOrg, pixelsImg, IMG_SIZE);

    smoothBlurNoFilter(backupOrg, pixelsImg);

    memcpy(image->data, pixelsImg, IMG_SIZE);

	free(pixelsImg);
	free(backupOrg);
}
void doConvolutionBlurWithFilter(Image *image) {

    pixel* pixelsImg = malloc(m*n*sizeof(pixel));
    pixel* backupOrg = malloc(m*n*sizeof(pixel));

    //coping image:
    memcpy(pixelsImg, image->data, IMG_SIZE);
    memcpy(backupOrg, pixelsImg, IMG_SIZE);

    smoothBlurWithFilter(backupOrg, pixelsImg);

    memcpy(image->data, pixelsImg, IMG_SIZE);

    free(pixelsImg);
    free(backupOrg);
}

void doConvolutionSharp(Image *image) {
    int picSize = m*n*sizeof(pixel);
    pixel* pixelsImg = malloc(picSize);
    pixel* backupOrg = malloc(picSize);

    //coping image:
    memcpy(pixelsImg, image->data, IMG_SIZE);
    memcpy(backupOrg, pixelsImg, IMG_SIZE);

    smoothSharpNoFilter(backupOrg, pixelsImg);

    memcpy(image->data, pixelsImg, IMG_SIZE);

    free(pixelsImg);
    free(backupOrg);
}

void myfunction(Image *image, char* srcImgpName, char* blurRsltImgName, char* sharpRsltImgName,
                char* filteredBlurRsltImgName, char* filteredSharpRsltImgName, char flag) {
    if (flag == '1') {
		// blur image
		doConvolutionBlurNoFilter(image);

		// write result image to file
		writeBMP(image, srcImgpName, blurRsltImgName);	

		// sharpen the resulting image
        doConvolutionSharp(image);
		
		// write result image to file
		writeBMP(image, srcImgpName, sharpRsltImgName);	
	} else {
		// apply extermum filtered kernel to blur image
        doConvolutionBlurWithFilter(image);

		// write result image to file
		writeBMP(image, srcImgpName, filteredBlurRsltImgName);

		// sharpen the resulting image
        doConvolutionSharp(image);

		// write result image to file
		writeBMP(image, srcImgpName, filteredSharpRsltImgName);	
	}
}

