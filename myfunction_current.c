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


int calcIndex(int i, int j, int n) {}

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
    register int dimLoop =m-1;
    register int i,j;
    for (i = 1 ; i < dimLoop; i++) {
        for (j =  1 ; j < dimLoop; j++) {
            register int location = LOCATION(i,j,m);
            register pixel_sum sum;
            sum.red = 0, sum.green = 0, sum.blue = 0;
            register int k;
            // sum all pixels around current pixel:
            for (k = i-1; k <= i+1; ++k) { // for a row in src
                register int currLoc = LOCATION(k,j,m);
                register pixel p1 =  src[currLoc-1];
                register pixel p2 =  src[currLoc];
                register pixel p3 =  src[currLoc+1];

                sum.red += p1.red;
                sum.green += p1.green;
                sum.blue += p1.blue;

                sum.red += p2.red;
                sum.green += p2.green;
                sum.blue += p2.blue;

                sum.red += p3.red;
                sum.green += p3.green;
                sum.blue += p3.blue;
            }
            dst[location].red = sum.red / 9;
            dst[location].green = sum.green / 9;
            dst[location].blue = sum.blue / 9;
        }
    }
}

void smoothBlurWithFilter(pixel *src, pixel *dst) {
    register int dimLoop =m-1; /// maybe in a different way
    register int i,j;
    for (i = 1 ; i < dimLoop; i++) {
        for (j =  1 ; j < dimLoop; j++) {
            register int location = LOCATION(i, j, m);
            pixel_sum sum;
            sum.red = 0, sum.green = 0, sum.blue = 0;
            register int maxi = 776, mini = -1;
            register int sumTotal, minRow, maxRow, minCol, maxCol,k;
            // sum all pixels around current pixel:
            for (k = i-1; k <= i+1; ++k) { // for a row in src
                register int currLoc = LOCATION(k,j,m);
                register pixel p1 =  src[currLoc-1];
                register pixel p2 =  src[currLoc];
                register pixel p3 =  src[currLoc+1];

                sum.red += p1.red;
                sum.green += p1.green;
                sum.blue += p1.blue;
                sumTotal = p1.red + p1.green + p1.blue;
                if (sumTotal <= maxi) {
                    maxi = sumTotal;
                    maxRow = k, maxCol = j-1;
                }
                if (sumTotal > mini) {
                    mini = sumTotal;
                    minRow = k, minCol = j-1;
                }

                sum.red += p2.red;
                sum.green += p2.green;
                sum.blue += p2.blue;
                sumTotal = p2.red + p2.green + p2.blue;
                if (sumTotal <= maxi) {
                    maxi = sumTotal;
                    maxRow = k, maxCol = j;
                }
                if (sumTotal > mini) {
                    mini = sumTotal;
                    minRow = k, minCol = j;
                }

                sum.red += p3.red;
                sum.green += p3.green;
                sum.blue += p3.blue;
                sumTotal = p3.red + p3.green + p3.blue;
                if (sumTotal <= maxi) {
                    maxi = sumTotal;
                    maxRow = k, maxCol = j+1;
                }
                if (sumTotal > mini) {
                    mini = sumTotal;
                    minRow = k, minCol = j+1;
                }
            }
            register int minLoc = minRow*m +minCol;
            register int maxLoc = maxRow*m +maxCol;
            sum.red -= src[minLoc].red;
            sum.green -= src[minLoc].green;
            sum.blue -= src[minLoc].blue;
            sum.red -= src[maxLoc].red;
            sum.green -= src[maxLoc].green;
            sum.blue -= src[maxLoc].blue;
            //divide by 7:
            dst[location].red = sum.red * 0.1428571429;
            dst[location].green = sum.green * 0.1428571429;
            dst[location].blue = sum.blue * 0.1428571429;
        }
    }
}

void smoothSharpNoFilter(pixel *src, pixel *dst) {
    register int dimLoop =m-1; /// maybe in a different way
    register int i,j,k;
    for (i = 1 ; i < dimLoop; i++) {
        for (j =  1 ; j < dimLoop; j++) {
            register int location = LOCATION(i, j, m);
            register pixel_sum sum;
            sum.red = 0, sum.green = 0, sum.blue = 0;
            // sum all pixels around current pixel:
            for (k = i-1; k <= i+1; k+=2) { // for a row in src-  2 rows- first and last (1,3)
                register int currLoc = LOCATION(k,j,m);
                register pixel p1 =  src[currLoc-1];
                register pixel p2 =  src[currLoc];
                register pixel p3 =  src[currLoc+1];
                sum.red -= p1.red;
                sum.green -= p1.green;
                sum.blue -=  p1.blue;

                sum.red -=  p2.red;
                sum.green -=  p2.green;
                sum.blue -=  p2.blue;

                sum.red -=  p3.red;
                sum.green -=  p3.green;
                sum.blue -=  p3.blue;
            }
            //add the middle row to sum:
            register int currLoc = LOCATION(i,j,m);
            pixel p1 =  src[currLoc-1];
            sum.red -= p1.red;
            sum.green -= p1.green;
            sum.blue -= p1.blue;

            sum.red += 9 * src[currLoc].red;
            sum.green += 9 * src[currLoc].green;
            sum.blue += 9 * src[currLoc].blue;

            pixel p2 = src[currLoc+1];
            sum.red -= p2.red;
            sum.green -= p2.green;
            sum.blue -= p2.blue;

            dst[location].red = min(max(sum.red, 0), 255);
            dst[location].green =min(max(sum.green, 0), 255);
            dst[location].blue = min(max(sum.blue, 0), 255);
        }
    }
}

void smooth(int dim, pixel *src, pixel *dst, int kernelSize, int kernel[kernelSize][kernelSize], int kernelScale, bool filter) {}
void charsToPixels(Image *charsImg, pixel* pixels) {}

void pixelsToChars(pixel* pixels, Image *charsImg) {}

void copyPixels(pixel* src, pixel* dst) {}

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

    pixel* pixelsImg = malloc(IMG_SIZE);
    pixel* backupOrg = malloc(IMG_SIZE);

    //coping image:
    memcpy(pixelsImg, image->data, IMG_SIZE);
    //copyPixels(pixelsImg, backupOrg);
    memcpy(backupOrg, pixelsImg, IMG_SIZE);

    smoothBlurNoFilter(backupOrg, pixelsImg);

    memcpy(image->data, pixelsImg, IMG_SIZE);

    free(pixelsImg);
    free(backupOrg);
}
void doConvolutionBlurWithFilter(Image *image) {

    pixel* pixelsImg = malloc(IMG_SIZE);
    pixel* backupOrg = malloc(IMG_SIZE);

    //coping image:
    memcpy(pixelsImg, image->data, IMG_SIZE);
    //copyPixels(pixelsImg, backupOrg);
    memcpy(backupOrg, pixelsImg, IMG_SIZE);

    smoothBlurWithFilter(backupOrg, pixelsImg);

    memcpy(image->data, pixelsImg, IMG_SIZE);

    free(pixelsImg);
    free(backupOrg);
}

void doConvolutionSharp(Image *image) {
    pixel* pixelsImg = malloc(IMG_SIZE);
    pixel* backupOrg = malloc(IMG_SIZE);

    //coping image:
    memcpy(pixelsImg, image->data, IMG_SIZE);
    //copyPixels(pixelsImg, backupOrg);
    memcpy(backupOrg, pixelsImg, IMG_SIZE);

    smoothSharpNoFilter(backupOrg, pixelsImg);

    memcpy(image->data, pixelsImg, IMG_SIZE);

    free(pixelsImg);
    free(backupOrg);
}

void myfunction(Image *image, char* srcImgpName, char* blurRsltImgName, char* sharpRsltImgName,
                char* filteredBlurRsltImgName, char* filteredSharpRsltImgName, char flag) {
    pixel* pixelsImg = malloc(IMG_SIZE);
    pixel* backupOrg = malloc(IMG_SIZE);
    //coping image:
    memcpy(pixelsImg, image->data, IMG_SIZE);
    memcpy(backupOrg, pixelsImg, IMG_SIZE);
    if (flag == '1') {
        // blur image
        //doConvolutionBlurNoFilter(image);
        smoothBlurNoFilter(backupOrg, pixelsImg);

        memcpy(image->data, pixelsImg, IMG_SIZE);


        // write result image to file
        writeBMP(image, srcImgpName, blurRsltImgName);

        // sharpen the resulting image
        //doConvolutionSharp(image);
        //coping image:
        memcpy(pixelsImg, image->data, IMG_SIZE);
        memcpy(backupOrg, pixelsImg, IMG_SIZE);

        smoothSharpNoFilter(backupOrg, pixelsImg);

        memcpy(image->data, pixelsImg, IMG_SIZE);

        // write result image to file
        writeBMP(image, srcImgpName, sharpRsltImgName);
    } else {
        // apply extermum filtered kernel to blur image
        //doConvolutionBlurWithFilter(image);

        smoothBlurWithFilter(backupOrg, pixelsImg);

        memcpy(image->data, pixelsImg, IMG_SIZE);

        // write result image to file
        writeBMP(image, srcImgpName, filteredBlurRsltImgName);

        // sharpen the resulting image
        //doConvolutionSharp(image);
        memcpy(pixelsImg, image->data, IMG_SIZE);
        memcpy(backupOrg, pixelsImg, IMG_SIZE);

        smoothSharpNoFilter(backupOrg, pixelsImg);

        memcpy(image->data, pixelsImg, IMG_SIZE);

        // write result image to file
        writeBMP(image, srcImgpName, filteredSharpRsltImgName);
    }
    free(pixelsImg);
    free(backupOrg);
}

