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

void initialize_pixel_sum(pixel_sum *sum) {}

static void assign_sum_to_pixel(pixel *current_pixel, pixel_sum sum, int kernelScale) {}

static void sum_pixels_by_weight(pixel_sum *sum, pixel p, int weight) {}

static pixel applyKernel(int dim, int i, int j, pixel *src, int kernelSize, int kernel[kernelSize][kernelSize], int kernelScale, bool filter) {}



/*
* The smooth function that blur the image if the command is 1- no filter.
 * Optimization explanations: (1) use "register" before int/pixel to tell the computer to put this member in register
 * and not to call it from memory every time.
 * (2) few less arguments to this function because the others are pixed (not changing trought all project in this state-
 * blur with no filter).
 * (3) not using help function- to save call in assembly.
*/
void smoothBlurNoFilter(pixel *src, pixel *dst) {
    register int dimLoop =m-1;
    register int i,j;
    // check every pixel in separate:
    for (i = 1 ; i < dimLoop; i++) {
        for (j =  1 ; j < dimLoop; j++) {
            register int location = LOCATION(i,j,m); // use define to do this calculte instead of every location.
            register pixel_sum sum;
            sum.red = 0, sum.green = 0, sum.blue = 0;
            register int k;
            // sum all pixels around current pixel:
            // for a row in src- one loop for rows- the columns checks inside (because there are only 3):
            for (k = i-1; k <= i+1; ++k) { // for a row in src
                // use define to do this calculte instead of every location:
                register int currLoc = LOCATION(k,j,m);
                // go local to these pixels instead of reaching from memory every time:
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
            //there is no need to fin min max because the result will always be that- every pixel is between 0 and 255.
            dst[location].red = sum.red / 9;
            dst[location].green = sum.green / 9;
            dst[location].blue = sum.blue / 9;
        }
    }
}

/*
* The smooth function that blur the image if the command is not 1- with filter.
 * Optimization explanations: (1) use "register" before int/pixel to tell the computer to put this member in register
 * and not to call it from memory every time.
 * (2) few less arguments to this function because the others are pixed (not changing trought all project in this state-
 * blur with no filter).
 * (3) not using help function- to save call in assembly.
*/
void smoothBlurWithFilter(pixel *src, pixel *dst) {
    register int dimLoop =m-1;
    register int i,j;
    // check every pixel in separate:
    for (i = 1 ; i < dimLoop; i++) {
        for (j =  1 ; j < dimLoop; j++) {
            register int location = LOCATION(i, j, m);  // use define to do this calculte instead of every location.
            pixel_sum sum;
            sum.red = 0, sum.green = 0, sum.blue = 0;
            register int maxi = 776, mini = -1;
            register int sumTotal, minRow, maxRow, minCol, maxCol,k;
            // sum all pixels around current pixel:
            // for a row in src- one loop for rows- the columns checks inside (because there are only 3):
            for (k = i-1; k <= i+1; ++k) {
                // use define to do this calculte instead of every location:
                register int currLoc = LOCATION(k,j,m);
                // go local to these pixels instead of reaching from memory every time:
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
            // instead of duplicate in (-1)- using minus ---> (-1)*x = -x:
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

/*
* The smooth function that sharp the image (filtered or not).
 * Optimization explanations: (1) use "register" before int/pixel to tell the computer to put this member in register
 * and not to call it from memory every time.
 * (2) few less arguments to this function because the others are pixed (not changing trought all project in this state-
 * blur with no filter).
 * (3) not using help function- to save call in assembly.
*/
void smoothSharpNoFilter(pixel *src, pixel *dst) {
    register int dimLoop =m-1; /// maybe in a different way
    register int i,j,k;
    // check every pixel in separate:
    for (i = 1 ; i < dimLoop; i++) {
        for (j =  1 ; j < dimLoop; j++) {
            register int location = LOCATION(i, j, m); // use define to do this calculte instead of every location.
            register pixel_sum sum;
            sum.red = 0, sum.green = 0, sum.blue = 0;
            // sum all pixels around current pixel:
            // for a row in src- one loop for 2 rows(the first and the third)- the columns checks inside
            // (because there are only 3)- the middle row in the matrix is checked after this loop:
            for (k = i-1; k <= i+1; k+=2) { // for a row in src-  2 rows- first and last (1,3)
                // use define to do this calculte instead of every location:
                register int currLoc = LOCATION(k,j,m);
                // go local to these pixels instead of reaching from memory every time:
                register pixel p1 =  src[currLoc-1];
                register pixel p2 =  src[currLoc];
                register pixel p3 =  src[currLoc+1];
                // instead of duplicate in (-1)- using minus from zero ---> (-1)*x = -x:
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
            register int currLoc = LOCATION(i,j,m); // use define to do this calculte instead of every location.

            //calculate the middle row in the matrix (because it's different from the other 2 row (in the above loop)):
            pixel p1 =  src[currLoc-1]; // go local to these pixels instead of reaching from memory every time.
            // instead of duplicate in (-1)- using minus from zero ---> (-1)*x = -x:
            sum.red -= p1.red;
            sum.green -= p1.green;
            sum.blue -= p1.blue;

            sum.red += 9 * src[currLoc].red;
            sum.green += 9 * src[currLoc].green;
            sum.blue += 9 * src[currLoc].blue;

            pixel p2 = src[currLoc+1]; // go local to these pixels instead of reaching from memory every time.
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

void doConvolution(Image *image, int kernelSize, int kernel[kernelSize][kernelSize], int kernelScale, bool filter) {}

/**
 * we dont need this function- debug:
 */
void doConvolutionBlurNoFilter(Image *image) {

//    pixel* pixelsImg = malloc(IMG_SIZE);
//    pixel* backupOrg = malloc(IMG_SIZE);
//
//    //coping image:
//    memcpy(pixelsImg, image->data, IMG_SIZE);
//    //copyPixels(pixelsImg, backupOrg);
//    memcpy(backupOrg, pixelsImg, IMG_SIZE);
//
//    smoothBlurNoFilter(backupOrg, pixelsImg);
//
//    memcpy(image->data, pixelsImg, IMG_SIZE);
//
//    free(pixelsImg);
//    free(backupOrg);
}
/**
 * we dont need this function- debug:
 */
void doConvolutionBlurWithFilter(Image *image) {

//    pixel* pixelsImg = malloc(IMG_SIZE);
//    pixel* backupOrg = malloc(IMG_SIZE);
//
//    //coping image:
//    memcpy(pixelsImg, image->data, IMG_SIZE);
//    //copyPixels(pixelsImg, backupOrg);
//    memcpy(backupOrg, pixelsImg, IMG_SIZE);
//
//    smoothBlurWithFilter(backupOrg, pixelsImg);
//
//    memcpy(image->data, pixelsImg, IMG_SIZE);
//
//    free(pixelsImg);
//    free(backupOrg);
}
/**
 * we dont need this function- debug:
 */
void doConvolutionSharp(Image *image) {
//    pixel* pixelsImg = malloc(IMG_SIZE);
//    pixel* backupOrg = malloc(IMG_SIZE);
//
//    //coping image:
//    memcpy(pixelsImg, image->data, IMG_SIZE);
//    //copyPixels(pixelsImg, backupOrg);
//    memcpy(backupOrg, pixelsImg, IMG_SIZE);
//
//    smoothSharpNoFilter(backupOrg, pixelsImg);
//
//    memcpy(image->data, pixelsImg, IMG_SIZE);
//
//    free(pixelsImg);
//    free(backupOrg);
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
        smoothBlurNoFilter(backupOrg, pixelsImg);

        memcpy(image->data, pixelsImg, IMG_SIZE);

        // write result image to file
        writeBMP(image, srcImgpName, blurRsltImgName);

        //coping image:
        memcpy(pixelsImg, image->data, IMG_SIZE);
        memcpy(backupOrg, pixelsImg, IMG_SIZE);

        // sharpen the resulting image
        smoothSharpNoFilter(backupOrg, pixelsImg);

        memcpy(image->data, pixelsImg, IMG_SIZE);

        // write result image to file
        writeBMP(image, srcImgpName, sharpRsltImgName);

    } else {
        //blur image
        smoothBlurWithFilter(backupOrg, pixelsImg);

        memcpy(image->data, pixelsImg, IMG_SIZE);

        // write result image to file
        writeBMP(image, srcImgpName, filteredBlurRsltImgName);

        //coping image:
        memcpy(pixelsImg, image->data, IMG_SIZE);
        memcpy(backupOrg, pixelsImg, IMG_SIZE);

        // sharpen the resulting image
        smoothSharpNoFilter(backupOrg, pixelsImg);

        memcpy(image->data, pixelsImg, IMG_SIZE);

        // write result image to file
        writeBMP(image, srcImgpName, filteredSharpRsltImgName);
    }
    free(pixelsImg);
    free(backupOrg);
}

