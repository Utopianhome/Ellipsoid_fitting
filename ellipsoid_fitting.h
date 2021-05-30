#ifndef __ELLIPSOID_FITTING_H__
#define __ELLIPSOID_FITTING_H__
/* SPDX-License-Identifier: GPL-2.0 */
/* Copyright (C) 2021, Jia Hu
 * ellipsoid_fitting.h 
 * /
 
 /* main API:
  * int ellipsoid_fitting(struct Matrix data, struct Matrix D, struct ellipsoid character);
  * @ return value --- int: indicate error, 0 represent no bugs happening.
  * @ struct Matrix data: format { entry: [x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4 ...], row = N, col = 3}
  * @ struct Matrix D:    format { entry: D = [x1^2, y1^2, z1^2, 2*x1*y1, 2*x1*z1, 2*y1*z1, 2x1, 2y1, 2z1, ....], row = N, col = 9}
  * @ struct ellipsoid character: format {offset: double[3], gain: double[3]}
  * 
  * attention: allocate the space for struct Matrix data, struct Matrix D and struct ellipsoid character outside the function.
  * recommend reserve enough space for matrix
  */

#include <stdio.h>
#include <math.h>

struct Matrix {
    double* entry;
    int row;
    int col;
    int reserve_space;
};

struct ellipsoid {
    double offset[3];
    double gain[3];
};

int ellipsoid_fitting(struct Matrix data, struct Matrix D, struct ellipsoid *character);
int ellipsoid_fitting_append(struct Matrix *D, struct ellipsoid *character, double appenddata[3]);
void printMatrix(const struct Matrix T);
void printEllipsoid(const struct ellipsoid value);

/* print function*/
void printMatrix(const struct Matrix T) {
    for(int i = 0; i < T.row; i++) {
        for(int j = 0; j < T.col; j++) {
            printf("%.8f\t", T.entry[i*T.col+j]);
        }
        printf("\n");
    }
}

void printEllipsoid(const struct ellipsoid value) {
    printf("offset:\t");
    for(int i = 0; i < 3; i++) {
        printf("%.6f\t", value.offset[i]);
    }
    printf("\n");
    printf("gain:\t");
    for(int i = 0; i < 3; i++) {
        printf("%.6f\t", value.gain[i]);
    }
    printf("\n");
}
/* construct D matrix D = [x1^2, y1^2, z1^2, 2*x1*y1, 2*x1*z1, 2*y1*z1, 2x1, 2y1, 2z1, ....] */
int constructD(struct Matrix result, struct Matrix data) {
    if(result.col != 9 && result.row != data.row) {
        return -1;
    }
    
    for(int i = 0; i < result.row; i++) {
        result.entry[i*result.col + 0] = data.entry[i*data.col + 0] * data.entry[i*data.col + 0];    //xx
        result.entry[i*result.col + 1] = data.entry[i*data.col + 1] * data.entry[i*data.col + 1];    //yy
        result.entry[i*result.col + 2] = data.entry[i*data.col + 2] * data.entry[i*data.col + 2];    //zz
        result.entry[i*result.col + 3] = 2 * data.entry[i*data.col + 0] * data.entry[i*data.col + 1];//2xy
        result.entry[i*result.col + 4] = 2 * data.entry[i*data.col + 0] * data.entry[i*data.col + 2];//2xz
        result.entry[i*result.col + 5] = 2 * data.entry[i*data.col + 1] * data.entry[i*data.col + 2];//2yz
        result.entry[i*result.col + 6] = 2 * data.entry[i*data.col + 0] ; // 2x
        result.entry[i*result.col + 7] = 2 * data.entry[i*data.col + 1] ; // 2y
        result.entry[i*result.col + 8] = 2 * data.entry[i*data.col + 2] ; // 2z
    }
    return 0;
}

/* construct A = inv(D)*D */
int transposeMultiply(struct Matrix result, struct Matrix T) { // inv(T)*T
    if(result.col != T.col && result.row != T.col) {
        return -1;
    }
    double tmp_sum = 0;
    for(int i = 0; i < result.row; i++) {
        for(int j = 0; j < result.col; j++) {
            for(int k = 0; k < T.row; k++) {
                tmp_sum += T.entry[i + k*T.col] * T.entry[j + k*T.col];
            }
            result.entry[i*result.col + j] = tmp_sum;
            tmp_sum = 0;
        }

    }
    return 0;
}

/* construct inv(D)*Y:  inv(D)*Dv = inv(D)*Y */
int transposeMultiplyByOnes(struct Matrix result, struct Matrix A) {
    if(result.row != A.col && result.col != 1 ) {
        return -1;
    }
    double tmp_sum = 0;
    for(int i = 0; i < A.col; i++) {
        for(int j = 0; j < A.row; j++) {
            tmp_sum += A.entry[i + j*A.col];
        }
        result.entry[i] = tmp_sum;
        tmp_sum = 0;
    }
    return 0;
}

/* LUP Decompostion A */
int LUPDecomposition(struct Matrix LU, struct Matrix P, struct Matrix A) {
    if(P.col != A.row && P.row != LU.row && LU.col != A.col) {
        return -1;
    }
    int row = A.row;
    int col = A.col;
    // assign the A to LU
    for(int i = 0; i < row; i++) {
        for(int j = 0; j < col; j++) {
            LU.entry[i*col+j] = A.entry[i*col+j];
            P.entry[i*col+j] = 0;
        }
    }
    for(int i = 0; i < row; i++) {
        P.entry[i*col + i] = 1;
    }

    //
    double max_value = 0;
    double tmp_value = 0;
    double current_value = 0;
    int max_index = -1;
    double beta = 0;

    for(int i = 0; i < row-1; i++) {
        max_index = i;
        for(int ir = i; ir < row; ir++) { // find the abs maximum
            current_value = LU.entry[ir*col+i];
            current_value = current_value < 0 ? -current_value : current_value;
            if(ir == i) {
                max_value = current_value;
                continue;
            }
            if(current_value > max_value) {
                max_value = current_value;
                max_index = ir;
            }
        }
        if(max_index != i) {
            P.entry[i*col + max_index] = 1;
            P.entry[i*col + i] = 0;
            P.entry[max_index*col+max_index] = 0;
            P.entry[max_index*col+i] = 1;
            for(int jj = 0; jj < col; jj++) { // exchange the max line with the head line
                tmp_value = LU.entry[i*col+jj];
                LU.entry[i*col+jj] = LU.entry[max_index*col+jj];
                LU.entry[max_index*col+jj] = tmp_value;
            }
        }
        for(int istart = i+1; istart < row; istart++) {
            beta = 0;
            if(LU.entry[i*col + i] == 0) {
                break;
            }
            beta = LU.entry[istart*col + i] / LU.entry[i*col + i];
                
            for(int js = i+1; js < col; js++) {
            LU.entry[istart*col+js] = LU.entry[istart*col+js] - beta * LU.entry[i*col+js];
            }
            LU.entry[istart*col + i] = beta;
        }
    }
    return 0;
}

/* solve equation PAv=PY, LUv = PY */
int solveLUP(struct Matrix result, struct Matrix LU, struct Matrix P, struct Matrix M) {
    // LUx=M , solution x is required. PA = LU, PAx = PM, LUx = PM, Ly = PM, Ux = y
    if(P.col != M.row && result.col != M.col && result.row != LU.col) {
        return -1;
    }
    //PA
    double sum = 0;
    for(int i = 0; i < M.row; i++) {
        sum = 0;
        for(int j = 0; j < P.col; j++) {
            sum += P.entry[i*P.col+j]*M.entry[j];
        }
        result.entry[i] = sum;
    }

    for(int i = 0; i < M.row; i++) {
        M.entry[i] = result.entry[i];
    }
    double tmp = 0; // calculate the y
    for(int i = 0; i < M.row; i++) {
        for(int j = 0; j < i; j++) {
            tmp += LU.entry[i*LU.col + j] * result.entry[j];
        }
        result.entry[i] = M.entry[i] - tmp;
        tmp = 0;
    }
    // calculate the x
    for(int i = result.row - 1; i >= 0; i--) {
        for(int j = i + 1; j < LU.col; j++) {
            tmp += LU.entry[i*LU.col + j] * result.entry[j];
        }
        if(LU.entry[i*LU.col + i] == 0) {
            return -2;
        }
        result.entry[i] = (result.entry[i] - tmp) / LU.entry[i*LU.col + i];
        tmp = 0;
    }
    return 0;
}

/* from the result calculate the offset and gain */
int getOffsetGain(struct Matrix v, struct ellipsoid *result) {  // result[9,1] A[4, 4]
    double g = 1 + (v.entry[6]*v.entry[6] / v.entry[0]
                + v.entry[7]*v.entry[7] / v.entry[1]
                + v.entry[8]*v.entry[8] / v.entry[2]);
    for(int i = 0; i < 3; i++) {
        result->offset[i] = -v.entry[6+i] / v.entry[0 + i];
        result->gain[i] = g / v.entry[i];
        if(result->gain[i] < 0) {
            return -1;
        }
        result->gain[i] = sqrt(result->gain[i]);
    }
    return 0;
}
/* main API */
int ellipsoid_fitting(struct Matrix data, struct Matrix D, struct ellipsoid *result) 
{
    double dataA[9][9];
    double dataY[9][1];
    double dataLU[9][9];
    double datav[9][1];
    double dataP[9][9];
    
    struct Matrix A = {.entry = &dataA[0][0], 9, 9};
    struct Matrix Y = {.entry = &dataY[0][0], 9, 1};
    struct Matrix LU = {.entry = &dataLU[0][0], 9, 9};
    struct Matrix P = {.entry = &dataP[0][0], 9, 9};
    struct Matrix v = {.entry = &datav[0][0], 9, 1};
    if(constructD(D, data) != 0) {
        //printf("constructD error\n");
        return -1;
    }
    if(transposeMultiply(A, D) != 0) {
        //printf("transposeMultiply error\n");
        return -1;
    }
    if(transposeMultiplyByOnes(Y, D) != 0) {
        //printf("transposeMultiplyByOnes error\n");
        return -1;
    } 
    if(LUPDecomposition(LU, P, A) != 0) {
        //printf("LUPDecompostion error\n");
        return -1;
    }
    if(solveLUP(v, LU, P, Y) != 0) {
        //printf("solveLUP error, code: %d\n", solveLUP(v, LU, P, Y));
        return -1;
    }
    if(getOffsetGain(v, result) != 0) {
        return -1;
    }
}

int appenddata(struct Matrix *D, double* data) {
    if((D->reserve_space - 9 * D->row) < 9) {
        printf("no enough space\n");
        return -1;
    }
    D->reserve_space = D->reserve_space - 9;
    int endrow = D->row;
    D->entry[9*endrow + 0] = data[0] * data[0];    //xx
    D->entry[9*endrow + 1] = data[1] * data[1];    //yy
    D->entry[9*endrow + 2] = data[2] * data[2];    //zz
    D->entry[9*endrow + 3] = 2 * data[0] * data[1];//2xy
    D->entry[9*endrow + 4] = 2 * data[0] * data[2];//2xz
    D->entry[9*endrow + 5] = 2 * data[1] * data[2];//2yz
    D->entry[9*endrow + 6] = 2 * data[0]; // 2x
    D->entry[9*endrow + 7] = 2 * data[1]; // 2y
    D->entry[9*endrow + 8] = 2 * data[2]; // 2z
    D->row++;
    return 0;
}

int ellipsoid_fitting_append(struct Matrix *D, struct ellipsoid *result, double *data) // we need chang the D data
{
    double dataA[9][9];
    double dataY[9][1];
    double dataLU[9][9];
    double datav[9][1];
    double dataP[9][9];
    
    struct Matrix A = {.entry = &dataA[0][0], 9, 9};
    struct Matrix Y = {.entry = &dataY[0][0], 9, 1};
    struct Matrix LU = {.entry = &dataLU[0][0], 9, 9};
    struct Matrix P = {.entry = &dataP[0][0], 9, 9};
    struct Matrix v = {.entry = &datav[0][0], 9, 1};
    if(appenddata(D, data) != 0) {
        return -1;
    }
    if(transposeMultiply(A, *D) != 0) {
        //printf("transposeMultiply error\n");
        return -1;
    }
    if(transposeMultiplyByOnes(Y, *D) != 0) {
        //printf("transposeMultiplyByOnes error\n");
        return -1;
    } 
    if(LUPDecomposition(LU, P, A) != 0) {
        //printf("LUPDecompostion error\n");
        return -1;
    }
    if(solveLUP(v, LU, P, Y) != 0) {
        //printf("solveLUP error, code: %d\n", solveLUP(v, LU, P, Y));
        return -1;
    }
    if(getOffsetGain(v, result) != 0) {
        return -1;
    }
}

#endif





