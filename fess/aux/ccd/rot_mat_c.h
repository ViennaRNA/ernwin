#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void rotation_matrix_c(double *axis, double theta, double *mat);
void get_closer_rotation_matrix_c(double *TH, double *point, double *M, double *F, double *out_rot_mat);
void ccd_c(double *moving, int len_moving, double *fixed, long *points, int len_points, int moving_end, int iterations, double rmsd_threshold);
