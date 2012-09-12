#include <math.h>

void rotation_matrix_c(double *axis, double theta, double *mat);
void get_closer_rotation_matrix_c(double *TH, double *point, double *M, double *F, double *out_rot_mat);
void ccd_c(double *moving, double *fixed, int len_moving, int iterations);
