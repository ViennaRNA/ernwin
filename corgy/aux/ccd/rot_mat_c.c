#include "rot_mat_c.h"

double c_dot(double *v1, double *v2) {
    return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

void rotation_matrix_c(double *axis1, double theta, double *mat1) {
    double *axis = (double *)axis1;
    double *mat = (double *)mat1;

    double x = sqrt(axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]);
    double a = cos(theta / 2.0);
    double b = -(axis[0] / x) * sin(theta / 2.0);
    double c = -(axis[1] / x) * sin(theta / 2.0);
    double d = -(axis[2] / x) * sin(theta / 2.0);

    mat[0] = a*a + b*b - c*c - d*d;
    mat[1] = 2 * (b*c - a*d);
    mat[2] = 2 * (b*d + a*c);

    mat[3*1 + 0] = 2*(b*c+a*d);
    mat[3*1 + 1] = a*a+c*c-b*b-d*d;
    mat[3*1 + 2] = 2*(c*d-a*b);

    mat[3*2 + 0] = 2*(b*d-a*c);
    mat[3*2 + 1] = 2*(c*d+a*b);
    mat[3*2 + 2] = a*a+d*d-b*b-c*c;
}


void get_closer_rotation_matrix_c(double *TH, double *point, double *M, double *F, double *out_rot_mat) {
    int i, j;

    double sum1 = 0;
    double TH_norm[3];
    double d1[3];

    double R[9], O[9];

    //TH_norm = normalize(TH)
    for (i = 0; i < 3; i++) 
        sum1 += TH[i] * TH[i];
    sum1 = sqrt(sum1);

    //TH_arr = array([TH_norm, TH_norm, TH_norm])
    for (i = 0; i < 3; i++)
        TH_norm[i] = TH[i] / sum1;

    for (i = 0; i < 3; i++) {
        sum1 = 0;

        for (j = 0; j < 3; j++) {
            sum1 += (M[i * 3 + j] - point[j]) * TH_norm[j];
        }
        
        d1[i] = sum1;
    }

    for ( i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            O[i*3 + j] = TH_norm[j] * d1[i] + point[j];
            R[i*3 + j] = M[i*3 + j] - O[i*3 + j];

        }
    }

    //double fs_sum = 0;
    //double fr_sum = 0;

    double temp_fs_sum = 0;
    double temp_fr_sum = 0;

    double n3 = 0, d3 = 0, a, r_row_sums[3], s_row_sums[3];
    double S[9], F1[9];

    for (i = 0; i < 3; i++) {
        // S = cross(R, TH)
        S[3*i + 0] = R[3*i + 1] * TH[2] - R[3*i + 2] * TH[1];
        S[3*i + 1] = R[3*i + 2] * TH[0] - R[3*i + 0] * TH[2];
        S[3*i + 2] = R[3*i + 0] * TH[1] - R[3*i + 1] * TH[0];

        s_row_sums[i] = 0;
        r_row_sums[i] = 0;

        //s_row_sums = np.sum(np.abs(S) ** 2, axis=-1) ** (1./2.)
        //r_row_sums = np.sum(np.abs(R) ** 2, axis=-1) ** (1./2.)
        //F = F - O

        temp_fs_sum = 0;
        temp_fr_sum = 0;

        for (j = 0; j < 3; j++)  {
            s_row_sums[i] += S[3*i + j] * S[3*i + j];
            r_row_sums[i] += R[3*i + j] * R[3*i + j];

            F1[3*i + j] = F[3*i + j] - O[3*i + j];
        }

        s_row_sums[i] = sqrt(s_row_sums[i]);
        r_row_sums[i] = sqrt(r_row_sums[i]);

        //np.sum(F * S, axis=-1)
        //np.sum(F * R, axis=-1)
        
        for (j = 0; j < 3; j++) {
            S[3*i + j] /= s_row_sums[i];
            R[3*i + j] /= r_row_sums[i];

            temp_fs_sum += F1[3*i + j] * S[3*i + j];
            temp_fr_sum += F1[3*i + j] * R[3*i + j];
        }

        n3 += r_row_sums[i] * temp_fs_sum;
        d3 += r_row_sums[i] * temp_fr_sum;
    }

    a = atan2(n3, d3);

    rotation_matrix_c(TH, a, out_rot_mat);
}

void mat_times_vec(double *mat, double *vec, double *out_vec, double *offset) {
    double temp[3];

    int i, j;

    for (i = 0; i < 3; i++) {
        temp[i] = 0;

        for (j = 0; j < 3; j++) 
            temp[i] += mat[i*3 + j] * (vec[j] - offset[j]);

        temp[i] += offset[i];
    }

    memcpy(out_vec, temp, 3 * sizeof(double));
}

double calc_rmsd_c(double *v1, double *v2, int n)
{
    double sum = 0;
    int i, j;

    for (i = 0; i < n; i++)
        for (j = 0; j < 3; j++)
            sum += (v1[j] - v2[j]) * (v1[j] - v2[j]);

    sum /= n;
    return sqrt(sum);
}

void print_rows(double *m, int n)
{
    int i, j;

    for (i = 0; i < n; i++) {
        printf("[");
        for (j = 0; j < 3; j++)
            printf("%f ", m[i*3 + j]);
        printf("]\n");
    }
}

void ccd_c(double *moving, int len_moving, double *fixed, long *points, int len_points, int moving_end, int iterations, double rmsd_threshold) {
    double rot_mat[9], rmsd;
    int k, i, prev_i, j;
    //print_rows(&moving[(len_moving-3) * 3], 3);

    for (k = 0; k < iterations; k++) {
        for (i = 0; i < len_points; i++) {
            double TH[3];
            int l = points[i];

            TH[0] = moving[l*3 + 0] - moving[(l-1) * 3 + 0];
            TH[1] = moving[l*3 + 1] - moving[(l-1) * 3 + 1];
            TH[2] = moving[l*3 + 2] - moving[(l-1) * 3 + 2];

            get_closer_rotation_matrix_c(TH, &moving[(l-1) * 3], &moving[(moving_end) * 3], fixed, rot_mat);

            for (j = l+1; j < len_moving; j++) 
                mat_times_vec(rot_mat, &moving[j*3], &moving[j*3], &moving[l*3]);

        }

        rmsd = calc_rmsd_c( &moving[(moving_end - 3) * 3], fixed, 3);
        
        if (rmsd < rmsd_threshold)
            return;
    }
}
