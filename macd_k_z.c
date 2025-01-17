#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#define N 16       // ñ�W�i�Ϊ��ס]16 �Ӵ����^
#define MAX_ITER 1000 // ���N����
#define alpha 0.1   // �ڤɧE���o�i�����u���Y��
#define Tc 1.0     // �����g��
#define p 0.2     // ���i�Ҿ��v
#define oversampling 2*N     // ���˲v�]�C�Ӵ��������q�ơ^

// �b�����i�� psi(t)
double psi(double t) {
    if (fabs(t) < Tc / 2.0) {
        return ccos(M_PI * t / Tc);
    }
    return 0.0;
}
// ��l��ñ�W�i�� c_k^(0)(t)
void initialize_signature_waveform(double complex *k_z, double *s_k, int user_index, int length) {
    int extended_length = 2 * length; // 2N
    double step = Tc / oversampling; // �ɶ��B��
    int total_samples = N * oversampling;          // ���( extended_length * oversampling )

	// ��l���A�G�H����� +1 �� -1
    s_k[0] = (rand() % 2 ? 1 : -1);
    // �ͦ����i�ҧǦC
    for (int i = 1; i < extended_length; i++) {
        double rand_prob = (double)rand() / RAND_MAX; // �H�����v [0, 1)
        if (s_k[i - 1] == 1) {
            s_k[i] = (rand_prob < p) ? -1 : 1;   // +1 -> -1 �����v p
        } else {
            s_k[i] = (rand_prob < p) ? 1 : -1;   // -1 -> +1 �����v p
        }
    }
	// c_k^(0)(t)
    for (int i = 0; i < total_samples; i++) {
        double t = i * step; // ��e�ɶ�
        double complex value = 0.0;

        // �ھڤ��� (6) �p��i��
        for (int n = 0; n < extended_length; n++) {
            double phase_shift = pow(-1, user_index) * n * M_PI / 2.0; // (-1)^k n ���ۦ찾��
            double complex phase = cexp(I * phase_shift); // j^((-1)^k n)
            value += phase * s_k[n] * psi(t - n * Tc / 2.0);
        }

        k_z[i] = value; // �O�s�i��
    }
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ���T���W��
void normalize(double complex *k_z, int length) {
	double max ;
	double k_z_amp[length] ;
	double k_z_arg[length];
    for (int i = 0; i < length; i++) {
        k_z_amp[i] = sqrt(creal(k_z[i]) * creal(k_z[i]) + cimag(k_z[i]) * cimag(k_z[i]));     //�s�� c_k���T
		k_z_arg[i] = atan2( creal(k_z[i]) , cimag(k_z[i] ));     //�s�� c_k�ۦ�
    }
	max = k_z_amp[0];
	for (int i = 0; i < length; i++){
		if(k_z_amp[i] >= max )
			max = k_z_amp[i];
	}
	for (int i = 0; i < length; i++)
		k_z_amp[i] = k_z_amp[i] / max;
	for (int i = 0; i < length; i++)
		k_z[i] = k_z_amp[i] * cexp(I * k_z_arg[i]);
}
// �ڤɧE���o�i���ͦ�
void generate_rrc_filter(double *filter, int length ) {
	double Ts = N*Tc ;
    for (int i = 0; i < length; i++) {
        double t = (i - length / 2) * Tc; // ���ߤƮɶ�
        if (fabs(t) < 1e-6) {
            filter[i] = 1.0; // �קK t = 0 ���������s
        } else if (fabs(fabs(t) - 1.0 / (2.0 * alpha)) < 1e-6) {
            filter[i] = (alpha / sqrt(2)) * ((1 + 2 / M_PI) * sin(M_PI / (4 * alpha)) +
                                             (1 - 2 / M_PI) * cos(M_PI / (4 * alpha)));
        } else {
            filter[i] = (sin(M_PI * t * (1 - alpha)) +
                         4 * alpha * t * cos(M_PI * t * (1 + alpha))) /
                        (M_PI * t * (1 - pow(4 * alpha * t, 2)));
        }
    }
}


// �D�y�{�A�ëO�s�S�w���N���ƪ��ƾ�
void generate_ce_waveform(double complex *k_z, int length, int save_iters[], int num_saves, const char *output_file ) {
    double filter[length];
    generate_rrc_filter(filter, length ); // �ͦ��ڤɧE���o�i��

    FILE *file = fopen(output_file, "w");
    if (!file) {
        perror("Failed to open output file");
        exit(EXIT_FAILURE);
    }
    fprintf(file, "Iteration,Time,Phase\n"); // �g�J���Y

    int save_index = 0;
    for (int iter = 0; iter <= MAX_ITER; iter++) {
        // �o�i�]²�ƪ��`�����n�^
        for (int i = 0; i < length; i++) {
            double complex temp = 0.0;
            for (int j = 0; j < length; j++) {
                temp += k_z[j] * filter[(i - j + length) % length];
            }
            k_z[i] = temp;
        }

        normalize(k_z, length); // ���W�ƴT��

        // �O�s�ƾڡ]�Ȧb���w���N���ơ^
        if (save_index < num_saves && iter == save_iters[save_index]) {
            printf("This is length : %d \n",length);
            for (int i = 0; i < length; i++) {
                printf("This is ( (float)(i * Tc) )/oversampling : %f \n",( (float)(i * Tc) )/oversampling);
                fprintf(file, "%d,%f,%f,%f\n", iter, ( (float)(i * Tc) )/oversampling, cimag(k_z[i]), creal(k_z[i]) ) ; // �p��ۦ�
            }
            save_index++;
        }
    }

    fclose(file);
}

int main() {
    int total_samples = N * oversampling; // �`�����I��       // ���( 2 * N * oversampling )
    double s_k[2 * N]; // ��l�����ǦC
    double complex k_z[total_samples]; // ñ�W�i��
    int save_iters[] = {0, 1, 10, 100, 1000}; // �n�O�s�����N����
    int num_saves = sizeof(save_iters) / sizeof(save_iters[0]);

    srand((unsigned int)time(NULL)); // ��l���H���ƺؤl

    initialize_signature_waveform(k_z,s_k, 1 ,N) ;// ��l�ƪi�Ρ]���i�ҹL�{�^
    generate_ce_waveform(k_z, total_samples, save_iters, num_saves,"ce_waveform_data.csv"); // �O�s�ƾ�

    printf("Data saved to ce_waveform_data.csv\n");
    return 0;
}
