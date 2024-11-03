#include <omp.h>
#include <unistd.h>

#include <cmath>
#include <vector>

#include "fitsio.h"
#include "stdio.h"

using namespace std;

typedef struct elipse {
    int o_x, o_y;
    double alpha, theta;
    int beta;
} elipse;

vector<int> get_borders(double *myimage, long *naxes, int hebras1);
vector<elipse> get_elipses(double *myimage, long *naxes, int hebras1, int hebras2, int numBetas, double relativeVotes, int alphaMin, vector<int> borders);

int main(int argc, char *argv[]) {
    int status;
    fitsfile *fptr;
    long fpixel = 1, naxis = 2, exposure;
    long naxes[2];

    int numBetas, alphaMin, hebras1, hebras2;
    double relativeVotes;
    char *input_file;

    int opt;
    while ((opt = getopt(argc, argv, "i:a:r:b:u:d:")) != -1) {
        switch (opt) {
            case 'i':
                input_file = optarg;
                break;
            case 'a':
                alphaMin = atoi(optarg);
                break;
            case 'r':
                relativeVotes = atof(optarg);
                break;
            case 'b':
                numBetas = atoi(optarg);
                break;
            case 'u':
                hebras1 = atoi(optarg);
                break;
            case 'd':
                hebras2 = atoi(optarg);
                break;
            default:
                fprintf(stderr, "Usage: [-i file_name] [-b numBetas] [-v relativeVotes] [-a alphaMin] [-h hebras1] [-H hebras2] <input_file>\n");
                exit(EXIT_FAILURE);
        }
    }

    status = 0; /* initialize status before calling fitsio routines */

    fits_open_file(&fptr, input_file, READONLY, &status);
    fits_get_img_size(fptr, 2, naxes, &status);

    double *myimage = (double *)malloc(naxes[0] * naxes[1] * sizeof(double));
    fits_read_img(fptr, TDOUBLE, fpixel, naxes[0] * naxes[1], NULL, myimage, NULL, &status);

    vector<int> borders = vector<int>();

    // 1 thread
    double start_1 = omp_get_wtime();
    borders = get_borders(myimage, naxes, 1);
    double start_h1 = omp_get_wtime();
    get_elipses(myimage, naxes, 1, 1, numBetas, relativeVotes, alphaMin, borders);
    double stop_h1 = omp_get_wtime();
    double stop_1 = omp_get_wtime();

    double duration_h1 = stop_h1 - start_h1;
    double duration_1 = stop_1 - start_1;

    // N threads
    double start_n = omp_get_wtime();
    borders = get_borders(myimage, naxes, hebras1);
    double start_hn = omp_get_wtime();
    vector<elipse> newElipses = get_elipses(myimage, naxes, hebras1, hebras2, numBetas, relativeVotes, alphaMin, borders);
    double stop_hn = omp_get_wtime();
    double stop_n = omp_get_wtime();

    double duration_hn = stop_hn - start_hn;
    double duration_n = stop_n - start_n;

    fits_close_file(fptr, &status);

    for (elipse e : newElipses) {
        printf("%d, %d, %f, %d, %f\n", e.o_x, e.o_y, e.alpha, e.beta, e.theta);
    }

    printf("%.2f\n", duration_h1);
    printf("%.2f\n", duration_hn);

    // porcentaje serial
    double porcentaje = 25.0 / 200.0;
    printf("%.2f%%\n", porcentaje * 100);

    double speedup_total = duration_1 / duration_n;
    printf("%.2f\n", speedup_total);

    printf("%.2f\n", duration_1);
    printf("%.2f\n", duration_n);

    double speedup_hough = duration_h1 / duration_hn;
    printf("%.2f\n", speedup_hough);

    return 0;
}

vector<int> get_borders(double *myimage, long *naxes, int hebras1) {
    vector<int> borders = vector<int>();
#pragma omp parallel num_threads(hebras1)
    {
        vector<int> local_borders = vector<int>();
#pragma omp for
        for (int i = 0; i < naxes[0] * naxes[1]; i++) {
            if (myimage[i] > 0.0) {
                local_borders.insert(local_borders.end(), i);
            }
        }

#pragma omp critical
        {
            for (int i = 0; i < local_borders.size(); i++) {
                borders.insert(borders.end(), local_borders[i]);
            }
        }
    }
    return borders;
}

vector<elipse> get_elipses(double *myimage, long *naxes, int hebras1, int hebras2, int numBetas, double relativeVotes, int alphaMin, vector<int> borders) {
    vector<elipse> newElipses = vector<elipse>();

    omp_set_nested(1);
#pragma omp parallel num_threads(hebras1)
    {
#pragma omp for
        for (int t : borders) {
            // calculate axis for point t
            int y_axis_t = t / naxes[0];
            int x_axis_t = t % naxes[0];
            for (int u : borders) {
                vector<int> votes = vector<int>(numBetas, 0);

                // calculate axis for point u
                int y_axis_u = u / naxes[0];
                int x_axis_u = u % naxes[0];

                // compute o_x, o_y, alpha, theta
                int o_y = (y_axis_t + y_axis_u) / 2;
                int o_x = (x_axis_t + x_axis_u) / 2;
                double alpha = sqrt(pow(x_axis_u - x_axis_t, 2) + pow(y_axis_u - y_axis_t, 2)) / 2;
                double theta = atan2(y_axis_u - y_axis_t, x_axis_u - x_axis_t);

                // check if alpha is greater than alphaMin
                if (alpha < alphaMin) continue;
#pragma omp parallel num_threads(hebras2)
                {
#pragma omp for
                    for (int k : borders) {
                        // check if k is t or u
                        if (k == t || k == u) continue;

                        // calculate axis for point u
                        int y_axis_k = k / naxes[0];
                        int x_axis_k = k % naxes[0];

                        // calculate distance from center to k
                        double delta = sqrt(pow(y_axis_k - o_y, 2) + pow(x_axis_k - o_x, 2));

                        if (delta > alpha) continue;

                        // compute beta, gamma
                        double gamma = sin(theta) * (y_axis_k - o_y) + cos(theta) * (x_axis_k - o_x);
                        double beta = sqrt(
                            (pow(alpha, 2) * pow(delta, 2) - pow(alpha, 2) * pow(gamma, 2)) /
                            (pow(alpha, 2) - pow(gamma, 2)));

                        // check if beta is valid and if it is less than alpha
                        if (isnan(beta) || beta < 0 || beta > alpha) continue;

                        // discretize beta and add vote
                        int beta_index = beta / ((double)max(naxes[0], naxes[1]) / (2 * numBetas));
#pragma omp atomic
                        votes[beta_index]++;
                    }
                }
#pragma omp critical
                {
                    // save elipse if it has the minimum number of votes
                    for (int i = 0; i < numBetas; i++) {
                        // calculate ellipse circumference
                        int CE = M_PI * (3 * (alpha + i) - sqrt((3 * alpha + i) * (alpha + 3 * i)));
                        if (votes[i] > CE * relativeVotes) {
                            theta = theta * (180 / M_PI);
                            elipse newElipse = {o_x, o_y, alpha, theta, i};
                            newElipses.insert(newElipses.end(), newElipse);
                        }
                    }
                }
            }
        }
    }

    return newElipses;
}