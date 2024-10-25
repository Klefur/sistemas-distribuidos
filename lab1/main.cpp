#include <omp.h>
#include <unistd.h>

#include <chrono>
#include <cmath>
#include <vector>

#include "fitsio.h"
#include "stdio.h"

using namespace std;
using namespace std::chrono;

typedef struct elipse {
    int o_x, o_y;
    double alpha, theta;
    int beta;
} elipse;

int main(int argc, char *argv[]) {
    int status;
    fitsfile *fptr; /* pointer to the FITS file; defined in fitsio.h */
    long fpixel = 1, naxis = 2, exposure;
    long naxes[2];

    int numBetas, relativeVotes, alphaMin, hebras1, hebras2;
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
                relativeVotes = atoi(optarg);
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

    int variable = 4;

    auto start = high_resolution_clock::now();
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

    vector<elipse> newElipses = vector<elipse>();

#pragma omp parallel num_threads(hebras1)
    {
#pragma omp for
        for (int t : borders) {
            // calculate axis for point t
            int x_axis_t = t / naxes[0];
            int y_axis_t = t % naxes[0];
            for (int u : borders) {
                vector<int> votes = vector<int>(numBetas, 0);

                // calculate axis for point u
                int x_axis_u = u / naxes[0];
                int y_axis_u = u % naxes[0];

                // compute o_x, o_y, alpha, theta
                int o_x = (x_axis_t + x_axis_u) / 2;
                int o_y = (y_axis_t + y_axis_u) / 2;
                double alpha = sqrt(pow(x_axis_u - x_axis_t, 2) + pow(y_axis_u - y_axis_t, 2)) / 2;
                double theta = atan2(y_axis_u - y_axis_t, x_axis_u - x_axis_t);

                if (alpha < alphaMin) continue;

                for (int k : borders) {
                    if (k == t && k == u) {
                        continue;
                    }

                    // calculate axis for point u
                    int x_axis_k = k / naxes[0];
                    int y_axis_k = k % naxes[0];

                    // calculate distance from center to k
                    double delta = sqrt(pow(y_axis_k - o_y, 2) + pow(x_axis_k - o_x, 2));

                    if (delta > alpha) continue;

                    // compute beta, gamma
                    double gamma = sin(theta) * (y_axis_k - o_y) - cos(theta) * (x_axis_k - o_x);
                    double beta = sqrt((pow(alpha, 2) * pow(delta, 2) - pow(alpha, 2) * pow(gamma, 2)) / (pow(alpha, 2) - pow(gamma, 2)));

                    if (isnan(beta) || beta < 0 || beta > (max(naxes[0], naxes[1])) / 2) continue;

                    // discretize beta and add vote
                    int beta_index = beta / ((double)max(naxes[0], naxes[1]) / (2 * numBetas));
                    votes[beta_index] += 1;
                }

#pragma omp critical
                {
                    // save elipse if it has the minimum number of votes
                    for (int i = 0; i < numBetas; i++) {
                        if (votes[i] > relativeVotes) {
                            elipse newElipse = {o_x, o_y, alpha, theta, i};
                            newElipses.insert(newElipses.end(), newElipse);
                        }
                    }
                }
            }
        }
    }
    // After function call
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);

    fits_close_file(fptr, &status);

    for (elipse e : newElipses) {
        printf("%d, %d, %f, %f, %d\n", e.o_x, e.o_y, e.alpha, e.theta, e.beta);
    }


    double seconds = (double)duration.count() / 1000000.00;
    printf("Time taken by function: %f seconds\n", seconds);
    return 0;
}
