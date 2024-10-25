#include <omp.h>

#include <vector>

#include "fitsio.h"
#include "stdio.h"

using namespace std;

typedef struct elipse {
    int o_x, o_y;
    double alpha, beta, theta;
} elipse;

int main(int argc, char *argv[]) {
    int status;
    fitsfile *fptr; /* pointer to the FITS file; defined in fitsio.h */
    long fpixel = 1, naxis = 2, exposure;
    long naxes[2];

    status = 0; /* initialize status before calling fitsio routines */

    fits_open_file(&fptr, argv[1], READONLY, &status);
    fits_get_img_size(fptr, 2, naxes, &status);

    double *myimage = (double *)malloc(naxes[0] * naxes[1] * sizeof(double));
    fits_read_img(fptr, TDOUBLE, fpixel, naxes[0] * naxes[1], NULL, myimage, NULL, &status);

    vector<int> borders = vector<int>();

    int variable = 4;

#pragma omp parallel
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

    for (int t: borders) {
        vector<elipse> foundElipses = vector<elipse>();

        int x_axis_t = t / naxes[0];
        int y_axis_t = t % naxes[0];

        for (int u: borders) {
            vector<int> votes = vector<int>();

            int x_axis_u = u / naxes[0];
            int y_axis_u = u % naxes[0];

            // compute o_x, o_y, alpha, theta
            int o_x, o_y;
            double alpha, theta;

            for (int k: borders) {
                if (k == t && k == u) {
                    continue;
                }

                int x_axis_k = k / naxes[0];
                int y_axis_k = k % naxes[0];

                // compute beta
                double beta = 0;

                // add vote
            }

            // delete votes with less than 3 votes
            for (int vote: votes) {
                if (vote > 3) {
                    elipse e;
                    e.o_x = o_x;
                    e.o_y = o_y;
                    e.alpha = alpha;
                    e.beta = vote[];

                    foundElipses.insert(foundElipses.end(), e);
                }
            }

            

            for (elipse e : foundElipses) {
                // compute distance between elipse and new elipse
                // if distance is less than 10, merge elipses

            }
    }

    printf("Total elements = %ld\n", borders.size());

    fits_close_file(fptr, &status);

    return 0;
}
