#ifndef WEIGHTED_COUNT_CALCS_H
#define WEIGHTED_COUNT_CALCS_H





int getWeightedCountCExt_main(uint8_t *x,
                   double *wcount, double *resp,
                   int wcount_dim0, int wcount_dim1,
                   int wcount_dim2, int x_dim0,
                   int x_dim1, int n_threads); 

void *getWeightedCountCExt_worker(uint8_t *x, double *resp,
        double *wcount, int wcount_dim1, int wcount_dim2,
        int x_dim0, int x_dim1, int startRow, int endRow);


int getWeightedCountCExt_single_thread(uint8_t *x, double *wcount, double *resp,
                   int wc_dim0, int wc_dim1, int wc_dim2,
                   int x_dim0, int x_dim1);

#endif
