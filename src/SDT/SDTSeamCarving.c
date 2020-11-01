#include "SDTSeamCarving.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <math.h>

#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "SDTCommon.h"
#include "SDTComplex.h"
#include "SDTFFT.h"

/*
 * HALF_WIDTH is a value use to compute the dimension of neighbourhood.
 * Set 1 to choose 3 neighbours; set 2 to choose 5 neighbours and so on.
 */
#define HALF_WIDTH 1

/*
 * Defining low thresholds in linear scale:
 *  - 80dB;
 *  - 60dB;
 */
#define dB80 pow(10,-80.0/20)
#define dB60 pow(10,-60.0/20)

typedef double realtype;

typedef unsigned int bool;
#define false       0
#define true        1

/*
 * Defining the support structure to describe a seam:
 * - start: index where the mode starts in the seam;
 * - peak: index of the max magnitude in the mode;
 * - end: index where the mode ends;
 * - decay: in an impact sound, the estimated decay time in seconds;
 * - magnitudes: vector with all magnitude values in the seam;
 * - seam: a vector of y-axis coordinates representing the seam through time (x-axis);
 */
typedef struct SDTSeam {
    long start;
    long peak;
    long end;
    long decay;
    gsl_vector *magnitudes, *phases;
    gsl_vector_int *seam;
} SDTSeam;

/*
 * Creating an SDTMode:
 * - nFrames: is the number of DFTs frames analyzed.
 */
SDTSeam *SDTSeam_new(long nFrames) {
    SDTSeam *x;

    x = (SDTSeam *) malloc(sizeof(SDTSeam));
    x->start = 0;
    x->peak = 0;
    x->end = 0;
    x->decay = 0;
    x->magnitudes = gsl_vector_alloc(nFrames);
    x->phases = gsl_vector_alloc(nFrames);
    x->seam = gsl_vector_int_alloc(nFrames);
    return x;
}

/*
 * Free the mode allocated memory.
 */
void SDTSeam_free(SDTSeam *x) {
    gsl_vector_int_free(x->seam);
    gsl_vector_free(x->magnitudes);
    gsl_vector_free(x->phases);
    free(x);
}

//-----------------------------------------------------------------------//

/*
 * This structure contains the follow:
 * - buffer: variable used to read samples;
 * - mags: variable used to store all mags of the spectrum;
 * - modes: the seams extracted;
 * - mX: a GSL view of the mags in a matrix shape (nFrames, fftSize);
 * - nFrames: the number of DFTs;
 * - nSamples: the number of samples;
 * - nModes: the number of mode (seams) to compute;
 * - winSize: window size for the STFT;
 * - fftSize: fft size for the STFT;
 * - hopSize: hop size for the STFT;
 */
struct SDTSeamCarving {
    double *buffer, *mags, *phases, magMax;
    SDTSeam **seams;
    gsl_matrix_view mX, pX;
    long nFrames, nSamples, nModes, bufferSize, winSize, fftSize, hopSize;
};

SDTSeamCarving *SDTSeamCarving_new(long nModes, long bufferSize, long winSize) {
    SDTSeamCarving *x;
    x = (SDTSeamCarving *) malloc(sizeof(SDTSeamCarving));
    x->buffer = (double *) calloc(bufferSize, sizeof(double));
    x->mags = NULL;
    x->phases = NULL;
    x->seams = (SDTSeam **) calloc(nModes, sizeof(SDTSeam *));
    x->nSamples = 0;
    x->magMax = 0;
    x->nFrames = 0;
    x->nModes = nModes;
    x->bufferSize = bufferSize;
    x->winSize = winSize;
    x->fftSize = winSize;           // Experimentally better
    x->hopSize = winSize / 4;       // ---------------------
    return x;
}

void SDTSeamCarving_free(SDTSeamCarving *x) {
    long i;

    free(x->buffer);
    free(x->mags);
    free(x->phases);

    for (i = 0; i < x->nModes; i++) {
        if (x->seams[i]) SDTSeam_free(x->seams[i]);
    }
    free(x->seams);

    free(x);
}

void SDTSeamCarving_setOverlap(SDTSeamCarving *x, double f) {
    x->hopSize = SDT_clip((1.0 - f) * x->winSize, 1, x->winSize);
}

long SDTSeamCarving_readSamples(SDTSeamCarving *x, double *in, long n) {
    long start, end, span;
    start = x->nSamples;
    end = SDT_clip(start + n, start, x->bufferSize);
    if (start == end) return 0;
    span = end - start;
    memcpy(&x->buffer[start], in, span * sizeof(double));
    x->nSamples = end;
    return span;
}

long SDTSeamCarving_clearSamples(SDTSeamCarving *x, long n) {
    long start, end, span;

    end = x->nSamples;
    start = SDT_clip(end - n, 0, end);
    if (start == end) return 0;
    span = end - start;
    memset(&x->buffer[start], 0, span * sizeof(double));
    x->nSamples = start;
    return span;
}

void SDTSeamCarving_STFT(SDTSeamCarving *x) {

    double *currMags, *currPhases, *window;
    long frameIndex, binIndex, nSamples;
    SDTComplex *fft;
    SDTFFT *fftPlan;

    nSamples = x->nSamples < x->bufferSize - x->winSize ? x->nSamples : x->bufferSize - x->winSize;
    x->nFrames = ceil(nSamples / x->hopSize);
    x->mags = (double *) calloc(x->nFrames * x->fftSize, sizeof(double));
    x->phases = (double *) calloc(x->nFrames * x->fftSize, sizeof(double));

    // reminder: winSize == fftSize, following experimental results.
    fft = (SDTComplex *) calloc(x->fftSize, sizeof(SDTComplex));
    fftPlan = SDTFFT_new(x->fftSize / 2);
    window = (double *) calloc(x->winSize, sizeof(double));

    // Compute the STFT
    for (frameIndex = 0; frameIndex < x->nFrames; frameIndex++) {
        memcpy(window, &x->buffer[frameIndex * x->hopSize], x->winSize * sizeof(double));
        SDT_hanning(window, x->winSize);
        SDTFFT_fftr(fftPlan, window, fft);
        currMags = &x->mags[frameIndex * x->fftSize];
        currPhases = &x->phases[frameIndex * x->fftSize];
        for (binIndex = 0; binIndex < x->fftSize; binIndex++) {
            currMags[binIndex] = SDTComplex_abs(fft[binIndex]);
            currPhases[binIndex] = SDTComplex_angle(fft[binIndex]);
        }

    }

    x->mX = gsl_matrix_view_array(x->mags, x->nFrames, x->fftSize);
    x->pX = gsl_matrix_view_array(x->phases, x->nFrames, x->fftSize);

    // Normalization between 0 and 1
    x->magMax = gsl_matrix_max(&x->mX.matrix);
    gsl_matrix_scale(&x->mX.matrix, 1 / x->magMax);


    free(window);
    free(fft);
    free(fftPlan);

}

// Extracts the k-th seam from the spectrogram.
void SDTSeamCarving_findVerticalSeam(SDTSeamCarving *x, long k) {

    /*
     * Allocate and initialize the support matrices
     * - Q: this is the pointer to the distTo matrix
     *      used to contains all costs, the matrix
     *      should be initialized to -inf, but we do
     *      not have negative values, so -1.0 will suffice.
     * - P: this is a pointer to a support matrix used
     *      to find the seam indices. 0 is the default
     *      straight column direction.
     */
    long rows = x->nFrames, cols = x->fftSize;
    gsl_matrix *Q = gsl_matrix_alloc(rows, cols);
    gsl_matrix_set_all(Q, -1.0);
    gsl_matrix_int *P = gsl_matrix_int_alloc(rows, cols);
    gsl_matrix_int_set_zero(P);

    /*
     * First we initialize the first row costs,
     * because we have no history so far.
     * This corresponds to:
     *     - Q[0][j] = mX[0][j];
     */
    for (int j = 0; j < cols; j++) {
        double mX_0j = gsl_matrix_get(&x->mX.matrix, 0, j);
        gsl_matrix_set(Q, 0, j, mX_0j);
    }

    /*
     * Now the algorithm will compute the longest path.
     *
     * Giving a parametric window of possible neighbours at each iteration,
     * it will find and eventually update the maximal cost.
     * We do not use gsl sub vectors here due to computational convenience.
     */
    for (long i = 1; i < rows; i++) {
        for (long j = 0; j < cols; j++) {
            // The neighbours elements goes from:
            for (int t = -HALF_WIDTH; t <= HALF_WIDTH; t++) {
                // Managing the borders of the matrix:
                if (j + t >= 0 && j + t < cols) {
                    /*
                     * Defining the possible descending cost from that branch
                     * as the previous row element at j+t column position plus
                     * the magnitude value, in linear scale, in the spectrogram
                     *     - cost = Q[i-1][j+t] + mX[i][j]
                     */
                    double cost = (gsl_matrix_get(Q, i - 1, j + t) + gsl_matrix_get(&x->mX.matrix, i, j));
                    if (cost > gsl_matrix_get(Q, i, j)) {
                        // Updating Q[i][j] and P with the new greater cost.
                        gsl_matrix_set(Q, i, j, cost);
                        gsl_matrix_int_set(P, i, j, t);
                    }
                }
            }
        }
    }

    /*
     * Now the backtracking must be performed.
     * Reading the P matrix, we follow the seam
     * indexes which maximize the overall cost.
     */

    gsl_vector_int *seam = gsl_vector_int_alloc(rows);

    /*
     * Find the max index in the last row.
     * When there are several equal maximum elements
     * then the lowest index is considered.
     *
     * The gsl_matrix_row is a view and points to the
     * same memory block. De-allocating Q will suffice
     * in order to clear also the view.
     */
    gsl_vector last_row = gsl_matrix_row(Q, rows - 1).vector;
    gsl_vector_int_set(seam, rows - 1, gsl_vector_max_index(&last_row));

    /*
     * Backtracking:
     * finding the column index for the i-th row.
     *
     * This corresponds to:
     *  -  seam[i] = seam[i + 1] + P[i + 1][seam[i + 1]];
     *  where P[i + 1][seam[i + 1]] is called P_index.
     */
    for (long i = rows - 2; i >= 0; i--) {
        const int P_index = gsl_matrix_int_get(P, i + 1, gsl_vector_int_get(seam, i + 1));
        gsl_vector_int_set(
                seam,
                i,
                gsl_vector_int_get(seam, i + 1) + P_index
        );
    }

    /*
     * Now we need only the computed seam positions.
     * De-allocating the support matrices.
     */
    gsl_matrix_free(Q);
    gsl_matrix_int_free(P);

    /*
     * Register the full seam
     */
    x->seams[k] = SDTSeam_new(rows);
    x->seams[k]->seam = seam;

    /*
     * Find properties
     */
    SDTSeamCarving_trimSeam(x, k);
    SDTSeamCarving_backwardIntegration(x, k);
    SDTSeamCarving_eraseSeam(x, k);

}

// Computes the magnitudes and indexes of the k-th seam
void SDTSeamCarving_trimSeam(SDTSeamCarving *x, long k) {


    long rows = x->nFrames;
    gsl_vector *magnitudes = gsl_vector_alloc(rows);
    gsl_vector *phases = gsl_vector_alloc(rows);

    x->seams[k]->magnitudes = magnitudes;
    x->seams[k]->phases = phases;

    /*
     * Compute the seam properties reading the magnitudes
     * values (seam_magnitudes) and identifying the first
     * point greater than 80dB (top).
     *
     * Each magnitude value is taken as mX[i][seam[i]];
     */
    long peak = 0, top = 0, bottom = rows - 1;
    double magnitude_max = 0;

    for (long i = 0; i < rows; i++) {
        double magnitude = gsl_matrix_get(&x->mX.matrix, i, gsl_vector_int_get(x->seams[k]->seam, i));
        double phase = gsl_matrix_get(&x->pX.matrix, i, gsl_vector_int_get(x->seams[k]->seam, i));
        gsl_vector_set(magnitudes, i, magnitude);
        gsl_vector_set(phases, i, phase);
        if (magnitude > magnitude_max) {
            magnitude_max = magnitude;
            peak = i;
            if (!top && magnitude > dB80)
                // The first value greater than the threshold
                top = i;
        }
    }

    if (magnitude_max < dB80) {
        // fputs("[Info] Seam values are too low.\n", stderr);
        return;
    }

    /*
     * Find the point where the magnitude decreases below
     * 60db of attenuation from the peak (bottom).
     *
     * There could be a case where peak = top = 0;
     */
    const double stop_condition = magnitude_max * dB60;

    for (long i = peak; i < rows; i++) {
        if (gsl_vector_get(magnitudes, i) <= stop_condition) {
            bottom = i > 0 ? i - 1 : 0;
            break;
        }
    }

    if (top >= bottom) {
        // fputs("[Info] No values for the current seam.\n", stderr);
        return;
    }

    x->seams[k]->start = top;
    x->seams[k]->peak = peak;
    x->seams[k]->end = bottom;

}

void SDTSeamCarving_eraseSeam(SDTSeamCarving *x, long k) {

    SDTSeam *s;
    long i, j, currMagYPos;
    double currMag;

    s = x->seams[k];

    if (!s) {
        return;
    }

    // Removing the seam from the spectrogram,
    // barely erasing the area around the seam position
    // until the magnitude is less than 80dB
    for (i = 0; i < x->nFrames; i++) {
        currMagYPos = gsl_vector_int_get(s->seam, i);

        currMag = gsl_matrix_get(&x->mX.matrix, i, currMagYPos);
        for (j = currMagYPos + 1; j < x->fftSize && currMag > dB60; j++) {
            currMag = gsl_matrix_get(&x->mX.matrix, i, j);
            gsl_matrix_set(&x->mX.matrix, i, j, 0.0);
        }

        currMag = gsl_matrix_get(&x->mX.matrix, i, currMagYPos);
        for (j = currMagYPos - 1; j > 0 && currMag > dB60; j--) {
            currMag = gsl_matrix_get(&x->mX.matrix, i, j);
            gsl_matrix_set(&x->mX.matrix, i, j, 0.0);
        }

        gsl_matrix_set(&x->mX.matrix, i, gsl_vector_int_get(s->seam, i), 0.0);
    }
}

void SDTSeamCarving_popMode(SDTSeamCarving *x, long k, double *decays, double *mags, double *freqs, double *phases) {

    SDTSeam *s;

    s = x->seams[k];
    mags[k] = 0.0;
    freqs[k] = 0.0;
    phases[k] = 0.0;
    decays[k] = 0.0;

    if (!s)
        return;

    mags[k] = gsl_vector_get(s->magnitudes, s->peak) * x->magMax;

    if (mags[k] < 1) {
        mags[k] = 0.0;
        return;
    }

    decays[k] = (float) (s->decay * x->hopSize) / SDT_sampleRate;
    phases[k] = SDT_wrap(gsl_vector_get(s->phases, s->peak));
    freqs[k] = gsl_vector_int_get(s->seam, s->peak) * SDT_sampleRate / (float) x->fftSize;

}

// f is the time
void SDTSeamCarving_popMode_t(SDTSeamCarving *x, long k, double time, double *mags, double *freqs, double *phases) {

    SDTSeam *s;
    double frameTime;

    s = x->seams[k];
    mags[k] = 0.0;
    freqs[k] = 0.0;
    phases[k] = 0.0;

    if (!s)
        return;

    frameTime = (int) floor(time / (x->hopSize * SDT_timeStep));

    if (frameTime >= s->start && frameTime < s->end) {
        mags[k] = gsl_vector_get(s->magnitudes, frameTime) * x->magMax;
        phases[k] = SDT_wrap(gsl_vector_get(s->phases, frameTime));
        freqs[k] = gsl_vector_int_get(s->seam, frameTime) * SDT_sampleRate / (float) x->fftSize;
    }
}

/** @brief Compute the moore-penrose pseudo-inverse of a matrix.
 * Source of this code: <a href="https://gist.github.com/turingbirds/5e99656e08dbe1324c99">gist</a>.
@param[in] A Input matrix; this matrix is destroyed but its memory still must to be free.
@param[in] rcond A real number specifying the singular value threshold for inclusion.
@return The new allocated A_pinv matrix containing the result. */
gsl_matrix *moore_penrose_pinv(gsl_matrix *A, const realtype rcond) {
    gsl_matrix *V, *Sigma_pinv, *U, *A_pinv;
    gsl_matrix *_tmp_mat = NULL;
    gsl_vector *_tmp_vec;
    gsl_vector *u;
    realtype x, cutoff;
    size_t i, j;
    unsigned int n = A->size1;
    unsigned int m = A->size2;
    bool was_swapped = false;


    if (m > n) {
        /* libgsl SVD can only handle the case m <= n - transpose matrix */
        was_swapped = true;
        _tmp_mat = gsl_matrix_alloc(m, n);
        gsl_matrix_transpose_memcpy(_tmp_mat, A);
        A = _tmp_mat;
        i = m;
        m = n;
        n = i;
    }

    /* do SVD */
    V = gsl_matrix_alloc(m, m);
    u = gsl_vector_alloc(m);
    _tmp_vec = gsl_vector_alloc(m);
    gsl_linalg_SV_decomp(A, V, u, _tmp_vec);
    gsl_vector_free(_tmp_vec);

    /* compute Σ⁻¹ */
    Sigma_pinv = gsl_matrix_alloc(m, n);
    gsl_matrix_set_zero(Sigma_pinv);
    cutoff = rcond * gsl_vector_max(u);

    for (i = 0; i < m; ++i) {
        if (gsl_vector_get(u, i) > cutoff) {
            x = 1. / gsl_vector_get(u, i);
        } else {
            x = 0.;
        }
        gsl_matrix_set(Sigma_pinv, i, i, x);
    }

    /* libgsl SVD yields "thin" SVD - pad to full matrix by adding zeros */
    U = gsl_matrix_alloc(n, n);
    gsl_matrix_set_zero(U);

    for (i = 0; i < n; ++i) {
        for (j = 0; j < m; ++j) {
            gsl_matrix_set(U, i, j, gsl_matrix_get(A, i, j));
        }
    }

    if (_tmp_mat != NULL) {
        gsl_matrix_free(_tmp_mat);
    }

    /* two dot products to obtain pseudoinverse */
    _tmp_mat = gsl_matrix_alloc(m, n);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1., V, Sigma_pinv, 0., _tmp_mat);

    if (was_swapped) {
        A_pinv = gsl_matrix_alloc(n, m);
        gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1., U, _tmp_mat, 0., A_pinv);
    } else {
        A_pinv = gsl_matrix_alloc(m, n);
        gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1., _tmp_mat, U, 0., A_pinv);
    }

    gsl_matrix_free(_tmp_mat);
    gsl_matrix_free(U);
    gsl_matrix_free(Sigma_pinv);
    gsl_vector_free(u);
    gsl_matrix_free(V);

    return A_pinv;

}

// Estimates the decay index in the x axis for the k-th seam.
void SDTSeamCarving_backwardIntegration(SDTSeamCarving *x, long k) {

    long rows = x->nFrames;

    /*
     * Prepare the support curve vector
     * of backward cumulated sum in dB scale
     */
    gsl_vector *y = gsl_vector_alloc(rows);

    /*
     * Compute the backward integration and converts values in dB.
     */
    double bwint = 0;
    for (long i = rows - 1; i >= 0; i--) {
        bwint += gsl_vector_get(x->seams[k]->magnitudes, rows - (i + 1));
        gsl_vector_set(y, i, 20 * log10(bwint));
    }

    /*
     *  Compute the linear regression.
     *
     *  First prepare a matrix as follow:
     *      - The first row must be all ones, to allow a non-zero intercept in the line equation;
     *      - The second row is the x axis.
     *
     *  This should appear like this:
     *      [1,1,1,1,1,...]
     *      [0,1,2,3,4,...] <- like x = np.arange(y) in Python
     */
    gsl_matrix *X = gsl_matrix_alloc(2, rows);

    /*
     * Using rows as cols, consider the transposed matrix.
     *     X[0][i]=1
     *     X[0][i]=i
     */
    for (int i = 0; i < rows; i++) {
        gsl_matrix_set(X, 0, i, 1);
        gsl_matrix_set(X, 1, i, i);
    }

    /*
     * Compute the (Moore-Penrose) pseudo-inverse and
     * using the normal equation, obtains the new y elements.
     * This is equivalent to the Python form:
     *     w = np.dot(pinv(X @ X.T) @ X, y.T)
     *     Xt = np.dot(w.T, X)
     */
    gsl_matrix *X_XT = gsl_matrix_alloc(2, 2);                          // (XX') ... (2,rows) @ (rows,2)
    gsl_blas_dgemm(CblasNoTrans, CblasTrans,
                   1.0, X, X,
                   0.0, X_XT);

    gsl_matrix *pinv = moore_penrose_pinv(X_XT, 1E-15);                   // (XX')^-1 ... (2,2)

    gsl_matrix *pinv_X = gsl_matrix_alloc(2, rows);                         // (XX')^-1 X ... (2,rows)
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,
                   1.0, pinv, X,
                   0.0, pinv_X);

    gsl_vector *w = gsl_vector_alloc(2);                                    // becomes a vector (2,1)
    gsl_blas_dgemv(CblasNoTrans, 1.0, pinv_X, y, 0.0, w);

    gsl_vector *Xt = gsl_vector_alloc(rows);
    gsl_blas_dgemv(CblasTrans, 1.0, X, w, 0.0, Xt);         // And finally returns a (rows,1)

    /*
     * computing the time intersection in the x axis:
     *     attenuated_curve = np.abs(Xt - Xt[0] + 60)
     *     decay = np.argmin(attenuated_curve)
     */
    const double Xt0 = -gsl_vector_get(Xt, 0) + 60;
    for (int i = 0; i < rows; i++)
        gsl_vector_set(Xt, i, fabs(gsl_vector_get(Xt, i) + Xt0));


    x->seams[k]->decay = gsl_vector_min_index(Xt);

    /*
     * Cleaning support vectors and matrices.
     */
    gsl_vector_free(y);
    gsl_matrix_free(X);
    gsl_matrix_free(X_XT);
    gsl_matrix_free(pinv);
    gsl_matrix_free(pinv_X);
    gsl_vector_free(w);
    gsl_vector_free(Xt);


};


