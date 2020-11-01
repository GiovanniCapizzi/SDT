/** @file SDTSeamCarving.h
@defgroup modaltracker Modal analyzer using seam carving
This object detects the prominent modal components
 in an audio segment and tracks their temporal evolution
 across time using the seam carving algorithm.
@{ */

#ifndef SDT_SEAMCARVING_H
#define SDT_SEAMCARVING_H

#ifdef __cplusplus
extern "C" {
#endif

/** @brief Opaque data structure for a modal analyzer. */
typedef struct SDTSeamCarving SDTSeamCarving;

/** @brief Instantiates a modal analyzer object.
@param[in] bufferSize Size of the audio buffer, in samples
@param[in] winSize Size of the analysis window, in samples
@return Pointer to the new instance */
SDTSeamCarving *SDTSeamCarving_new(long nModes, long bufferSize, long winSize);

/** @brief Destroys a modal analyzer instance.
@param[in] x Pointer to the instance to destroy */
extern void SDTSeamCarving_free(SDTSeamCarving *x);

/** @brief Sets the analysis window overlapping ratio.
Accepted values go from 0.0 to 1.0, with 0.0 meaning no overlap
and 1.0 meaning total overlap.
@param[in] x Pointer to the instance
@param[in] f Overlap ratio [0.0, 1.0] */
extern void SDTSeamCarving_setOverlap(SDTSeamCarving *x, double f);

/** @brief Reads incoming sound samples and stores them into the object's audio buffer.
@param[in] x Pointer to the instance
@param[in] in Samples to read
@param[in] n Number of samples to read
@return Number of samples actually stored in the buffer */
extern long SDTSeamCarving_readSamples(SDTSeamCarving *x, double *in, long n);

/** @brief Clears n samples from the object's audio buffer.
@param[in] x Pointer to the instance
@param[in] n Number of samples to clear
@return Number of samples actually cleared from the buffer */
extern long SDTSeamCarving_clearSamples(SDTSeamCarving *x, long n);

/** @brief Compute the STFT for the seam extraction.
@param[in] x Pointer to the instance */
extern void SDTSeamCarving_STFT(SDTSeamCarving *x);

/** @brief Extracts the k-th seam from the spectrogram.
@param[in] x Pointer to the instance
@param[in] k The seam index to compute */
extern void SDTSeamCarving_findVerticalSeam(SDTSeamCarving *x, long k);

/** @brief Computes the magnitudes and indexes of the k-th seam.
@param[in] x Pointer to the instance
@param[in] k The seam index to trim */
extern void SDTSeamCarving_trimSeam(SDTSeamCarving *x, long k);

/** @brief Estimates the decay index in the x axis for the k-th seam.
@param[in] x Pointer to the instance
@param[in] k The seam index to analyze */
extern void SDTSeamCarving_backwardIntegration(SDTSeamCarving *x, long k);

void SDTSeamCarving_eraseSeam(SDTSeamCarving *x, long k);

void SDTSeamCarving_popMode(SDTSeamCarving *x, long k, double *decays, double *mags, double *freqs, double *phases);

void SDTSeamCarving_popMode_t(SDTSeamCarving *x, long k, double time, double *mags, double *freqs, double *phases);

#ifdef __cplusplus
};
#endif


#endif //SDT_SEAMCARVING_H

/** @} */
