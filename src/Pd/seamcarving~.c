#include "m_pd.h"
#include "SDT/SDTCommon.h"
#include "SDT/SDTSeamCarving.h"
#include <stdio.h>
#ifdef NT
#pragma warning( disable : 4244 )
#pragma warning( disable : 4305 )
#endif

static t_class *seamcarving_class;

typedef struct _seamcarving {
  t_object obj;
  SDTSeamCarving *seamcarving;
  t_float f;
  t_outlet *out0, *out1;
  int nModes, nSamples, pickup, isRecording, analized;
} t_seamcarving;

void seamcarving_overlap(t_seamcarving *x, t_float f) {
  SDTSeamCarving_setOverlap(x->seamcarving, f);
}

void seamcarving_clear(t_seamcarving *x) {
  SDTSeamCarving_clearSamples(x->seamcarving, x->nSamples);
  x->analized = 0;
}

void seamcarving_start(t_seamcarving *x) {
  seamcarving_clear(x);
  x->isRecording = 1;
}

void seamcarving_stop(t_seamcarving *x) {
  SDTSeamCarving_STFT(x->seamcarving);
  x->isRecording = 0;
  outlet_bang(x->out1);
}

void seamcarving_bang(t_seamcarving *x) {
  int nModes = x->nModes;
  double mags[nModes], decays[nModes], freqs[nModes], phases[nModes];
  t_atom magAtoms[nModes+1], freqAtoms[nModes], decayAtoms[nModes], phaseAtoms[nModes];
  int i;

  if (x->isRecording) seamcarving_stop(x);
  SETFLOAT(&magAtoms[0], x->pickup);
  for(i = 0; i < nModes; i++){
    if(!x->analized)
      SDTSeamCarving_findVerticalSeam(x->seamcarving, i);
    SDTSeamCarving_popMode(x->seamcarving, i, decays, mags, freqs, phases);
    
    post("phases %f", phases[i]);
    post("mags %f", mags[i]);
    post("freqs %f", freqs[i]);
    post("decays %f", decays[i]);
    
    SETFLOAT(&magAtoms[i+1], mags[i]);
    SETFLOAT(&phaseAtoms[i], phases[i]);
    SETFLOAT(&freqAtoms[i], freqs[i]);
    SETFLOAT(&decayAtoms[i], decays[i]);
  }
  x->analized = 1;

  outlet_anything(x->out0, gensym("pickup"), nModes + 1, magAtoms);
  outlet_anything(x->out0, gensym("freqs"), nModes, freqAtoms);
  outlet_anything(x->out0, gensym("decays"), nModes, decayAtoms);
}

void seamcarving_float(t_seamcarving *x, t_float f) {
  int nModes = x->nModes;
  double mags[nModes], freqs[nModes], phases[nModes];
  t_atom magAtoms[nModes+1], freqAtoms[nModes], decayAtoms[nModes], phaseAtoms[nModes];
  int i;

  if (x->isRecording) seamcarving_stop(x);
  SETFLOAT(&magAtoms[0], x->pickup);
  for(i = 0; i < nModes; i++){
    if(!x->analized)
      SDTSeamCarving_findVerticalSeam(x->seamcarving, i);
    SDTSeamCarving_popMode_t(x->seamcarving, i, f, mags, freqs, phases);
    SETFLOAT(&magAtoms[i+1], mags[i]);
    SETFLOAT(&phaseAtoms[i], phases[i]);
    SETFLOAT(&freqAtoms[i], freqs[i]);
    SETFLOAT(&decayAtoms[i], 0.0);
  }
  x->analized=1;

  outlet_anything(x->out0, gensym("pickup"), nModes + 1, magAtoms);
  outlet_anything(x->out0, gensym("freqs"), nModes, freqAtoms);
  outlet_anything(x->out0, gensym("decays"), nModes, decayAtoms);
}




t_int *seamcarving_perform(t_int *w) {
  t_seamcarving *x = (t_seamcarving *)(w[1]);
  t_float *in = (t_float *)(w[2]);
  int n = (int)w[3];
  double samples[n];
  long i, samplesRead;

  if (x->isRecording) {
    for (i = 0; i < n; i++) samples[i] = in[i];
    samplesRead = SDTSeamCarving_readSamples(x->seamcarving, samples, n);
    if (!samplesRead) {
      seamcarving_stop(x);
    }
  }
  return w + 4;
}

void seamcarving_dsp(t_seamcarving *x, t_signal **sp) {
  SDT_setSampleRate(sp[0]->s_sr);
  dsp_add(seamcarving_perform, 3, x, sp[0]->s_vec, sp[0]->s_n);
}

void seamcarving_pickup(t_seamcarving *x, t_float f) {
  x->pickup = SDT_clip(f, 0, x->nModes);
}

void *seamcarving_new(t_symbol *s, long argc, t_atom *argv) {
  unsigned int tmpSize, windowSize;
  t_float nModes, nSamples;

  t_seamcarving *x = (t_seamcarving *)pd_new(seamcarving_class);
  nModes = argc > 0 && argv[0].a_type == A_FLOAT ? atom_getfloat(&argv[0]) : 8;
  nSamples = argc > 1 && argv[1].a_type == A_FLOAT ? atom_getfloat(&argv[1]) : 441000;
  if (argc > 2 && argv[2].a_type == A_FLOAT) {
    tmpSize = atom_getfloat(&argv[2]);
    windowSize = SDT_nextPow2(tmpSize);
    if (tmpSize != windowSize) {
      post("sdt.seamcarving~: Window size must be a power of 2, setting it to %d", windowSize);
    }
  }
  else {
    windowSize = 1024;
  }
  x->seamcarving = SDTSeamCarving_new(nModes, nSamples, windowSize);
  x->out0 = outlet_new(&x->obj, NULL);
  x->out1 = outlet_new(&x->obj, gensym("bang"));
  x->nModes = nModes;
  x->nSamples = nSamples;
  x->isRecording = 0;
  return (x);
}

void seamcarving_free(t_seamcarving *x) {
  outlet_free(x->out0);
  outlet_free(x->out1);
  SDTSeamCarving_free(x->seamcarving);
}


void seamcarving_tilde_setup(void) {
  seamcarving_class = class_new(gensym("seamcarving~"), (t_newmethod)seamcarving_new, (t_method)seamcarving_free, sizeof(t_seamcarving), CLASS_DEFAULT, A_GIMME, 0);
  CLASS_MAINSIGNALIN(seamcarving_class, t_seamcarving , f);
  class_addmethod(seamcarving_class, (t_method)seamcarving_overlap, gensym("overlap"), A_FLOAT, 0);
  class_addmethod(seamcarving_class, (t_method)seamcarving_clear, gensym("clear"), 0);
  class_addmethod(seamcarving_class, (t_method)seamcarving_pickup, gensym("pickup"), A_FLOAT, 0);
  class_addmethod(seamcarving_class, (t_method)seamcarving_float, gensym("float"), A_FLOAT, 0);
  class_addmethod(seamcarving_class, (t_method)seamcarving_start, gensym("start"), 0);
  class_addmethod(seamcarving_class, (t_method)seamcarving_stop, gensym("stop"), 0);
  class_addmethod(seamcarving_class, (t_method)seamcarving_bang, gensym("bang"), 0);
  class_addmethod(seamcarving_class, (t_method)seamcarving_dsp, gensym("dsp"), 0);
}