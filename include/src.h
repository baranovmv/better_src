#ifndef SRC_SRC_H
#define SRC_SRC_H

#if defined(__GNUC__)
#define SRC_API __attribute__((visibility("default")))
#else /* !__GNUC__ */
#error "unsupported compiler"
#endif /* __GNUC__ */

#ifdef __cplusplus
extern "C" {
#endif

typedef struct src_t src_t;

typedef enum
{
    SRC_PROFILE_DEFAULT,
} SrcProfile_t;

typedef enum
{
    MONO,
    STEREO
} SrcNumChannels_t;

typedef enum
{
    S8 =  8,
    S16 = 16,
    S20 = 20,
    S24 = 24,
    F32 = 0x80+32
} SrcSampleType;

SRC_API src_t * src_open(SrcProfile_t profile, SrcNumChannels_t nchannels, SrcSampleType sample_type,
    unsigned int in_fs, unsigned int out_fs);
SRC_API void src_set_scale(src_t * src, float coeff);
SRC_API int src_push_samples(src_t *src, float *samples, unsigned int nsamples);
SRC_API unsigned int src_pop_samples(src_t *src, float *out, unsigned int max_out);
SRC_API float src_left_to_process(const src_t * src);
SRC_API void src_close(src_t * src);

#ifdef __cplusplus
} /* extern "C" */
#endif


#endif //SRC_SRC_H
