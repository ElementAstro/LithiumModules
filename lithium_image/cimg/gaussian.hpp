#include <gsl/gsl_multifit_nlin.h>
#include <array>
#include "type.hpp"

class GaussianFitTraitsT
{
private:
public:
    struct CurveParamsT
    {
        // b = base, p = peak, c = center in x, w = mean width (FWHM)
        enum TypeE
        {
            B_IDX = 0,
            P_IDX,
            C_IDX,
            W_IDX,
            _Count
        };
        struct TypeT : public std::array<float, TypeE::_Count>
        {
            TypeT(const gsl_vector *inVec);

            TypeT() = default;
        };
    };

    /* Makes a guess for b, p, c and w based on the supplied data */
    static void makeGuess(const GslMultiFitParmsT &inData, gsl_vector *guess);

    /* y = b + p * exp(-0.5f * ((t - c) / w) * ((t - c) / w)) */
    static float fx(float x, const CurveParamsT::TypeT &inParms);

    /* Calculates f(x) = b + p * e^[0.5*((x-c)/w)] for each data point. */
    static int gslFx(const gsl_vector *x, void *inGslParams, gsl_vector *outResultVec);

    /* Calculates the Jacobian (derivative) matrix of f(x) = b + p * e^[0.5*((x-c)/w)^2] for each data point */
    static int gslDfx(const gsl_vector *x, void *params, gsl_matrix *J);

    /* Invokes f(x) and f'(x) */
    static int gslFdfx(const gsl_vector *x, void *params, gsl_vector *f, gsl_matrix *J);
};