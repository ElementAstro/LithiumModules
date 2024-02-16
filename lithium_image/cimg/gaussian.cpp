#include "gaussian.hpp"
#include "type.hpp"

GaussianFitTraitsT::CurveParamsT::TypeT::TypeT(const gsl_vector *inVec = 0)
{
    for (size_t i = 0; i < TypeE::_Count; ++i)
    {
        TypeE idx = static_cast<TypeE>(i);
        (*this)[i] = (inVec ? gsl_vector_get(inVec, idx) : 0);
    }
}

/* Makes a guess for b, p, c and w based on the supplied data */
void GaussianFitTraitsT::makeGuess(const GslMultiFitParmsT &inData, gsl_vector *guess)
{
    size_t numDataPoints = inData.size();
    float y_mean = 0;
    float y_max = inData.at(0).pt.y;
    float c = inData.at(0).pt.x;

    for (size_t i = 0; i < numDataPoints; ++i)
    {
        const DataPointT &dataPoint = inData.at(i).pt;

        y_mean += dataPoint.y;

        if (y_max < dataPoint.y)
        {
            y_max = dataPoint.y;
            c = dataPoint.x;
        }
    }

    y_mean /= (float)numDataPoints;
    float w = (inData.at(numDataPoints - 1).pt.x - inData.at(0).pt.x) / 10.0;

    gsl_vector_set(guess, CurveParamsT::B_IDX, y_mean);
    gsl_vector_set(guess, CurveParamsT::P_IDX, y_max);
    gsl_vector_set(guess, CurveParamsT::C_IDX, c);
    gsl_vector_set(guess, CurveParamsT::W_IDX, w);
}

/* y = b + p * exp(-0.5f * ((t - c) / w) * ((t - c) / w)) */
float GaussianFitTraitsT::fx(float x, const CurveParamsT::TypeT &inParms)
{
    float b = inParms[CurveParamsT::B_IDX];
    float p = inParms[CurveParamsT::P_IDX];
    float c = inParms[CurveParamsT::C_IDX];
    float w = inParms[CurveParamsT::W_IDX];
    float t = ((x - c) / w);
    t *= t;
    return (b + p * exp(-0.5f * t));
}

/* Calculates f(x) = b + p * e^[0.5*((x-c)/w)] for each data point. */
int GaussianFitTraitsT::gslFx(const gsl_vector *x, void *inGslParams, gsl_vector *outResultVec)
{
    CurveParamsT::TypeT curveParams(x);                                      // Store the current coefficient values
    const GslMultiFitParmsT *gslParams = ((GslMultiFitParmsT *)inGslParams); // Store parameter values

    // Execute Levenberg-Marquart on f(x)
    for (size_t i = 0; i < gslParams->size(); ++i)
    {
        const GslMultiFitDataT &gslData = gslParams->at(i);
        float yi = GaussianFitTraitsT::fx((float)gslData.pt.x, curveParams);
        gsl_vector_set(outResultVec, i, (yi - gslData.y) / gslData.sigma);
    }
    return GSL_SUCCESS;
}

/* Calculates the Jacobian (derivative) matrix of f(x) = b + p * e^[0.5*((x-c)/w)^2] for each data point */
int GaussianFitTraitsT::gslDfx(const gsl_vector *x, void *params, gsl_matrix *J)
{

    // Store parameter values
    const GslMultiFitParmsT *gslParams = ((GslMultiFitParmsT *)params);

    // Store current coefficients
    float p = gsl_vector_get(x, CurveParamsT::P_IDX);
    float c = gsl_vector_get(x, CurveParamsT::C_IDX);
    float w = gsl_vector_get(x, CurveParamsT::W_IDX);

    // Store non-changing calculations
    float w2 = w * w;
    float w3 = w2 * w;

    for (size_t i = 0; i < gslParams->size(); ++i)
    {
        const GslMultiFitDataT &gslData = gslParams->at(i);
        float x_minus_c = (gslData.pt.x - c);
        float e = exp(-0.5f * (x_minus_c / w) * (x_minus_c / w));

        gsl_matrix_set(J, i, CurveParamsT::B_IDX, 1 / gslData.sigma);
        gsl_matrix_set(J, i, CurveParamsT::P_IDX, e / gslData.sigma);
        gsl_matrix_set(J, i, CurveParamsT::C_IDX, (p * e * x_minus_c) / (gslData.sigma * w2));
        gsl_matrix_set(J, i, CurveParamsT::W_IDX, (p * e * x_minus_c * x_minus_c) / (gslData.sigma * w3));
    }
    return GSL_SUCCESS;
}

/* Invokes f(x) and f'(x) */
int GaussianFitTraitsT::gslFdfx(const gsl_vector *x, void *params, gsl_vector *f, gsl_matrix *J)
{
    gslFx(x, params, f);
    gslDfx(x, params, J);

    return GSL_SUCCESS;
}