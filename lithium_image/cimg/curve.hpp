#include "type.hpp"

template <class FitTraitsT>
class CurveFitTmplT
{
public:
    typedef typename FitTraitsT::CurveParamsT CurveParamsT;

    /**
     * DataAccessor allows specifying how x,y data is accessed.
     * See http://en.wikipedia.org/wiki/Approximation_error for expl. of rel and abs errors.
     */
    template <typename DataAccessorT>
    static int
    fitGslLevenbergMarquart(const typename DataAccessorT::TypeT &inData, typename CurveParamsT::TypeT *outResults,
                            double inEpsAbs, double inEpsRel, size_t inNumMaxIter = 500)
    {
        GslMultiFitParmsT gslMultiFitParms(inData.size());

        // Fill in the parameters
        for (typename DataAccessorT::TypeT::const_iterator it = inData.begin(); it != inData.end(); ++it)
        {
            size_t idx = std::distance(inData.begin(), it);
            const DataPointT &dataPoint = DataAccessorT::getDataPoint(idx, it);
            gslMultiFitParms[idx].y = dataPoint.y;
            gslMultiFitParms[idx].sigma = 0.1f;
            gslMultiFitParms[idx].pt = dataPoint;
        }

        // Fill in function info
        gsl_multifit_function_fdf f;
        f.f = FitTraitsT::gslFx;
        f.df = FitTraitsT::gslDfx;
        f.fdf = FitTraitsT::gslFdfx;
        f.n = inData.size();
        f.p = FitTraitsT::CurveParamsT::_Count;
        f.params = &gslMultiFitParms;

        gsl_vector *guess = gsl_vector_alloc(FitTraitsT::CurveParamsT::_Count); // Allocate the guess vector

        FitTraitsT::makeGuess(gslMultiFitParms, guess); // Make initial guesses based on the data

        // Create a Levenberg-Marquardt solver with n data points and m parameters
        gsl_multifit_fdfsolver *solver = gsl_multifit_fdfsolver_alloc(gsl_multifit_fdfsolver_lmsder,
                                                                      inData.size(), FitTraitsT::CurveParamsT::_Count);
        gsl_multifit_fdfsolver_set(solver, &f, guess); // Initialize the solver

        int status, i = 0;

        // Iterate to to find a result
        do
        {
            i++;
            status = gsl_multifit_fdfsolver_iterate(solver); // returns 0 in case of success
            if (status)
            {
                break;
            }
            status = gsl_multifit_test_delta(solver->dx, solver->x, inEpsAbs, inEpsRel);
        } while (status == GSL_CONTINUE && i < inNumMaxIter);

        // Store the results to be returned to the user (copy from gsl_vector to result structure)
        for (size_t i = 0; i < FitTraitsT::CurveParamsT::_Count; ++i)
        {
            typename FitTraitsT::CurveParamsT::TypeE idx = static_cast<typename FitTraitsT::CurveParamsT::TypeE>(i);
            (*outResults)[idx] = gsl_vector_get(solver->x, idx);
        }

        // Free GSL memory
        gsl_multifit_fdfsolver_free(solver);
        gsl_vector_free(guess);

        return status;
    }
};