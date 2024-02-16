#include <iostream>
#include <assert.h>

#include "cimg/CImg.h"
#include <CCfits/CCfits>

#include <functional>
#include <array>
#include <vector>

#include "type.hpp"

#include "gaussian.hpp"

#include "atom/log/loguru.hpp"
#include "atom/type/json.hpp"
using json = nlohmann::json;

using namespace std;
using namespace cimg_library;
using namespace CCfits;

void readFile(CImg<float> &inImg, const std::string &inFilename, long *outBitPix);

bool insideCircle(float inX /*pos of x*/, float inY /*pos of y*/, float inCenterX, float inCenterY, float inRadius);

void thresholdOtsu(const CImg<float> &inImg, long inBitPix, CImg<float> *outBinImg);

void getAndRemoveNeighbours(PixelPosT inCurPixelPos, PixelPosSetT *inoutWhitePixels, PixelPosListT *inoutPixelsToBeProcessed);

template <typename T>
void clusterStars(const CImg<T> &inImg, StarInfoListT *outStarInfos);

float calcIx2(const CImg<float> &img, int x);

float calcJy2(const CImg<float> &img, int y);

void calcIntensityWeightedCenter(const CImg<float> &inImg, float *outX, float *outY);

void calcSubPixelCenter(const CImg<float> &inImg, float *outX, float *outY, size_t inNumIter /*num iterations*/);

void calcCentroid(const CImg<float> &inImg, const FrameT &inFrame, PixSubPosT *outPixelPos, PixSubPosT *outSubPixelPos, size_t inNumIterations);

float calcHfd(const CImg<float> &inImage, unsigned int inOuterDiameter);

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

FrameT rectify(const FrameT &inFrame);

bool drawStarImage(const std::string &filename, int &star_index, double &total_hfd, bool is_save = false, bool is_plot = false)
{
    const unsigned int outerHfdDiameter = 21;
    StarInfoListT starInfos;
    std::vector<std::list<StarInfoT *>> starBuckets;
    CImg<float> img;
    long bitPix = 0;

    json ret;

    // Read file to CImg
    try
    {
        DLOG_F(INFO, "Opening file {} ...", filename);
        readFile(img, filename, &bitPix);
    }
    catch (const CCfits::FitsException &e)
    {
        LOG_F(ERROR, "Read FITS failed: {}", e.message());
        return false;
    }

    CImg<unsigned char> rgbImg(img.width(), img.height(), 1 /*depth*/, 3 /*3 channels - RGB*/);
    float min = img.min(), mm = img.max() - min;

    cimg_forXY(img, x, y)
    {
        int value = 255.0 * (img(x, y) - min) / mm;
        rgbImg(x, y, 0 /*red*/) = value;
        rgbImg(x, y, 1 /*green*/) = value;
        rgbImg(x, y, 2 /*blue*/) = value;
    }

    // AD noise reduction --> In: Loaded image, Out: Noise reduced image
    // NOTE: This step takes a while for big images... too long for usage in a loop ->
    //       Should only be used on image segments, later...
    //
    // http://cimg.sourceforge.net/reference/structcimg__library_1_1CImg.html
    CImg<float> &aiImg = img.blur_anisotropic(30.0f, /*amplitude*/
                                              0.7f,  /*sharpness*/
                                              0.3f,  /*anisotropy*/
                                              0.6f,  /*alpha*/
                                              1.1f,  /*sigma*/
                                              0.8f,  /*dl*/
                                              30,    /*da*/
                                              2,     /*gauss_prec*/
                                              0,     /*interpolation_type*/
                                              false  /*fast_approx*/
    );

    // Thresholding (Otsu) --> In: Noise reduced image, Out: binary image
    CImg<float> binImg;
    thresholdOtsu(aiImg, bitPix, &binImg);

    // Clustering --> In: binary image from thresholding, Out: List of detected stars, subimg-boundaries (x1,y1,x2,y2) for each star
    clusterStars(binImg, &starInfos);

    DLOG_F(INFO, "Found {} stars", starInfos.size());

    ret["star_index"] = starInfos.size();

    // Calc brightness boundaries for possible focusing stars
    float maxPossiblePixValue = std::pow(2.0, bitPix) - 1;

    // For each star
    for (StarInfoListT::iterator it = starInfos.begin(); it != starInfos.end(); ++it)
    {
        const FrameT &frame = it->clusterFrame;
        FrameT &cogFrame = it->cogFrame;
        FrameT &hfdFrame = it->hfdFrame;
        PixSubPosT &cogCentroid = it->cogCentroid;
        PixSubPosT &subPixelInterpCentroid = it->subPixelInterpCentroid;
        float &hfd = it->hfd;
        float &fwhmHorz = it->fwhmHorz;
        float &fwhmVert = it->fwhmVert;
        float &maxPixValue = it->maxPixValue;
        bool &saturated = it->saturated;

        FrameT squareFrame = rectify(frame);

        // Centroid calculation --> In: Handle to full noise reduced image, subimg-boundaries (x1,y1,x2,y2), Out: (x,y) - abs. centroid coordinates
        calcCentroid(aiImg, squareFrame, &cogCentroid, &subPixelInterpCentroid, 10 /* num iterations */);
        std::get<0>(cogCentroid) += std::get<0>(squareFrame);
        std::get<1>(cogCentroid) += std::get<1>(squareFrame);
        std::get<0>(subPixelInterpCentroid) += std::get<0>(squareFrame);
        std::get<1>(subPixelInterpCentroid) += std::get<1>(squareFrame);

        // Calculate cog boundaries
        float maxClusterEdge = std::max(std::fabs(std::get<0>(frame) - std::get<2>(frame)), std::fabs(std::get<1>(frame) - std::get<3>(frame)));
        float cogHalfEdge = std::ceil(maxClusterEdge / 2.0);
        float cogX = std::get<0>(cogCentroid);
        float cogY = std::get<1>(cogCentroid);
        std::get<0>(cogFrame) = cogX - cogHalfEdge - 1;
        std::get<1>(cogFrame) = cogY - cogHalfEdge - 1;
        std::get<2>(cogFrame) = cogX + cogHalfEdge + 1;
        std::get<3>(cogFrame) = cogY + cogHalfEdge + 1;

        // HFD calculation --> In: image, Out: HFD value
        // Subtract mean value from image which is required for HFD calculation
        size_t hfdRectDist = std::floor(outerHfdDiameter / 2.0);
        std::get<0>(hfdFrame) = cogX - hfdRectDist;
        std::get<1>(hfdFrame) = cogY - hfdRectDist;
        std::get<2>(hfdFrame) = cogX + hfdRectDist;
        std::get<3>(hfdFrame) = cogY + hfdRectDist;

        CImg<float> hfdSubImg = aiImg.get_crop(std::get<0>(hfdFrame), std::get<1>(hfdFrame), std::get<2>(hfdFrame), std::get<3>(hfdFrame));
        maxPixValue = hfdSubImg.max();
        // saturated = (maxPixValue > lowerBound && maxPixValue < upperBound);
        saturated = (maxPixValue == maxPossiblePixValue);

        CImg<float> imgHfdSubMean(hfdSubImg);
        double mean = hfdSubImg.mean();

        cimg_forXY(hfdSubImg, x, y)
        {
            imgHfdSubMean(x, y) = (hfdSubImg(x, y) < mean ? 0 : hfdSubImg(x, y) - mean);
        }

        // Calc the HFD
        hfd = calcHfd(imgHfdSubMean, outerHfdDiameter /*outer diameter in px*/);

        total_hfd += hfd;

        // FWHM calculation --> In: Handle to full noise reduced image, abs. centroid coordinates, Out: FWHM value
        MyDataContainerT vertDataPoints, horzDataPoints;

        cimg_forX(imgHfdSubMean, x)
        {
            horzDataPoints.push_back(make_pair(x, imgHfdSubMean(x, std::floor(imgHfdSubMean.height() / 2.0 + 0.5))));
        }
        cimg_forY(imgHfdSubMean, y)
        {
            vertDataPoints.push_back(make_pair(y, imgHfdSubMean(std::floor(imgHfdSubMean.width() / 2.0 + 0.5), y)));
        }

        // Do the LM fit
        typedef CurveFitTmplT<GaussianFitTraitsT> GaussMatcherT;
        typedef GaussMatcherT::CurveParamsT CurveParamsT;
        CurveParamsT::TypeT gaussCurveParmsHorz, gaussCurveParmsVert;

        GaussMatcherT::fitGslLevenbergMarquart<MyDataAccessorT>(horzDataPoints, &gaussCurveParmsHorz, 0.1f /*EpsAbs*/, 0.1f /*EpsRel*/);
        fwhmHorz = gaussCurveParmsHorz[CurveParamsT::W_IDX];

        GaussMatcherT::fitGslLevenbergMarquart<MyDataAccessorT>(vertDataPoints, &gaussCurveParmsVert, 0.1f /*EpsAbs*/, 0.1f /*EpsRel*/);
        fwhmVert = gaussCurveParmsVert[CurveParamsT::W_IDX];
    }

    ret["hfd"] = total_hfd / starInfos.size();

    DLOG_F(INFO, "Finished processing {} stars", starInfos.size());
    DLOG_F(INFO, "{}", ret.dump(4));

    if (is_plot)
    {
        const int factor = 4;
        CImg<unsigned char> &rgbResized = rgbImg.resize(factor * rgbImg.width(), factor * rgbImg.height(),
                                                        -100 /*size_z*/, -100 /*size_c*/, 1 /*interpolation_type*/);

        // Draw cluster boundaries and square cluster boundaries
        const unsigned char red[3] = {255, 0, 0}, green[3] = {0, 255, 0}, yellow[3] = {255, 255, 0};
        const unsigned char black[3] = {0, 0, 0}, blue[3] = {0, 0, 255}, white[3] = {255, 255, 255};
        const size_t cCrossSize = 3;

        // Mark all stars in RGB image
        for (StarInfoListT::iterator it = starInfos.begin(); it != starInfos.end(); ++it)
        {
            StarInfoT *curStarInfo = &(*it);
            PixSubPosT &cogCentroid = curStarInfo->cogCentroid;
            float &hfd = curStarInfo->hfd;
            float &fwhmHorz = curStarInfo->fwhmHorz;
            float &fwhmVert = curStarInfo->fwhmVert;
            float &maxPixValue = curStarInfo->maxPixValue;

            std::string crood = fmt::format("{},{}", std::get<0>(curStarInfo->cogCentroid), std::get<1>(curStarInfo->cogCentroid));

            ret[crood] = {
                {"maxPix", maxPixValue}, {"sat", curStarInfo->saturated}, {"hfd", hfd}, {"fwhmHorz", fwhmHorz}, {"fwhmVert", fwhmVert}};

            const FrameT &frame = curStarInfo->clusterFrame;
            FrameT squareFrame(rectify(frame));
            rgbResized.draw_rectangle(std::floor(factor * (std::get<0>(frame) - 1) + 0.5), std::floor(factor * (std::get<1>(frame) - 1) + 0.5),
                                      std::floor(factor * (std::get<2>(frame) + 1) + 0.5), std::floor(factor * (std::get<3>(frame) + 1) + 0.5),
                                      red, 1 /*opacity*/, ~0 /*pattern*/);

            rgbResized.draw_rectangle(std::floor(factor * (std::get<0>(squareFrame) - 1) + 0.5), std::floor(factor * (std::get<1>(squareFrame) - 1) + 0.5),
                                      std::floor(factor * (std::get<2>(squareFrame) + 1) + 0.5), std::floor(factor * (std::get<3>(squareFrame) + 1) + 0.5),
                                      blue, 1 /*opacity*/, ~0 /*pattern*/);

            // Draw centroid crosses and centroid boundaries
            const PixSubPosT &subPos = curStarInfo->cogCentroid;
            const FrameT &cogFrame = curStarInfo->cogFrame;
            const FrameT &hfdFrame = curStarInfo->hfdFrame;

            rgbResized.draw_line(std::floor(factor * (std::get<0>(subPos) - cCrossSize) + 0.5), std::floor(factor * std::get<1>(subPos) + 0.5),
                                 std::floor(factor * (std::get<0>(subPos) + cCrossSize) + 0.5), std::floor(factor * std::get<1>(subPos) + 0.5), green, 1 /*opacity*/);

            rgbResized.draw_line(std::floor(factor * std::get<0>(subPos) + 0.5), std::floor(factor * (std::get<1>(subPos) - cCrossSize) + 0.5),
                                 std::floor(factor * std::get<0>(subPos) + 0.5), std::floor(factor * (std::get<1>(subPos) + cCrossSize) + 0.5), green, 1 /*opacity*/);

            rgbResized.draw_rectangle(std::floor(factor * std::get<0>(cogFrame) + 0.5), std::floor(factor * std::get<1>(cogFrame) + 0.5),
                                      std::floor(factor * std::get<2>(cogFrame) + 0.5), std::floor(factor * std::get<3>(cogFrame) + 0.5),
                                      green, 1 /*opacity*/, ~0 /*pattern*/);

            // Draw HFD
            rgbResized.draw_rectangle(std::floor(factor * std::get<0>(hfdFrame) + 0.5), std::floor(factor * std::get<1>(hfdFrame) + 0.5),
                                      std::floor(factor * std::get<2>(hfdFrame) + 0.5), std::floor(factor * std::get<3>(hfdFrame) + 0.5),
                                      yellow, 1 /*opacity*/, ~0 /*pattern*/);

            rgbImg.draw_circle(std::floor(factor * std::get<0>(subPos) + 0.5), std::floor(factor * std::get<1>(subPos) + 0.5), factor * outerHfdDiameter / 2, yellow, 1 /*pattern*/, 1 /*opacity*/);
            rgbImg.draw_circle(std::floor(factor * std::get<0>(subPos) + 0.5), std::floor(factor * std::get<1>(subPos) + 0.5), factor * hfd / 2, yellow, 1 /*pattern*/, 1 /*opacity*/);

            // Draw text
            const bool &saturated = curStarInfo->saturated;

            std::ostringstream oss;
            oss.precision(4);
            oss << "HFD=" << hfd << endl
                << "FWHM H=" << fwhmHorz << endl
                << "FWHM V=" << fwhmVert << endl
                << "MAX=" << (int)maxPixValue << endl
                << "SAT=" << (saturated ? "Y" : "N");

            rgbImg.draw_text(std::floor(factor * std::get<0>(subPos) + 0.5), std::floor(factor * std::get<1>(subPos) + 0.5), oss.str().c_str(), white /*fg color*/, black /*bg color*/, 0.7 /*opacity*/, 9 /*font-size*/);
        }
    }

    if (is_save)
    {
        rgbImg.save("out.png");
    }
    return true;
}

void readFile(CImg<float> &inImg, const std::string &inFilename, long *outBitPix = 0)
{
    std::unique_ptr<FITS> pInfile = std::make_unique<FITS>(inFilename, Read, true);
    PHDU &image = pInfile->pHDU();

    if (outBitPix)
    {
        *outBitPix = image.bitpix();
    }

    inImg.resize(image.axis(0) /*x*/, image.axis(1) /*y*/, 1 /*z*/, 1 /*1 color*/);

    // NOTE: At this point we assume that there is only 1 layer.
    // TODO: Color Image support
    std::valarray<unsigned long> imgData;
    image.read(imgData);
    cimg_forXY(inImg, x, y) { inImg(x, inImg.height() - y - 1) = imgData[inImg.offset(x, y)]; }
}

void writeFile(const CImg<float> &inImg, const string &inFilename)
{
    long naxes[2] = {inImg.width(), inImg.height()};
    auto pFits = make_unique<FITS>("!" + inFilename, USHORT_IMG, 2, naxes);

    long nelements = accumulate(begin(naxes), end(naxes), 1L, multiplies<long>());
    std::valarray<int> array(nelements);

    cimg_forXY(inImg, x, y)
    {
        array[inImg.offset(x, y)] = static_cast<int>(inImg(x, inImg.height() - y - 1));
    }

    long fpixel(1);
    pFits->pHDU().write(fpixel, nelements, array);
}

/**
 * Get all pixels inside a radius: http://stackoverflow.com/questions/14487322/get-all-pixel-array-inside-circle
 * Algorithm: http://en.wikipedia.org/wiki/Midpoint_circle_algorithm
 */
bool insideCircle(float inX /*pos of x*/, float inY /*pos of y*/, float inCenterX, float inCenterY, float inRadius)
{
    return (pow(inX - inCenterX, 2.0) + pow(inY - inCenterY, 2.0) <= pow(inRadius, 2.0));
}

void thresholdOtsu(const CImg<float> &inImg, long inBitPix, CImg<float> *outBinImg)
{
    CImg<> hist = inImg.get_histogram(std::pow(2.0, inBitPix));

    float sum = 0;
    cimg_forX(hist, pos) { sum += pos * hist[pos]; }

    float numPixels = inImg.width() * inImg.height();
    float sumB = 0, wB = 0, max = 0.0;
    float threshold1 = 0.0, threshold2 = 0.0;

    cimg_forX(hist, i)
    {
        wB += hist[i];

        if (!wB)
        {
            continue;
        }

        float wF = numPixels - wB;

        if (!wF)
        {
            break;
        }

        sumB += i * hist[i];

        float mF = (sum - sumB) / wF;
        float mB = sumB / wB;
        float diff = mB - mF;
        float bw = wB * wF * std::pow(diff, 2.0);

        if (bw >= max)
        {
            threshold1 = i;
            if (bw > max)
            {
                threshold2 = i;
            }
            max = bw;
        }
    } // end loop

    float th = (threshold1 + threshold2) / 2.0;

    *outBinImg = inImg; // Create a copy
    outBinImg->threshold(th);
}

/**
 * Removes all white neighbours arond pixel from whitePixels
 * if they exist and adds them to pixelsToBeProcessed.
 */
void getAndRemoveNeighbours(PixelPosT inCurPixelPos, PixelPosSetT *inoutWhitePixels, PixelPosListT *inoutPixelsToBeProcessed)
{
    const size_t _numPixels = 8, _x = 0, _y = 1;
    const int offsets[_numPixels][2] = {{-1, -1}, {0, -1}, {1, -1}, {-1, 0}, {1, 0}, {-1, 1}, {0, 1}, {1, 1}};

    for (size_t p = 0; p < _numPixels; ++p)
    {
        PixelPosT curPixPos(std::get<0>(inCurPixelPos) + offsets[p][_x], std::get<1>(inCurPixelPos) + offsets[p][_y]);
        PixelPosSetT::iterator itPixPos = inoutWhitePixels->find(curPixPos);

        if (itPixPos != inoutWhitePixels->end())
        {
            const PixelPosT &curPixPos = *itPixPos;
            inoutPixelsToBeProcessed->push_back(curPixPos);
            inoutWhitePixels->erase(itPixPos); // Remove white pixel from "white set" since it has been now processed
        }
    }
    return;
}

template <typename T>
void clusterStars(const CImg<T> &inImg, StarInfoListT *outStarInfos)
{
    PixelPosSetT whitePixels;

    cimg_forXY(inImg, x, y)
    {
        if (inImg(x, y))
        {
            whitePixels.insert(whitePixels.end(), PixelPosT(x, y));
        }
    }

    // Iterate over white pixels as long as set is not empty
    while (whitePixels.size())
    {
        PixelPosListT pixelsToBeProcessed;

        PixelPosSetT::iterator itWhitePixPos = whitePixels.begin();
        pixelsToBeProcessed.push_back(*itWhitePixPos);
        whitePixels.erase(itWhitePixPos);

        FrameT frame(inImg.width(), inImg.height(), 0, 0);

        while (!pixelsToBeProcessed.empty())
        {
            PixelPosT curPixelPos = pixelsToBeProcessed.front();

            // Determine boundaries (min max in x and y directions)
            if (std::get<0>(curPixelPos) /*x*/ < std::get<0>(frame) /*x1*/)
            {
                std::get<0>(frame) = std::get<0>(curPixelPos);
            }
            if (std::get<0>(curPixelPos) /*x*/ > std::get<2>(frame) /*x2*/)
            {
                std::get<2>(frame) = std::get<0>(curPixelPos);
            }
            if (std::get<1>(curPixelPos) /*y*/ < std::get<1>(frame) /*y1*/)
            {
                std::get<1>(frame) = std::get<1>(curPixelPos);
            }
            if (std::get<1>(curPixelPos) /*y*/ > std::get<3>(frame) /*y2*/)
            {
                std::get<3>(frame) = std::get<1>(curPixelPos);
            }

            getAndRemoveNeighbours(curPixelPos, &whitePixels, &pixelsToBeProcessed);
            pixelsToBeProcessed.pop_front();
        }

        // Create new star-info and set cluster-frame.
        // NOTE: we may use new to avoid copy of StarInfoT...
        StarInfoT starInfo;
        starInfo.clusterFrame = frame;
        outStarInfos->push_back(starInfo);
    }
}

float calcIx2(const CImg<float> &img, int x)
{
    float Ix = 0;
    cimg_forY(img, y) { Ix += std::pow(img(x, y), 2.0) * (float)x; }
    return Ix;
}

float calcJy2(const CImg<float> &img, int y)
{
    float Iy = 0;
    cimg_forX(img, x) { Iy += pow(img(x, y), 2.0) * (float)y; }
    return Iy;
}

// Calculate Intensity Weighted Center (IWC)
void calcIntensityWeightedCenter(const CImg<float> &inImg, float *outX, float *outY)
{
    assert(outX && outY);

    // Determine weighted centroid - See http://cdn.intechopen.com/pdfs-wm/26716.pdf
    float Imean2 = 0, Jmean2 = 0, Ixy2 = 0;

    for (size_t i = 0; i < inImg.width(); ++i)
    {
        Imean2 += calcIx2(inImg, i);
        cimg_forY(inImg, y) { Ixy2 += std::pow(inImg(i, y), 2.0); }
    }

    for (size_t i = 0; i < inImg.height(); ++i)
    {
        Jmean2 += calcJy2(inImg, i);
    }

    *outX = Imean2 / Ixy2;
    *outY = Jmean2 / Ixy2;
}

void calcSubPixelCenter(const CImg<float> &inImg, float *outX, float *outY, size_t inNumIter = 10 /*num iterations*/)
{
    // Sub pixel interpolation
    float c, a1, a2, a3, a4, b1, b2, b3, b4;
    float a1n, a2n, a3n, a4n, b1n, b2n, b3n, b4n;

    assert(inImg.width() == 3 && inImg.height() == 3);

    b1 = inImg(0, 0);
    a2 = inImg(1, 0);
    b2 = inImg(2, 0);
    a1 = inImg(0, 1);
    c = inImg(1, 1);
    a3 = inImg(2, 1);
    b4 = inImg(0, 2);
    a4 = inImg(1, 2);
    b3 = inImg(2, 2);

    for (size_t i = 0; i < inNumIter; ++i)
    {
        float c2 = 2 * c;
        float sp1 = (a1 + a2 + c2) / 4;
        float sp2 = (a2 + a3 + c2) / 4;
        float sp3 = (a3 + a4 + c2) / 4;
        float sp4 = (a4 + a1 + c2) / 4;

        // New maximum is center
        float newC = std::max({sp1, sp2, sp3, sp4});

        // Calc position of new center
        float ad = pow(2.0, -((float)i + 1));

        if (newC == sp1)
        {
            *outX = *outX - ad; // to the left
            *outY = *outY - ad; // to the top

            // Calculate new sub pixel values
            b1n = (a1 + a2 + 2 * b1) / 4;
            b2n = (c + b2 + 2 * a2) / 4;
            b3n = sp3;
            b4n = (b4 + c + 2 * a1) / 4;
            a1n = (b1n + c + 2 * a1) / 4;
            a2n = (b1n + c + 2 * a2) / 4;
            a3n = sp2;
            a4n = sp4;
        }
        else if (newC == sp2)
        {
            *outX = *outX + ad; // to the right
            *outY = *outY - ad; // to the top

            // Calculate new sub pixel values
            b1n = (2 * a2 + b1 + c) / 4;
            b2n = (2 * b2 + a3 + a2) / 4;
            b3n = (2 * a3 + b3 + c) / 4;
            b4n = sp4;
            a1n = sp1;
            a2n = (b2n + c + 2 * a2) / 4;
            a3n = (b2n + c + 2 * a3) / 4;
            a4n = sp3;
        }
        else if (newC == sp3)
        {
            *outX = *outX + ad; // to the right
            *outY = *outY + ad; // to the bottom

            // Calculate new sub pixel values
            b1n = sp1;
            b2n = (b2 + 2 * a3 + c) / 4;
            b3n = (2 * b3 + a3 + a4) / 4;
            b4n = (2 * a4 + b4 + c) / 4;
            a1n = sp4;
            a2n = sp2;
            a3n = (b3n + 2 * a3 + c) / 4;
            a4n = (b3n + 2 * a4 + c) / 4;
        }
        else
        {
            *outX = *outX - ad; // to the left
            *outY = *outY + ad; // to the bottom

            // Calculate new sub pixel values
            b1n = (2 * a1 + b1 + c) / 4;
            b2n = sp2;
            b3n = (c + b3 + 2 * a4) / 4;
            b4n = (2 * b4 + a1 + a4) / 4;
            a1n = (b4n + 2 * a1 + c) / 4;
            a2n = sp1;
            a3n = sp3;
            a4n = (b4n + 2 * a4 + c) / 4;
        }

        c = newC; // Oi = Oi+1

        a1 = a1n;
        a2 = a2n;
        a3 = a3n;
        a4 = a4n;

        b1 = b1n;
        b2 = b2n;
        b3 = b3n;
        b4 = b4n;
    }
}

void calcCentroid(const CImg<float> &inImg, const FrameT &inFrame, PixSubPosT *outPixelPos, PixSubPosT *outSubPixelPos = 0, size_t inNumIterations = 10)
{
    // Get frame sub img
    CImg<float> subImg = inImg.get_crop(std::get<0>(inFrame), std::get<1>(inFrame), std::get<2>(inFrame), std::get<3>(inFrame));

    float &xc = std::get<0>(*outPixelPos);
    float &yc = std::get<1>(*outPixelPos);

    // 1. Calculate the IWC
    calcIntensityWeightedCenter(subImg, &xc, &yc);

    if (outSubPixelPos)
    {
        // 2. Round to nearest integer and then iteratively improve.
        int xi = std::floor(xc + 0.5);
        int yi = std::floor(yc + 0.5);

        CImg<float> img3x3 = inImg.get_crop(xi - 1 /*x0*/, yi - 1 /*y0*/, xi + 1 /*x1*/, yi + 1 /*y1*/);

        // 3. Interpolate using sub-pixel algorithm
        float xsc = xi, ysc = yi;
        calcSubPixelCenter(img3x3, &xsc, &ysc, inNumIterations);

        std::get<0>(*outSubPixelPos) = xsc;
        std::get<1>(*outSubPixelPos) = ysc;
    }
}

float calcHfd(const CImg<float> &inImage, unsigned int inOuterDiameter)
{
    // Sum up all pixel values in whole circle
    float outerRadius = inOuterDiameter / 2;
    float sum = 0, sumDist = 0;
    int centerX = std::ceil(inImage.width() / 2.0);
    int centerY = std::ceil(inImage.height() / 2.0);

    cimg_forXY(inImage, x, y)
    {
        if (insideCircle(x, y, centerX, centerY, outerRadius))
        {
            sum += inImage(x, y);
            sumDist += inImage(x, y) * std::sqrt(std::pow((float)x - (float)centerX, 2.0f) + std::pow((float)y - (float)centerY, 2.0f));
        }
    }
    // NOTE: Multiplying with 2 is required since actually just the HFR is calculated above
    return (sum ? 2.0 * sumDist / sum : std::sqrt(2.0) * outerRadius);
}

FrameT
rectify(const FrameT &inFrame)
{
    float border = 3;
    float border2 = 2.0 * border;
    float width = std::fabs(std::get<0>(inFrame) - std::get<2>(inFrame)) + border2;
    float height = std::fabs(std::get<1>(inFrame) - std::get<3>(inFrame)) + border2;
    float L = std::max(width, height);
    float x0 = std::get<0>(inFrame) - (std::fabs(width - L) / 2.0) - border;
    float y0 = std::get<1>(inFrame) - (std::fabs(height - L) / 2.0) - border;
    return FrameT(x0, y0, x0 + L, y0 + L);
}
