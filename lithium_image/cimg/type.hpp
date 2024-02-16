#pragma once

#include <list>
#include <set>
#include <list>
#include <vector>

typedef std::tuple<int /*x*/, int /*y*/> PixelPosT;
typedef std::set<PixelPosT> PixelPosSetT;
typedef std::list<PixelPosT> PixelPosListT;

typedef std::tuple<float, float> PixSubPosT;
typedef std::tuple<float /*x1*/, float /*y1*/, float /*x2*/, float /*y2*/> FrameT;

struct StarInfoT
{
  FrameT clusterFrame;
  FrameT cogFrame;
  FrameT hfdFrame;
  PixSubPosT cogCentroid;
  PixSubPosT subPixelInterpCentroid;
  float hfd;
  float fwhmHorz;
  float fwhmVert;
  float maxPixValue;
  bool saturated;
};
typedef std::list<StarInfoT> StarInfoListT;

struct DataPointT
{
    float x;
    float y;
    DataPointT(float inX = 0, float inY = 0) : x(inX), y(inY) {}
};

typedef std::vector<DataPointT> DataPointsT;

struct GslMultiFitDataT
{
    float y;
    float sigma;
    DataPointT pt;
};

typedef std::vector<GslMultiFitDataT> GslMultiFitParmsT;

typedef std::list<PixSubPosT> MyDataContainerT;

class MyDataAccessorT
{
public:
    typedef MyDataContainerT TypeT;
    static DataPointT getDataPoint(size_t inIdx, TypeT::const_iterator inIt)
    {
        const PixSubPosT &pos = *inIt;
        DataPointT dp(std::get<0>(pos) /*inIdx*/, std::get<1>(pos) /*y*/);
        return dp;
    }
};