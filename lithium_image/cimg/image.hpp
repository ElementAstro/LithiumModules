/*
 * image.hpp
 *
 * Copyright (C) 2023-2024 Max Qian <lightapt.com>
 */

/*************************************************

Date: 2023-12-1

Description: Image processing plugin for Lithium

**************************************************/

#ifndef LITHIUM_IMAGE_CIMG_IMAGE_HPP
#define LITHIUM_IMAGE_CIMG_IMAGE_HPP

#include "atom/components/templates/shared_component.hpp"

#include "atom/search/cache.hpp"

#include "atom/type/json.hpp"
using json = nlohmann::json;

namespace cimg_library {

  // Declare the four classes of the CImg Library.
  template<typename T=float> struct CImg;
  template<typename T=float> struct CImgList;
  struct CImgDisplay;
  struct CImgException;
}

class ImageProcessingPlugin : public SharedComponent
{
public:
    ImageProcessingPlugin();

    ~ImageProcessingPlugin();

    static std::shared_ptr<ImageProcessingPlugin> createShared();

    json _drawStarImage(const json &params);

private:

    void readFile(const std::string &key,const std::string &filename, long *outBitPix);

    void writeFile(const std::string &key, const std::string &filename);

    bool readImage(const std::string &key, const std::string &filename);

    bool readColorImage(const std::string &key, const std::string &filename);

    bool saveImage(const std::string &key, const std::string &filename);

    std::string imageToBase64(const std::string &key);

    void base64ToImage(const std::string &key, const std::string &img);

    bool calcDarkNoise(const std::string &key, float &average_dark, float &sigma_dark, float &sigma_readout);

    bool detectStars(const std::string &key, int threshold, int max_radius);

    bool whiteBalance(const std::string &key);

    void blur(const std::string &key, const std::string& mode, float param);

    void crop(const std::string &key, const int &x, const int &y, const int &w, const int &h);

    void compressImage(const std::string& key, int compress_ratio);

private:
    // Image cache in memory
    mutable std::unique_ptr<ResourceCache<std::shared_ptr<cimg_library::CImg<float>>>> m_cache;
};

#endif
