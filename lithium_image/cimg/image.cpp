/*
 * image->cpp
 *
 * Copyright (C) 2023-2024 Max Qian <lightapt.com>
 */

/*************************************************

Date: 2023-12-1

Description: Image processing plugin for Lithium

**************************************************/

#include "image.hpp"
#include "draw.hpp"

#include <fstream>
#include <sstream>

#include <CCfits/CCfits>
#include "cimg/CImg.h"

#include "atom/algorithm/base64.hpp"
#include "atom/log/loguru.hpp"

ImageProcessingPlugin::ImageProcessingPlugin()
{
    m_cache = std::make_unique<ResourceCache<std::shared_ptr<cimg_library::CImg<float>>>>();
    LOG_F(INFO, "ImageProcessingPlugin()");
    RegisterFunc("drawStarImage", &ImageProcessingPlugin::_drawStarImage, this);
}

ImageProcessingPlugin::~ImageProcessingPlugin()
{
    // Do nothing
}

json ImageProcessingPlugin::_drawStarImage(const json &params)
{
    DLOG_F(INFO, "_drawStarImage(): {}", params.dump());
    json result;
    if (!params.contains("filename"))
    {
        LOG_F(ERROR, "Missing parameter: filename");
        result["error"] = "Missing parameter: filename";
        return result;
    }
    if (!params.contains("threshold"))
    {
        LOG_F(ERROR, "Missing parameter: threshold");
        result["error"] = "Missing parameter: threshold";
        return result;
    }
    std::string filename = params["filename"].get<std::string>();
    int threshold = params["threshold"].get<int>();
    int max_radius = params["max_radius"].get<int>();
    bool is_save = params["is_save"].get<bool>();
    bool is_plot = params["is_plot"].get<bool>();

    int star_index = 0;
    double total_hfd = 0;

    if (!drawStarImage(filename, star_index, total_hfd, is_save, is_plot))
    {
        LOG_F(ERROR, "Failed to draw star image");
        return;
    }
    LOG_F(INFO, "Star index: {}, HFD: {}", star_index, total_hfd);
}

void ImageProcessingPlugin::readFile(const std::string &key, const std::string &filename, long *outBitPix = 0)
{
    std::unique_ptr<CCfits::FITS> pInfile = std::make_unique<CCfits::FITS>(filename, CCfits::Read, true);
    CCfits::PHDU &image = pInfile->pHDU();

    if (outBitPix)
    {
        *outBitPix = image.bitpix();
    }

    std::shared_ptr<cimg_library::CImg<float>> inImg = std::make_shared<cimg_library::CImg<float>>();

    inImg->resize(image.axis(0) /*x*/, image.axis(1) /*y*/, 1 /*z*/, 1 /*1 color*/);

    // NOTE: At this point we assume that there is only 1 layer.
    // TODO: Color Image support
    std::valarray<unsigned long> imgData;
    image.read(imgData);
    cimg_forXY((*inImg), x, y) { (*inImg)(x, inImg->height() - y - 1) = imgData[inImg->offset(x, y)]; }

    // Load into cache
    // The image should be processed successfully within 120 seconds
    m_cache.insert(key, inImg, std::chrono::seconds(120));
}->
void ImageProcessingPlugin::writeFile(const std::string &key, const string &filename)
{
    if (key.empty() || filename.empty())
        return;

    std::shared_ptr<cimg_library::CImg<float>> inImg = m_cache->get(key);

    if (!inImg)
    {
        return;
    }
    long naxes[2] = {inImg->width(), inImg->height()};
    auto pFits = std::make_unique<CCfits::FITS>(filename, USHORT_IMG, 2, naxes);

    long nelements = std::accumulate(std::begin(naxes), std::end(naxes), 1L, std::multiplies<long>());
    std::valarray<int> array(nelements);

    cimg_forXY((*inImg), x, y)
    {
        array[inImg->offset(x, y)] = static_cast<int>((*inImg)(x, inImg->height() - y - 1));
    }

    long fpixel(1);
    pFits->pHDU().write(fpixel, nelements, array);
}

bool ImageProcessingPlugin::readImage(const std::string &key, const std::string &filename)
{
    if (key.empty() || filename.empty())
        return;
    if (std::ifstream(filename))
    {
        std::shared_ptr<cimg_library::CImg<float>> img = std::make_shared<cimg_library::CImg<float>>();
        try
        {
            img->load(filename.c_str());
        }
        catch (const cimg_library::CImgIOException &e)
        {
            LOG_F(ERROR, "Error reading image file: {}", e.what());
            return false;
        }
        if (img)
        {
            m_cache->insert(key, img, std::chrono::seconds(120));
        }
        else
        {
        }
    }
    else
    {
        LOG_F(ERROR, "Error: image file does not exist");
        return false;
    }
    return true;
}

bool ImageProcessingPlugin::readColorImage(const std::string &key, const std::string &filename)
{
    if (key.empty() || filename.empty())
        return;
    if (std::ifstream(filename))
    {
        std::shared_ptr<cimg_library::CImg<float>> img = std::make_shared<cimg_library::CImg<float>>();
        try
        {
            img->load(filename.c_str());
            if (img->spectrum() == 1)
            {
                return false;
            }
        }
        catch (const cimg_library::CImgIOException &e)
        {
            LOG_F(ERROR, "Error reading image file: {}", e.what());
            return false;
        }
        if (img)
        {
            m_cache->insert(key, img, std::chrono::seconds(120));
        }
        else
        {
        }
    }
    else
    {
        LOG_F(ERROR, "Error: image file does not exist");
        return false;
    }
    return true;
}

bool ImageProcessingPlugin::saveImage(const std::string &key, const std::string &filename)
{
    if (key.empty() || filename.empty())
        return;
    try
    {
        std::shared_ptr<cimg_library::CImg<float>> img = m_cache->get(key);
        if (!img)
        {
            return false;
        }
        img->save(filename.c_str());
    }
    catch (const cimg_library::CImgIOException &e)
    {
        LOG_F(ERROR, "Error writing image file: {}", e.what());
        return false;
    }
    return true;
}

std::string ImageProcessingPlugin::imageToBase64(const std::string &key)
{
    // 读取图片
    try
    {
        std::shared_ptr<cimg_library::CImg<float>> img = m_cache->get(key);
        if (!img)
        {
            return "";
        }
        // 将图片数据转为字符串
        std::ostringstream oss;
        // image.save(oss, "jpg");
        std::string imageData = oss.str();

        // 将图片数据进行Base64编码
    }
    catch (const cimg_library::CImgIOException &e)
    {
        LOG_F(ERROR, "Error convert image to base64: {}", e.what());
        return "";
    }
}

void ImageProcessingPlugin::base64ToImage(const std::string &key, const std::string &img)
{
    // std::string imageData = Atom::Utils::base64Decode(img);

    // 创建CImg对象并保存为图片文件
    // std::istringstream iss(imageData);
    // CImg<unsigned char> image;
    // image.load(iss);
}

bool ImageProcessingPlugin::calcDarkNoise(const std::string &key, float &average_dark, float &sigma_dark, float &sigma_readout)
{
    std::shared_ptr<cimg_library::CImg<float>> img = m_cache->get(key);
    cimg_library::CImg<float> dark_img = *img;
    if (!img)
    {
        return false;
    }
    // 计算均值
    average_dark = 0;
    cimg_forXY(dark_img, x, y)
    {
        average_dark += dark_img(x, y);
    }
    average_dark /= (dark_img.width() * dark_img.height());

    // 计算标准差
    sigma_dark = 0;
    sigma_readout = 0;
    cimg_forXY(dark_img, x, y)
    {
        sigma_dark += pow(dark_img(x, y) - average_dark, 2);
        if (x < dark_img.width() - 1)
        {
            sigma_readout += pow(dark_img(x, y) - dark_img(x + 1, y), 2);
        }
        if (y < dark_img.height() - 1)
        {
            sigma_readout += pow(dark_img(x, y) - dark_img(x, y + 1), 2);
        }
    }
    sigma_dark = sqrt(sigma_dark / (dark_img.width() * dark_img.height()));
    sigma_readout = sqrt(sigma_readout / (2 * (dark_img.width() - 1) * dark_img.height() + 2 * dark_img.width() * (dark_img.height() - 1)));
}

bool ImageProcessingPlugin::detectStars(const std::string &key, int threshold, int max_radius)
{
    // 读取图像并记录日志
    DLOG_F(INFO, "Loading image: {}", key);
    std::shared_ptr<cimg_library::CImg<float>> img = m_cache->get(key);
    if (!img)
    {
        return false;
    }

    // 转换成灰度图像
    cimg_library::CImg<unsigned char> gray = img->get_RGBtoYCbCr().get_channel(0);

    // 阈值处理
    cimg_library::CImg<unsigned char> binary = gray.threshold(threshold);

    // 星点检测
    cimg_library::CImg<unsigned char> stars(binary.width(), binary.height(), 1, 1, 0);
    int count = 0;
    cimg_forXY(binary, x, y)
    {
        if (binary(x, y) == 0)
        {
            for (int r = 1; r <= max_radius; r++)
            {
                bool is_star = true;
                for (int t = 0; t < stars.spectrum(); t++)
                {
                    int tx = x + r * std::cos(t * cimg_library::cimg::PI / 4);
                    int ty = y + r * std::sin(t * cimg_library::cimg::PI / 4);
                    if (tx < 0 || tx >= binary.width() || ty < 0 || ty >= binary.height() || binary(tx, ty) != 0)
                    {
                        is_star = false;
                        break;
                    }
                }
                if (is_star)
                {
                    const unsigned char red[] = {255, 0, 0};
                    stars.draw_circle(x, y, r, red, 1);
                    count++;
                }
            }
        }
    }
    DLOG_F(INFO, "Finished detecting {} stars in image: {}", count, key);
    return true;
}

bool ImageProcessingPlugin::whiteBalance(const std::string &key)
{
    std::shared_ptr<cimg_library::CImg<float>> img = m_cache->get(key);
    if (!img)
    {
        return false;
    }

    // 计算每个通道的平均值
    double r = 0, g = 0, b = 0;
    cimg_forXY(*img, x, y)
    {
        r += (*img)(x, y, 0);
        g += (*img)(x, y, 1);
        b += (*img)(x, y, 2);
    }
    int size = img->width() * img->height();
    r /= size;
    g /= size;
    b /= size;

    // 调整每个通道的白平衡
    cimg_forXY((*img), x, y)
    {
        double factor_r = r / (*img)(x, y, 0);
        double factor_g = g / (*img)(x, y, 1);
        double factor_b = b / (*img)(x, y, 2);
        (*img)(x, y, 0) = cimg_library::cimg::cut((*img)(x, y, 0) * factor_r, 0, 255);
        (*img)(x, y, 1) = cimg_library::cimg::cut((*img)(x, y, 1) * factor_g, 0, 255);
        (*img)(x, y, 2) = cimg_library::cimg::cut((*img)(x, y, 2) * factor_b, 0, 255);
    }

    m_cache->insert(key, img , std::chrono::seconds(120));
    return true;
}

void ImageProcessingPlugin::blur(const std::string &key, const std::string& mode, float param)
{
    std::shared_ptr<cimg_library::CImg<float>> img = m_cache->get(key);
    if (!img)
    {
        return;
    }

    // 模糊处理
    if (mode == "gaussian")
    {
        img->get_blur(param);
    }
    else if (mode == "median")
    {
        img->blur_median(param);
    }
    else if (mode == "mean")
    {
        img->get_blur(static_cast<int>(param));
    }
    m_cache->insert(key, img , std::chrono::seconds(120));
}

void ImageProcessingPlugin::crop(const std::string &key, const int &x, const int &y, const int &w, const int &h)
{
    std::shared_ptr<cimg_library::CImg<float>> img = m_cache->get(key);
    if (!img)
    {
        return;
    }

    // 裁剪图像
    img->crop(x, y, w, h);
    m_cache->insert(key, img , std::chrono::seconds(120));
}

void ImageProcessingPlugin::compressImage(const std::string& key, int compress_ratio)
{
    // 计算压缩后的图像宽度和高度，并创建新图像
    cimg_library::CImg<unsigned char> img;
    int new_width = img.width() / compress_ratio, new_height = img.height() / compress_ratio;
    cimg_library::CImg<unsigned char> new_img(new_width, new_height, 1, img.spectrum());

    // 对每个新像素点，计算与原图像对应的多个像素的平均值，并输出调试日志
    DLOG_F(INFO, "Compress the image with ratio {}.", compress_ratio);
    cimg_forXY(new_img, x, y)
    {
        int sum_r = 0, sum_g = 0, sum_b = 0, count = 0;
        for (int i = 0; i < compress_ratio; i++)
        {
            for (int j = 0; j < compress_ratio; j++)
            {
                int px = x * compress_ratio + i, py = y * compress_ratio + j;
                if (px < img.width() && py < img.height())
                {
                    sum_r += img(px, py, 0, 0);
                    if (img.spectrum() > 1)
                    {
                        sum_g += img(px, py, 0, 1);
                        sum_b += img(px, py, 0, 2);
                    }
                    count++;
                }
            }
        }
        new_img(x, y, 0, 0) = sum_r / count;
        if (img.spectrum() > 1)
        {
            new_img(x, y, 0, 1) = sum_g / count;
            new_img(x, y, 0, 2) = sum_b / count;
        }
    }

    img = new_img; // 将压缩后的图像覆盖原始图像
}