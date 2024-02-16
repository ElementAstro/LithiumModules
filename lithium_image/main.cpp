/*
 * main.cpp
 *
 * Copyright (C) 2023-2024 Max Qian <lightapt.com>
 */

/*************************************************

Date: 2023-12-1

Description: Image processing plugin main

**************************************************/

#include "config.h"

#if ENABLE_CIMG
#include "cimg/image.hpp"
#endif
#if ENABLE_OPENCV
#include "opencv/image.hpp"
#endif
#include <memory>

extern "C"
{
    std::shared_ptr<ImageProcessingPlugin> GetInstance()
    {
        return std::make_shared<ImageProcessingPlugin>();
    }

    json GetInfo()
    {
        json config;
        config["name"] = "lithium_image";
        config["version"] = "1.0.0";
        config["description"] = "Image processing plugin";
        config["author"] = "Max Qian";
        config["email"] = "astro_air@126.com";
        config["url"] = "lightapt.com";
        config["license"] = "GPLv3";
        config["copyright"] = "2023 Max Qian. All rights reserved";
        return config;
    }
}