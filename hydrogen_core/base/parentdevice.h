#pragma once

#include "basedevice.h"

/**
 * @class HYDROGEN::ParentDevice
 * @brief The class is used to create device instances.
 * Class copying is not allowed. When an object is destroyed,
 * the property list (HYDROGEN::Property) is cleared to prevent a circular reference along with the properties.
 * The base HYDROGEN::BaseDevice class and its HYDROGEN::Properties exist as long as they are used by other objects.
 */

namespace HYDROGEN
{

    class ParentDevicePrivate;
    class ParentDevice : public BaseDevice
    {
        DECLARE_PRIVATE(ParentDevice)

        ParentDevice(const ParentDevice &) = delete;
        ParentDevice &operator=(const ParentDevice &) = delete;

    public:
        enum Type
        {
            Valid,
            Invalid
        };

    public:
        explicit ParentDevice(Type type);
        ~ParentDevice();

    public:
        ParentDevice(ParentDevice &&other) = default;
        ParentDevice &operator=(ParentDevice &&other) = default;

    protected:
        ParentDevice(const std::shared_ptr<ParentDevicePrivate> &dd);
    };

}
