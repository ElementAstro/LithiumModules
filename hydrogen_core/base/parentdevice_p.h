#pragma once

#include "basedevice_p.h"
#include <atomic>

namespace HYDROGEN
{

    class ParentDevicePrivate : public BaseDevicePrivate
    {
    public:
        std::atomic_int ref{0};
    };

}
