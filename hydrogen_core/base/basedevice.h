/*******************************************************************************
  Copyright(c) 2011 Jasem Mutlaq. All rights reserved.

 This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Library General Public
 License version 2 as published by the Free Software Foundation.

 This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Library General Public License for more details.

 You should have received a copy of the GNU Library General Public License
 along with this library; see the file COPYING.LIB.  If not, write to
 the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 Boston, MA 02110-1301, USA.
*******************************************************************************/

#pragma once

#include "hydrogenbase.h"
#include "hydrogenproperty.h"
#include "hydrogenproperties.h"

#include <string>
#include <vector>
#include <cstdint>

#include "hydrogenpropertytext.h"
#include "hydrogenpropertynumber.h"
#include "hydrogenpropertyswitch.h"
#include "hydrogenpropertylight.h"
#include "hydrogenpropertyblob.h"

// #define MAXRBUF 2048 // #PS: defined in indibase.h

/** @class HYDROGEN::BaseDevice
 *  @brief Class to provide basic INDI device functionality.
 *
 *  HYDROGEN::BaseDevice is the base device for all INDI devices and contains a list of all properties defined by the device either explicity or via a skeleton file.
 *  You don't need to subclass HYDROGEN::BaseDevice class directly, it is inheritied by HYDROGEN::DefaultDevice which takes care of building a standard INDI device. Moreover, HYDROGEN::BaseClient
 *  maintains a list of HYDROGEN::BaseDevice objects as they get defined from the INDI server, and those objects may be accessed to retrieve information on the object properties or message log.
 *
 *  @author Jasem Mutlaq
 */
namespace HYDROGEN
{

    class LilXmlElement;
    class BaseDevicePrivate;
    class BaseDevice
    {
        DECLARE_PRIVATE(BaseDevice)
    public:
        typedef HYDROGEN::Properties Properties;
        // typedef std::vector<HYDROGEN::Property*> Properties;

        /*! INDI error codes. */
        enum HYDROGEN_ERROR
        {
            HYDROGEN_DEVICE_NOT_FOUND = -1,    /*!< INDI Device was not found. */
            HYDROGEN_PROPERTY_INVALID = -2,    /*!< Property has an invalid syntax or attribute. */
            HYDROGEN_PROPERTY_DUPLICATED = -3, /*!< INDI Device was not found. */
            HYDROGEN_DISPATCH_ERROR = -4       /*!< Dispatching command to driver failed. */
        };

        /*! Used for switch Enabled/Disabled or On/Off type properties */
        enum
        {
            HYDROGEN_ENABLED,
            HYDROGEN_DISABLED
        };

        /*! Used for watchProperty callback method. */
        enum WATCH
        {
            WATCH_NEW = 0,      /*!< Applies to discovered properties only. */
            WATCH_UPDATE,       /*!< Applies to updated properties only. */
            WATCH_NEW_OR_UPDATE /*!< Applies when a property appears or is updated, i.e. both of the above. */
        };

        /** @brief The DRIVER_INTERFACE enum defines the class of devices the driver implements. A driver may implement one or more interfaces. */
        enum DRIVER_INTERFACE
        {
            GENERAL_INTERFACE = 0,              /**< Default interface for all INDI devices */
            TELESCOPE_INTERFACE = (1 << 0),     /**< Telescope interface, must subclass HYDROGEN::Telescope */
            CCD_INTERFACE = (1 << 1),           /**< CCD interface, must subclass HYDROGEN::CCD */
            GUIDER_INTERFACE = (1 << 2),        /**< Guider interface, must subclass HYDROGEN::GuiderInterface */
            FOCUSER_INTERFACE = (1 << 3),       /**< Focuser interface, must subclass HYDROGEN::FocuserInterface */
            FILTER_INTERFACE = (1 << 4),        /**< Filter interface, must subclass HYDROGEN::FilterInterface */
            DOME_INTERFACE = (1 << 5),          /**< Dome interface, must subclass HYDROGEN::Dome */
            GPS_INTERFACE = (1 << 6),           /**< GPS interface, must subclass HYDROGEN::GPS */
            WEATHER_INTERFACE = (1 << 7),       /**< Weather interface, must subclass HYDROGEN::Weather */
            AO_INTERFACE = (1 << 8),            /**< Adaptive Optics Interface */
            DUSTCAP_INTERFACE = (1 << 9),       /**< Dust Cap Interface */
            LIGHTBOX_INTERFACE = (1 << 10),     /**< Light Box Interface */
            DETECTOR_INTERFACE = (1 << 11),     /**< Detector interface, must subclass HYDROGEN::Detector */
            ROTATOR_INTERFACE = (1 << 12),      /**< Rotator interface, must subclass HYDROGEN::RotatorInterface */
            SPECTROGRAPH_INTERFACE = (1 << 13), /**< Spectrograph interface */
            CORRELATOR_INTERFACE = (1 << 14),   /**< Correlators (interferometers) interface */
            AUX_INTERFACE = (1 << 15),          /**< Auxiliary interface */

            SENSOR_INTERFACE = SPECTROGRAPH_INTERFACE | DETECTOR_INTERFACE | CORRELATOR_INTERFACE
        };

    public:
        BaseDevice();
        virtual ~BaseDevice();

    public: // property
        /** @brief Register the property to be able to observe and update.
         *  @param property any property from the HYDROGEN::PropertyXXX family.
         */
        void registerProperty(const HYDROGEN::Property &property);
        void registerProperty(const HYDROGEN::Property &property,
                              HYDROGEN_PROPERTY_TYPE type); // backward compatiblity (PentaxCCD, PkTriggerCordCCD)

        /** @brief Remove a property
         *  @param name name of property to be removed. Pass NULL to remove the whole device.
         *  @param errmsg buffer to store error message.
         *  @return 0 if successul, -1 otherwise.
         */
        int removeProperty(const char *name, char *errmsg);

        /** @brief Call the callback function if property is available.
         *  @param name of property.
         *  @param callback as an argument of the function you can use HYDROGEN::PropertyNumber, HYDROGEN::PropertySwitch etc.
         *  @param watch you can decide whether the callback should be executed only once (WATCH_NEW) on discovery of the property or
         *  also on every change of the value (WATCH_UPDATE) or both (WATCH_NEW_OR_UPDATE)
         *  @note the behavior is analogous to BaseMediator::newProperty/updateProperty
         */
        void watchProperty(const char *name, const std::function<void(HYDROGEN::Property)> &callback, WATCH watch = WATCH_NEW);

        /** @brief Return a property and its type given its name.
         *  @param name of property to be found.
         *  @param type of property found.
         *  @return If property is found, it is returned. To be used you must use static_cast with given the type of property
         *  returned.
         */
        Property getProperty(const char *name, HYDROGEN_PROPERTY_TYPE type = HYDROGEN_UNKNOWN) const;

        /** @brief Return a list of all properties in the device. */
        Properties getProperties();
        const Properties getProperties() const;

    public:
        /** @return Return vector number property given its name */
        HYDROGEN::PropertyNumber getNumber(const char *name) const;
        /** @return Return vector text property given its name */
        HYDROGEN::PropertyText getText(const char *name) const;
        /** @return Return vector switch property given its name */
        HYDROGEN::PropertySwitch getSwitch(const char *name) const;
        /** @return Return vector light property given its name */
        HYDROGEN::PropertyLight getLight(const char *name) const;
        /** @return Return vector BLOB property given its name */
        HYDROGEN::PropertyBlob getBLOB(const char *name) const;

    public: // deprecated
        /** @return Return property state */
        IPState getPropertyState(const char *name) const;
        /** @return Return property permission */
        IPerm getPropertyPermission(const char *name) const;

        /** @brief Return a property and its type given its name.
         *  @param name of property to be found.
         *  @param type of property found.
         *  @return If property is found, the raw void * pointer to the IXXXVectorProperty is returned. To be used you must use static_cast with given the type of property
         *  returned. For example, INumberVectorProperty *num = static_cast<INumberVectorProperty> getRawProperty("FOO", HYDROGEN_NUMBER);
         *
         *  @note This is a low-level function and should not be called directly unless necessary. Use getXXX instead where XXX
         *  is the property type (Number, Text, Switch..etc).
         */
        void *getRawProperty(const char *name, HYDROGEN_PROPERTY_TYPE type = HYDROGEN_UNKNOWN) const;

    public:
        /** @brief Add message to the driver's message queue.
         *  @param msg Message to add.
         */
        void addMessage(const std::string &msg);
        void checkMessage(XMLEle *root);
        void doMessage(XMLEle *msg);

        /** @return Returns a specific message. */
        const std::string &messageQueue(size_t index) const;

        /** @return Returns last message message. */
        const std::string &lastMessage() const;

    public:
        /** @return True if the device exists */
        bool isValid() const;

        /** @return True if the device is connected (CONNECT=ON), False otherwise */
        bool isConnected() const;

        /** @brief indicates that the device is ready
         *  @note the BaseMediator::newDevice method will be called
         */
        void attach();

        /** @brief indicates that the device is being removed
         *  @note the BaseMediator::removeDevice method will be called
         */
        void detach();

        /** @brief Set the driver's mediator to receive notification of news devices and updated property values. */
        void setMediator(HYDROGEN::BaseMediator *mediator);

        /** @returns Get the meditator assigned to this driver */
        HYDROGEN::BaseMediator *getMediator() const;

        /** @brief Set the device name
         *  @param dev new device name
         */
        void setDeviceName(const char *dev);

        /** @return Returns the device name */
        const char *getDeviceName() const;

        /** @brief Check that the device name matches the argument. **/
        bool isDeviceNameMatch(const char *otherName) const;

        /** @brief Check that the device name matches the argument. **/
        bool isDeviceNameMatch(const std::string &otherName) const;

        /** @return driver name
         *  @note This can only be valid if DRIVER_INFO is defined by the driver.
         */
        const char *getDriverName() const;

        /** @return driver executable name
         *  @note This can only be valid if DRIVER_INFO is defined by the driver.
         */
        const char *getDriverExec() const;

        /** @return driver version
         *  @note This can only be valid if DRIVER_INFO is defined by the driver.
         */
        const char *getDriverVersion() const;

        /** @brief getDriverInterface returns ORed values of @ref HYDROGEN::BaseDevice::DRIVER_INTERFACE "DRIVER_INTERFACE". It presents the device classes supported by the driver.
         *  @return driver device interface descriptor.
         *  @note For example, to know if the driver supports CCD interface, check the retruned value:
         *  @code{.cpp}
         *  if (device.getDriverInterface() & CCD_INTERFACE)
         *       cout << "We received a camera!" << endl;
         *  @endcode
         */
        uint16_t getDriverInterface() const;

    public:
        /** @brief Build driver properties from a skeleton file.
         *  @param filename full path name of the file.
         *  @return true if successful, false otherwise.
         *
         *  A skeloton file defines the properties supported by this driver. It is a list of defXXX elements enclosed by @<INDIDriver>@
         *  and @</INDIDriver>@ opening and closing tags. After the properties are created, they can be rerieved, manipulated, and defined
         *  to other clients.
         *  @see An example skeleton file can be found under examples/tutorial_four_sk.xml
         */
        bool buildSkeleton(const char *filename);

        /** @brief Build a property given the supplied XML element (defXXX)
         *  @param root XML element to parse and build.
         *  @param errmsg buffer to store error message in parsing fails.
         *  @param isDynamic set to true if property is loaded from an XML file.
         *  @return 0 if parsing is successful, -1 otherwise and errmsg is set
         */
        int buildProp(const HYDROGEN::LilXmlElement &root, char *errmsg, bool isDynamic = false);

        /** @brief handle SetXXX commands from client */
        int setValue(const HYDROGEN::LilXmlElement &root, char *errmsg);

        /** @brief Return full path to shared INDI file (typically installed to /usr/share/indi on Linux)
         *  @example On Linux, calling getSharedFilePath("drivers.xml") will return /usr/share/indi/drivers.xml
         */
        static std::string getSharedFilePath(std::string fileName);

    public:
        HYDROGEN_DEPRECATED("Do not use BaseDevice as pointer.")
        operator BaseDevice *();

        HYDROGEN_DEPRECATED("Do not use BaseDevice as pointer.")
        BaseDevice *operator->();

        HYDROGEN_DEPRECATED("Use comparison to true.")
        bool operator!=(std::nullptr_t) const
        {
            return isValid();
        }

        HYDROGEN_DEPRECATED("Use comparison to false.")
        bool operator==(std::nullptr_t) const
        {
            return !isValid();
        }

        operator bool() const
        {
            return isValid();
        }
        operator bool()
        {
            return isValid();
        }

    protected:
        friend class AbstractBaseClientPrivate;
        std::shared_ptr<BaseDevicePrivate> d_ptr;
        BaseDevice(BaseDevicePrivate &dd);
        BaseDevice(const std::shared_ptr<BaseDevicePrivate> &dd);
    };

}
