/*******************************************************************************
  Copyright(c) 2016 Jasem Mutlaq. All rights reserved.
  Copyright(c) 2022 Pawel Soja <kernel32.pl@gmail.com>

 HYDROGEN Qt Client

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

#include "abstractbaseclient_p.h"
#include "hydrogenlilxml.h"

#include <QTcpSocket>

namespace HYDROGEN
{

class BaseClientQtPrivate : public AbstractBaseClientPrivate
{
    public:
        BaseClientQtPrivate(BaseClientQt *parent);
        ~BaseClientQtPrivate() = default;

    public:
        ssize_t sendData(const void *data, size_t size) override;

    public:
        void listenHYDROGEN();

    public:
        QTcpSocket clientSocket;
        LilXmlParser xmlParser;

};
}
