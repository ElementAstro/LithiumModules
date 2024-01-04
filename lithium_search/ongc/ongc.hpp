/*
 * ongc.hpp
 *
 * Copyright (C) 2023-2024 Max Qian <lightapt.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

/*************************************************

Date: 2023-7-13

Description: C++ Version of PyOngc

**************************************************/

#pragma once

#include "dso.hpp"
#include "core/plugin/plugin.hpp"
#include "atom/search/sqlite.hpp"

class OpenNGC : public Plugin
{
public:
    OpenNGC(const std::string &path, const std::string &version, const std::string &author, const std::string &description);
    ~OpenNGC();

    void loadDatabase(const json &params);
    void unloadDatabase(const json &params);

private:

    void recognizeName(const std::string &name, std::string &catalog, std::string &objectname);

    std::unique_ptr<SqliteDB> db;

    bool is_database_loaded;
};