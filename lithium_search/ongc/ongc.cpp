/*
 * ongc.cpp
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

#include "ongc.hpp"
#include "ongc_utils.hpp"

#include <cmath>
#include <regex>

#include <sqlite3.h>

#if ENABLE_FASTHASH
#include "emhash/hash_table8.hpp"
#else
#include <unordered_map>
#endif
#include "atom/type/json.hpp"
#include "atom/log/loguru.hpp"

using json = nlohmann::json;

#define DBPATH "ognc.db"

template <int i>
std::vector<std::string> get_identifiers_helper(const std::tuple<std::string, std::vector<std::string>, std::vector<std::string>, std::vector<std::string>, std::vector<std::string>> &identifiers)
{
    return std::get<i>(identifiers);
}

template <>
std::vector<std::string> get_identifiers_helper<0>(const std::tuple<std::string, std::vector<std::string>, std::vector<std::string>, std::vector<std::string>, std::vector<std::string>> &identifiers)
{
    return {std::get<0>(identifiers)};
}

OpenNGC::OpenNGC(const std::string &path, const std::string &version, const std::string &author, const std::string &description)
    : Plugin(path, version, author, description)
{
    loadDatabase({{"path", ""}, {"name", "ongc.db"}});
}

OpenNGC::~OpenNGC()
{
}

void OpenNGC::loadDatabase(const json &params)
{
    if (!is_database_loaded)
    {
        try
        {
            db = std::make_unique<SqliteDB>(params["path"].empty() ? DBPATH : params["path"].get<std::string>().c_str());
            is_database_loaded = true;
        }
        catch (const std::exception &e)
        {
            LOG_F(ERROR, "Failed to load database: {}", e.what());
        }
    }
    else
    {
        LOG_F(WARNING, "Database is already loaded");
    }
}

void OpenNGC::unloadDatabase(const json &params)
{
    if (!is_database_loaded)
    {
        LOG_F(WARNING, "Database is not loaded");
    }
    else
    {
        db.reset();
        is_database_loaded = false;
    }
}