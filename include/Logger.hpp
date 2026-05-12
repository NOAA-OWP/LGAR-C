#ifndef LASAM_LOGGER_HPP
#define LASAM_LOGGER_HPP

#include "ewts/module_constants.hpp"
#include "ewts/logger.hpp"
#include "ewts/log_levels.hpp"

#define LOG(...) ::ewts::GetLogger(::ewts::modules::EWTS_ID_LASAM).Log(__VA_ARGS__)
#define GetLogLevel() ::ewts::GetLogger(::ewts::modules::EWTS_ID_LASAM).GetLogLevel()
#define IsLoggingEnabled() ::ewts::GetLogger(::ewts::modules::EWTS_ID_LASAM).IsLoggingEnabled()

using ewts::EwtsInit;
using ewts::LogLevel;

inline constexpr const char* LASAM_MODULE_ID = ewts::modules::EWTS_ID_LASAM;

#endif /* LASAM_LOGGER_HPP */
