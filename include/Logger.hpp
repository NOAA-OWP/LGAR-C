#ifndef LOGGER_HPP
#define LOGGER_HPP

#include <stdlib.h>
#include <string>
#include <iostream>
#include <memory>
#include <sstream>
#include <fstream>
#include <ctime>
#include <stdarg.h> // for variable args: va_list

#define LOG (Logger::GetInstance())->Log

enum class LogLevel {
	NONE = 0,
	DEBUG = 1,
	INFO = 2,
	ERROR = 3,
	WARN = 4,
	FATAL = 5,
};

#undef NGEN

enum class LoggingModule {
	NGEN = 0,
	NOAHOWP, 
	SNOW17, 
	UEB, 
	CFE, 
	SACSMA, 
	LASAM, 
	SMP, 
	SFT, 
	TROUTE, 
	SCHISM, 
	SFINCS, 
	GC2D, 
	TOPOFLOW,
	MODULE_COUNT
};

/**
* Logger Class Used to Output Details of Current Application Flow
*/
class Logger {
  public:
	static std::shared_ptr<Logger> GetInstance();
	void SetLogPreferences(LogLevel level);
	void Log(std::string message, LogLevel messageLevel);
	LogLevel GetLogLevel(const std::string& logLevel);
	std::string createTimestamp();
	static void setup_logger(void);
	std::string getLogFilePath();
	static void debug_log(const char* message, ...);

  private:
	LogLevel logLevel;
	std::fstream logFile;
	static std::shared_ptr<Logger> loggerInstance;
	std::string logFilePath;

};



#endif
