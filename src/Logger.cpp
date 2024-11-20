#include "Logger.hpp"
#include <cstring>
#include <string>
#include <chrono>
#include <thread>
#include <cassert>

std::shared_ptr<Logger> Logger::loggerInstance;

using namespace std;

std::string module_name[static_cast<int>(LoggingModule::MODULE_COUNT)] 
{
	"NGEN    ",
	"NOAHOWP ", 
	"SNOW17  ", 
	"UEB     ", 
	"CFE     ", 
	"SACSMA  ", 
	"LASAM   ", 
	"SMP     ", 
	"SFT     ", 
	"TROUTE  ", 
	"SCHISM  ", 
	"SFINCS  ", 
	"GC2D    ", 
	"TOPOFLOW",
};

/**
* Configure Logger Preferences
* @param logFile
* @param level: LogLevel::ERROR by Default
* @param output: LogOutput::CONSOLE by Default
* @return void
*/
void Logger::SetLogPreferences(LogLevel level = LogLevel::ERROR) {
	std::stringstream ss("");

	// set the logging level
 	logLevel = level;
    
	// get the log file path
	ss << getenv("NGEN_LOG_FILE_PATH");
	logFilePath = ss.str();
	ss.str("");

	logFile.open(logFilePath, std::ios::app);
	if (!logFile.good()) {
		std::cerr << "Warning: Can't Open shared Log File referenced from NGEN_LOG_FILE_PATH env. variable" << std::endl;
    	// create a local log file for LASAM module instead
		std::string fwd_slash = "/";
    	std::string logFileName = "lasam_log.txt";
    	std::string logFileDir = "./run-logs/ngen_" + Logger::createTimestamp() + fwd_slash;
   		int status;
		std::string mkdir_cmd = "mkdir -p " + logFileDir;
		const char *cstr = mkdir_cmd.c_str();
   		status = system(cstr);
   		if (status == -1)
   		   std::cerr << "Error(" << (errno) << ") creating log file directory: " << logFileDir << std::endl;
   		else
   		   std::cout << "Log directory: " << logFileDir <<std::endl;
		// create a local log file for LASAM module
		logFilePath = logFileDir+logFileName;
		logFile.open(logFilePath, ios::out | ios::app);
		if (!logFile.good()) {
			std::cerr << "Can't Open local directory Log File:" << logFilePath <<std::endl;			
		}
		else {
			std::cout << "Logging instead into: " << logFilePath << std::endl;
		}
			
	}
	else {
		std::cout << "Log File Path:" << logFilePath << std::endl;
	}

}

/**
* Get Single Logger Instance or Create new Object if Not Created
* @return std::shared_ptr<Logger>
*/
std::shared_ptr<Logger> Logger::GetInstance() {
	if (loggerInstance == nullptr) {
		loggerInstance = std::shared_ptr<Logger>(new Logger());
	}

	return loggerInstance;
}

/**
* Log given message with defined parameters and generate message to pass on Console or File
* @param message: Log Message
* @param messageLevel: Log Level, LogLevel::DEBUG by default
*/
void Logger::Log(std::string message, LogLevel messageLevel = LogLevel::DEBUG) {
	LoggingModule module=LoggingModule::LASAM;

	// don't log if messageLevel < logLevel 
	if (messageLevel >= logLevel) {
		std::string logType;
		//Set Log Level Name
		switch (messageLevel) {
		case LogLevel::FATAL:
			logType = "FATAL ";
			break;
		case LogLevel::DEBUG:
			logType = "DEBUG ";
			break;
		case LogLevel::INFO:
			logType = "INFO  ";
			break;
		case LogLevel::WARN:
			logType = "WARN  ";
			break;
		case LogLevel::ERROR:
			logType = "ERROR ";
			break;
		default:
			logType = "NONE  ";
			break;
		}

		std::string final_message;
		std::string mod_name;
		mod_name = module_name[static_cast<int>(module)];
		std::string separator = " ";
		// log the message while handling multiline cases
		final_message = createTimestamp() + separator + mod_name + separator + logType + message;
		if (!logFile.bad()) {
			logFile << final_message;
			std::cout << final_message;
			logFile.flush();
		}

	}
}

/**
* Convert String Representation of Log Level to LogLevel Type
* @param logLevel : String log level
* @return LogLevel
*/
LogLevel Logger::GetLogLevel(const std::string& logLevel) {
	if (logLevel == "DEBUG") {
		return LogLevel::DEBUG;
	}
	else if (logLevel == "INFO") {
		return LogLevel::INFO;
	}
	else if (logLevel == "WARN") {
		return LogLevel::ERROR;
	}
	else if (logLevel == "ERROR") {
		return LogLevel::ERROR;
	}
	else if (logLevel == "FATAL") {
		return LogLevel::ERROR;
	}

	return LogLevel::NONE;
}

using std::chrono::system_clock;

std::string Logger::createTimestamp() {
    std::chrono::_V2::system_clock::time_point currentTime = std::chrono::system_clock::now();
    char buffer1[120];
    char buffer2[120];
    std::stringstream ss;

    long transformed = currentTime.time_since_epoch().count() / 1000000;
    
    long millis = transformed % 1000;
    
    std::time_t tt;
    tt = system_clock::to_time_t ( currentTime );
    tm *timeinfo = gmtime (&tt);
    strftime (buffer1,100,"%FT%H:%M:%S",timeinfo);
    sprintf(buffer2, ":%03d", (int)millis);
	ss << buffer1 << buffer2;
    
    return ss.str();
}

void Logger::setup_logger(void) {
    // One time log preferences
    (Logger::GetInstance())->SetLogPreferences(LogLevel::INFO);
}

std::string Logger::getLogFilePath() {
	return logFilePath;
}