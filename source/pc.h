#pragma once

#pragma warning(push)
#pragma warning( disable : 4245) 
#pragma warning( disable : 4701)
#include <boost/crc.hpp>
#pragma warning(pop)
#include <boost/system/error_code.hpp>

#ifndef __INTELLISENSE__
	#include <spdlog/spdlog.h>
	#include <spdlog/sinks/basic_file_sink.h>
	#include <spdlog/fmt/ostr.h>
#endif
#include "../../framework/source/log/Logging.h"
#include "../../framework/source/StringUtilities.h"
#include "../../framework/source/io/File.h"

