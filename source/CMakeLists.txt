#export CXX=g++-8
cmake_minimum_required(VERSION 3.16)
project(Jde.Mat)
set( CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/../../Framework/cmake )
include(SetUpConfigurations)
#Custom######################################################
#message( "REPO_DIR=$ENV{REPO_DIR}" )
include_directories( "$ENV{REPO_DIR}/matio/src" )
include_directories( "$ENV{REPO_DIR}/eigen" )
#############################################################
add_library( Jde.Mat SHARED MatLabFile.cpp MatLabVariable.cpp )
#############################################################
target_precompile_headers(Jde.Mat
  PRIVATE
	<boost/crc.hpp>
	<boost/system/error_code.hpp>
	<spdlog/spdlog.h>
	<spdlog/sinks/basic_file_sink.h>
	<spdlog/fmt/ostr.h>
	"../../Framework/source/log/Logging.h"
	"../../Framework/source/StringUtilities.h"
	"../../Framework/source/io/File.h"
)
