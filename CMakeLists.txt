CMAKE_MINIMUM_REQUIRED(VERSION 3.12)

# Name of the project
PROJECT(2PSP)

set(SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR}/src)

file(GLOB_RECURSE HEADERS ${SOURCE_DIR}/*.h)
file(GLOB_RECURSE SOURCES ${SOURCE_DIR}/*.cpp)

find_package(coal REQUIRED)
# Set the executable.
ADD_EXECUTABLE(${CMAKE_PROJECT_NAME} ${SOURCES} ${HEADERS})

# Get the EIGEN environment variable. Since EIGEN is a header-only library, we
# just need to add it to the include directory.
SET(EIGEN3_INCLUDE_DIR "$ENV{EIGEN3_INCLUDE_DIR}")
IF(NOT EIGEN3_INCLUDE_DIR)
	# The environment variable was not set
	SET(ERR_MSG "Please point the environment variable EIGEN3_INCLUDE_DIR to the root directory of your EIGEN installation.")
	IF(WIN32)
		# On Windows, try the default location
		MESSAGE(STATUS "Looking for EIGEN in ${DEF_DIR_EIGEN}")
		IF(IS_DIRECTORY ${DEF_DIR_EIGEN})
			MESSAGE(STATUS "Found!")
			SET(EIGEN3_INCLUDE_DIR ${DEF_DIR_EIGEN})
		ELSE()
			MESSAGE(FATAL_ERROR ${ERR_MSG})
		ENDIF()
	ELSE()
		MESSAGE(FATAL_ERROR ${ERR_MSG})
	ENDIF()
ENDIF()
INCLUDE_DIRECTORIES(${EIGEN3_INCLUDE_DIR})

# Use c++17
SET_TARGET_PROPERTIES(${CMAKE_PROJECT_NAME} PROPERTIES CXX_STANDARD 17)

target_include_directories(${CMAKE_PROJECT_NAME} PRIVATE ${SOURCE_DIR}/ ${COAL_INCLUDE_DIRS})
target_link_libraries(${CMAKE_PROJECT_NAME} PRIVATE coal::coal)