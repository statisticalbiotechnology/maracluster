# - Try to find Cassandra cpp-driver package
# This module defines
# CASSANDRA_INCLUDE_DIR, where to find cassandra.h, etc.
# CASSANDRA_LIBRARY, the library to link against to use cassandra.
# CASSANDRA_FOUND, If false, don't try to use cassandra.

find_library(CASSANDRA_LIBRARY NAMES cassandra)
find_library(LIBUV_LIBRARY
  NAMES uv libuv)

find_path(CASSANDRA_INCLUDE_DIR NAMES cassandra.h )
find_path(LIBUV_INCLUDE_DIR
  NAMES uv.h)
  
set(CASSANDRA_LIBRARIES ${CASSANDRA_LIBRARY} ${LIBUV_LIBRARY} )
set(CASSANDRA_INCLUDE_DIRS ${CASSANDRA_INCLUDE_DIR} ${LIBUV_INCLUDE_DIR} )
# if the include and the library are found then we have it
if (CASSANDRA_INCLUDE_DIRS AND CASSANDRA_LIBRARIES)
  SET( CASSANDRA_FOUND 1 )
endif()

MARK_AS_ADVANCED( 
  CASSANDRA_INCLUDE_DIRS
  CASSANDRA_LIBRARIES
)
