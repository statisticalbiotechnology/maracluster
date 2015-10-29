

###############################################################################
# FIND EXTERNAL LIBRARIES REQUIRED BY MERGESPEC
###############################################################################

# INCLUDE HEADER FILES IN SRC DIRECTORY
include_directories(${CMAKE_CURRENT_BINARY_DIR})

MESSAGE( STATUS "Checking packages:" )

#set(Boost_INCLUDE_DIRS ${BOOST_INCLUDEDIR})
#set(Boost_LIBRARIES boost_filesystem-mt-s boost_iostreams-mt-s boost_regex-mt-s boost_date_time-mt-s boost_thread-mt-s boost_serialization-mt-s boost_system-mt-s)
find_package(Boost ${BOOST_MIN_VERSION} COMPONENTS filesystem iostreams regex thread serialization system chrono REQUIRED)
if(Boost_FOUND)
  message(STATUS "Package Boost found at: ${Boost_INCLUDE_DIRS}")
else(Boost_FOUND)
  message(FATAL_ERROR "Boost has not been found")
endif(Boost_FOUND)
include_directories(${Boost_INCLUDE_DIRS})

find_package(ZLIB REQUIRED)
if(ZLIB_FOUND)
  message(STATUS "Zlib found")
else(ZLIB_FOUND)
  message(FATAL_ERROR "The package Zlib has not been found")
endif(ZLIB_FOUND)
include_directories(${ZLIB_INCLUDE_DIR})

find_package(ProteoWizard REQUIRED)
if(PWIZ_FOUND)
  message(STATUS "ProteoWizard found at: ${PWIZ_INCLUDE_DIR}")
else(PWIZ_FOUND)
  message(FATAL_ERROR "The package ProteoWizard has not been found")
endif(PWIZ_FOUND)
include_directories(${PWIZ_INCLUDE_DIR})

find_package(Threads REQUIRED)

#############################################################################
# COMPILE MARACLUSTER
#############################################################################

if(FINGERPRINT_FILTER)
  set(BFM_SRC "BinaryFingerprintMethods.cpp")
else(FINGERPRINT_FILTER)
  set(BFM_SRC "")
endif(FINGERPRINT_FILTER)

add_library(maraclusterlibrary STATIC SparseClustering.cpp MatrixLoader.cpp PvalueCalculator.cpp PvalueFilterAndSort.cpp PeakDistribution.cpp PercolatorInterface.cpp BinSpectra.cpp BinAndRank.cpp InterpolationMerge.cpp RankMerge.cpp ClusterMerge.cpp PeakCounts.cpp ScanMergeInfo.cpp SpectrumFileList.cpp SpectrumHandler.cpp MSFileHandler.cpp MSFileExtractor.cpp MSFileMerger.cpp MZIntensityPair.cpp MSClusterMerge.cpp Option.cpp MyException.cpp  ScanId.cpp PvalueTriplet.cpp ${BFM_SRC})

add_library(batchlibrary STATIC BatchGlobals.cpp BatchPvalues.cpp BatchPvalueVectors.cpp BatchSpectra.cpp BatchSpectrumClusters.cpp BatchSpectrumFiles.cpp)

#add_executable(extractspec extractSpectra.cpp)
#add_executable(msgffixmzml msgfFixMzML.cpp)

add_executable(maracluster maracluster.cpp)

if(APPLE)
  set_property(TARGET maracluster PROPERTY LINK_SEARCH_START_STATIC FALSE)
  set_property(TARGET maracluster PROPERTY LINK_SEARCH_END_STATIC TRUE)
#  set_property(TARGET mergespec PROPERTY MACOSX_BUNDLE ON)
elseif(MINGW)
#  set_target_properties(maracluster PROPERTIES LINK_FLAGS "-Wl,-Bstatic -static-libgcc -lcurl -liconv")
  set_target_properties(maracluster PROPERTIES LINK_FLAGS "-static-libgcc")
endif()

if(PROFILING AND UNIX)
  set_target_properties(maracluster PROPERTIES COMPILE_FLAGS -pg LINK_FLAGS -pg)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wno-unused-result -Wno-overflow -Wall -g -pg -fprofile-arcs -ftest-coverage")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -Wno-unused-result -Wno-overflow -Wall -g -pg -fprofile-arcs -ftest-coverage")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -Wno-unused-result -Wno-overflow -Wall -g -pg -fprofile-arcs -ftest-coverage")
  set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O3 -Wno-unused-result -Wno-overflow -Wall -g -pg -fprofile-arcs -ftest-coverage")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -g -DNDEBUG -Wno-unused-result -Wno-overflow -Wall -g -pg -fprofile-arcs -ftest-coverage")
  set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O0 -g -DNDEBUG -Wno-unused-result -Wno-overflow -Wall -g -pg -fprofile-arcs -ftest-coverage")
endif() 

set(COMMON_LIBRARIES ${PWIZ_LIBRARIES} ${Boost_LIBRARIES} ${ZLIB_LIBRARY} ${CMAKE_THREAD_LIBS_INIT})

# SET LIBRARIES FOR MARACLUSTER
#target_link_libraries(extractspec maraclusterlibrary ${COMMON_LIBRARIES})
#target_link_libraries(msgffixmzml maraclusterlibrary ${COMMON_LIBRARIES})

target_link_libraries(maracluster batchlibrary maraclusterlibrary ${COMMON_LIBRARIES})

# INSTALL MARACLUSTER
#install(TARGETS extractspec msgffixmzml maracluster EXPORT MARACLUSTER DESTINATION bin) # Important to use relative path here (used by CPack)!
install(TARGETS maracluster EXPORT MARACLUSTER DESTINATION bin) # Important to use relative path here (used by CPack)!

