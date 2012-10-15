## generate libdivsufsort.pc ##
set(prefix "${CMAKE_INSTALL_PREFIX}")
set(W64BIT "")
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/libdivsufsort.pc.cmake" "${CMAKE_CURRENT_BINARY_DIR}/libdivsufsort.pc" @ONLY)
if(BUILD_DIVSUFSORT64)
  set(W64BIT "64")
  configure_file("${CMAKE_CURRENT_SOURCE_DIR}/libdivsufsort.pc.cmake" "${CMAKE_CURRENT_BINARY_DIR}/libdivsufsort64.pc" @ONLY)
endif(BUILD_DIVSUFSORT64)
