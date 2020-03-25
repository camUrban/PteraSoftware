#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "osqp::osqpstatic" for configuration "Release"
set_property(TARGET osqp::osqpstatic APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(osqp::osqpstatic PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "C"
  IMPORTED_LOCATION_RELEASE "/home/travis/build/casadi/binaries/casadi/python_install/casadi/libosqp.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS osqp::osqpstatic )
list(APPEND _IMPORT_CHECK_FILES_FOR_osqp::osqpstatic "/home/travis/build/casadi/binaries/casadi/python_install/casadi/libosqp.a" )

# Import target "osqp::osqp" for configuration "Release"
set_property(TARGET osqp::osqp APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(osqp::osqp PROPERTIES
  IMPORTED_IMPLIB_RELEASE "/home/travis/build/casadi/binaries/casadi/python_install/casadi/libosqp.dll.a"
  IMPORTED_LOCATION_RELEASE "/home/travis/build/casadi/binaries/casadi/python_install/casadi/libosqp.dll"
  )

list(APPEND _IMPORT_CHECK_TARGETS osqp::osqp )
list(APPEND _IMPORT_CHECK_FILES_FOR_osqp::osqp "/home/travis/build/casadi/binaries/casadi/python_install/casadi/libosqp.dll.a" "/home/travis/build/casadi/binaries/casadi/python_install/casadi/libosqp.dll" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
