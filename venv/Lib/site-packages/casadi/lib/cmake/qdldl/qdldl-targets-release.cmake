#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "qdldl::qdldlstatic" for configuration "Release"
set_property(TARGET qdldl::qdldlstatic APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(qdldl::qdldlstatic PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "C"
  IMPORTED_LOCATION_RELEASE "/home/travis/build/casadi/binaries/casadi/python_install/casadi/libqdldl.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS qdldl::qdldlstatic )
list(APPEND _IMPORT_CHECK_FILES_FOR_qdldl::qdldlstatic "/home/travis/build/casadi/binaries/casadi/python_install/casadi/libqdldl.a" )

# Import target "qdldl::qdldl" for configuration "Release"
set_property(TARGET qdldl::qdldl APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(qdldl::qdldl PROPERTIES
  IMPORTED_IMPLIB_RELEASE "/home/travis/build/casadi/binaries/casadi/python_install/casadi/libqdldl.dll.a"
  IMPORTED_LOCATION_RELEASE "/home/travis/build/casadi/binaries/casadi/python_install/casadi/libqdldl.dll"
  )

list(APPEND _IMPORT_CHECK_TARGETS qdldl::qdldl )
list(APPEND _IMPORT_CHECK_FILES_FOR_qdldl::qdldl "/home/travis/build/casadi/binaries/casadi/python_install/casadi/libqdldl.dll.a" "/home/travis/build/casadi/binaries/casadi/python_install/casadi/libqdldl.dll" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
