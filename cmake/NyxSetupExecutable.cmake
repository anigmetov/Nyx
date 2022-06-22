#
# Macro to setup a nyx executables
#
macro (nyx_setup_executable _srcs _inputs)
   cmake_parse_arguments( "" ""
      "BASE_NAME;RUNTIME_SUBDIR;EXTRA_DEFINITIONS;MAIN" "PROPERTIES" ${ARGN} )

   if (_BASE_NAME)
      set(_base_name ${_BASE_NAME})
   else ()
      string(REGEX REPLACE ".*Exec/" "" _base_name ${CMAKE_CURRENT_LIST_DIR})
      string(REPLACE "/" "_" _base_name ${_base_name})
   endif ()

   # Prepend "nyx_" to base name
   set(_base_name "nyx_${_base_name}")

   if (_RUNTIME_SUBDIR)
      set(_exe_dir ${CMAKE_CURRENT_BINARY_DIR}/${_RUNTIME_SUBDIR})
   else ()
      set(_exe_dir ${CMAKE_CURRENT_BINARY_DIR})
   endif ()

   set( _exe_name  ${_base_name} )

   add_executable( ${_exe_name} )
   target_sources( ${_exe_name} PRIVATE ${${_srcs}} )
   target_link_libraries(${_exe_name} nyxcore)

   set(EXE_NAME ${_exe_name})

   if (_MAIN)
      target_sources( ${_exe_name} PRIVATE ${_MAIN} )
      set(DEFAULT_MAIN 0)
   else ()
      target_sources( ${_exe_name} PRIVATE ${NYX_DEFAULT_MAIN} )
   endif ()

   set_target_properties(${_exe_name}
      PROPERTIES
      RUNTIME_OUTPUT_DIRECTORY ${_exe_dir} )

   if (_EXTRA_DEFINITIONS)
      target_compile_definitions(${_exe_name} PRIVATE ${_EXTRA_DEFINITIONS})
   endif ()

   # Find out which include directory is needed
   set(_includes ${${_srcs}})
   list(FILTER _includes INCLUDE REGEX "\\.H$")
   foreach(_item IN LISTS _includes)
      get_filename_component( _include_dir ${_item} DIRECTORY)
      if (NOT _include_dir)
         set(_include_dir ${CMAKE_CURRENT_LIST_DIR})
      endif ()
      target_include_directories( ${_exe_name} PRIVATE  ${_include_dir} )
   endforeach()

   #target_include_directories( ${_exe_name} PRIVATE  ${LOW_FIVE_DIR}/include ${LOW_FIVE_DIR}/ext/fmt/include ${LOW_FIVE_DIR}/ext/diy/include )
   #target_link_libraries( ${_exe_name} ${LOW_FIVE_BUILD_DIR}/src/liblowfive.so
   #                                            ${LOW_FIVE_BUILD_DIR}/src/liblowfive-dist.a)

   if (Nyx_GPU_BACKEND STREQUAL "CUDA")
      setup_target_for_cuda_compilation( ${_exe_name} )
   endif ()

   if (${_inputs})
      file( COPY ${${_inputs}} DESTINATION ${_exe_dir} )
   endif ()

   # Add shared library with same function for LyA only, to test lowfive with henson
   if (${_exe_name} STREQUAL "nyx_LyA")

      set( _lib_name  "henson_${_base_name}" )

      add_library( ${_lib_name} SHARED)
      set_target_properties(${_lib_name} PROPERTIES PREFIX "")
      set_target_properties(${_lib_name} PROPERTIES SUFFIX ".hx")

      target_sources( ${_lib_name} PRIVATE ${${_srcs}} )
      target_sources( ${_lib_name} PRIVATE ${NYX_DEFAULT_MAIN} )

      target_link_libraries(${_lib_name} ${HENSON_BUILD_DIR}/libhenson-pmpi.so)
      target_link_libraries(${_lib_name} ${HENSON_BUILD_DIR}/libhenson.a)
      target_link_libraries(${_lib_name} nyxcore)

      set_target_properties(${_lib_name} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${_exe_dir} )

      if (_EXTRA_DEFINITIONS)
         target_compile_definitions(${_lib_name} PRIVATE ${_EXTRA_DEFINITIONS})
      endif ()

      # Find out which include directory is needed
      set(_includes ${${_srcs}})
      list(FILTER _includes INCLUDE REGEX "\\.H$")
      foreach(_item IN LISTS _includes)
         get_filename_component( _include_dir ${_item} DIRECTORY)
         if (NOT _include_dir)
            set(_include_dir ${CMAKE_CURRENT_LIST_DIR})
         endif ()
         target_include_directories( ${_lib_name} PRIVATE  ${_include_dir} )
      endforeach()

      # Henson-specific
      set(linker_flags "-Wl,--export-dynamic")
      set(linker_flags "${linker_flags} -Wl,-u,henson_set_contexts,-u,henson_set_namemap")
      set_target_properties(${_lib_name} PROPERTIES LINK_FLAGS ${linker_flags})
   endif()
endmacro ()
