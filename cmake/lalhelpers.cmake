


include(GNUInstallDirs)

function(_process_header_files component hdrs_in hdrs_out)
    set(_tmp "")
    foreach(_file IN LISTS hdrs_in)
        set(_extended_path "${LAL_INC_DIR}/${component}/${_file}")
        if (NOT EXISTS "${_extended_path}")
            message(SEND_ERROR "File does not exist:\n${_extended_path}")
        else()
            list(APPEND _tmp "${_extended_path}")
        endif()
    endforeach()
    set(${hdrs_out} "${_tmp}" PARENT_SCOPE)
endfunction()


function(lal_component)

    cmake_parse_arguments(LAL_LIB
            ""                          # options
            "NAME;TYPE"                 # one value args
            "SOURCES;HEADERS;PRIVATE_DEPS;PUBLIC_DEPS;DEFN"        # multivalue args
            ${ARGN}
            )
    message(STATUS "Adding component lal::${LAL_LIB_NAME}")

    set(_lib_name "libalgebra_lite_${LAL_LIB_NAME}")

    _process_header_files("${LAL_LIB_NAME}" "${LAL_LIB_HEADERS}" _hdrs)


    if (LAL_LIB_TYPE STREQUAL "STATIC" OR LAL_LIB_TYPE STREQUAL "SHARED")
        set(_srcs "")
        foreach(_file IN LISTS LAL_LIB_SOURCES)
            if (IS_ABSOLUTE "${_file}")
                list(APPEND "${_file}")
            else()
                list(APPEND "${CMAKE_CURRENT_LIST_DIR}/${_file}")
            endif()
        endforeach()


        add_library(${_lib_name} ${LAL_LIB_TYPE})
        add_library(lal::${LAL_LIB_NAME} ALIAS ${_lib_name})

        target_sources(${_lib_name} PRIVATE "${_srcs}" "${_hdrs}")

        set_target_properties(${_lib_name} PROPERTIES
                POSITION_INDEPENDENT_CODE ON
                LINKER_LANGUAGE CXX
                )

        target_link_libraries(${_lib_name}
                PUBLIC "${LAL_LIB_PUBLIC_DEPS}"
                PRIVATE "${LAL_LIB_PRIVATE_DEPS}")

        target_include_directories(${_lib_name} PRIVATE
                "${LAL_INC_DIR}/${LAL_LIB_NAME}")

        if (LAL_LiB_DEFN)
            target_compile_definitions(${_lib_name}
                    PRIVATE "${LAL_LIB_DEFN}")
        endif()


        install(TARGETS ${_lib_name}
                RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
                LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
                ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        )
    else()

        add_library(${_lib_name} INTERFACE)
        add_library(lal::${LAL_LIB_NAME} ALIAS ${_lib_name})

        target_sources(${_lib_name} INTERFACE "${_hdrs}")

        target_link_libraries(${_lib_name} INTERFACE
                "${LAL_LIB_PUBLIC_DEPS}"
                )
        target_compile_definitions(${_lib_name} INTERFACE
                "${LAL_LIB_DEFN}")

    endif()


    install(FILES ${_hdrs} DESTINATION
            ${CMAKE_INSTALL_INCLUDEDIR}/libalgebra_lite/${LAL_LIB_NAME})

endfunction()
