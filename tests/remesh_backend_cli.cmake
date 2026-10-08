if(NOT DEFINED EXE OR NOT DEFINED TEST_DIR)
    message(FATAL_ERROR "EXE and TEST_DIR are required")
endif()

file(MAKE_DIRECTORY "${TEST_DIR}")
set(input "${TEST_DIR}/octahedron.off")
file(WRITE "${input}" "OFF\n6 8 0\n1 0 0\n-1 0 0\n0 1 0\n0 -1 0\n0 0 1\n0 0 -1\n3 4 0 2\n3 4 2 1\n3 4 1 3\n3 4 3 0\n3 5 2 0\n3 5 1 2\n3 5 3 1\n3 5 0 3\n")

set(fields rar csf cgal-adaptive cgal-adaptive-radius)
if(WITH_IPOPT)
    list(APPEND fields rar-chen cgal-adaptive-chen cgal-adaptive-radius-chen)
endif()
foreach(field IN LISTS fields)
    foreach(backend cgal custom)
        set(output "${TEST_DIR}/${field}-${backend}.off")
        execute_process(
            COMMAND "${EXE}" "${input}" "${output}"
                --field "${field}" --remesh-backend "${backend}"
                --epsilon 0.04 --min-edge 0.15 --max-edge 0.7
                --no-export-field
            RESULT_VARIABLE result OUTPUT_VARIABLE stdout ERROR_VARIABLE stderr
            TIMEOUT 60)
        if(NOT "${result}" STREQUAL "0" OR NOT stdout MATCHES "remesh_backend=${backend}")
            message(FATAL_ERROR "${field}/${backend} failed: ${result}\n${stdout}\n${stderr}")
        endif()
        file(SHA256 "${output}" hash_${backend})
    endforeach()
    if(NOT hash_cgal STREQUAL hash_custom)
        message(FATAL_ERROR "Output differs between backends for ${field}")
    endif()
    message(STATUS "CLI parity: ${field}")
endforeach()

# Exercise automatic output naming and the default backend.
foreach(backend default custom)
    set(backend_args)
    set(expected cgal)
    if(backend STREQUAL "custom")
        set(backend_args --remesh-backend custom)
        set(expected custom)
    endif()
    execute_process(
        COMMAND "${EXE}" "${input}" --field rar ${backend_args}
            --epsilon 0.04 --min-edge 0.15 --max-edge 0.7 --no-export-field
        RESULT_VARIABLE result OUTPUT_VARIABLE stdout ERROR_VARIABLE stderr
        TIMEOUT 60)
    if(NOT "${result}" STREQUAL "0" OR NOT stdout MATCHES "__backend-${expected}\\.off")
        message(FATAL_ERROR "Automatic naming failed: ${stdout}\n${stderr}")
    endif()
endforeach()

execute_process(COMMAND "${EXE}" "${input}" --remesh-backend invalid
    RESULT_VARIABLE result ERROR_VARIABLE stderr)
if("${result}" STREQUAL "0" OR NOT stderr MATCHES "Unknown remesh backend")
    message(FATAL_ERROR "Unknown backend was not rejected: ${stderr}")
endif()
execute_process(COMMAND "${EXE}" "${input}" --remesh-backend
    RESULT_VARIABLE result ERROR_VARIABLE stderr)
if("${result}" STREQUAL "0" OR NOT stderr MATCHES "Missing value")
    message(FATAL_ERROR "Missing backend was not rejected: ${stderr}")
endif()
