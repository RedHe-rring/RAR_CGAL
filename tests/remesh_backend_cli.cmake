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

# Diagnostics use the field directory, including an explicit prefix.
set(prefix "${TEST_DIR}/shared-field/field")
set(diagnostic_output "${TEST_DIR}/diagnostic-result.off")
foreach(level stage operation)
    execute_process(
        COMMAND "${EXE}" "${input}" "${diagnostic_output}"
            --field rar --remesh-backend custom --remesh-diagnostics --diagnostic-level "${level}"
            --epsilon 0.04 --min-edge 0.15 --max-edge 0.7 --export-field "${prefix}"
        RESULT_VARIABLE result OUTPUT_VARIABLE stdout ERROR_VARIABLE stderr TIMEOUT 60)
    if(NOT "${result}" STREQUAL "0")
        message(FATAL_ERROR "Diagnostics failed: ${stdout}\n${stderr}")
    endif()
    string(REGEX MATCH "Remesh diagnostics: ([^\r\n]+)" match "${stdout}")
    set(run "${CMAKE_MATCH_1}")
    file(TO_CMAKE_PATH "${run}" run)
    string(FIND "${run}" "${TEST_DIR}/shared-field/remesh_diagnostics/run_" prefix_position)
    if(NOT prefix_position EQUAL 0 OR NOT EXISTS "${prefix}.csv" OR NOT EXISTS "${prefix}.ply")
        message(FATAL_ERROR "Diagnostics do not share the field directory: ${stdout}")
    endif()
    file(SHA256 "${diagnostic_output}" diagnostic_hash)
    file(SHA256 "${TEST_DIR}/rar-custom.off" reference_hash)
    if(NOT diagnostic_hash STREQUAL reference_hash)
        message(FATAL_ERROR "Enabling diagnostics changed the final mesh")
    endif()
    file(READ "${run}/manifest.json" manifest)
    string(JSON actual_sha1 GET "${manifest}" input_sha1)
    file(SHA1 "${input}" expected_sha1)
    if(NOT actual_sha1 STREQUAL expected_sha1)
        message(FATAL_ERROR "Incorrect diagnostic input checksum: ${actual_sha1}")
    endif()
    file(READ "${run}/status.json" status)
    string(JSON state GET "${status}" state)
    if(NOT state STREQUAL "complete" OR NOT EXISTS "${run}/iter_005/result.off")
        message(FATAL_ERROR "Incomplete diagnostic run")
    endif()
endforeach()

execute_process(COMMAND "${EXE}" "${input}" --remesh-diagnostics
    RESULT_VARIABLE result ERROR_VARIABLE stderr)
if("${result}" STREQUAL "0" OR NOT stderr MATCHES "require --remesh-backend custom")
    message(FATAL_ERROR "CGAL backend incorrectly accepted diagnostics: ${stderr}")
endif()
execute_process(COMMAND "${EXE}" "${input}" --diagnostic-level invalid
    RESULT_VARIABLE result ERROR_VARIABLE stderr)
if("${result}" STREQUAL "0" OR NOT stderr MATCHES "Unknown diagnostic level")
    message(FATAL_ERROR "Invalid diagnostic level was not rejected: ${stderr}")
endif()
execute_process(COMMAND "${EXE}" "${input}" --diagnostic-level
    RESULT_VARIABLE result ERROR_VARIABLE stderr)
if("${result}" STREQUAL "0" OR NOT stderr MATCHES "Missing value")
    message(FATAL_ERROR "Missing diagnostic level was not rejected: ${stderr}")
endif()

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
