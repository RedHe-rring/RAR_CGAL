file(MAKE_DIRECTORY "${TEST_DIR}")
set(input "${TEST_DIR}/octahedron.off")
file(WRITE "${input}" "OFF\n6 8 0\n1 0 0\n-1 0 0\n0 1 0\n0 -1 0\n0 0 1\n0 0 -1\n3 4 0 2\n3 4 2 1\n3 4 1 3\n3 4 3 0\n3 5 2 0\n3 5 1 2\n3 5 3 1\n3 5 0 3\n")
set(args --field cgal-adaptive-sdf-thickness --epsilon 0.04 --min-edge 0.15 --max-edge 0.7 --iterations 1)
foreach(backend cgal custom)
    set(prefix "${TEST_DIR}/${backend}/field")
    execute_process(COMMAND "${EXE}" "${input}" "${TEST_DIR}/${backend}.off" ${args}
        --remesh-backend "${backend}" --export-field "${prefix}"
        RESULT_VARIABLE result OUTPUT_VARIABLE stdout ERROR_VARIABLE stderr TIMEOUT 60)
    if(NOT "${result}" STREQUAL "0" OR NOT stdout MATCHES "Chen disabled")
        message(FATAL_ERROR "SDF ${backend} failed: ${stdout}\n${stderr}")
    endif()
    file(SHA256 "${TEST_DIR}/${backend}.off" hash_${backend})
    file(READ "${prefix}.csv" csv)
    if(NOT csv MATCHES "thickness,vertex_area,epsilon_local" OR NOT csv MATCHES "min_edge_clipped")
        message(FATAL_ERROR "Missing initial SDF diagnostics")
    endif()
    file(STRINGS "${prefix}.csv" rows)
    list(LENGTH rows count)
    if(NOT count EQUAL 7)
        message(FATAL_ERROR "Field is not exported on the 6 original vertices")
    endif()
    file(READ "${prefix}_sdf.json" metadata)
    string(JSON postprocess GET "${metadata}" sdf_postprocess)
    string(JSON chen GET "${metadata}" chen_applied)
    string(JSON self_intersections GET "${metadata}" input_self_intersections)
    string(JSON thickness_ratio GET "${metadata}" thickness_epsilon_ratio)
    if(postprocess OR chen OR self_intersections OR NOT thickness_ratio EQUAL 0.25 OR NOT EXISTS "${prefix}_sdf_faces.csv")
        message(FATAL_ERROR "Incorrect SDF metadata")
    endif()
endforeach()
if(NOT hash_cgal STREQUAL hash_custom)
    message(FATAL_ERROR "SDF backend parity failed")
endif()

# Chen mode preserves the SDF initial target and exports the corrected target.
set(chen_args --field cgal-adaptive-sdf-thickness-chen --epsilon 0.04
    --min-edge 0.15 --max-edge 0.7 --sdf-rays 31 --sdf-cone-angle 90
    --beta 1.05 --iterations 1)
if(WITH_IPOPT)
    foreach(backend cgal custom)
        set(prefix "${TEST_DIR}/chen-${backend}/field")
        execute_process(COMMAND "${EXE}" "${input}" "${TEST_DIR}/chen-${backend}.off"
            ${chen_args} --remesh-backend "${backend}" --export-field "${prefix}"
            RESULT_VARIABLE result OUTPUT_VARIABLE stdout ERROR_VARIABLE stderr TIMEOUT 60)
        if(NOT "${result}" STREQUAL "0" OR NOT stdout MATCHES "Chen correction:")
            message(FATAL_ERROR "SDF+Chen ${backend} failed: ${stdout}\n${stderr}")
        endif()
        file(SHA256 "${TEST_DIR}/chen-${backend}.off" chen_hash_${backend})
        file(STRINGS "${prefix}.csv" chen_rows)
        list(GET chen_rows 0 header)
        if(NOT header MATCHES "curvature,raw_target_length,target_length,thickness")
            message(FATAL_ERROR "SDF+Chen raw/corrected targets missing")
        endif()
        file(READ "${prefix}_sdf.json" chen_metadata)
        string(JSON chen_enabled GET "${chen_metadata}" chen_applied)
        string(JSON chen_beta GET "${chen_metadata}" chen_beta)
        string(JSON chen_field GET "${chen_metadata}" field)
        string(JSON chen_rays GET "${chen_metadata}" sdf_rays)
        string(JSON chen_angle GET "${chen_metadata}" sdf_cone_angle_degrees)
        if(NOT chen_enabled OR NOT chen_beta EQUAL 1.05 OR
           NOT chen_field STREQUAL "cgal-adaptive-sdf-thickness-chen" OR
           NOT chen_rays EQUAL 31 OR NOT chen_angle EQUAL 90)
            message(FATAL_ERROR "Incorrect SDF+Chen metadata")
        endif()
        list(REMOVE_AT chen_rows 0)
        foreach(row IN LISTS chen_rows)
            string(REPLACE "," ";" values "${row}")
            list(GET values 5 raw_target)
            list(GET values 6 corrected_target)
            if(corrected_target GREATER raw_target OR corrected_target LESS 0.15)
                message(FATAL_ERROR "SDF+Chen target outside [min_edge, raw_target]")
            endif()
        endforeach()
    endforeach()
    if(NOT chen_hash_cgal STREQUAL chen_hash_custom)
        message(FATAL_ERROR "SDF+Chen backend parity failed")
    endif()

    # An asymmetric shape must actually enter the IPOPT correction path.
    set(irregular "${TEST_DIR}/irregular.off")
    file(READ "${input}" irregular_off)
    string(REPLACE "1 0 0" "3 0 0" irregular_off "${irregular_off}")
    string(REPLACE "0 0 1" "0 0 0.2" irregular_off "${irregular_off}")
    file(WRITE "${irregular}" "${irregular_off}")
    set(prefix "${TEST_DIR}/irregular/field")
    execute_process(COMMAND "${EXE}" "${irregular}" "${TEST_DIR}/irregular-result.off"
        --field cgal-adaptive-sdf-thickness-chen --epsilon 0.1
        --min-edge 0.01 --max-edge 1 --beta 1.01 --iterations 0
        --export-field "${prefix}"
        RESULT_VARIABLE result OUTPUT_VARIABLE stdout ERROR_VARIABLE stderr TIMEOUT 60)
    if(NOT "${result}" STREQUAL "0")
        message(FATAL_ERROR "SDF+Chen active correction failed: ${stdout}\n${stderr}")
    endif()
    file(READ "${prefix}_sdf.json" active_metadata)
    string(JSON changed GET "${active_metadata}" chen_changed_vertex_count)
    string(JSON before GET "${active_metadata}" chen_max_gradient_before)
    string(JSON after GET "${active_metadata}" chen_max_gradient_after)
    if(NOT changed GREATER 0 OR NOT before GREATER after OR after GREATER 0.0101)
        message(FATAL_ERROR "SDF+Chen did not enforce the beta=1.01 gradient limit")
    endif()

    execute_process(COMMAND "${EXE}" "${input}" --field cgal-adaptive-sdf-thickness-chen
        --epsilon 0.04 --beta 1.05 --iterations 0 --no-export-field
        RESULT_VARIABLE result OUTPUT_VARIABLE stdout ERROR_VARIABLE stderr TIMEOUT 60)
    if(NOT "${result}" STREQUAL "0" OR
       NOT stdout MATCHES "__sdf-h0-chen__e-0p04__mn-0p001__mx-0p5__r25__a120__i0__q3__p1__f0__beta-1p05__b-cgal")
        message(FATAL_ERROR "SDF+Chen compact automatic naming failed: ${stdout}\n${stderr}")
    endif()
else()
    execute_process(COMMAND "${EXE}" "${input}" ${chen_args}
        RESULT_VARIABLE result OUTPUT_VARIABLE stdout ERROR_VARIABLE stderr TIMEOUT 60)
    if("${result}" STREQUAL "0" OR NOT stderr MATCHES "requires IPOPT")
        message(FATAL_ERROR "SDF+Chen must report missing IPOPT: ${stdout}\n${stderr}")
    endif()
endif()

execute_process(COMMAND "${EXE}" "${input}" ${chen_args} --beta 1
    RESULT_VARIABLE result OUTPUT_VARIABLE stdout ERROR_VARIABLE stderr TIMEOUT 60)
if("${result}" STREQUAL "0" OR NOT stderr MATCHES "beta must be")
    message(FATAL_ERROR "SDF+Chen accepted invalid beta: ${stdout}\n${stderr}")
endif()

execute_process(COMMAND "${EXE}" "${input}" "${TEST_DIR}/observed.off" ${args}
    --remesh-backend custom --remesh-diagnostics --no-export-field
    RESULT_VARIABLE result OUTPUT_VARIABLE stdout ERROR_VARIABLE stderr TIMEOUT 60)
if(NOT "${result}" STREQUAL "0")
    message(FATAL_ERROR "SDF observer failed: ${stdout}\n${stderr}")
endif()
file(SHA256 "${TEST_DIR}/observed.off" hash_observed)
if(NOT hash_observed STREQUAL hash_custom)
    message(FATAL_ERROR "SDF observer changed the mesh")
endif()
string(REGEX MATCH "Remesh diagnostics: ([^\r\n]+)" match "${stdout}")
set(run "${CMAKE_MATCH_1}")
file(READ "${run}/manifest.json" manifest)
string(JSON rays GET "${manifest}" sdf_rays)
string(JSON angle GET "${manifest}" sdf_cone_angle_degrees)
if(NOT rays EQUAL 25 OR NOT angle EQUAL 120)
    message(FATAL_ERROR "SDF parameters missing from remesh manifest")
endif()

# Default automatic epsilon, parameterized output names, and field-only mode.
execute_process(COMMAND "${EXE}" "${input}" --field cgal-adaptive-sdf-thickness
    --sdf-rays 31 --sdf-cone-angle 90 --iterations 1 --no-export-field
    --remesh-backend custom --remesh-diagnostics
    RESULT_VARIABLE result OUTPUT_VARIABLE stdout ERROR_VARIABLE stderr TIMEOUT 60)
if(NOT "${result}" STREQUAL "0" OR
   NOT stdout MATCHES "__sdf-h0__e-[^,]+__mn-0p001__mx-0p5__r31__a90__i1__q3__p1__f0__b-custom")
    message(FATAL_ERROR "SDF automatic epsilon/naming failed: ${stdout}\n${stderr}")
endif()
string(REGEX MATCH "Remesh diagnostics: ([^\r\n]+)" match "${stdout}")
set(auto_run "${CMAKE_MATCH_1}")
file(TO_CMAKE_PATH "${auto_run}" auto_run)
if(auto_run STREQUAL "" OR NOT EXISTS "${auto_run}/iter_001/result.off")
    message(FATAL_ERROR "Compact automatic path did not hold diagnostics: ${stdout}\n${stderr}")
endif()

function(expect_failure pattern)
    execute_process(COMMAND "${EXE}" "${input}" ${args} ${ARGN}
        RESULT_VARIABLE result OUTPUT_VARIABLE stdout ERROR_VARIABLE stderr TIMEOUT 60)
    if("${result}" STREQUAL "0" OR NOT stderr MATCHES "${pattern}")
        message(FATAL_ERROR "Expected ${pattern}: ${stdout}\n${stderr}")
    endif()
endfunction()
foreach(value -1 0 2.5 4294967296)
    expect_failure("sdf-rays" --sdf-rays "${value}")
endforeach()
foreach(value 0 180 nan 120junk)
    expect_failure("sdf-cone-angle" --sdf-cone-angle "${value}")
endforeach()
expect_failure("SDF options require" --field cgal-adaptive --sdf-rays 25)
expect_failure("Missing value" --sdf-rays)
expect_failure("finite" --epsilon nan)

set(open "${TEST_DIR}/open.off")
file(WRITE "${open}" "OFF\n3 1 0\n0 0 0\n1 0 0\n0 1 0\n3 0 1 2\n")
execute_process(COMMAND "${EXE}" "${open}" --field cgal-adaptive-sdf-thickness --iterations 0
    RESULT_VARIABLE result OUTPUT_VARIABLE stdout ERROR_VARIABLE stderr TIMEOUT 60)
if("${result}" STREQUAL "0" OR NOT stderr MATCHES "closed mesh")
    message(FATAL_ERROR "Open SDF input was not rejected before auto curvature: ${stderr}")
endif()
