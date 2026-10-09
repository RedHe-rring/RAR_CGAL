# Run on every build; configure_file only changes the header when metadata changes.
execute_process(COMMAND git -C "${SOURCE_DIR}" rev-parse HEAD
    RESULT_VARIABLE git_result OUTPUT_VARIABLE source_commit OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET)
set(dirty_json null)
if(git_result EQUAL 0)
    execute_process(COMMAND git -C "${SOURCE_DIR}" status --porcelain
        RESULT_VARIABLE status_result OUTPUT_VARIABLE source_status ERROR_QUIET)
    if(status_result EQUAL 0)
        set(dirty_json false)
        if(NOT source_status STREQUAL "")
            set(dirty_json true)
        endif()
    endif()
else()
    set(source_commit unknown)
endif()
configure_file("${SOURCE_DIR}/cmake/BuildInfo.h.in" "${OUTPUT}" @ONLY)
