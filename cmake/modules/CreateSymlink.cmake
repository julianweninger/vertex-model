# This function creates a symlink
# Call it with two arguments: The target path and the link path
function(create_symlink)
    if (NOT ${ARGC} EQUAL 2)
        message(SEND_ERROR "Called 'create_symlink' with wrong number of "
                           "arguments")
    endif ()

    set(TARGET ${ARGV0})
    set(LINK ${ARGV1})

    if (CMAKE_VERSION VERSION_GREATER_EQUAL 3.14)
        file (CREATE_LINK ${TARGET} ${LINK} SYMBOLIC)
    else ()
        execute_process(
            COMMAND ${CMAKE_COMMAND} -E create_symlink ${TARGET} ${LINK}
            RESULT_VARIABLE RETURN_VALUE
            OUTPUT_QUIET
        )
        if (NOT RETURN_VALUE STREQUAL 0)
            message(SEND_ERROR "Error creating a symlink: ${RETURN_VALUE}")
        endif ()
    endif ()
endfunction()
