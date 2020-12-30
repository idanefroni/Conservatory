# CMake script for the CoMap package.
# Authors:
#   Julien Dutheil
#   Francois gindraud (2017)
# Created: 22/08/2009

# Build manpages.
# In practice, they are just compressed from the text files using COMPRESS_PROGRAM
# Manpages are built and installed as part of "all" if a COMPRESS_PROGRAM is found.

# Take all manpages files in the directory
file (GLOB manpage_files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.1)

if (NOT COMPRESS_BIN)
  # Just install manpages from source
  foreach (manpage_file ${manpage_files})
    install (FILES ${CMAKE_CURRENT_SOURCE_DIR}/${manpage_file} DESTINATION ${CMAKE_INSTALL_MANDIR}/man1)
  endforeach (manpage_file)
else ()
  # Create a list of manpage targets
  set (manpage-targets)

  foreach (manpage_file ${manpage_files})
    # Compress manpage, install, add to manpage target list
    set (input ${CMAKE_CURRENT_SOURCE_DIR}/${manpage_file})
    set (output ${CMAKE_CURRENT_BINARY_DIR}/${manpage_file}.${COMPRESS_EXT})
    add_custom_command (
      OUTPUT ${output}
      COMMAND ${COMPRESS_BIN} ${COMPRESS_ARGS} ${input} > ${output}
      DEPENDS ${input}
      COMMENT "Compressing manpage ${manpage_file}"
      VERBATIM
      )
    install (FILES ${output} DESTINATION ${CMAKE_INSTALL_MANDIR}/man1)
    list (APPEND manpage-targets ${output})
    unset (input)
    unset (output)
  endforeach (manpage_file)

  # Add target "man", built with "all" (needed because install will fail if not built).
  add_custom_target (man ALL DEPENDS ${manpage-targets})
endif ()
