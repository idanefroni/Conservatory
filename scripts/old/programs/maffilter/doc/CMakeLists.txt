# CMake script for MafFilter
# Authors:
#   Julien Dutheil
#   Francois Gindraud (2017)
# Created: 22/08/2009

# Builds info, html, pdf doc.
# Info doc is built and install as part of "all" if makeinfo is found.
# Html doc is proposed as a "html" optional target if makeinfo is found.
# Pdf doc is proposed as a "pdf" optional target if makeinfo AND texi2dvi are found.

find_program (MAKEINFO NAMES makeinfo texi2any DOC "makeinfo doc generator program")
if (NOT MAKEINFO)
  message (STATUS "makeinfo program not found: 'info' and 'html' target disabled (builds info/html doc)")
else (NOT MAKEINFO)
  message (STATUS "Found ${MAKEINFO}: 'info' and 'html' target enabled (builds info/html doc)")

  set (input ${CMAKE_CURRENT_SOURCE_DIR}/maffilter.texi)

  # Build and install info page
  set (output ${CMAKE_CURRENT_BINARY_DIR}/maffilter.info)
  add_custom_command (
    OUTPUT ${output}
    COMMAND ${MAKEINFO} --no-split -o ${output} ${input}
    DEPENDS ${input}
    COMMENT "Generating info page"
    VERBATIM
    )

  # Install, and have "info" built with "all" (install needs the file to be built)
  if (NOT COMPRESS_BIN)
    # Install uncompressed info page
    install (FILES ${output} DESTINATION ${CMAKE_INSTALL_INFODIR})
    add_custom_target (info ALL DEPENDS ${output})
  else ()
    # Compress and install compressed file
    set (compressed_ouput ${output}.${COMPRESS_EXT})
    add_custom_command (
      OUTPUT ${compressed_ouput}
      COMMAND ${COMPRESS_BIN} ${COMPRESS_ARGS} ${output} > ${compressed_ouput}
      DEPENDS ${output}
      COMMENT "Compressing info page"
      VERBATIM
      )
    install (FILES ${compressed_ouput} DESTINATION ${CMAKE_INSTALL_INFODIR})
    add_custom_target (info ALL DEPENDS ${compressed_ouput})
  endif ()

  # Also provide a "html" target that builds html doc (not installed, and not part of "all").
  set (output ${CMAKE_CURRENT_BINARY_DIR}/maffilter.html)
  set (makeinfo-css "http://www.w3.org/StyleSheets/Core/Steely")
  add_custom_command (
    OUTPUT ${output}
    COMMAND ${MAKEINFO} --html --css-ref=${makeinfo-css} --no-split -o ${output} ${input}
    DEPENDS ${input}
    COMMENT "Generating html doc"
    VERBATIM
    )
  add_custom_target (html DEPENDS ${output})

  # Also provide a "html-multipages" target that builds html doc (not installed, and not part of "all").
  set (output ${CMAKE_CURRENT_BINARY_DIR}/html/)
  set (makeinfo-css "http://www.w3.org/StyleSheets/Core/Steely")
  add_custom_command (
    OUTPUT ${output}
    COMMAND ${MAKEINFO} --html --no-headers --css-ref=${makeinfo-css} -o ${output} ${input}
    DEPENDS ${input}
    COMMENT "Generating html (multipages) doc"
    VERBATIM
    )
  add_custom_target (html-multipages DEPENDS ${output})


  # Provide a "pdf" target that builds pdf doc (not installed, not part of "all").
  find_program (TEXIDVI NAMES texi2dvi)
  if (TEXIDVI)
    message (STATUS "Found texi2dvi: 'pdf' target enabled (builds pdf doc)")
    set (output ${CMAKE_CURRENT_BINARY_DIR}/maffilter.pdf)
    add_custom_command (
      OUTPUT ${output}
      COMMAND ${MAKEINFO} --pdf --Xopt=--clean -o ${output} ${input}
      DEPENDS ${input}
      COMMENT "Generating pdf doc"
      VERBATIM
      )
    add_custom_target (pdf DEPENDS ${output})
  endif ()
endif ()
