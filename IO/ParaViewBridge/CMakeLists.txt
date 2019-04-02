cmake_minimum_required(VERSION 2.8)

project(ParaViewBridge)
set(PluginName ParaViewBridge )
set(ParaViewBridge_VERSION_MAJOR 1 )
set(ParaViewBridge_VERSION_MINOR 7 )
set(ParaViewBridge_VERSION_PATCH 0 )


function(ESCAPEXML pythontoescape inputfile outputfile)
  add_custom_command(OUTPUT  ${CMAKE_CURRENT_BINARY_DIR}/Generated/${outputfile}
  		  COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/Generated
		  COMMAND python escapexmlstring.py -f ${pythontoescape}  -r ${inputfile} -o ${CMAKE_CURRENT_BINARY_DIR}/Generated/${outputfile}
		  DEPENDS escapexmlstring.py ${CMAKE_SOURCE_DIR}/${inputfile} ${pythontoescape}
		  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
		  VERBATIM)
endfunction()

function(ESCAPEXMLII pythontoescapeI pythontoescapeII inputfile outputfile)
  add_custom_command(OUTPUT  ${CMAKE_CURRENT_BINARY_DIR}/Generated/${outputfile}
  		  COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/Generated
		  COMMAND python escapexmlstring.py -f ${pythontoescapeI} ${pythontoescapeII} -r ${inputfile} -o ${CMAKE_CURRENT_BINARY_DIR}/Generated/${outputfile}
		  DEPENDS escapexmlstring.py ${CMAKE_SOURCE_DIR}/${inputfile} ${pythontoescapeI} ${pythontoescapeII}
		  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
		  VERBATIM)
endfunction()

function(ESCAPEXMLIII pythontoescapeI pythontoescapeII pythontoescapeIII inputfile outputfile)
  add_custom_command(OUTPUT  ${CMAKE_CURRENT_BINARY_DIR}/Generated/${outputfile}
  		  COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/Generated
		  COMMAND python escapexmlstring.py -f ${pythontoescapeI} ${pythontoescapeII} ${pythontoescapeIII} -r ${inputfile} -o ${CMAKE_CURRENT_BINARY_DIR}/Generated/${outputfile}
		  DEPENDS escapexmlstring.py ${CMAKE_SOURCE_DIR}/${inputfile} ${pythontoescapeI} ${pythontoescapeII}
		  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
		  VERBATIM)
endfunction()


INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})
# if compiling inside of outside paraview
IF (ParaView_SOURCE_DIR)
  INCLUDE_DIRECTORIES(${VTK_INCLUDE_DIRS})
  set(ParaView_OK TRUE)
ELSE (ParaView_SOURCE_DIR)
  FIND_PACKAGE(ParaView  QUIET PATHS ${ParaView_DIR} )
  INCLUDE(${PARAVIEW_USE_FILE})
endif(ParaView_SOURCE_DIR)

set(SERVER_MANAGER_SOURCES_files )

######## Python readers and filters ##########
ESCAPEXMLIII(UniversalReader.py  UniversalReaderExt.txt UniversalReaderRequest.py UniversalReader.xml.in UniversalReader.xml)
ESCAPEXMLII(UniversalWriter.py  UniversalWriterExt.txt UniversalWriter.xml.in UniversalWriter.xml)

##############################################
set(SERVER_MANAGER_XML_files
${CMAKE_CURRENT_BINARY_DIR}/Generated/UniversalReader.xml
${CMAKE_CURRENT_BINARY_DIR}/Generated/UniversalWriter.xml
 )

#######################################################################################

ADD_PARAVIEW_PLUGIN(
  ${PluginName}
  "${ParaViewBridge_VERSION_MAJOR}.${ParaViewBridge_VERSION_MINOR}.${ParaViewBridge_VERSION_PATCH}"
  GUI_INTERFACES  ${IFACES}
  SOURCES
      ${MOC_SRCS}
      ${IFACE_SRCS}
      ${SOURCES_files}
  SERVER_MANAGER_XML
      ${SERVER_MANAGER_XML_files}
  SERVER_MANAGER_SOURCES
      ${SERVER_MANAGER_SOURCES_files}
  REQUIRED_ON_SERVER
 )

install(TARGETS ${PluginName} LIBRARY DESTINATION lib/paraview-${PARAVIEW_VERSION_MAJOR}.${PARAVIEW_VERSION_MINOR}/plugins)

