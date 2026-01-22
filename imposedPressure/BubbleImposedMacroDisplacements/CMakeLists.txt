project(OperaHPCBubbleBis)

add_custom_command(
  OUTPUT "src/libBehaviour.so"
  COMMAND mfront --obuild --interface=generic ${PROJECT_BINARY_DIR}/Elasticity.mfront
  DEPENDS ${PROJECT_BINARY_DIR}/Elasticity.mfront
  COMMENT "Performing mfront code generation")

add_custom_target(lib-elasticity ALL
  DEPENDS "src/libBehaviour.so")

# build executable

file( GLOB_RECURSE SRCS *.cxx )

add_executable(BubbleImposedMacroDisplacements
  ${SRCS}
)

add_dependencies(BubbleImposedMacroDisplacements lib-elasticity)

target_link_libraries(BubbleImposedMacroDisplacements
 PRIVATE TFELMath TFELException
 PUBLIC mfem-mgis::MFEMMGIS)

target_include_directories(BubbleImposedMacroDisplacements
   PUBLIC
   ${TFEL_INCLUDE_PATH} 
   ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_directories(BubbleImposedMacroDisplacements
   PUBLIC
   ${TFEL_LIBRARY_PATH}) 


# copy and install

file(COPY Elasticity.mfront DESTINATION ${PROJECT_BINARY_DIR})

install(TARGETS BubbleImposedMacroDisplacements DESTINATION bin)
