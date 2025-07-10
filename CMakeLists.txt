##### Library ######
add_library(
  SolarSystemTools STATIC
  src/BinnedExposureSun.cxx
  src/CosineBinner2D.cxx
  src/ExposureCubeSun.cxx
  src/ExposureSun.cxx
  src/FitsSolarProfile.cxx
  src/HealpixExposureSun.cxx
  src/LikeExposureSun.cxx
  src/SolarProfile.cxx
  src/SolarTemplate.cxx
)
target_link_libraries(
  SolarSystemTools
  PUBLIC astro
  tip
  map_tools
  healpix
  st_facilities
  facilities
  st_app
  irfLoader
  evtbin
  Likelihood
)
target_include_directories(
  SolarSystemTools PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:>
)

###### Executables ######
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

add_executable(gtltcubesun src/makeSolarExposureCube/makeSolarExposureCube.cxx)
add_executable(gtltsumsun src/gtltsumsun/gtltsumsun.cxx)
add_executable(gtexpcubesun src/gtexpcubesun/gtexpcubesun.cxx)
add_executable(gtexphpsun src/gtexphpsun/gtexphpsun.cxx)
add_executable(gtsuntemp src/gtsuntemp/gtsuntemp.cxx)

target_link_libraries(gtltcubesun PRIVATE SolarSystemTools)
target_link_libraries(gtltsumsun PRIVATE SolarSystemTools)
target_link_libraries(gtexpcubesun PRIVATE SolarSystemTools CLHEP::VectorS)
target_link_libraries(gtexphpsun PRIVATE SolarSystemTools)
target_link_libraries(gtsuntemp PRIVATE SolarSystemTools)

###############################################################
# Installation
###############################################################
install(DIRECTORY SolarSystemTools DESTINATION ${FERMI_INSTALL_INCLUDEDIR})
install(DIRECTORY data/ DESTINATION ${FERMI_INSTALL_DATADIR})
install(DIRECTORY pfiles/ DESTINATION ${FERMI_INSTALL_PFILESDIR})

install(
  TARGETS SolarSystemTools gtltcubesun gtltsumsun gtexpcubesun gtexphpsun gtsuntemp
  # EXPORT fermiTargets
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
