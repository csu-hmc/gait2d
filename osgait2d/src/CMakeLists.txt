add_executable(gait2d Gait2D.cpp)
target_link_libraries(gait2d ${OPENSIMSIMBODY_LIBRARIES} yaml-cpp)

add_executable(control GainScheduledController.cpp)
target_link_libraries(control ${OPENSIMSIMBODY_LIBRARIES})
