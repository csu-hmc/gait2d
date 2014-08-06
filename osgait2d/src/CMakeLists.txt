set(FILES Gait2D.cpp)
add_executable(gait2d ${FILES})
target_link_libraries(gait2d ${OPENSIMSIMBODY_LIBRARIES} yaml-cpp)
