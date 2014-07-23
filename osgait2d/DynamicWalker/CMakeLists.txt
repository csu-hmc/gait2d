set(FILES DynamicWalkerBuildModel.cpp)
add_executable(osimDynamicWalker ${FILES})
target_link_libraries(osimDynamicWalker ${OPENSIMSIMBODY_LIBRARIES} yaml-cpp)
