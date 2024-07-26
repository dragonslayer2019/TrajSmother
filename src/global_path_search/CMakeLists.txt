# cmake_minimum_required(VERSION 3.0.0)
# project(global_path_search VERSION 0.1.0 LANGUAGES C CXX)

# set(CMAKE_CXX_STANDARD 14)

# # 手动包含头文件路径
# include_directories("/usr/include")

# find_package(nlohmann_json REQUIRED)

# include(FetchContent)
# FetchContent_Declare(
#   googletest
#   URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
# )

# set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
# FetchContent_MakeAvailable(googletest)

# set(SOURCES
#     global_path_search.cpp
#     src/oct_cube.cpp
#     src/a_star.cpp
#     src/octree.cpp
#     src/vertice.cpp
#     src/zyb_octree.cpp
# )

# add_library(global_path_search ${SOURCES})

# target_include_directories(global_path_search
#     PRIVATE 
#         ${PROJECT_SOURCE_DIR}/include
# )

# include(CTest)
# enable_testing()

# set (TESTS
#     global_path_search_test.cpp
#     src/oct_cube_test.cpp
#     src/a_star_test.cpp
# )

# add_executable(unit_tests ${TESTS})

# target_include_directories(unit_tests
#     PRIVATE
#         ${PROJECT_SOURCE_DIR}/include
# )

# target_link_libraries(
#     unit_tests
#     global_path_search
#     GTest::gtest_main
# )

# include(GoogleTest)
# gtest_discover_tests(unit_tests)


# set(CPACK_PROJECT_NAME ${PROJECT_NAME})
# set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
# include(CPack)


cmake_minimum_required(VERSION 3.0.0)
project(global_path_search VERSION 0.1.0 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 14)

######################
## 找到catkin的依赖
find_package(catkin REQUIRED COMPONENTS
  roscpp
  std_msgs
  geometry_msgs
)

## 声明catkin的包
catkin_package(
  INCLUDE_DIRS include
  LIBRARIES global_path_search
  CATKIN_DEPENDS roscpp std_msgs geometry_msgs
)

## 指定头文件路径
include_directories(
  include
  ${catkin_INCLUDE_DIRS}
)

# 手动包含头文件路径
include_directories("/usr/include/nlohmann")

# 源文件列表
set(SOURCES
    global_path_search.cpp
    src/oct_cube.cpp
    src/a_star.cpp
    src/octree.cpp
    src/vertice.cpp
    src/zyb_octree.cpp
)

# 创建库
add_library(global_path_search ${SOURCES})

target_include_directories(global_path_search
    PRIVATE 
        ${PROJECT_SOURCE_DIR}/include
)

# 增加全局路径搜索节点
add_executable(global_path_search_node global_path_search_node.cpp)
target_link_libraries(global_path_search_node global_path_search ${catkin_LIBRARIES})

# 测试文件列表
set(TESTS
    global_path_search_test.cpp
    src/oct_cube_test.cpp
    src/a_star_test.cpp
)

# 创建测试可执行文件
add_executable(unit_tests ${TESTS})

target_include_directories(unit_tests
    PRIVATE
        ${PROJECT_SOURCE_DIR}/include
)

target_link_libraries(
    unit_tests
    global_path_search
    ${catkin_LIBRARIES}
    gtest
    gtest_main
)

include(GoogleTest)
gtest_discover_tests(unit_tests)

# 安装目标
install(TARGETS global_path_search_node
  RUNTIME DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
)

install(TARGETS global_path_search
  LIBRARY DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
  ARCHIVE DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
)

install(DIRECTORY include/${PROJECT_NAME}/
  DESTINATION ${CATKIN_PACKAGE_INCLUDE_DESTINATION}
)

install(DIRECTORY launch
  DESTINATION ${CATKIN_PACKAGE_SHARE_DESTINATION}
)
