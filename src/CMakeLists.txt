find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
find_package(Eigen3 REQUIRED)
add_library(eigen INTERFACE IMPORTED)

set(SOURCE_FILES sources/main.cpp
		 sources/smoke_sim.cpp
		 headers/smoke_sim.h
		 sources/mac_grid.cpp
		 headers/mac_grid.h
		 sources/vec.cpp
		 headers/vec.h
		 headers/matrix.h
		 sources/grid_data.cpp
		 headers/grid_data.h
		 headers/grid_data_matrix.h
		 sources/camera.cpp
		 headers/camera.h
		 headers/custom_output.h
		 headers/clock.h
		 headers/timer.h
		 sources/fps.cpp
		 headers/fps.h
		 sources/constants.cpp
		 headers/constants.h
		 headers/basic_math.h
		 headers/open_gl_headers.h)

add_SMOKE_executable(smoke ${SOURCE_FILES})
include_directories( ${OPENGL_INCLUDE_DIR}  ${GLUT_INCLUDE_DIRS} )
target_link_libraries(smoke ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES} )
target_include_directories(smoke SYSTEM PUBLIC ${EIGEN3_INCLUDE_DIR})
target_link_libraries(smoke eigen)
target_link_libraries(smoke partio)
