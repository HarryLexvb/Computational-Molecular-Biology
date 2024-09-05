cmake_minimum_required(VERSION 3.28)
project(test)

set(CMAKE_CXX_STANDARD 26)

add_executable(test main.cpp)

# Buscar y vincular OpenMP
find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
    target_link_libraries(test PRIVATE OpenMP::OpenMP_CXX)
    # Asegúrate de que OpenMP esté habilitado en el compilador
    set_target_properties(test PROPERTIES
            COMPILE_OPTIONS "${OpenMP_CXX_FLAGS}"
            LINK_FLAGS "${OpenMP_CXX_FLAGS}"
    )
endif()
