add_library(libsources ${stokessources} ${matredsources})
target_link_libraries(libsources NeoPZ::pz)