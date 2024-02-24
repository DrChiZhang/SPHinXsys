find_package(pybind11 CONFIG REQUIRED)

pybind11_add_module(pysphinxsys
  main.cpp
)

# All the same postfix, otherwise the module name ist not the same as the exported name and python gets confused
set_target_properties(pysphinxsys PROPERTIES
        DEBUG_POSTFIX ""
        RELWITHDEBINFO_POSTFIX ""
        MINSIZEREL_POSTFIX "")

target_link_libraries(pysphinxsys PRIVATE sphinxsys_2d)