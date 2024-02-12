# SDSP
Work-In-Progress header only DSP helpers, specifically for the JUCE framework. 


# Documentation
[todo]

# Features 

# Usage 
The library itself is header only, but can also be used as a JUCE Module.
## Adding to your project:
```
mkdir modules
cd modules
git submodule add https://github.com/MeijisIrlnd/SDSP.git
```
If you're using the juce cmake api, after your call to `juce_add_target(...)`, adding 
`juce_add_module(${CMAKE_CURRENT_SOURCE_DIR}/modules/SDSP)` will make the module available. To actually use the module, add `SDSP` to your target link libraries call: 

```
cmake_minimum_required(VERSION ...)
juce_add_plugin(Foo ....)
juce_add_module(${CMAKE_CURRENT_SOURCE_DIR}/modules/SDSP)
target_sources(Foo .....)
target_link_libraries(Foo PRIVATE SDSP ...)
```



