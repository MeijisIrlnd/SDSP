# SDSP
A perenially-work-in-progress header only DSP helper lib, specifically for the JUCE framework. SDSP is where oddities and occasionally useful bits of code I write inevitably end up, for reuse across multiple projects, or for people to either use themselves, or scavenge nice parts from for their own DSP Libraries or projects. The code is in the slow process of being cleaned up and refactored (ideally in a non-breaking way), so expect some sharp edges and bear-traps along the way.


# Documentation
[todo: Learn to use doxygen, this would actually be really nice to have]

# Features 
## Filters
- (Cascadable) Biquad Filters
- A helper class for interpolating Biquad coefficients
- Implementations of the RBJ Cookbook filters
- First and second order Allpass filters, and a ModulatedAPF class for modulating the allpass' coefficients (useful in phasers for example)
- An implementation of the Cytomic TPT State Variable Filter
- A Bessel filter
- A Single Pole LPF for smoothing / very rough filtering
- A Tilt EQ (low and high shelves)
## Utility 
- Helpers for loading audio from either binary or disk, with optional on-load resampling.
- A fast Hadamard and Householder matrix implementation
- A curiosity found in the depths of the internet called an Ornstein Ulenbech process implemetation
- A few windowing functions [todo, more]
## Misc
- Shape configurable oscillators, implemented as both trig-based, and wavetable, with polyblep.
- A clamp utility function called protectYourEars, to avoid hearing damage and other fun beartraps.
- A simple (non thread-safe) FIFO
- Various configurations of circular buffers
- Some buffer helpers (BufferTransfer)

# Usage 
The library itself is header only, but can also be used as a JUCE Module.
## Adding to your project:
```
mkdir modules
cd modules
git submodule add https://github.com/MeijisIrlnd/SDSP.git
```
## Juce Module
### CMake
If you're using the juce cmake api, after your call to `juce_add_target(...)`, adding 
`juce_add_module(${CMAKE_CURRENT_SOURCE_DIR}/modules/SDSP)` will make the module available. To actually use the module, add `SDSP` to your target link libraries call - A pseudocode example of what your CMakeLists.txt might look like is included below:

```
cmake_minimum_required(VERSION ...)
juce_add_plugin(Foo ....)
juce_add_module(${CMAKE_CURRENT_SOURCE_DIR}/modules/SDSP)
target_sources(Foo .....)
target_link_libraries(Foo PRIVATE SDSP ...)
```
### Projucer
If you're using Projucer, open your project's jucer file, and navigate to the exporters page. Clicking the plus icon in that panel will bring up a dialog, select "add a module from a specified folder". In the file explorer, point it to the root directory of SDSP, and it will become available in your project. You can also have SDSP in your User Modules global path specified in the Projucer, and it should then show up in the list of available modules in the "add a module" dialog.

## Include
As SDSP is a header only library, you can also skip the juce module setup entirely, and just ensure that SDSP is on your project's include path. If you're using cmake, that should be as simple as `target_include_directories(Foo PRIVATE /Path/To/SDSP)`. 
If you're using Projucer, you'd just add the path to SDSP in the "Header Search Paths" field, in the exporter and target you want to use SDSP with. 

# Tests
The test coverage is a work-in-progress as well, but there *are* some in the `Tests` directory. The CMakeLists.txt in here simply sets an `SDSP_TEST_SOURCE` variable to the test sources, and expects the parent project to define the test target itself [todo obviously, not this]. Tests are written with Catch2!

# TODO: 
- SDSP_HEADER_ONLY as an option, allowing you to build as a static lib to help with compilation times as an alternative
- FetchContent interface
- Tests as a part of this repo as opposed to dumping responsibility on the consumer
- More test coverage in general
  



