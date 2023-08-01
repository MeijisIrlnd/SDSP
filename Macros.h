#pragma once 
#ifndef SDSP_UNUSED 
#define SDSP_UNUSED [[maybe_unused]]
#endif

#ifndef SDSP_INLINE 
#ifdef _MSC_VER
#define SDSP_INLINE __forceinline
#else 
#define SDSP_INLINE __attribute__((always_inline))
#endif
#endif

#ifndef SDSP_NODISCARD
#define SDSP_NODISCARD [[nodiscard]]
#endif

#ifndef SDSP_NOINLINE
#ifdef _MSC_VER
#define SDSP_NOINLINE __declspec(noinline)
#else
#define SDSP_NOINLINE __attribute__((noinline))
#endif
#endif
