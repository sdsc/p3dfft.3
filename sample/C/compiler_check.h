#if (defined(__GNUC__) || defined(__GNUG__)) && !(defined(__clang__) || defined(__INTEL_COMPILER))
#	define COMPILER_DETECTED "GCC/G++"
#	define COMPILER_V_DETECTED (__VERSION__)

#elif defined(__clang__)
#	define COMPILER_DETECTED "Clang"
#	define COMPILER_V_DETECTED (__VERSION__)

#elif defined(__ICC) || defined(__INTEL_COMPILER)
#	define COMPILER_DETECTED (__VERSION__)
#	define COMPILER_V_DETECTED (__INTEL_COMPILER)

#elif defined(__HP_cc)
#	define COMPILER_DETECTED "HP-C"
#	define COMPILER_V_DETECTED (__HP_cc)

#elif defined(__HP_aCC)
#	define COMPILER_DETECTED "HP-aC"
#	define COMPILER_V_DETECTED (__HP_aCC)

#elif defined(__IBMC__) || defined(__IBMCPP__)
#	define COMPILER_DETECTED "IBM XL"
# define COMPILER_V_DETECTED (__XlC__)

#elif defined(_MSC_VER)
#	define COMPILER_DETECTED "Microsoft Visual Studio"
#	define COMPILER_V_DETECTED (_MSC_VER)

#elif defined(__PGI)
#	define COMPILER_DETECTED "PGI"
# define COMPILER_V_DETECTED (__PGIC__ * 10000 \
														 + __PGIC_MINOR__ * 100 \
														 + __PGIC_PATCHLEVEL__)

#elif defined(__SUNPRO_C)
#	define COMPILER_DETECTED "Oracle Solaris Studio"
#	define COMPILER_V_DETECTED (__SUNPRO_C)

#elif defined(__SUNPRO_CC)
#	define COMPILER_DETECTED "Oracle Solaris Studio"
#	define COMPILER_V_DETECTED (__SUNPRO_CC)

#elif defined(_CRAYC)
#	define COMPILER_DETECTED "Cray"
#	define COMPILER_V_DETECTED (_RELEASE * 100 \
														 + _RELEASE_MINOR)

#else
#	define COMPILER_DETECTED "Unknown"
# define COMPILER_V_DETECTED "Unknown"

#endif
