#pragma once
#ifdef JDE_Matlab
	#ifdef _MSC_VER
		#define JDE_MATLAB_VISIBILITY __declspec( dllexport )
	#else
		#define JDE_MATLAB_VISIBILITY __attribute__((visibility("default")))
	#endif
#else 
	#ifdef _MSC_VER
		#define JDE_MATLAB_VISIBILITY __declspec( dllimport )
		#if NDEBUG
			#pragma comment(lib, "Jde.Math.IO.MatLab.lib")
		#else
			#pragma comment(lib, "Jde.Math.IO.MatLab.lib")
		#endif
	#else
		#define JDE_MATLAB_VISIBILITY
	#endif
#endif 

