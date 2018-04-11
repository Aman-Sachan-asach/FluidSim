#pragma once
//OpenGl Headers
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>

#include <cmath>
#include <cstdarg>
#include <iostream>
#include <sstream>
#include "../external/vec.h"

#ifdef _DEBUG
const int theDim[3] = {4, 4, 1};
#else
const int theDim[3] = {12, 12, 1};
#endif

////////////////////////////////////////////////////////////

// Simulation Controls
const int    millisecondsPerFrame       = 10;
const double dt                         = 0.04;

const double gridCellSize               = 0.5;
const double One_By_gridCellSize        = 1.0/gridCellSize;

const double fluidDensity               = 1.0;
const double buoyancyAlpha              = 0.08; // Gravity's effect on the smoke particles.
const double buoyancyBeta               = 0.37; // Buoyancy's effect due to temperature difference.	
const double buoyancyAmbientTemperature = 0.0;  // Ambient temperature.
const double vorticityEpsilon 		    = 0.10;

const double solidBoundaryConstant = (fluidDensity * gridCellSize) / dt; // why not squared instead of just gridCellSize

// Used for math operations and constants that aren't included in the C++ standard library.
const double PI =					3.1415926535897932384626422832795028841971;
const double ONE_OVER_PI =			0.3183098861837906715377675267450287240689;
const double TWO_PI =				6.2831853071795864769252867665590057683943;
const double FOUR_PI =				12.566370614359172953850573533118011536788;
const double ONE_OVER_FOUR_PI =		0.0795774715459476678844418816862571810172;
const double E =					2.7182818284590452353602874713526624977572;
const double DOUBLE_EPSILON =       0.00000000000000000000000000001; //std::numeric_limits<double>::epsilon()

/////////////////////TODO///////////////////////////////////////

// Print Utilities
template <typename T>
inline void Print (T const& X) 
{
	std::ostringstream stream;
	stream << X;
	std::cout << stream.str();
	fflush(stdout);
}

template <typename T>
inline void PrintLine (T const& X) 
{
	std::ostringstream stream;
	stream << X << std::endl;
	std::cout << stream.str();
	fflush(stdout); 
}

////////////////////////////////////////////////////////////
//Some compilers dont supply things like std::isNaN or std::log2
namespace Globals
{
	inline vec4 Lerp(vec4 a, vec4 b, double t)
	{
		return (1-t)*a + t*b;
	}

	// Proper int mod with an always-positive result (unlike %).
	inline unsigned int mod(int x, int y)
	{
		int result = x % y;
		if (result < 0) result += y;
		return result;
	}

	// Proper double mod with an always-positive result (unlike fmod).
	inline double mod(double x, double y)
	{
		return x - y * std::floor(x / y);
	}

	inline double radiansToDegrees(double radians)
	{
		return radians * 180.0 / PI;
	}

	inline double degreesToRadians(double degrees)
	{
		return degrees / 180.0 * PI;
	}

	inline double average(double n1, double n2)
	{
		return ((n1 + n2) / 2);
	}
	
	inline double round(double n)
	{
		return std::floor(n + 0.5);
	}

	inline double log2(double n)
	{
		return std::log(n) / std::log(2.0);
	}

	inline bool isNaN(double n)
	{
		return (n != n);
	}
}