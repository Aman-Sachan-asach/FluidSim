#pragma once

#ifdef __linux__
    #include <GL/gl.h>
    #include <GL/glu.h>
    #include <GL/glut.h>
#elif _WIN32
    #include "../../GL/GL.h"
    #include "../../GL/GLU.h"
    #include "../../GL/glut.h"
#elif __APPLE__
    #include <OpenGL/gl.h>
    #include <OpenGL/glu.h>
    #include <GLUT/glut.h>
#endif

#include "custom_output.h"
#include "../external/vec.h"
#include <cmath>
#include <cstdarg>
#include <iostream>

#ifdef _DEBUG
const int theDim[3] = {4, 4, 1};
#else
const int theDim[3] = {32, 32, 1};
#endif

// Simulation Controls
const int theMillisecondsPerFrame          = 10;
const double theCellSize                   = 0.5;
const double theAirDensity                 = 1.0;
const double theBuoyancyAlpha              = 0.08; // Gravity's effect on the smoke particles.
const double theBuoyancyBeta               = 0.37; // Buoyancy's effect due to temperature difference.	
const double theBuoyancyAmbientTemperature = 0.0;  // Ambient temperature.
const double theVorticityEpsilon = 0.10;


// Used for math operations and constants that aren't included in the C++ standard library.
const double PI =					3.1415926535897932384626422832795028841971;
const double ONE_OVER_PI =			0.3183098861837906715377675267450287240689;
const double TWO_PI =				6.2831853071795864769252867665590057683943;
const double FOUR_PI =				12.566370614359172953850573533118011536788;
const double ONE_OVER_FOUR_PI =		0.0795774715459476678844418816862571810172;
const double E =					2.7182818284590452353602874713526624977572;

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

	// std::log2 is not avaiable in all compilers.
	inline double log2(double n)
	{
		return std::log(n) / std::log(2.0);
	}

	// std::isnan is not avaiable in all compilers.
	inline bool isNaN(double n)
	{
		return (n != n);
	}
}