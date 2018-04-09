// Used for math operations and constants that aren't included in the C++ standard library.

#pragma once
#include "custom_output.h"
#include <cmath>
#include <cstdarg>
#include <iostream>

const double PI =					3.1415926535897932384626422832795028841971;
const double ONE_OVER_PI =			0.3183098861837906715377675267450287240689;
const double TWO_PI =				6.2831853071795864769252867665590057683943;
const double FOUR_PI =				12.566370614359172953850573533118011536788;
const double ONE_OVER_FOUR_PI =		0.0795774715459476678844418816862571810172;
const double E =					2.7182818284590452353602874713526624977572;

namespace BasicMath 
{
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