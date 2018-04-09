// See Bridson's Fluids Notes, section 4.3.1 Putting Them In Matrix-Vector Form (on page 31, or PDF page 43)
// for an explanation of this symmetric sparse matrix data structure.

#pragma once
#include "grid_data.h"

class GridDataMatrix 
{
public:
	GridDataMatrix() 
	{
		diag.initialize();
		plusI.initialize();
		plusJ.initialize();
		plusK.initialize();
	}
	GridData diag;
	GridData plusI;
	GridData plusJ;
	GridData plusK;
};