#include "../headers/mac_grid.h"

// Globals
MACGrid target;

// NOTE: x -> cols, z -> rows, y -> stacks
MACGrid::RenderMode MACGrid::theRenderMode = SHEETS;
bool MACGrid::theDisplayVel = false;//true

#define FOR_EACH_CELL \
	for(int k = 0; k < gridDim[MACGrid::Z]; k++)  \
		for(int j = 0; j < gridDim[MACGrid::Y]; j++) \
			for(int i = 0; i < gridDim[MACGrid::X]; i++) 

#define FOR_EACH_CELL_REVERSE \
	for(int k = gridDim[MACGrid::Z] - 1; k >= 0; k--)  \
		for(int j = gridDim[MACGrid::Y] - 1; j >= 0; j--) \
			for(int i = gridDim[MACGrid::X] - 1; i >= 0; i--) 

#define FOR_EACH_FACE \
	for(int k = 0; k < gridDim[MACGrid::Z]+1; k++) \
		for(int j = 0; j < gridDim[MACGrid::Y]+1; j++) \
			for(int i = 0; i < gridDim[MACGrid::X]+1; i++) 

MACGrid::MACGrid()
{
	reset();
}

MACGrid::MACGrid(const MACGrid& orig): mU(orig.mU), mV(orig.mV), mW(orig.mW), mP(orig.mP), mD(orig.mD), mT(orig.mT)
{
	reset();
}

MACGrid& MACGrid::operator=(const MACGrid& orig)
{
	if (&orig == this)
	{
	  return *this;
	}
	mU = orig.mU;
	mV = orig.mV;
	mW = orig.mW;
	mP = orig.mP;
	mD = orig.mD;
	mT = orig.mT;   

	return *this;
}

void MACGrid::reset()
{
	mU.initialize();
	mV.initialize();
	mW.initialize();
	mP.initialize();
	mD.initialize();
	mT.initialize(0.0);

	calculateAMatrix();
	calculatePreconditioner(AMatrix);
}

void MACGrid::updateSources()
{
	// Set initial values for density, temperature, velocity
	// Sources shooting smoke Horizontally
	// for(int sourceIndex=0; sourceIndex<2; sourceIndex++)
	// {
	// 	for(int i=sourcePosMin[sourceIndex][0]; i<sourcePosMax[sourceIndex][0]; i++)
	// 	{
	// 		for(int j=sourcePosMin[sourceIndex][1]; j<sourcePosMax[sourceIndex][1]; j++)
	// 		{
	// 			for(int k=sourcePosMin[sourceIndex][2]; k<sourcePosMax[sourceIndex][2]; k++)
	// 			{
	// 				mU(i+(sign[sourceIndex]*1),j,k) = sourceVelocity[sourceIndex][0];
	// 				mU(i+(sign[sourceIndex]*2),j,k) = sourceVelocity[sourceIndex][0];

	// 				// mV(i,j+1,k) = sourceVelocity[sourceIndex][1];
	// 				// mV(i,j+2,k) = sourceVelocity[sourceIndex][1];

	// 				mD(i,j,k) = fluidDensity;
	// 				mT(i,j,k) = 1.0;
	// 			}
	// 		}
	// 	}
	// }

	// Sources shooting smoke Vertically
	for(int i=sourcePosMin[2][0]; i<sourcePosMax[2][0]; i++)
	{
		for(int j=sourcePosMin[2][1]; j<sourcePosMax[2][1]; j++)
		{
			for(int k=sourcePosMin[2][2]; k<sourcePosMax[2][2]; k++)
			{
				mV(i,j-1,k) = sourceVelocity[2][1];
				mV(i,j-2,k) = sourceVelocity[2][1];

				mD(i,j  ,k) = fluidDensity;
				mT(i,j  ,k) = 1.0;
			}
		}
	}


	// Refresh particles in source.
	for(int sourceIndex=2; sourceIndex<numSources; sourceIndex++)
	{
		for(int i=sourcePosMin[sourceIndex][0]; i<sourcePosMax[sourceIndex][0]; i++)
		{
			for(int j=sourcePosMin[sourceIndex][1]; j<sourcePosMax[sourceIndex][1]; j++)
			{
				for (int k=sourcePosMin[sourceIndex][2]; k<=sourcePosMax[sourceIndex][2]; k++) 
				{
					vec3 cell_center(gridCellSize*(i+0.5), gridCellSize*(j+0.5), gridCellSize*(k+0.5));
					for(int p=0; p<10; p++)
					{
						double a = ((float) rand() / RAND_MAX - 0.5) * gridCellSize;
						double b = ((float) rand() / RAND_MAX - 0.5) * gridCellSize;
						double c = ((float) rand() / RAND_MAX - 0.5) * gridCellSize;
						vec3 shift(a, b, c);
						vec3 xp = cell_center + shift;
						rendering_particles.push_back(xp);
					}
				}
			}
		}
	}
}

//----------------------------------------------
// Advection, External Forces, and Projection //
//----------------------------------------------

/* 
	Advection is being calculated in a semi-lagrangian fashion:
		We are trying to calculate the value of some quantity 'q' at some position 'currentPosition'
		-- Determine the currentPosition

		Treat the quantity q as being stored in some imaginary particle (particles are a Lagrangian concept). 
		Advection just says we are going to move the quantity q using some imaginary particle from its previous 
		position to its new position (think diffusion).
		-- "currState = prevState"
	
	Because Advection deals with imaginary particles and not actual particles, when we retrieve the 'prevState'
	we get an averaged out value because every timestep we average quantities during the projection face when we
	calculate the velocity at the center of the cell
*/ 

void MACGrid::advectVelocity()
{
	// CHECK Section
	// Calculate new velocities and store in target and also in mU, mV, mW

	FOR_EACH_FACE 
	{
		if (isValidFace(MACGrid::X, i, j, k)) 
		{
			vec3 currentPosition = getFacePosition(MACGrid::X, i, j, k);
			vec3 rewoundPosition = getRewoundPosition(currentPosition);
			vec3 newVelocity = getVelocity(rewoundPosition);
			target.mU(i,j,k) = newVelocity[0];
		}
		if (isValidFace(MACGrid::Y, i, j, k)) 
		{
			vec3 currentPosition = getFacePosition(MACGrid::Y, i, j, k);
			vec3 rewoundPosition = getRewoundPosition(currentPosition);
			vec3 newVelocity = getVelocity(rewoundPosition);
			target.mV(i,j,k) = newVelocity[1];
		}
		if (isValidFace(MACGrid::Z, i, j, k)) 
		{
			vec3 currentPosition = getFacePosition(MACGrid::Z, i, j, k);
			vec3 rewoundPosition = getRewoundPosition(currentPosition);
			vec3 newVelocity = getVelocity(rewoundPosition);
			target.mW(i,j,k) = newVelocity[2];
		}
	}

	mU = target.mU;
	mV = target.mV;
	mW = target.mW;
}

void MACGrid::advectTemperature()
{
	// CHECK Section
	// Calculate new temp and store in target and mT

	FOR_EACH_CELL 
	{
		vec3 currentPosition = getCenter(i,j,k);
		vec3 rewoundPosition = getRewoundPosition(currentPosition);
		double newTemperature = getTemperature(rewoundPosition);
		target.mT(i,j,k) = newTemperature;
	}

	mT = target.mT;
}

void MACGrid::advectDensity()
{
	// CHECK Section
	// Calculate new densities and store in target and mD

	FOR_EACH_CELL 
	{
		vec3 currentPosition = getCenter(i,j,k);
		vec3 rewoundPosition = getRewoundPosition(currentPosition);
		double newDensity = getDensity(rewoundPosition);
		target.mD(i,j,k) = newDensity;
	}

	mD = target.mD;
}

void MACGrid::advectRenderingParticles() 
{
	rendering_particles_vel.resize(rendering_particles.size());
	for (size_t p = 0; p < rendering_particles.size(); p++) 
	{
		vec3 currentPosition = rendering_particles[p];
		vec3 currentVelocity = getVelocity(currentPosition);
		vec3 nextPosition = currentPosition + currentVelocity * dt;
		vec3 clippedNextPosition = clipToGrid(nextPosition, currentPosition);
		// Keep going...
		vec3 nextVelocity = getVelocity(clippedNextPosition);
		vec3 averageVelocity = (currentVelocity + nextVelocity) / 2.0;
		vec3 betterNextPosition = currentPosition + averageVelocity * dt;
		vec3 clippedBetterNextPosition = clipToGrid(betterNextPosition, currentPosition);
		rendering_particles[p] = clippedBetterNextPosition;
		rendering_particles_vel[p] = averageVelocity;
	}
}

void MACGrid::addExternalForces()
{
	computeBouyancy();
	computeVorticityConfinement();
}

void MACGrid::setPressureHighLow( int& i, int& j, int& k, const GridData& p, vec3& pLow, vec3& pHigh )
{
	if (isValidFace(MACGrid::X, i, j, k)) 
	{
		if (i-1 >= 0) 
		{
			pLow[0] = p(i-1,j,k);
		}

		if (i < gridDim[MACGrid::X]) 
		{
			pHigh[0] = p(i,j,k);
		}

		if (i-1 < 0) 
		{
			pLow[0] = pHigh[0] - solidBoundaryConstant * (mU(i,j,k) - 0);
		}

		if (i >= gridDim[MACGrid::X]) 
		{
			pHigh[0] = pLow[0] + solidBoundaryConstant * (mU(i,j,k) - 0);
		}

	}
	if (isValidFace(MACGrid::Y, i, j, k)) 
	{
		if (j-1 >= 0) 
		{
			pLow[1] = p(i,j-1,k);
		}

		if (j < gridDim[MACGrid::Y]) 
		{
			pHigh[1] = p(i,j,k);
		}

		if (j-1 < 0) 
		{
			pLow[1] = pHigh[1] - solidBoundaryConstant * (mV(i,j,k) - 0);
		}

		if (j >= gridDim[MACGrid::Y]) 
		{
			pHigh[1] = pLow[1] + solidBoundaryConstant * (mV(i,j,k) - 0);
		}
	}
	if (isValidFace(MACGrid::Z, i, j, k)) 
	{
		if (k-1 >= 0) 
		{
			pLow[2] = p(i,j,k-1);
		}

		if (k < gridDim[MACGrid::Z]) 
		{
			pHigh[2] = p(i,j,k);
		}

		if (k-1 < 0) 
		{
			pLow[2] = pHigh[2] - solidBoundaryConstant * (mW(i,j,k) - 0);
		}

		if (k >= gridDim[MACGrid::Z]) 
		{
			pHigh[2] = pLow[2] + solidBoundaryConstant * (mW(i,j,k) - 0);
		}
	}
}

double MACGrid::calcDivergence( int& i, int& j, int& k )
{
	// Construct the vector of divergences d:
	double vel_LowX  = mU(i,j,k);
	double vel_HighX = mU(i+1,j,k);
	double vel_LowY  = mV(i,j,k);
	double vel_HighY = mV(i,j+1,k);
	double vel_LowZ  = mW(i,j,k);
	double vel_HighZ = mW(i,j,k+1);

	// Use 0 for solid boundary velocities:
	if (i == 0)
	{
		vel_LowX = 0.0;
	}
	if (j == 0)
	{
		vel_LowY = 0.0;
	}
	if (k == 0) 
	{
		vel_LowZ = 0.0;
	}

	if (i+1 == gridDim[MACGrid::X])
	{
		vel_HighX = 0.0;
	} 
	if (j+1 == gridDim[MACGrid::Y])
	{
		vel_HighY = 0.0;
	}
	if (k+1 == gridDim[MACGrid::Z]) 
	{
		vel_HighZ = 0.0;
	}

	return ((vel_HighX - vel_LowX) + (vel_HighY - vel_LowY) + (vel_HighZ - vel_LowZ)) / gridCellSize;
}

void MACGrid::project()
{
	/*
		Resource: Section 4.3 on page 29 of Bridson's 2007 SIGGRAPH fluid course notes.
		Resource: Fig 4.1 on page 34 of Bridson's 2007 SIGGRAPH fluid course notes.

		The projection step is used to calculate the new velocity at every grid cell (u^(n+1)) and 
		store an interpolated value at every face of the grid cell.
		This velocity that we compute for every grid cell is supposed to be divergence free to maintain
		the incompressibility of the fluid (which is an assumption we use to heavily simply the navier 
		stokes equations).

		The new velocity u^(n+1) is calculated using 'pressure projection' which we solve using the 
		Conjugate Gradient Algorithm:
			Solve Ap = d for pressure, where p is pressure, 
											 A is the matrix of the divergence of the gradient,
										 and d is the divergence

			MAp = Md, where M is very close to the inverse of A and so MA is basically an identity matrix
			Therefore, p = Md; (This is the PCG (preconditionedConjugateGradient) method)

			Construct d
			Solve for p

		This pressure is used to get the gradient of p which is used to calculate the new velocity u^(n+1)
			Subtract pressure from our velocity and save in target
	*/

	// CHECK Section
	const double constant = (fluidDensity * (gridCellSize * gridCellSize))/dt; // Why square not cube

	GridData p = GridData();
	GridData d = GridData();

	// Construct d
	FOR_EACH_CELL 
	{
		d(i,j,k) = -calcDivergence(i,j,k);
	}

	// Use PCG method to calculate p
	preconditionedConjugateGradient(AMatrix, p, d, 200, 0.01);

	FOR_EACH_CELL 
	{
		p(i,j,k) *= constant;
		// Save pressure values
		target.mP(i,j,k) = p(i,j,k);
	}

	// Calculate the pressure gradient, use it to calculate the new velocity u^(n+1) and store 
	// this new velocity at the grid cell faces
	FOR_EACH_FACE 
	{
		// Initialize the pressure values to 0.
		vec3 pLow  = vec3(0.0, 0.0, 0.0);
		vec3 pHigh = vec3(0.0, 0.0, 0.0);

		const double deltaT_By_Density = dt / fluidDensity; // Bottom of page 27 in course notes
		setPressureHighLow( i, j, k, p, pLow, pHigh );

		// Update the velocities:
		if (isValidFace(MACGrid::X, i, j, k)) 
		{
			target.mU(i,j,k) = mU(i,j,k) - deltaT_By_Density * (pHigh[0] - pLow[0]) / gridCellSize;
		}
		if (isValidFace(MACGrid::Y, i, j, k)) 
		{
			target.mV(i,j,k) = mV(i,j,k) - deltaT_By_Density * (pHigh[1] - pLow[1]) / gridCellSize;
		}
		if (isValidFace(MACGrid::Z, i, j, k)) 
		{
			target.mW(i,j,k) = mW(i,j,k) - deltaT_By_Density * (pHigh[2] - pLow[2]) / gridCellSize;
		}
	}

	#ifdef _DEBUG
		checkBorderVelocities();
	#endif

	// Then save the result to our object
	mP = target.mP; 
	mU = target.mU;
	mV = target.mV;
	mW = target.mW;

	#ifdef _DEBUG
		// IMPLEMENT THIS AS A SANITY CHECK: assert (checkDivergence());
		std::ostringstream num_str1;
		std::ostringstream num_str2;

		FOR_EACH_CELL 
		{
			double divergence = calcDivergence(i,j,k);
			if (abs(divergence) > 0.02 ) 
			{
				PrintLine("WARNING: Divergent! ");
				num_str1 << "Divergence: " << divergence << "\n";
				PrintLine( num_str1.str() );
				num_str2 << "Cell: " << i << ", " << j << ", " << k;
				PrintLine( num_str2.str() );
			}
		}
	#endif
}

//--------------
// Simulation //
//--------------

void MACGrid::computeBouyancy()
{
	// CHECK Section 
	// Resource: At the end of section 2, in the Paper titled Visual Simulation of Smoke.
	// (http://physbam.stanford.edu/~fedkiw/papers/stanford2001-01.pdf)
	// Calculate bouyancy and store in target

	FOR_EACH_FACE 
	{
		if (isValidFace(MACGrid::Y, i, j, k)) 
		{
			vec3 position = getFacePosition(MACGrid::Y, i, j, k);
			double deltaTemperature = getTemperature(position) - ambientTemperature;
			double buoyancyForce = -buoyancyAlpha * getDensity(position) + buoyancyBeta * deltaTemperature;
			target.mV(i,j,k) = mV(i,j,k) + buoyancyForce;
		}
	}

	mV = target.mV;
}

void MACGrid::computeVorticityConfinement()
{
	// CHECK Section 
	// Calculate vorticity confinement forces and apply the forces to the current velocity

	// Save these values
	target.mU = mU;
	target.mV = mV;
	target.mW = mW;

	GridData tempX = GridData();
	GridData tempY = GridData();
	GridData tempZ = GridData();
	GridData tempM = GridData();

	const double One_By_twoDeltaCellSize = 1.0/(2.0 * gridCellSize);

	FOR_EACH_CELL 
	{
		const double a = (mW(i,j+1,k) - mW(i,j-1,k))  -  (mV(i,j,k+1) - mV(i,j,k-1));
		const double b = (mU(i,j,k+1) - mU(i,j,k-1))  -  (mW(i+1,j,k) - mW(i-1,j,k));
		const double c = (mV(i+1,j,k) - mV(i-1,j,k))  -  (mU(i,j+1,k) - mU(i,j-1,k));

		double vorticityX = a * One_By_twoDeltaCellSize;
		double vorticityY = b * One_By_twoDeltaCellSize;
		double vorticityZ = c * One_By_twoDeltaCellSize;
		vec3 vorticity(vorticityX, vorticityY, vorticityZ);

		// Store the vorticity (as separate components) and also store the magnitude of the vorticity:
		tempX(i,j,k) = vorticityX;
		tempY(i,j,k) = vorticityY;
		tempZ(i,j,k) = vorticityZ;
		tempM(i,j,k) = vorticity.Length();
	}

	FOR_EACH_CELL 
	{
		const double gradientX = (tempM(i+1,j,k) - tempM(i-1,j,k)) * One_By_twoDeltaCellSize;
		const double gradientY = (tempM(i,j+1,k) - tempM(i,j-1,k)) * One_By_twoDeltaCellSize;
		const double gradientZ = (tempM(i,j,k+1) - tempM(i,j,k-1)) * One_By_twoDeltaCellSize;
		vec3 gradient(gradientX, gradientY, gradientZ);
		gradient.Normalize();
		
		// Get the stored vorticity:
		vec3 vorticity(tempX(i,j,k), tempY(i,j,k), tempZ(i,j,k));

		// Calculate the confinement force:
		vec3 fConf = vorticityEpsilon * gridCellSize * (gradient.Cross(vorticity));

		// Spread fConf to the surrounding faces:
		if (isValidFace(0, i,j,k))   target.mU(i,j,k)   += fConf[0] * 0.5;
		if (isValidFace(0, i+1,j,k)) target.mU(i+1,j,k) += fConf[0] * 0.5;
		if (isValidFace(1, i,j,k))   target.mV(i,j,k)   += fConf[1] * 0.5;
		if (isValidFace(1, i,j+1,k)) target.mV(i,j+1,k) += fConf[1] * 0.5;
		if (isValidFace(2, i,j,k))   target.mW(i,j,k)   += fConf[2] * 0.5;
		if (isValidFace(2, i,j,k+1)) target.mW(i,j,k+1) += fConf[2] * 0.5;
	}

	mU = target.mU;
	mV = target.mV;
	mW = target.mW;
}

//-------------------------
// Miscellaneous Helpers //
//-------------------------

vec3 MACGrid::getRewoundPosition(const vec3 & currentPosition) 
{
	/*
	// EULER (RK1):
	vec3 currentVelocity = getVelocity(currentPosition);
	vec3 rewoundPosition = currentPosition - currentVelocity * dt;
	vec3 clippedRewoundPosition = clipToGrid(rewoundPosition, currentPosition);
	return clippedRewoundPosition;
	*/

	// HEUN / MODIFIED EULER (RK2):
	vec3 currentVelocity = getVelocity(currentPosition);
	vec3 rewoundPosition = currentPosition - currentVelocity * dt;
	vec3 clippedRewoundPosition = clipToGrid(rewoundPosition, currentPosition);
	// Keep going...
	vec3 rewoundVelocity = getVelocity(clippedRewoundPosition);
	vec3 averageVelocity = (currentVelocity + rewoundVelocity) / 2.0;
	vec3 betterRewoundPosition = currentPosition - averageVelocity * dt;
	vec3 clippedBetterRewoundPosition = clipToGrid(betterRewoundPosition, currentPosition);
	return clippedBetterRewoundPosition;
}

vec3 MACGrid::clipToGrid(const vec3& outsidePoint, const vec3& insidePoint) 
{
	vec3 clippedPoint = outsidePoint;

	for (int i = 0; i < 3; i++) 
	{
		if (clippedPoint[i] < 0) 
		{
			vec3 distance = clippedPoint - insidePoint;
			double newDistanceI = 0 - insidePoint[i];
			double ratio = newDistanceI / distance[i];
			clippedPoint = insidePoint + distance * ratio;
		}
		if (clippedPoint[i] > getSize(i)) 
		{
			vec3 distance = clippedPoint - insidePoint;
			double newDistanceI = getSize(i) - insidePoint[i];
			double ratio = newDistanceI / distance[i];
			clippedPoint = insidePoint + distance * ratio;
		}
	}

#ifdef _DEBUG
	// Make sure the point is now in the grid:
	if (clippedPoint[0] < 0 || clippedPoint[0] > getSize(0) || 
		clippedPoint[1] < 0 || clippedPoint[1] > getSize(1) || 
		clippedPoint[2] < 0 || clippedPoint[2] > getSize(2)) 
	{
		PrintLine("WARNING: Clipped point is outside grid!");
	}
#endif

	return clippedPoint;
}

double MACGrid::getSize(int dimension) 
{
	return gridDim[dimension] * gridCellSize;
}

int MACGrid::getCellIndex(int i, int j, int k)
{
	return i + j * gridDim[MACGrid::X] + k * gridDim[MACGrid::Y] * gridDim[MACGrid::X];
}

int MACGrid::getNumberOfCells()
{
	return gridDim[MACGrid::X] * gridDim[MACGrid::Y] * gridDim[MACGrid::Z];
}

vec3 MACGrid::getFacePosition(int dimension, int i, int j, int k)
{
	if (dimension == 0) 
	{
		return vec3(i * gridCellSize, (j + 0.5) * gridCellSize, (k + 0.5) * gridCellSize);
	} 
	else if (dimension == 1) 
	{
		return vec3((i + 0.5) * gridCellSize, j * gridCellSize, (k + 0.5) * gridCellSize);
	} 
	else if (dimension == 2) 
	{
		return vec3((i + 0.5) * gridCellSize, (j + 0.5) * gridCellSize, k * gridCellSize);
	}

	return vec3(0,0,0);
}

//------------------
// Perconditioner //
//------------------

void MACGrid::calculateAMatrix() 
{
	FOR_EACH_CELL 
	{
		int numFluidNeighbors = 0;
		if (i-1 >= 0) 
		{
			AMatrix.plusI(i-1,j,k) = -1;
			numFluidNeighbors++;
		}
		if (i+1 < gridDim[MACGrid::X]) 
		{
			AMatrix.plusI(i,j,k) = -1;
			numFluidNeighbors++;
		}
		if (j-1 >= 0) 
		{
			AMatrix.plusJ(i,j-1,k) = -1;
			numFluidNeighbors++;
		}
		if (j+1 < gridDim[MACGrid::Y]) 
		{
			AMatrix.plusJ(i,j,k) = -1;
			numFluidNeighbors++;
		}
		if (k-1 >= 0) 
		{
			AMatrix.plusK(i,j,k-1) = -1;
			numFluidNeighbors++;
		}
		if (k+1 < gridDim[MACGrid::Z]) 
		{
			AMatrix.plusK(i,j,k) = -1;
			numFluidNeighbors++;
		}
		// Set the diagonal:
		AMatrix.diag(i,j,k) = numFluidNeighbors;
	}
}

bool MACGrid::preconditionedConjugateGradient(const GridDataMatrix & A, GridData & p, const GridData & d, 
											  int maxIterations, double tolerance) 
{
	// Solves Ap = d for p.

	FOR_EACH_CELL 
	{
		p(i,j,k) = 0.0; // Initial guess p = 0.	
	}

	GridData r = d; // Residual vector.

	/*
	PrintLine("r: ");
	FOR_EACH_CELL 
	{
		PrintLine(r(i,j,k));
	}
	*/
	GridData z = GridData();
	applyPreconditioner(r, A, z); // Auxillary vector.
	/*
	PrintLine("z: ");
	FOR_EACH_CELL 
	{
		PrintLine(z(i,j,k));
	}
	*/

	GridData s = z; // Search vector;

	double sigma = dotProduct(z, r);

	for (int iteration = 0; iteration < maxIterations; iteration++) 
	{
		double rho = sigma;

		apply(A, s, z); // z = applyA(s);

		double alpha = rho/dotProduct(z, s);

		GridData alphaTimesS = GridData();
		multiply(alpha, s, alphaTimesS);
		add(p, alphaTimesS, p);

		GridData alphaTimesZ = GridData();
		multiply(alpha, z, alphaTimesZ);
		subtract(r, alphaTimesZ, r);

		if (maxMagnitude(r) <= tolerance) 
		{
			//PrintLine("PCG converged in " << (iteration + 1) << " iterations.");
			return true;
		}

		applyPreconditioner(r, A, z); // z = applyPreconditioner(r);

		double sigmaNew = dotProduct(z, r);

		double beta = sigmaNew / rho;

		GridData betaTimesS = GridData();
		multiply(beta, s, betaTimesS);
		add(z, betaTimesS, s);

		sigma = sigmaNew;
	}

	PrintLine( "BAD pressure value generated because the PCG algorithm didn't converge." );
	return false;
}

void MACGrid::calculatePreconditioner(const GridDataMatrix & A) 
{
	precon.initialize();

	// CHECK Section
	// Build the modified incomplete Cholesky preconditioner, essentially calculate precon(i,j,k) for all cells
	// Resource: Fig 4.2 on page 36 of Bridson's 2007 SIGGRAPH fluid course notes.

	const double tau = 0.97; // Tuning constant.
	FOR_EACH_CELL 
	{
		//if (A.diag(i,j,k) != 0.0)  // If cell is a fluid... put in check later for more complex scenes
		{
			const double a = A.plusI(i-1,j,k) * precon(i-1,j,k);
			const double b = A.plusJ(i,j-1,k) * precon(i,j-1,k);
			const double c = A.plusK(i,j,k-1) * precon(i,j,k-1);
			const double d = a*a + b*b + c*c;

			const double f = a * (A.plusJ(i-1,j,k) + A.plusK(i-1,j,k)) * precon(i-1,j,k);
			const double g = b * (A.plusI(i,j-1,k) + A.plusK(i,j-1,k)) * precon(i,j-1,k);
			const double h = c * (A.plusI(i,j,k-1) + A.plusJ(i,j,k-1)) * precon(i,j,k-1);
			const double m = d + tau * (f + g + h);

			double e = A.diag(i,j,k) - m;
			precon(i,j,k) = 1.0 / sqrt(e + DOUBLE_EPSILON);
		}
	}
}

void MACGrid::applyPreconditioner(const GridData & r, const GridDataMatrix & A, GridData & z) 
{
	// Resource: Fig 4.3 on page 37 of Bridson's 2007 SIGGRAPH fluid course notes.
	// Solve Lq = r for q:
	GridData q = GridData();

	FOR_EACH_CELL 
	{
		//if (A.diag(i,j,k) != 0.0) // If cell is a fluid.
		{
			const double a = A.plusI(i - 1, j, k) * precon(i - 1, j, k) * q(i - 1, j, k);
			const double b = A.plusJ(i, j - 1, k) * precon(i, j - 1, k) * q(i, j - 1, k);
			const double c = A.plusK(i, j, k - 1) * precon(i, j, k - 1) * q(i, j, k - 1);

			double t = r(i, j, k) - (a+b+c);
			q(i, j, k) = t * precon(i, j, k);
		}
	}

	// Solve L^Tz = q for z:
	FOR_EACH_CELL_REVERSE 
	{
		//if (A.diag(i,j,k) != 0.0) // If cell is a fluid.
		{
			const double a = A.plusI(i, j, k) * precon(i, j, k) * z(i + 1, j, k);
			const double b = A.plusJ(i, j, k) * precon(i, j, k) * z(i, j + 1, k);
			const double c = A.plusK(i, j, k) * precon(i, j, k) * z(i, j, k + 1);

			double t = q(i, j, k) - (a+b+c);
			z(i, j, k) = t * precon(i, j, k);
		}
	}
}

//----------------------------
// Per Cell Math Operations //
//----------------------------

void MACGrid::apply(const GridDataMatrix & matrix, const GridData & vector, GridData & result) 
{
	FOR_EACH_CELL 
	{ 
		// For each row of the matrix.
		double diag   = 0.0;
		double plusI  = 0.0;
		double plusJ  = 0.0;
		double plusK  = 0.0;
		double minusI = 0.0;
		double minusJ = 0.0;
		double minusK = 0.0;

		diag = matrix.diag(i,j,k) * vector(i,j,k);
		if (isValidCell(i+1,j,k)) plusI = matrix.plusI(i,j,k) * vector(i+1,j,k);
		if (isValidCell(i,j+1,k)) plusJ = matrix.plusJ(i,j,k) * vector(i,j+1,k);
		if (isValidCell(i,j,k+1)) plusK = matrix.plusK(i,j,k) * vector(i,j,k+1);
		if (isValidCell(i-1,j,k)) minusI = matrix.plusI(i-1,j,k) * vector(i-1,j,k);
		if (isValidCell(i,j-1,k)) minusJ = matrix.plusJ(i,j-1,k) * vector(i,j-1,k);
		if (isValidCell(i,j,k-1)) minusK = matrix.plusK(i,j,k-1) * vector(i,j,k-1);

		result(i,j,k) = diag + plusI + plusJ + plusK + minusI + minusJ + minusK;
	}
}

double MACGrid::dotProduct(const GridData & vector1, const GridData & vector2) 
{
	double result = 0.0;

	FOR_EACH_CELL 
	{
		result += vector1(i,j,k) * vector2(i,j,k);
	}

	return result;
}

void MACGrid::add(const GridData & vector1, const GridData & vector2, GridData & result) 
{
	FOR_EACH_CELL 
	{
		result(i,j,k) = vector1(i,j,k) + vector2(i,j,k);
	}
}

void MACGrid::subtract(const GridData & vector1, const GridData & vector2, GridData & result) 
{
	FOR_EACH_CELL 
	{
		result(i,j,k) = vector1(i,j,k) - vector2(i,j,k);
	}
}

void MACGrid::multiply(const double scalar, const GridData & vector, GridData & result) 
{
	FOR_EACH_CELL 
	{
		result(i,j,k) = scalar * vector(i,j,k);
	}
}

double MACGrid::maxMagnitude(const GridData & vector) 
{
	double result = 0.0;

	FOR_EACH_CELL 
	{
		if (abs(vector(i,j,k)) > result) result = abs(vector(i,j,k));
	}

	return result;
}

//------------------
// Save Functions //
//------------------

void MACGrid::saveSmoke(const char* fileName) 
{
	std::ofstream fileOut(fileName);
	if (fileOut.is_open()) 
	{
		FOR_EACH_CELL 
		{
			fileOut << mD(i,j,k) << std::endl;
		}
		fileOut.close();
	}
}

void MACGrid::saveParticle(std::string filename)
{
	Partio::ParticlesDataMutable *parts = Partio::create();
	Partio::ParticleAttribute posH, vH;
	posH = parts->addAttribute("position", Partio::VECTOR, 3);
	vH = parts->addAttribute("v", Partio::VECTOR, 3);
	for (unsigned int i = 0; i < rendering_particles.size(); i++)
	{
		int idx = parts->addParticle();
		float *p = parts->dataWrite<float>(posH, idx);
		float *v = parts->dataWrite<float>(vH, idx);
		for (int k = 0; k < 3; k++)
		{
			p[k] = rendering_particles[i][k];
			v[k] = rendering_particles_vel[i][k];
		}
	}
	
	Partio::write(filename.c_str(), *parts);
	parts->release();
}

void MACGrid::saveDensity(std::string filename)
{
	Partio::ParticlesDataMutable *density_field = Partio::create();
	Partio::ParticleAttribute posH, rhoH;
	posH = density_field->addAttribute("position", Partio::VECTOR, 3);
	rhoH = density_field->addAttribute("density", Partio::VECTOR, 1);
	
	FOR_EACH_CELL
	{
		int idx = density_field->addParticle();
		float *p = density_field->dataWrite<float>(posH, idx);
		float *rho = density_field->dataWrite<float>(rhoH, idx);
		vec3 cellCenter = getCenter(i, j, k);
		for (int l = 0; l < 3; l++)
		{
			p[l] = cellCenter[l];
		}
		rho[0] = getDensity(cellCenter);
	}
	Partio::write(filename.c_str(), *density_field);
	density_field->release();
}

//---------------
// RenderColor //
//---------------

vec4 MACGrid::getRenderColor(int i, int j, int k)
{
	double value = mD(i, j, k); 
	vec4 coldColor(0.5, 0.5, 1.0, value);
	vec4 hotColor(1.0, 0.5, 0.5, value);
	return Globals::Lerp(coldColor, hotColor, mT(i, j, k));
}

vec4 MACGrid::getRenderColor(const vec3& pt)
{
	double value = getDensity(pt);
	vec4 coldColor(0.5, 0.5, 1.0, value);
	vec4 hotColor(1.0, 0.5, 0.5, value);
	return Globals::Lerp(coldColor, hotColor, getTemperature(pt));
}

//---------------------
// Drawing Functions //
//---------------------

void MACGrid::draw(const Camera& c)
{
	drawWireGrid();
	if (theDisplayVel) drawVelocities();   
	if (theRenderMode == CUBES) drawSmokeCubes(c);
	else drawSmoke(c);
}

void MACGrid::drawVelocities()
{
	// Draw line at each center
	//glColor4f(0.0, 1.0, 0.0, 1.0);
	glBegin(GL_LINES);
	FOR_EACH_CELL
	{
		vec3 pos = getCenter(i,j,k);
		vec3 vel = getVelocity(pos);
		if (vel.Length() > 0.0001)
		{
			//vel.Normalize(); 
			vel *= gridCellSize/2.0;
			vel += pos;
			glColor4f(1.0, 1.0, 0.0, 1.0);
			glVertex3dv(pos.n);
			glColor4f(0.0, 1.0, 0.0, 1.0);
			glVertex3dv(vel.n);
		}
	}
	glEnd();
}

void MACGrid::drawSmoke(const Camera& c)
{
	vec3 eyeDir = c.getBackward();
	double zresult = fabs(Dot(eyeDir, vec3(1,0,0)));
	double xresult = fabs(Dot(eyeDir, vec3(0,0,1)));
	//double yresult = fabs(Dot(eyeDir, vec3(0,1,0)));

	if (zresult < xresult)
	{
		drawZSheets(c.getPosition()[2] < 0);
	}
	else 
	{
		drawXSheets(c.getPosition()[0] < 0);
	}
}

void MACGrid::drawSmokeCubes(const Camera& c)
{
	std::multimap<double, MACGrid::Cube, std::greater<double> > cubes;
	FOR_EACH_CELL
	{
		MACGrid::Cube cube;
		cube.color = getRenderColor(i,j,k);
		cube.pos = getCenter(i,j,k);
		cube.dist = DistanceSqr(cube.pos, c.getPosition());
		cubes.insert(make_pair(cube.dist, cube));
	} 

	// Draw cubes from back to front
	std::multimap<double, MACGrid::Cube, std::greater<double> >::const_iterator it;
	for (it = cubes.begin(); it != cubes.end(); ++it)
	{
		drawCube(it->second);
	}
}

//-------------------------
// Draw Helper Functions //
//-------------------------

void MACGrid::drawZSheets(bool backToFront)
{
	// Draw K Sheets from back to front
	double back =  (gridDim[2])*gridCellSize;
	double top  =  (gridDim[1])*gridCellSize;
	double right = (gridDim[0])*gridCellSize;

	double stepsize = gridCellSize*0.25;

	double startk = back - stepsize;
	double endk = 0;
	double stepk = -gridCellSize;

	if (!backToFront)
	{
		startk = 0;
		endk = back;   
		stepk = gridCellSize;
	}

	for (double k = startk; backToFront? k > endk : k < endk; k += stepk)
	{
		for (double j = 0.0; j < top; )
		{
			glBegin(GL_QUAD_STRIP);
			for (double i = 0.0; i <= right; i += stepsize)
			{
				vec3 pos1 = vec3(i,j,k); 
				vec3 pos2 = vec3(i, j+stepsize, k); 

				vec4 color1 = getRenderColor(pos1);
				vec4 color2 = getRenderColor(pos2);

				glColor4dv(color1.n);
				glVertex3dv(pos1.n);

				glColor4dv(color2.n);
				glVertex3dv(pos2.n);
			} 
			glEnd();
			j+=stepsize;

			glBegin(GL_QUAD_STRIP);
			for (double i = right; i >= 0.0; i -= stepsize)
			{
				vec3 pos1 = vec3(i,j,k); 
				vec3 pos2 = vec3(i, j+stepsize, k); 

				vec4 color1 = getRenderColor(pos1);
				vec4 color2 = getRenderColor(pos2);

				glColor4dv(color1.n);
				glVertex3dv(pos1.n);

				glColor4dv(color2.n);
				glVertex3dv(pos2.n);
			} 
			glEnd();
			j+=stepsize;
		}
	}
}

void MACGrid::drawXSheets(bool backToFront)
{
	// Draw K Sheets from back to front
	double back =  (gridDim[2])*gridCellSize;
	double top  =  (gridDim[1])*gridCellSize;
	double right = (gridDim[0])*gridCellSize;

	double stepsize = gridCellSize*0.25;

	double starti = right - stepsize;
	double endi = 0;
	double stepi = -gridCellSize;

	if (!backToFront)
	{
		starti = 0;
		endi = right;   
		stepi = gridCellSize;
	}

	for (double i = starti; backToFront? i > endi : i < endi; i += stepi)
	{
		for (double j = 0.0; j < top; )
		{
			glBegin(GL_QUAD_STRIP);
			for (double k = 0.0; k <= back; k += stepsize)
			{
				vec3 pos1 = vec3(i,j,k); 
				vec3 pos2 = vec3(i, j+stepsize, k); 

				vec4 color1 = getRenderColor(pos1);
				vec4 color2 = getRenderColor(pos2);

				glColor4dv(color1.n);
				glVertex3dv(pos1.n);

				glColor4dv(color2.n);
				glVertex3dv(pos2.n);
			} 
			glEnd();
			j+=stepsize;

			glBegin(GL_QUAD_STRIP);
			for (double k = back; k >= 0.0; k -= stepsize)
			{
				vec3 pos1 = vec3(i,j,k); 
				vec3 pos2 = vec3(i, j+stepsize, k); 

				vec4 color1 = getRenderColor(pos1);
				vec4 color2 = getRenderColor(pos2);

				glColor4dv(color1.n);
				glVertex3dv(pos1.n);

				glColor4dv(color2.n);
				glVertex3dv(pos2.n);
			} 
			glEnd();
			j+=stepsize;
		}
	}
}

void MACGrid::drawWireGrid()
{
	// Display grid in light grey, draw top & bottom
	double xstart = 0.0;
	double ystart = 0.0;
	double zstart = 0.0;
	double xend = gridDim[0]*gridCellSize;
	double yend = gridDim[1]*gridCellSize;
	double zend = gridDim[2]*gridCellSize;

	glPushAttrib(GL_LIGHTING_BIT | GL_LINE_BIT);
	glDisable(GL_LIGHTING);
	glColor3f(0.25, 0.25, 0.25);

	glBegin(GL_LINES);
	for (int i = 0; i <= gridDim[0]; i++)
	{
		double x = xstart + i*gridCellSize;
		glVertex3d(x, ystart, zstart);
		glVertex3d(x, ystart, zend);

		glVertex3d(x, yend, zstart);
		glVertex3d(x, yend, zend);
	}

	for (int i = 0; i <= gridDim[2]; i++)
	{
		double z = zstart + i*gridCellSize;
		glVertex3d(xstart, ystart, z);
		glVertex3d(xend, ystart, z);

		glVertex3d(xstart, yend, z);
		glVertex3d(xend, yend, z);
	}

	glVertex3d(xstart, ystart, zstart);
	glVertex3d(xstart, yend, zstart);

	glVertex3d(xend, ystart, zstart);
	glVertex3d(xend, yend, zstart);

	glVertex3d(xstart, ystart, zend);
	glVertex3d(xstart, yend, zend);

	glVertex3d(xend, ystart, zend);
	glVertex3d(xend, yend, zend);
	glEnd();
	glPopAttrib();

	glEnd();
}

void MACGrid::drawFace(const MACGrid::Cube& cube)
{
	glColor4dv(cube.color.n);
	glPushMatrix();
		glTranslated(cube.pos[0], cube.pos[1], cube.pos[2]);      
		glScaled(gridCellSize, gridCellSize, gridCellSize);
		glBegin(GL_QUADS);
			glNormal3d( 0.0,  0.0, 1.0);
			glVertex3d(-LEN, -LEN, LEN);
			glVertex3d(-LEN,  LEN, LEN);
			glVertex3d( LEN,  LEN, LEN);
			glVertex3d( LEN, -LEN, LEN);
		glEnd();
	glPopMatrix();
}

void MACGrid::drawCube(const MACGrid::Cube& cube)
{
	glColor4dv(cube.color.n);
	glPushMatrix();
		glTranslated(cube.pos[0], cube.pos[1], cube.pos[2]);      
		glScaled(gridCellSize, gridCellSize, gridCellSize);
		glBegin(GL_QUADS);
			glNormal3d( 0.0, -1.0,  0.0);
			glVertex3d(-LEN, -LEN, -LEN);
			glVertex3d(-LEN, -LEN,  LEN);
			glVertex3d( LEN, -LEN,  LEN);
			glVertex3d( LEN, -LEN, -LEN);

			glNormal3d( 0.0,  0.0, -0.0);
			glVertex3d(-LEN, -LEN, -LEN);
			glVertex3d(-LEN,  LEN, -LEN);
			glVertex3d( LEN,  LEN, -LEN);
			glVertex3d( LEN, -LEN, -LEN);

			glNormal3d(-1.0,  0.0,  0.0);
			glVertex3d(-LEN, -LEN, -LEN);
			glVertex3d(-LEN, -LEN,  LEN);
			glVertex3d(-LEN,  LEN,  LEN);
			glVertex3d(-LEN,  LEN, -LEN);

			glNormal3d( 0.0, 1.0,  0.0);
			glVertex3d(-LEN, LEN, -LEN);
			glVertex3d(-LEN, LEN,  LEN);
			glVertex3d( LEN, LEN,  LEN);
			glVertex3d( LEN, LEN, -LEN);

			glNormal3d( 0.0,  0.0, 1.0);
			glVertex3d(-LEN, -LEN, LEN);
			glVertex3d(-LEN,  LEN, LEN);
			glVertex3d( LEN,  LEN, LEN);
			glVertex3d( LEN, -LEN, LEN);

			glNormal3d(1.0,  0.0,  0.0);
			glVertex3d(LEN, -LEN, -LEN);
			glVertex3d(LEN, -LEN,  LEN);
			glVertex3d(LEN,  LEN,  LEN);
			glVertex3d(LEN,  LEN, -LEN);
		glEnd();
	glPopMatrix();
}

//---------------------------
// Getter/Setter Functions //
//---------------------------

vec3 MACGrid::getVelocity(const vec3& pt)
{
	vec3 vel( getVelocityX(pt), getVelocityY(pt), getVelocityZ(pt) );
	return vel;
}

double MACGrid::getVelocityX(const vec3& pt)
{
	return mU.interpolate(pt);
}

double MACGrid::getVelocityY(const vec3& pt)
{
	return mV.interpolate(pt);
}

double MACGrid::getVelocityZ(const vec3& pt)
{
	return mW.interpolate(pt);
}

double MACGrid::getTemperature(const vec3& pt)
{
	return mT.interpolate(pt);
}

double MACGrid::getDensity(const vec3& pt)
{
	return mD.interpolate(pt);
}

vec3 MACGrid::getCenter(int i, int j, int k)
{
	double xstart = gridCellSize/2.0;
	double ystart = gridCellSize/2.0;
	double zstart = gridCellSize/2.0;

	double x = xstart + i*gridCellSize;
	double y = ystart + j*gridCellSize;
	double z = zstart + k*gridCellSize;
	return vec3(x, y, z);
}

//--------------
// Validation //
//--------------

bool MACGrid::isValidCell(int i, int j, int k)
{
	if (i >= gridDim[MACGrid::X] || j >= gridDim[MACGrid::Y] || k >= gridDim[MACGrid::Z]) {
		return false;
	}

	if (i < 0 || j < 0 || k < 0) {
		return false;
	}

	return true;
}

bool MACGrid::isValidFace(int dimension, int i, int j, int k)
{
	if (dimension == 0) 
	{
		if (i > gridDim[MACGrid::X] || j >= gridDim[MACGrid::Y] || k >= gridDim[MACGrid::Z]) 
		{
			return false;
		}
	} 
	else if (dimension == 1) 
	{
		if (i >= gridDim[MACGrid::X] || j > gridDim[MACGrid::Y] || k >= gridDim[MACGrid::Z]) 
		{
			return false;
		}
	} 
	else if (dimension == 2) 
	{
		if (i >= gridDim[MACGrid::X] || j >= gridDim[MACGrid::Y] || k > gridDim[MACGrid::Z]) 
		{
			return false;
		}
	}

	if (i < 0 || j < 0 || k < 0) 
	{
		return false;
	}

	return true;
}

void MACGrid::checkBorderVelocities()
{
	std::ostringstream num_str;

	FOR_EACH_FACE 
	{
		if (isValidFace(MACGrid::X, i, j, k)) 
		{
			if (i == 0) 
			{
				if (abs(target.mU(i,j,k)) > 0.0000001) 
				{
					num_str << "LOW X:  " << target.mU(i,j,k);
					PrintLine( num_str.str() );
					//target.mU(i,j,k) = 0;
				}
			}

			if (i == gridDim[MACGrid::X]) 
			{
				if (abs(target.mU(i,j,k)) > 0.0000001) 
				{
					num_str << "HIGH X: " << target.mU(i,j,k);
					PrintLine( num_str.str() );
					//target.mU(i,j,k) = 0;
				}
			}
		}
		if (isValidFace(MACGrid::Y, i, j, k)) 
		{
			if (j == 0) 
			{
				if (abs(target.mV(i,j,k)) > 0.0000001) 
				{
					num_str << "LOW Y:  " << target.mV(i,j,k);
					PrintLine( num_str.str() );
					//target.mV(i,j,k) = 0;
				}
			}

			if (j == gridDim[MACGrid::Y]) 
			{
				if (abs(target.mV(i,j,k)) > 0.0000001) 
				{
					num_str << "HIGH Y: " << target.mV(i,j,k);
					PrintLine( num_str.str() );
					//target.mV(i,j,k) = 0;
				}
			}
		}
		if (isValidFace(MACGrid::Z, i, j, k)) 
		{			
			if (k == 0) 
			{
				if (abs(target.mW(i,j,k)) > 0.0000001) 
				{
					num_str << "LOW Z:  " << target.mV(i,j,k);
					PrintLine( num_str.str() );
					//target.mW(i,j,k) = 0;
				}
			}

			if (k == gridDim[MACGrid::Z]) 
			{
				if (abs(target.mW(i,j,k)) > 0.0000001) 
				{
					num_str << "HIGH Z: " << target.mV(i,j,k);
					PrintLine( num_str.str() );
					//target.mW(i,j,k) = 0;
				}
			}
		}
	}
}