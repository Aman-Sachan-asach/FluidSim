#pragma once

#include <Partio.h>
#include <math.h>
#include <map>
#include <stdio.h>
#include <cstdlib>
#include <fstream> 

#include "../external/vec.h"
#include "globals.h"

#include "grid_data.h"
#include "camera.h"

#undef max
#undef min 
#define LEN 0.5

class Camera;

class MACGrid
{
public:
	enum Direction { X, Y, Z };

	MACGrid();	
	MACGrid(const MACGrid& orig);
	~MACGrid() {};
	
	MACGrid& operator=(const MACGrid& orig);

	void reset();

	void draw(const Camera& c);
	void updateSources();

	// Advection	
	void advectVelocity(double dt);
	void advectTemperature(double dt);
	void advectDensity(double dt);
	void advectRenderingParticles(double dt);

	// External Forces and Projection
	void addExternalForces(double dt);
	void project(double dt);

	// Save Functions
	void saveSmoke(const char* fileName);
	void saveParticle(std::string filename);
	void saveDensity(std::string filename);

protected:
	// Simulation
	void computeBouyancy(double dt);
	void computeVorticityConfinement(double dt);

	// Helper Functions
	vec3 getRewoundPosition(const vec3 & currentPosition, const double dt);
	vec3 clipToGrid(const vec3& outsidePoint, const vec3& insidePoint);

	double getSize(int dimension);
	int getCellIndex(int i, int j, int k);
	int getNumberOfCells();
	vec3 getFacePosition(int dimension, int i, int j, int k);

	//Preconditioner
	void calculateAMatrix();
	bool preconditionedConjugateGradient(const GridDataMatrix & A, GridData & p, const GridData & d, int maxIterations, double tolerance);
	void calculatePreconditioner(const GridDataMatrix & A);
	void applyPreconditioner(const GridData & r, const GridDataMatrix & A, GridData & z);

	// Per Cell Math Operations
	void apply(const GridDataMatrix & matrix, const GridData & vector, GridData & result);	
	double dotProduct(const GridData & vector1, const GridData & vector2);
	void add(const GridData & vector1, const GridData & vector2, GridData & result);
	void subtract(const GridData & vector1, const GridData & vector2, GridData & result);
	void multiply(const double scalar, const GridData & vector, GridData & result);
	double maxMagnitude(const GridData & vector);

	// Rendering
	struct Cube { vec3 pos; vec4 color; double dist; };
	
	vec4 getRenderColor(int i, int j, int k);
	vec4 getRenderColor(const vec3& pt);

	void drawVelocities();
	void drawSmokeCubes(const Camera& c);
	void drawSmoke(const Camera& c);

	void drawWireGrid();
	void drawCube(const MACGrid::Cube& c);
	void drawFace(const MACGrid::Cube& c);
	void drawZSheets(bool backToFront);
	void drawXSheets(bool backToFront);

	// Getter/Setter Functions
	vec3 getVelocity(const vec3& pt);
	double getVelocityX(const vec3& pt);
	double getVelocityY(const vec3& pt);
	double getVelocityZ(const vec3& pt);
	double getTemperature(const vec3& pt);
	double getDensity(const vec3& pt);
	vec3 getCenter(int i, int j, int k);

	// Validation
	bool isValidCell(int i, int j, int k);
	bool isValidFace(int dimension, int i, int j, int k);

protected:
	GridDataX mU; // X component of velocity, stored on X faces, size is (dimX+1)*dimY*dimZ
	GridDataY mV; // Y component of velocity, stored on Y faces, size is dimX*(dimY+1)*dimZ
	GridDataZ mW; // W component of velocity, stored on Z faces, size is dimX*dimY*(dimZ+1)
	GridData mP;  // Pressure, stored at grid centers, size is dimX*dimY*dimZ
	GridData mD;  // Density, stored at grid centers, size is dimX*dimY*dimZ
	GridData mT;  // Temperature, stored at grid centers, size is dimX*dimY*dimZ

	GridDataMatrix AMatrix;
	GridData precon;

public:
	// rendering particles
	std::vector<vec3> rendering_particles;
	std::vector<vec3> rendering_particles_vel;

	enum RenderMode { CUBES, SHEETS };
	static RenderMode theRenderMode;
	static bool theDisplayVel;
};