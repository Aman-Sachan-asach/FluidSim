#include "../headers/grid_data.h"

GridData::GridData() : mDfltValue(0.0), mMax(0.0,0.0,0.0) 
{
	initialize(mDfltValue);
}

GridData& GridData::operator=(const GridData& orig)
{
	if (this == &orig)
	{
		return *this;
	}

	mDfltValue = orig.mDfltValue;
	mData = orig.mData;
	mMax = orig.mMax;
	return *this;
}

double& GridData::operator()(int i, int j, int k)
{
	static double dflt = mDfltValue;  // HACK: Protect against setting the default value

	if( i<0 || i > gridDim[0]-1 || 
		j<0 || j > gridDim[1]-1 || 
		k<0 || k > gridDim[2]-1 ) 
	{
		return dflt;
	}

	int col   = i;
	int row   = k*gridDim[0];
	int stack = j*gridDim[0]*gridDim[2];

	return mData[col+row+stack];
}

const double GridData::operator()(int i, int j, int k) const
{
	static double dflt = mDfltValue;  // HACK: Protect against setting the default value

	if( i<0 || i > gridDim[0]-1 || 
		j<0 || j > gridDim[1]-1 || 
		k<0 || k > gridDim[2]-1 ) 
	{
		return dflt;
	}

	int col   = i;
	int row   = k*gridDim[0];
	int stack = j*gridDim[0]*gridDim[2];

	return mData[col+row+stack];
}

void GridData::initialize(double dfltValue)
{
	mDfltValue = dfltValue;
	mMax[0] = gridCellSize*gridDim[0];
	mMax[1] = gridCellSize*gridDim[1];
	mMax[2] = gridCellSize*gridDim[2];
	mData.resize(gridDim[0]*gridDim[1]*gridDim[2], false);
	std::fill(mData.begin(), mData.end(), mDfltValue);
}

std::vector<double>& GridData::data()
{
	return mData;
}

void GridData::getCell(const vec3& pt, int& i, int& j, int& k)
{
	vec3 pos = worldToSelf(pt); 
	i = (int) (pos[0]/gridCellSize);
	j = (int) (pos[1]/gridCellSize);
	k = (int) (pos[2]/gridCellSize);   
}

double GridData::interpolate(const vec3& pt)
{
	// SHARPER CUBIC INTERPOLATION:
	vec3 pos = worldToSelf(pt);

	int i = (int) (pos[0]/gridCellSize);
	int j = (int) (pos[1]/gridCellSize);
	int k = (int) (pos[2]/gridCellSize);

	double scale = 1.0/gridCellSize;  
	double fractx = scale*(pos[0] - i*gridCellSize);
	double fracty = scale*(pos[1] - j*gridCellSize);
	double fractz = scale*(pos[2] - k*gridCellSize);

#ifdef _DEBUG
	assert (fractx < 1.0 && fractx >= 0);
	assert (fracty < 1.0 && fracty >= 0);
	assert (fractz < 1.0 && fractz >= 0);
#endif

  

	double t[4][4];
	double u[4];
	double f;
#define ONE 1
#define EVAL(a,b,c) (*this)(a,b,c)
	for (int x = -1; x <= 2; x++) 
	{
		for (int y = -1; y <= 2; y++) 
		{
			t[x+ONE][y+ONE] = CINT( EVAL(i+x,j+y,k-1), EVAL(i+x,j+y,k+0), EVAL(i+x,j+y,k+1), EVAL(i+x,j+y,k+2), fractz );
		}
	}
#undef EVAL
	for (int x = -1; x <= 2; x++) 
	{
		u[x+ONE] = CINT( t[x+ONE][-1+ONE], t[x+ONE][0+ONE], t[x+ONE][1+ONE], t[x+ONE][2+ONE], fracty );
	}
	f = CINT( u[-1+ONE], u[0+ONE], u[1+ONE], u[2+ONE], fractx );
#undef ONE
	return f;
}

double GridData::CINT(double q_i_minus_1, double q_i, double q_i_plus_1, double q_i_plus_2, double x) const 
{
	// The slopes:
	double d_i = (q_i_plus_1 - q_i_minus_1) / 2.0;
	double d_i_plus_1 = (q_i_plus_2 - q_i) / 2.0;

	// Delta q:
	double delta_q = q_i_plus_1 - q_i;

	// Restrict the slopes:
	if (delta_q > 0) 
	{
		if (d_i < 0)        d_i = 0;
		if (d_i_plus_1 < 0) d_i_plus_1 = 0;
	} 
	else if (delta_q < 0) 
	{
		if (d_i > 0)        d_i = 0;
		if (d_i_plus_1 > 0) d_i_plus_1 = 0;
	}

	// The Hermite cubic:
	double q_x = q_i + d_i * x + (3.0 * delta_q - 2.0 * d_i - d_i_plus_1) * (x * x) + (-2.0 * delta_q + d_i + d_i_plus_1) * (x * x * x);

	return q_x;
}

vec3 GridData::worldToSelf(const vec3& pt) const
{
	vec3 out;
	out[0] = min(max(0.0, pt[0] - gridCellSize*0.5), mMax[0]);
	out[1] = min(max(0.0, pt[1] - gridCellSize*0.5), mMax[1]);
	out[2] = min(max(0.0, pt[2] - gridCellSize*0.5), mMax[2]);
	return out;
}

//-------------
// GridDataX //
//-------------

double& GridDataX::operator()(int i, int j, int k)
{
	static double dflt = 0;
	dflt = mDfltValue;  // Protect against setting the default value

	if (i < 0 || i > gridDim[0]) return dflt;

	if (j < 0) j = 0;
	if (j > gridDim[1]-1) j = gridDim[1]-1;
	if (k < 0) k = 0;
	if (k > gridDim[2]-1) k = gridDim[2]-1;

	int col = i;
	int row = k*(gridDim[0]+1);
	int stack = j*(gridDim[0]+1)*gridDim[2];
	return mData[stack + row + col];
}

const double GridDataX::operator()(int i, int j, int k) const
{
	static double dflt = 0;
	dflt = mDfltValue;  // Protect against setting the default value

	if (i < 0 || i > gridDim[0]) return dflt;

	if (j < 0) j = 0;
	if (j > gridDim[1]-1) j = gridDim[1]-1;
	if (k < 0) k = 0;
	if (k > gridDim[2]-1) k = gridDim[2]-1;

	int col = i;
	int row = k*(gridDim[0]+1);
	int stack = j*(gridDim[0]+1)*gridDim[2];
	return mData[stack + row + col];
}

void GridDataX::initialize(double dfltValue)
{
	GridData::initialize(dfltValue);
	mMax[0] = gridCellSize*(gridDim[0]+1);
	mMax[1] = gridCellSize*gridDim[1];
	mMax[2] = gridCellSize*gridDim[2];
	mData.resize((gridDim[0]+1)*gridDim[1]*gridDim[2], false);
	std::fill(mData.begin(), mData.end(), mDfltValue);
}

vec3 GridDataX::worldToSelf(const vec3& pt) const
{   
	vec3 out;
	out[0] = min(max(0.0, pt[0]), mMax[0]);
	out[1] = min(max(0.0, pt[1]-gridCellSize*0.5), mMax[1]);
	out[2] = min(max(0.0, pt[2]-gridCellSize*0.5), mMax[2]);
	return out;
}

//-------------
// GridDataY //
//-------------

double& GridDataY::operator()(int i, int j, int k)
{
	static double dflt = 0;
	dflt = mDfltValue;  // Protect against setting the default value

	if (j < 0 || j > gridDim[1]) return dflt;

	if (i < 0) i = 0;
	if (i > gridDim[0]-1) i = gridDim[0]-1;
	if (k < 0) k = 0;
	if (k > gridDim[2]-1) k = gridDim[2]-1;

	int col = i;
	int row = k*gridDim[0];
	int stack = j*gridDim[0]*gridDim[2];
	return mData[stack + row + col];
}

const double GridDataY::operator()(int i, int j, int k) const
{
	static double dflt = 0;
	dflt = mDfltValue;  // Protect against setting the default value

	if (j < 0 || j > gridDim[1]) return dflt;

	if (i < 0) i = 0;
	if (i > gridDim[0]-1) i = gridDim[0]-1;
	if (k < 0) k = 0;
	if (k > gridDim[2]-1) k = gridDim[2]-1;

	int col = i;
	int row = k*gridDim[0];
	int stack = j*gridDim[0]*gridDim[2];
	return mData[stack + row + col];
}

void GridDataY::initialize(double dfltValue)
{
	GridData::initialize(dfltValue);
	mMax[0] = gridCellSize*gridDim[0];
	mMax[1] = gridCellSize*(gridDim[1]+1);
	mMax[2] = gridCellSize*gridDim[2];
	mData.resize(gridDim[0]*(gridDim[1]+1)*gridDim[2], false);
	std::fill(mData.begin(), mData.end(), mDfltValue);
}

vec3 GridDataY::worldToSelf(const vec3& pt) const
{
	vec3 out;
	out[0] = min(max(0.0, pt[0]-gridCellSize*0.5), mMax[0]);
	out[1] = min(max(0.0, pt[1]), mMax[1]);
	out[2] = min(max(0.0, pt[2]-gridCellSize*0.5), mMax[2]);
	return out;
}

//-------------
// GridDataZ //
//-------------

double& GridDataZ::operator()(int i, int j, int k)
{
	static double dflt = 0;
	dflt = mDfltValue;  // Protect against setting the default value

	if (k < 0 || k > gridDim[2]) return dflt;

	if (i < 0) i = 0;
	if (i > gridDim[0]-1) i = gridDim[0]-1;
	if (j < 0) j = 0;
	if (j > gridDim[1]-1) j = gridDim[1]-1;

	int col = i;
	int row = k*gridDim[0];
	int stack = j*gridDim[0]*(gridDim[2]+1);

	return mData[stack + row + col];
}

const double GridDataZ::operator()(int i, int j, int k) const
{
	static double dflt = 0;
	dflt = mDfltValue;  // Protect against setting the default value

	if (k < 0 || k > gridDim[2]) return dflt;

	if (i < 0) i = 0;
	if (i > gridDim[0]-1) i = gridDim[0]-1;
	if (j < 0) j = 0;
	if (j > gridDim[1]-1) j = gridDim[1]-1;

	int col = i;
	int row = k*gridDim[0];
	int stack = j*gridDim[0]*(gridDim[2]+1);

	return mData[stack + row + col];
}

void GridDataZ::initialize(double dfltValue)
{
	GridData::initialize(dfltValue);
	mMax[0] = gridCellSize*gridDim[0];
	mMax[1] = gridCellSize*gridDim[1];
	mMax[2] = gridCellSize*(gridDim[2]+1);
	mData.resize(gridDim[0]*gridDim[1]*(gridDim[2]+1), false);
	std::fill(mData.begin(), mData.end(), mDfltValue);
}

vec3 GridDataZ::worldToSelf(const vec3& pt) const
{
	vec3 out;
	out[0] = min(max(0.0, pt[0]-gridCellSize*0.5), mMax[0]);
	out[1] = min(max(0.0, pt[1]-gridCellSize*0.5), mMax[1]);
	out[2] = min(max(0.0, pt[2]), mMax[2]);
	return out;
}
