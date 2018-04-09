#pragma once

#include <fstream>
#include <Partio.h>
#include "mac_grid.h"
#include "custom_output.h" 
#include "globals.h"

class Camera;

class SmokeSim
{
public:
   SmokeSim();
   virtual ~SmokeSim();

   virtual void reset();
   virtual void step();
   virtual void draw(const Camera& c);
   //virtual void setGridDimensions(int x, int y, int z); 
   virtual void setRecording(bool on, int width, int height);
   virtual bool isRecording();
	
	
	int getTotalFrames();

protected:
   virtual void drawAxes();
   virtual void grabScreen();

protected:
	MACGrid mGrid;
	
	int mFrameNum;
	int mTotalFrameNum;
	bool mRecordEnabled;
	
	int recordWidth;
	int recordHeight;
};