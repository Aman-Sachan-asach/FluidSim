#pragma once

#include <time.h>
#include <sys/time.h>
typedef double MMC_TIME_TYPE;

//------------------------------------------------------------------------
// NOTE: Use the Clock and FPS classes as interfaces to the timer class //
//------------------------------------------------------------------------

// brief Class for timing application events.
// Only one of these is needed per application.
class Timer
{
public:
    // Get a new timer, initially stopped
    inline Timer() : start_(0.0), last_(0.0), now_(0.0), elapsed_(0.0), freq_(1.0), invFreq_(1.0)
    {}

    // Start the timer
    inline void start ()
    {
        struct timeval tv;
        struct timezone tz;
        gettimeofday(&tv, &tz);
        start_ = (double) (tv.tv_sec + 1e-6 * tv.tv_usec);
        last_ = now_ = start_;
    }

    // Call this once per frame to advance internal timer state
    inline void inc ()
    {
        last_ = now_;
        struct timeval tv;
        struct timezone tz;
        gettimeofday(&tv, &tz);
        now_ = (double) (tv.tv_sec + 1e-6 * tv.tv_usec);
    }

    // Returns number of ms elapsed since start() was called. Use this to determine how long the timer has been running
    inline MMC_TIME_TYPE queryElapsed () const { return now_ - start_; }

    // Returns number of ms elapsed between the two previous calls to \p inc(). Use this to determine the time between frames
    inline MMC_TIME_TYPE queryInc () const { return now_ - last_; }

    // Returns (1.0 / frequency) of the internal timing clock actually being used
    inline double getInvFreq () const { return invFreq_; }

private:
    MMC_TIME_TYPE start_;
    MMC_TIME_TYPE last_;
    MMC_TIME_TYPE now_;
    MMC_TIME_TYPE elapsed_;
    MMC_TIME_TYPE freq_;
    double invFreq_;
};

/////////////////////////////////////////////////////////////////////

// High precision clock for timing events. 
// This class is implemented on top of the Timer class. Use as an interface to the Timer class
class Clock
{
public:
    inline Clock (): paused_(false)
    {
        timer_ = new Timer;
        invFreq_ = timer_->getInvFreq();
        timer_->start();
        reset();
    }

    // Reset the clock (0 time elapsed)
    inline void reset ()
    {
        start_ = timer_->queryElapsed();
        inc_ = curTime_ = 0;
    }

    // Call once per frame to update the internal clock state
    inline void inc ()
    {
        timer_->inc();
        if (!paused_)
        {
            inc_ = timer_->queryInc();
            curTime_ += inc_;
        }
        else
        {
            inc_ = 0;
        }
    }

    // All Query functions convert internal clock time to long in milliseconds.
    // Returns the amount of time (in ms) elapsed between last two calls to \p inc()
    inline long queryInc () const { return (long) (1000.0 * inc_ * invFreq_); }

    // Returns the amount of time (in ms) elapsed since clock creation or reset() was called
    inline long queryTime () const { return (long) (1000.0 * curTime_ * invFreq_); }

    // Pause the clock. */
    inline void pauseToggle ()
    {
        if (paused_)
        {
            paused_ = false;
        }
        else
        {
            paused_ = true;
            inc_ = 0;
        }
    }

private:
    Timer *timer_;
    MMC_TIME_TYPE inc_;
    MMC_TIME_TYPE curTime_;
    MMC_TIME_TYPE start_;
    bool paused_;
    double invFreq_;
};

/////////////////////////////////////////////////////////////////////

// Utility class for measuring framerate
class FpsTracker
{
public:
    FpsTracker(int smoothSteps = 4) : clock_(new Clock()), steps_(smoothSteps), nSnaps_(0), snaps_(new long[steps_])
    {
        // smoothSteps is the window size over which to average the time measurements
        memset(snaps_, 0, sizeof(long) * steps_);    
    }

    ~FpsTracker ()
    {
        delete snaps_;
        delete clock_;
    }

    // Specify the window size
    void setNumSteps (int smoothSteps)
    {
        // smoothSteps over which the time measurements will be averaged
        steps_ = smoothSteps;
        delete snaps_;
        snaps_ = new long[steps_];
        memset(snaps_, 0, sizeof(long) * steps_);
        nSnaps_ = 0;
    }

    // Makes a timestamp; measures time intervals between successful calls
    void timestamp ()
    {
        clock_->inc();
        snaps_[nSnaps_ % steps_] = clock_->queryInc();
        nSnaps_++;
    }

    // Get the average FPS (averaged over 'smoothSteps' interval)
    float fpsAverage () const
    {
        int count = (nSnaps_ < steps_) ? nSnaps_ : steps_;
        long sum = 0L;
        for (int i = 0; i < count; ++i)
        {
            sum += snaps_[i];
        }
        if (!sum) sum = 1L; // prevent div-by-zero
        return 1000.0f * count / (float) sum;
    }

    // Get the instantaneous FPS (estimated from the last frame only)
    float fpsInstant () const
    {
        long inc = clock_->queryInc();
        if (!inc) inc = 1L; // prevent div-by-zero
        return 1000.0f / (float) inc;
    }

private:
    Clock *clock_;
    int steps_;
    int nSnaps_;
    long *snaps_;
};

