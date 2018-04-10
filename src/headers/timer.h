#pragma once

#include <time.h>
#include <sys/time.h>
typedef double MMC_TIME_TYPE;

namespace mmc
{
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
} // namespace mmc

