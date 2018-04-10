#pragma once

#include <vector>
#include <string.h>
#include "clock.h"

namespace mmc
{
    class Clock;

    // brief Utility class for measuring framerate
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
} // namespace mmc