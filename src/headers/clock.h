#pragma once
#include "timer.h"

namespace mmc
{
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
} // namespace mmc