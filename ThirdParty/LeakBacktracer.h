//
// Created by Syl on 12/04/2023.
//
#include <juce_core/juce_core.h>
namespace SDSP
{
    class AdvancedLeakDetector {
    public:
        AdvancedLeakDetector() {
            getBackTraceHash().set((void*)this, juce::SystemStats::getStackBacktrace());
        }
        ~AdvancedLeakDetector() {
            getBackTraceHash().remove((void*)this);
        }
    private:
        typedef juce::HashMap<void*, juce::String> BackTraceHash;
        struct HashHolder {
            ~HashHolder() {
                if (traces.size() > 0)
                {
                    /* Memory leak info. */
                    DBG("Found " + juce::String(traces.size()) + " possible leaks");
                    for (BackTraceHash::Iterator i(traces); i.next();)
                    {
                        DBG("-----");
                        DBG(i.getValue());
                    }
                    jassertfalse;
                }
            }
            BackTraceHash traces;
        };
        BackTraceHash& getBackTraceHash() {
            static HashHolder holder;
            return holder.traces;
        }
    };
}
