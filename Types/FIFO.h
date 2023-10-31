//
// Created by Syl on 14/07/2023.
//

#ifndef SDSP_FIFO_H
#define SDSP_FIFO_H
#include "../Macros.h"
#include <vector>
namespace SDSP {
    class FIFO {
    public:
        explicit FIFO(size_t initialSize) {
            m_buffer.reserve(initialSize);
        }

        SDSP_INLINE void push(float toPush) noexcept {
            m_buffer.emplace_back(toPush);
        }

        [[nodiscard]] SDSP_INLINE float pop(std::ptrdiff_t offset = 0, bool erase = true) noexcept {
            if(m_buffer.empty()) return 0.0f;
            auto it = m_buffer.begin() + offset;
            if(it == m_buffer.end()) return 0.0f;
            auto out = *it;
            if(erase) {
                m_buffer.erase(m_buffer.begin() + offset);
            }
            return out;
        }

        SDSP_INLINE void clear() noexcept {
            m_buffer.clear();
        }
    private:
        std::vector<float> m_buffer;

    };
}
#endif //KALIDEPROTOTYPES2_FIFO_H
