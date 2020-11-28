// Copyright(c) 2015-present, Gabi Melman & spdlog contributors.
// Distributed under the MIT License (http://opensource.org/licenses/MIT)

#ifndef UTILS_LOGGER_SPDLOG_DETAILS_LOG_MSG_BUFFER_HPP_
#define UTILS_LOGGER_SPDLOG_DETAILS_LOG_MSG_BUFFER_HPP_

#include "./log_msg.hpp"

namespace spdlog {
namespace details {

// Extend log_msg with internal buffer to store its payload.
// This is needed since log_msg holds string_views that points to stack data.

class log_msg_buffer : public log_msg
{
    memory_buf_t buffer;
    void update_string_views();

public:
    log_msg_buffer() = default;
    explicit log_msg_buffer(const log_msg &orig_msg);
    log_msg_buffer(const log_msg_buffer &other);
    log_msg_buffer(log_msg_buffer &&other) SPDLOG_NOEXCEPT;
    log_msg_buffer &operator=(const log_msg_buffer &other);
    log_msg_buffer &operator=(log_msg_buffer &&other) SPDLOG_NOEXCEPT;
};

} // namespace details
} // namespace spdlog

#ifdef SPDLOG_HEADER_ONLY
#include "log_msg_buffer-inl.hpp"
#endif
#endif  // UTILS_LOGGER_SPDLOG_DETAILS_LOG_MSG_BUFFER_HPP_
