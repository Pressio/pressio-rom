// Copyright(c) 2015-present, Gabi Melman & spdlog contributors.
// Distributed under the MIT License (http://opensource.org/licenses/MIT)

#ifndef UTILS_LOGGER_SPDLOG_SINKS_SINK_HPP_
#define UTILS_LOGGER_SPDLOG_SINKS_SINK_HPP_

//#include "../details/log_msg.hpp"
// #include "../formatter.hpp"
// #include "../common.hpp"

namespace spdlog { namespace sinks {

class sink
{
  public:
  virtual ~sink() = default;
  virtual void log(const details::log_msg &msg) = 0;
  virtual void flush() = 0;
  virtual void set_pattern(const std::string &pattern) = 0;
  virtual void set_formatter(std::unique_ptr<spdlog::formatter> sink_formatter) = 0;

  void set_level(level::level_enum log_level)
  {
    level_.store(log_level, std::memory_order_relaxed);
  }

  level::level_enum level() const
  {
    return static_cast<spdlog::level::level_enum>(level_.load(std::memory_order_relaxed));
  }

  bool should_log(level::level_enum msg_level) const
  {
    return msg_level >= level_.load(std::memory_order_relaxed);
  }

  void setMpiRank(int rank) { mpiRank_ = rank; }

protected:
  int mpiRank_ = 0;

  // sink log level - default is all
  level_t level_{level::trace};
  };

} // namespace sinks
} // namespace spdlog

// #ifdef SPDLOG_HEADER_ONLY
// #include "sink-inl.hpp"
// #endif
#endif  // UTILS_LOGGER_SPDLOG_SINKS_SINK_HPP_
