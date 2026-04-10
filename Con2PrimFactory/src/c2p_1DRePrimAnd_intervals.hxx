#ifndef C2P_1DREPRIMAND_INTERVALS_HXX
#define C2P_1DREPRIMAND_INTERVALS_HXX

namespace Con2PrimFactory {

template <class T> class interval {
  T lo_, hi_;

public:
  using value_t = T;

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE interval()
      : lo_(T(0)), hi_(T(0)) {}

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE interval(T lo, T hi)
      : lo_(lo), hi_(hi) {}

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE T min() const {
    return lo_;
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE T max() const {
    return hi_;
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE T &min() { return lo_; }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE T &max() { return hi_; }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE T width() const {
    return hi_ - lo_;
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE T mid() const {
    return (lo_ + hi_) * T(0.5);
  }
};

} // namespace Con2PrimFactory
#endif
