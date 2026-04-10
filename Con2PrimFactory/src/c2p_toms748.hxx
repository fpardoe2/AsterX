#ifndef C2P_TOMS748_HXX
#define C2P_TOMS748_HXX

#include <algorithm>
#include <cmath>
#include <limits>
#include <type_traits>
#include <utility>

#include "c2p_1DRePrimAnd_intervals.hxx"
#include "c2p_utils.hxx"

namespace Con2PrimFactory {

namespace toms748_detail {

template <class T>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline bool same_sign(T a,
                                                                        T b) {
  return (a > T(0) && b > T(0)) || (a < T(0) && b < T(0));
}

template <class T>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline bool
in_open_interval(T x, T a, T b) {
  return (x > a) && (x < b);
}

template <class T>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline T midpoint(T a,
                                                                     T b) {
  return (a + b) * T(0.5);
}

template <class T>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline T
secant_step(T a, T fa, T b, T fb) {
  const T denom = fb - fa;
  if (denom == T(0))
    return std::numeric_limits<T>::quiet_NaN();
  return (a * fb - b * fa) / denom;
}

template <class T>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline T
iqi_step(T a, T fa, T b, T fb, T c, T fc) {
  const T dab = (fa - fb);
  const T dac = (fa - fc);
  const T dba = (fb - fa);
  const T dbc = (fb - fc);
  const T dca = (fc - fa);
  const T dcb = (fc - fb);

  if (dab == T(0) || dac == T(0) || dba == T(0) || dbc == T(0) ||
      dca == T(0) || dcb == T(0))
    return std::numeric_limits<T>::quiet_NaN();

  return a * fb * fc / (dab * dac) + b * fa * fc / (dba * dbc) +
         c * fa * fb / (dca * dcb);
}

template <class T>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline interval<T>
invalid_interval() {
  return {std::numeric_limits<T>::lowest(), std::numeric_limits<T>::max()};
}

template <class F, class T> struct has_stopif {
private:
  template <class U>
  static auto test(int)
      -> decltype(std::declval<U &>().stopif(std::declval<T>(),
                                             std::declval<T>(),
                                             std::declval<T>()),
                  std::true_type{});
  template <class> static std::false_type test(...);

public:
  static constexpr bool value = decltype(test<F>(0))::value;
};

template <class F, class T>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline bool
call_stopif_impl(F &f, T left, T right, T tol, std::true_type) {
  // Preserve RePrimAnd contract: stopif(left, right-left, tol)
  return f.stopif(left, right - left, tol);
}

template <class F, class T>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline bool
call_stopif_impl(F &, T left, T right, T tol, std::false_type) {
  const T width = std::abs(right - left);
  const T scale = std::max(T(1), std::max(std::abs(left), std::abs(right)));
  return width <= tol * scale;
}

template <class F, class T>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline bool
call_stopif(F &f, T left, T right, T tol) {
  using has_stopif_t = std::integral_constant<bool, has_stopif<F, T>::value>;
  return call_stopif_impl(f, left, right, tol, has_stopif_t{});
}

template <class T>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
update_bracket(T x, T fx, T &a, T &fa, T &b, T &fb) {
  if (same_sign(fa, fx)) {
    a = x;
    fa = fx;
  } else {
    b = x;
    fb = fx;
  }
}

} // namespace toms748_detail

template <class F, class T = typename F::value_t>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline interval<T>
toms748_solve(F &f, interval<T> bracket, T tol, unsigned int max_calls,
              ROOTSTAT &errs) {
  using namespace toms748_detail;

  T a = bracket.min();
  T b = bracket.max();
  if (b < a) {
    const T tmp = a;
    a = b;
    b = tmp;
  }

  if (!std::isfinite(a) || !std::isfinite(b)) {
    errs = ROOTSTAT::NOT_BRACKETED;
    return invalid_interval<T>();
  }

  T fa = f(a);
  T fb = f(b);
  unsigned int calls = 2u;

  // Match RePrimAnd endpoint precedence: right endpoint first.
  if (fb == T(0)) {
    errs = ROOTSTAT::SUCCESS;
    return {b, b};
  }
  if (fa == T(0)) {
    errs = ROOTSTAT::SUCCESS;
    return {a, a};
  }

  if (!std::isfinite(fa) || !std::isfinite(fb) || same_sign(fa, fb)) {
    errs = ROOTSTAT::NOT_BRACKETED;
    return invalid_interval<T>();
  }

  if (max_calls <= 2u) {
    errs = ROOTSTAT::NOT_CONVERGED;
    return {a, b};
  }

  if (call_stopif(f, a, b, tol)) {
    errs = ROOTSTAT::SUCCESS;
    return {a, b};
  }

  T prev_width = b - a;

  while (calls < max_calls) {
    const T a0 = a, b0 = b;
    const T fa0 = fa, fb0 = fb;
    const T width0 = b - a;

    T x1 = secant_step(a, fa, b, fb);
    if (!std::isfinite(x1) || !in_open_interval(x1, a, b))
      x1 = midpoint(a, b);

    T fx1 = f(x1);
    ++calls;

    if (!std::isfinite(fx1)) {
      if (calls >= max_calls)
        break;
      x1 = midpoint(a, b);
      fx1 = f(x1);
      ++calls;
      if (!std::isfinite(fx1)) {
        errs = ROOTSTAT::NOT_CONVERGED;
        return {a, b};
      }
    }

    if (fx1 == T(0)) {
      errs = ROOTSTAT::SUCCESS;
      return {x1, x1};
    }

    update_bracket(x1, fx1, a, fa, b, fb);
    if (call_stopif(f, a, b, tol)) {
      errs = ROOTSTAT::SUCCESS;
      return {a, b};
    }
    if (calls >= max_calls)
      break;

    T x2 = iqi_step(a0, fa0, b0, fb0, x1, fx1);
    if (!std::isfinite(x2) || !in_open_interval(x2, a, b))
      x2 = midpoint(a, b);

    const T width1 = b - a;
    if ((x2 - a) < T(0.125) * width1 || (b - x2) < T(0.125) * width1)
      x2 = midpoint(a, b);

    T fx2 = f(x2);
    ++calls;

    if (!std::isfinite(fx2)) {
      if (calls >= max_calls)
        break;
      x2 = midpoint(a, b);
      fx2 = f(x2);
      ++calls;
      if (!std::isfinite(fx2)) {
        errs = ROOTSTAT::NOT_CONVERGED;
        return {a, b};
      }
    }

    if (fx2 == T(0)) {
      errs = ROOTSTAT::SUCCESS;
      return {x2, x2};
    }

    update_bracket(x2, fx2, a, fa, b, fb);
    if (call_stopif(f, a, b, tol)) {
      errs = ROOTSTAT::SUCCESS;
      return {a, b};
    }

    const T new_width = b - a;
    if (calls < max_calls &&
        (new_width > T(0.5) * width0 || new_width > T(0.9) * prev_width)) {
      const T xm = midpoint(a, b);
      const T fxm = f(xm);
      ++calls;

      if (!std::isfinite(fxm)) {
        errs = ROOTSTAT::NOT_CONVERGED;
        return {a, b};
      }
      if (fxm == T(0)) {
        errs = ROOTSTAT::SUCCESS;
        return {xm, xm};
      }

      update_bracket(xm, fxm, a, fa, b, fb);
      if (call_stopif(f, a, b, tol)) {
        errs = ROOTSTAT::SUCCESS;
        return {a, b};
      }
    }

    prev_width = b - a;
  }

  errs = ROOTSTAT::NOT_CONVERGED;
  return {a, b};
}

} // namespace Con2PrimFactory

#endif
