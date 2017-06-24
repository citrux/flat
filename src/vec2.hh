#pragma once
#include <cmath>

struct Vec2 {
  double x;
  double y;
  Vec2(double _x, double _y) : x(_x), y(_y){};
  Vec2() : x(0), y(0){};
  double dot(Vec2 const &v) const { return x * v.x + y * v.y; }
  double len() const { return std::sqrt(this->dot(*this)); }
  Vec2 rotate(double angle) {
    double s = std::sin(angle);
    double c = std::cos(angle);
    return {x * c - y * s, x * s + y * c};
  }
};

inline bool operator==(Vec2 const &lhs, Vec2 const &rhs) {
  return (lhs.x == rhs.x) && (lhs.y == rhs.y);
}

inline Vec2 &operator+=(Vec2 &lhs, Vec2 const &rhs) {
  lhs.x += rhs.x;
  lhs.y += rhs.y;
  return lhs;
}
inline Vec2 &operator-=(Vec2 &lhs, Vec2 const &rhs) {
  lhs.x -= rhs.x;
  lhs.y -= rhs.y;
  return lhs;
}
inline Vec2 &operator*=(Vec2 &lhs, Vec2 const &rhs) {
  lhs.x *= rhs.x;
  lhs.y *= rhs.y;
  return lhs;
}
inline Vec2 &operator/=(Vec2 &lhs, Vec2 const &rhs) {
  lhs.x /= rhs.x;
  lhs.y /= rhs.y;
  return lhs;
}

template <typename T> inline Vec2 &operator+=(Vec2 &lhs, T rhs) {
  lhs.x += (double)rhs;
  lhs.y += (double)rhs;
  return lhs;
}
template <typename T> inline Vec2 &operator-=(Vec2 &lhs, T rhs) {
  lhs.x -= (double)rhs;
  lhs.y -= (double)rhs;
  return lhs;
}
template <typename T> inline Vec2 &operator*=(Vec2 &lhs, T rhs) {
  lhs.x *= (double)rhs;
  lhs.y *= (double)rhs;
  return lhs;
}
template <typename T> inline Vec2 &operator/=(Vec2 &lhs, T rhs) {
  lhs.x /= (double)rhs;
  lhs.y /= (double)rhs;
  return lhs;
}

template <typename T> inline Vec2 operator+(Vec2 lhs, T rhs) {
  return lhs += rhs;
}
template <typename T> inline Vec2 operator-(Vec2 lhs, T rhs) {
  return lhs -= rhs;
}
template <typename T> inline Vec2 operator*(Vec2 lhs, T rhs) {
  return lhs *= rhs;
}
template <typename T> inline Vec2 operator/(Vec2 lhs, T rhs) {
  return lhs /= rhs;
}

inline Vec2 operator*(double lhs, Vec2 rhs) { return rhs *= lhs; }
