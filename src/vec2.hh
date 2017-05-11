#pragma once
#include <cmath>

struct Vec2 {
    float x;
    float y;
    Vec2(float _x, float _y): x(_x), y(_y) {};
    Vec2(): x(0), y(0) {};
    float dot(Vec2 const & v) const {
        return x * v.x + y * v.y;
    }
    float len() const {
        return std::sqrt(this->dot(*this));
    }
    Vec2 rotate(float angle) {
        float s = std::sin(angle);
        float c = std::cos(angle);
        return {x * c - y * s, x * s + y * c};
    }
};

inline bool operator==(Vec2 const & lhs, Vec2 const & rhs) {
    return (lhs.x == rhs.x) && (lhs.y == rhs.y);
}

inline Vec2 & operator+=(Vec2 & lhs, Vec2 const & rhs) {
    lhs.x += rhs.x;
    lhs.y += rhs.y;
    return lhs;
}
inline Vec2 & operator-=(Vec2 & lhs, Vec2 const & rhs) {
    lhs.x -= rhs.x;
    lhs.y -= rhs.y;
    return lhs;
}
inline Vec2 & operator*=(Vec2 & lhs, Vec2 const & rhs) {
    lhs.x *= rhs.x;
    lhs.y *= rhs.y;
    return lhs;
}
inline Vec2 & operator/=(Vec2 & lhs, Vec2 const & rhs) {
    lhs.x /= rhs.x;
    lhs.y /= rhs.y;
    return lhs;
}

template <typename T>
inline Vec2 & operator+=(Vec2 & lhs, T rhs) {
    lhs.x += (float) rhs;
    lhs.y += (float) rhs;
    return lhs;
}
template <typename T>
inline Vec2 & operator-=(Vec2 & lhs, T rhs) {
    lhs.x -= (float) rhs;
    lhs.y -= (float) rhs;
    return lhs;
}
template <typename T>
inline Vec2 & operator*=(Vec2 & lhs, T rhs) {
    lhs.x *= (float) rhs;
    lhs.y *= (float) rhs;
    return lhs;
}
template <typename T>
inline Vec2 & operator/=(Vec2 & lhs, T rhs) {
    lhs.x /= (float) rhs;
    lhs.y /= (float) rhs;
    return lhs;
}

template <typename T>
inline Vec2 operator+(Vec2 lhs, T rhs) {return lhs += rhs;}
template <typename T>
inline Vec2 operator-(Vec2 lhs, T rhs) {return lhs -= rhs;}
template <typename T>
inline Vec2 operator*(Vec2 lhs, T rhs) {return lhs *= rhs;}
template <typename T>
inline Vec2 operator/(Vec2 lhs, T rhs) {return lhs /= rhs;}

inline Vec2 operator*(float lhs, Vec2 rhs) {return rhs *= lhs;}
