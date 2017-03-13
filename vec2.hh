#pragma once
#include <cmath>

struct Vec2 {
    float x;
    float y;
    Vec2(float _x, float _y): x(_x), y(_y) {};
    float dot(Vec2 const & v) {
        return x * v.x + y * v.y;
    }
    float len() {
        return sqrt(this.dot(this));
    }
};

bool operator==(Vec2 const & lhs, Vec2 const & rhs) {
    return (lhs.x == rhs.x) && (lhs.y == rhs.y);
}

Vec2 & operator+=(Vec2 & lhs, Vec2 const & rhs) {
    lhs.x += rhs.x;
    lhs.y += rhs.y;
    return lhs;
}
Vec2 & operator-=(Vec2 & lhs, Vec2 const & rhs) {
    lhs.x -= rhs.x;
    lhs.y -= rhs.y;
    return lhs;
}
Vec2 & operator*=(Vec2 & lhs, Vec2 const & rhs) {
    lhs.x *= rhs.x;
    lhs.y *= rhs.y;
    return lhs;
}
Vec2 & operator/=(Vec2 & lhs, Vec2 const & rhs) {
    lhs.x /= rhs.x;
    lhs.y /= rhs.y;
    return lhs;
}

template <typename T>
Vec2 & operator+=(Vec2 & lhs, T rhs) {
    lhs.x += (float) rhs;
    lhs.y += (float) rhs;
    return lhs;
}
template <typename T>
Vec2 & operator-=(Vec2 & lhs, T rhs) {
    lhs.x -= (float) rhs;
    lhs.y -= (float) rhs;
    return lhs;
}
template <typename T>
Vec2 & operator*=(Vec2 & lhs, T rhs) {
    lhs.x *= (float) rhs;
    lhs.y *= (float) rhs;
    return lhs;
}
template <typename T>
Vec2 & operator/=(Vec2 & lhs, T rhs) {
    lhs.x /= (float) rhs;
    lhs.y /= (float) rhs;
    return lhs;
}

template <typename T>
Vec2 operator+(Vec2 lhs, T rhs) {return lhs += rhs;}
template <typename T>
Vec2 operator-(Vec2 lhs, T rhs) {return lhs -= rhs;}
template <typename T>
Vec2 operator*(Vec2 lhs, T rhs) {return lhs *= rhs;}
template <typename T>
Vec2 operator/(Vec2 lhs, T rhs) {return lhs /= rhs;}