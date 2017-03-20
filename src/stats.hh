#include "vec2.hh"
/*
* ## Structure for results
*/
struct Data
{
    /* average values */
    Vec2 v;
    float power;

    int acoustic_phonon_scattering_count;
    int optical_phonon_scattering_count;

    float tau;
    
    Data() : v(Vec2(0, 0)),
        power(0),
        tau(0),
        acoustic_phonon_scattering_count(0),
        optical_phonon_scattering_count(0) {};
};


inline Data & operator+=(Data & lhs, const Data & rhs) {
    lhs.v += rhs.v;
    lhs.tau += rhs.tau;
    lhs.power += rhs.power;
    lhs.acoustic_phonon_scattering_count += rhs.acoustic_phonon_scattering_count;
    lhs.optical_phonon_scattering_count += rhs.optical_phonon_scattering_count;
    return lhs;
}

inline Data & operator-=(Data & lhs, const Data & rhs) {
    lhs.v -= rhs.v;
    lhs.tau -= rhs.tau;
    lhs.power -= rhs.power;
    lhs.acoustic_phonon_scattering_count -= rhs.acoustic_phonon_scattering_count;
    lhs.optical_phonon_scattering_count -= rhs.optical_phonon_scattering_count;
    return lhs;
}

inline Data & operator*=(Data & lhs, const Data & rhs) {
    lhs.v *= rhs.v;
    lhs.tau *= rhs.tau;
    lhs.power *= rhs.power;
    lhs.acoustic_phonon_scattering_count *= rhs.acoustic_phonon_scattering_count;
    lhs.optical_phonon_scattering_count *= rhs.optical_phonon_scattering_count;
    return lhs;
}

inline Data & operator/=(Data & lhs, int rhs) {
    lhs.v /= rhs;
    lhs.tau /= rhs;
    lhs.power /= rhs;
    lhs.acoustic_phonon_scattering_count /= rhs;
    lhs.optical_phonon_scattering_count /= rhs;
    return lhs;
}

inline Data sqrt(Data d) {
    d.v = (d.v.x > 0 && d.v.y > 0) ? Vec2(sqrt(d.v.x), sqrt(d.v.y)) : Vec2(0, 0);
    d.power = (d.power > 0) ? sqrt(d.power) : 0;
    d.tau = (d.tau > 0) ? sqrt(d.tau) : 0;
    d.acoustic_phonon_scattering_count = (d.acoustic_phonon_scattering_count > 0) ? sqrt(d.acoustic_phonon_scattering_count) : 0;
    d.optical_phonon_scattering_count = (d.optical_phonon_scattering_count > 0) ? sqrt(d.optical_phonon_scattering_count) : 0;
    return d;
}

inline Data operator+(Data lhs, Data const & rhs) {return lhs += rhs;}
inline Data operator-(Data lhs, Data const & rhs) {return lhs -= rhs;}
inline Data operator*(Data lhs, Data const & rhs) {return lhs *= rhs;}
inline Data operator/(Data lhs, int rhs) {return lhs /= rhs;}

template <typename T>
inline T mean(std::vector<T> const & values) {
    T res;
    for (T x: values) {
        res += x;
    }
    return res / values.size();
}

template <typename T>
inline T disp(std::vector<T> const & values) {
    T res;
    for (T x: values) {
        res += x * x;
    }
    return res / values.size() - mean(values) * mean(values);
}

template <typename T>
inline T stdev(std::vector<T> const & values) {
    return sqrt(disp(values));
}
