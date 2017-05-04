#include "vec2.hh"
/*
* ## Structure for results
* TODO: find better name for it
*/
struct Data
{
    /* average values */
    Vec2 v;
    std::vector<float> power;

    int acoustic_phonon_scattering_count;
    int optical_phonon_scattering_count;
    int vertical_transitions_count;

    float tau;
    std::vector<float> population;

    Data() : v(Vec2(0, 0)),
        power(0),
        tau(0),
        acoustic_phonon_scattering_count(0),
        optical_phonon_scattering_count(0),
        vertical_transitions_count(0) {};
};

std::vector<float> & operator+=(std::vector<float> & lhs, std::vector<float> const & rhs) {
    int len = std::max(lhs.size(), rhs.size());
    for (int i = 0; i < len; i++) {
        lhs[i] += rhs[i];
    }
    return lhs;
}

std::vector<float> & operator-=(std::vector<float> & lhs, std::vector<float> const & rhs) {
    int len = std::max(lhs.size(), rhs.size());
    for (int i = 0; i < len; i++) {
        lhs[i] -= rhs[i];
    }
    return lhs;
}

std::vector<float> & operator*=(std::vector<float> & lhs, std::vector<float> const & rhs) {
    int len = std::max(lhs.size(), rhs.size());
    for (int i = 0; i < len; i++) {
        lhs[i] *= rhs[i];
    }
    return lhs;
}

std::vector<float> & operator*=(std::vector<float> & lhs, float rhs) {
    int len = lhs.size();
    for (int i = 0; i < len; i++) {
        lhs[i] *= rhs;
    }
    return lhs;
}

std::vector<float> & operator/=(std::vector<float> & lhs, float rhs) {
    int len = lhs.size();
    for (int i = 0; i < len; i++) {
        lhs[i] /= rhs;
    }
    return lhs;
}

bool operator>(std::vector<float> & lhs, float rhs) {
    int len = lhs.size();
    for (int i = 0; i < len; i++) {
        if (lhs[i] <= rhs)
            return false;
    }
    return true;
}

std::vector<float> sqrt(std::vector<float> lhs) {
    for (auto & x: lhs) {
        x = sqrt(x);
    }
    return lhs;
}

inline Data & operator+=(Data & lhs, const Data & rhs) {
    lhs.v += rhs.v;
    lhs.tau += rhs.tau;
    lhs.power += rhs.power;
    lhs.acoustic_phonon_scattering_count += rhs.acoustic_phonon_scattering_count;
    lhs.optical_phonon_scattering_count += rhs.optical_phonon_scattering_count;
    lhs.vertical_transitions_count += rhs.vertical_transitions_count;
    lhs.population += rhs.population;
    return lhs;
}

inline Data & operator-=(Data & lhs, const Data & rhs) {
    lhs.v -= rhs.v;
    lhs.tau -= rhs.tau;
    lhs.power -= rhs.power;
    lhs.acoustic_phonon_scattering_count -= rhs.acoustic_phonon_scattering_count;
    lhs.optical_phonon_scattering_count -= rhs.optical_phonon_scattering_count;
    lhs.vertical_transitions_count -= rhs.vertical_transitions_count;
    lhs.population -= rhs.population;
    return lhs;
}

inline Data & operator*=(Data & lhs, const Data & rhs) {
    lhs.v *= rhs.v;
    lhs.tau *= rhs.tau;
    lhs.power *= rhs.power;
    lhs.acoustic_phonon_scattering_count *= rhs.acoustic_phonon_scattering_count;
    lhs.optical_phonon_scattering_count *= rhs.optical_phonon_scattering_count;
    lhs.vertical_transitions_count *= rhs.vertical_transitions_count;
    lhs.population *= rhs.population;
    return lhs;
}

template <typename T>
inline Data & operator*=(Data & lhs, T rhs) {
    lhs.v *= rhs;
    lhs.tau *= rhs;
    lhs.power *= rhs;
    lhs.acoustic_phonon_scattering_count *= rhs;
    lhs.optical_phonon_scattering_count *= rhs;
    lhs.vertical_transitions_count *= rhs;
    lhs.population *= rhs;
    return lhs;
}

template <typename T>
inline Data & operator/=(Data & lhs, T rhs) {
    lhs.v /= rhs;
    lhs.tau /= rhs;
    lhs.power /= rhs;
    lhs.acoustic_phonon_scattering_count /= rhs;
    lhs.optical_phonon_scattering_count /= rhs;
    lhs.vertical_transitions_count /= rhs;
    lhs.population /= rhs;
    return lhs;
}

inline Data sqrt(Data d) {
    d.v = (d.v.x > 0 && d.v.y > 0) ? Vec2(sqrt(d.v.x), sqrt(d.v.y)) : Vec2(0, 0);
    d.power = (d.power > 0) ? sqrt(d.power) : std::vector<float>(d.power.size(), 0);
    d.tau = (d.tau > 0) ? sqrt(d.tau) : 0;
    d.acoustic_phonon_scattering_count = (d.acoustic_phonon_scattering_count > 0) ? sqrt(d.acoustic_phonon_scattering_count) : 0;
    d.optical_phonon_scattering_count = (d.optical_phonon_scattering_count > 0) ? sqrt(d.optical_phonon_scattering_count) : 0;
    d.vertical_transitions_count = (d.vertical_transitions_count > 0) ? sqrt(d.vertical_transitions_count) : 0;
    d.population = (d.population > 0) ? sqrt(d.population) : std::vector<float>(d.population.size(), 0);
    return d;
}

inline Data operator+(Data lhs, Data const & rhs) {return lhs += rhs;}
inline Data operator-(Data lhs, Data const & rhs) {return lhs -= rhs;}
template <typename T>
inline Data operator*(Data lhs, T rhs) {return lhs *= rhs;}
template <typename T>
inline Data operator/(Data lhs, T rhs) {return lhs /= rhs;}

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
