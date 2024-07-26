// created by yuqing.wu on 06/06/2024

#include <algorithm>
#include <cmath>
#include <vector>
#include <string>
#include <memory>

#ifndef GEOM_H
#define GEOM_H

// namespace global_path_search {

typedef unsigned long long ull;

struct Point3D{
    double x, y, z;
    Point3D(double X, double Y, double Z) : x(X), y(Y), z(Z) {}
    Point3D() : Point3D(0., 0., 0.) {}
    Point3D(std::vector<double> XYZ) {
        if(XYZ.size() >= 3) {
            x = XYZ[0];
            y = XYZ[1];
            z = XYZ[2];
        } else {
            Point3D();
        }
    }
    Point3D operator[] (int bias) const {
        return {((bias & 1) ? x : -x), ((bias & 2) ? y : -y), ((bias & 4) ? z : -z)};
    }
    Point3D operator+(const Point3D &rhs) const { return {x + rhs.x, y + rhs.y, z + rhs.z}; }
    Point3D operator-(const Point3D &rhs) const { return {x - rhs.x, y - rhs.y, z - rhs.z}; }
    Point3D operator*(double k) const { return {x * k, y * k, z * k}; }
    Point3D operator* (int v) const {
        return {x * v, y * v, z * v};
    }
    bool operator!= (Point3D point) const {
        return !(point.x == x && point.y == y && point.z == z);
    }
    bool operator== (Point3D point) const {
        return (point.x == x && point.y == y && point.z == z);
    }
    friend bool operator < (const Point3D &x, const Point3D &y) {
        if (x.x != y.x) return x.x < y.x;
        if (x.y != y.y) return x.y < y.y;
        return x.z < y.z;
    }
    std::string ToString() const { return "(" + std::to_string(x) + ", " + std::to_string(y) + ", " + std::to_string(z) + ")"; }
};

struct Point3U {
    ull x, y, z;
    Point3U (ull X, ull Y, ull Z) : x(X), y(Y), z(Z) {}
    Point3U () : Point3U(0, 0, 0) {}
    Point3U operator[] (int bias) const {
        if (bias == 0) return {x, 0, 0};
        if (bias == 5) return {-x, 0, 0};
        if (bias == 1) return {0, y, 0};
        if (bias == 4) return {0, -y, 0};
        if (bias == 2) return {0, 0, z};
        if (bias == 3) return {0, 0, -z};
        return {0, 0, 0};
    }
    Point3U operator* (int v) const {
        return {x * v, y * v, z * v};
    }
    Point3U operator+ (Point3U point) const {
        return {x + point.x, y + point.y, z + point.z};
    }
    std::string ToString() const { return "(" + std::to_string(x) + ", " + std::to_string(y) + ", " + std::to_string(z) + ")"; }
};

inline double GetDis(Point3D a, Point3D b) {
    return std::sqrt(std::pow(a.x - b.x, 2.) + std::pow(a.y - b.y, 2.) + std::pow(a.z - b.z, 2.));
}

inline double fd(Point3D x, Point3D y) {
    return sqrt((x.x - y.x) * (x.x - y.x) + (x.y - y.y) * (x.y - y.y) + (x.z - y.z) * (x.z - y.z));
}

inline Point3D Interpolation(Point3D x, Point3D y, double ratio) {
    if (ratio <= 0.) {
        return x;
    } else if (ratio >= 1) {
        return y;
    }
    return x * ratio + y * (1. - ratio);
}

// } // namespace global_path_search

#endif // GEOM_H