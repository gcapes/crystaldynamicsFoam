/*---------------------------------------------------------------------------*\
    ebsdTools :: orientation.H

    Rotation utilities (plain C++11).  g = Z(phi2) X(Phi) Z(phi1) is the passive Bunge matrix,
    v_crystal = g v_sample; R = Rz(p3) Ry(p2) Rx(p1) = g^T is the active crystaldynamicsFoam
    matrix (StressUpdate.H), v_sample = R v_crystal.  Quaternions always represent R.

    Authors: Tom Flint and co-workers (crystaldynamicsFoam)
\*---------------------------------------------------------------------------*/

#ifndef ebsdTools_orientation_H
#define ebsdTools_orientation_H

#include <cmath>
#include <vector>
#include <string>
#include <stdexcept>
#include <algorithm>

namespace ebsd
{

const double PI = 3.14159265358979323846;

inline double deg2rad(double d) { return d*PI/180.0; }
inline double rad2deg(double r) { return r*180.0/PI; }

// --- Vec3
struct Vec3
{
    double x, y, z;
    Vec3() : x(0), y(0), z(0) {}
    Vec3(double X, double Y, double Z) : x(X), y(Y), z(Z) {}
    double& operator[](int i) { return i == 0 ? x : (i == 1 ? y : z); }
    double  operator[](int i) const { return i == 0 ? x : (i == 1 ? y : z); }
    Vec3 operator+(const Vec3& b) const { return Vec3(x+b.x, y+b.y, z+b.z); }
    Vec3 operator-(const Vec3& b) const { return Vec3(x-b.x, y-b.y, z-b.z); }
    Vec3 operator*(double s) const { return Vec3(x*s, y*s, z*s); }
    Vec3 operator/(double s) const { return Vec3(x/s, y/s, z/s); }
    Vec3 operator-() const { return Vec3(-x, -y, -z); }
    double dot(const Vec3& b) const { return x*b.x + y*b.y + z*b.z; }
    Vec3 cross(const Vec3& b) const
    {
        return Vec3(y*b.z - z*b.y, z*b.x - x*b.z, x*b.y - y*b.x);
    }
    double norm() const { return std::sqrt(dot(*this)); }
    Vec3 normalised() const
    {
        double n = norm();
        if (n < 1e-300) throw std::runtime_error("Vec3::normalised: zero vector");
        return (*this)/n;
    }
};
inline Vec3 operator*(double s, const Vec3& v) { return v*s; }


// --- Mat3
struct Mat3
{
    double m[3][3];

    Mat3() { for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) m[i][j] = 0.0; }

    static Mat3 identity()
    {
        Mat3 I; I.m[0][0] = I.m[1][1] = I.m[2][2] = 1.0; return I;
    }

    double& operator()(int i, int j) { return m[i][j]; }
    double  operator()(int i, int j) const { return m[i][j]; }

    Mat3 T() const
    {
        Mat3 t;
        for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) t.m[i][j] = m[j][i];
        return t;
    }

    Mat3 operator*(const Mat3& b) const
    {
        Mat3 c;
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
            {
                double s = 0.0;
                for (int k = 0; k < 3; k++) s += m[i][k]*b.m[k][j];
                c.m[i][j] = s;
            }
        return c;
    }

    Vec3 operator*(const Vec3& v) const
    {
        return Vec3
        (
            m[0][0]*v.x + m[0][1]*v.y + m[0][2]*v.z,
            m[1][0]*v.x + m[1][1]*v.y + m[1][2]*v.z,
            m[2][0]*v.x + m[2][1]*v.y + m[2][2]*v.z
        );
    }

    double trace() const { return m[0][0] + m[1][1] + m[2][2]; }

    double det() const
    {
        return m[0][0]*(m[1][1]*m[2][2] - m[1][2]*m[2][1])
             - m[0][1]*(m[1][0]*m[2][2] - m[1][2]*m[2][0])
             + m[0][2]*(m[1][0]*m[2][1] - m[1][1]*m[2][0]);
    }

    // Active rotation by angle (rad) about a unit axis (Rodrigues)
    static Mat3 axisAngle(const Vec3& axisIn, double angle)
    {
        Vec3 a = axisIn.normalised();
        double c = std::cos(angle), s = std::sin(angle), C = 1.0 - c;
        Mat3 R;
        R.m[0][0] = c + a.x*a.x*C;      R.m[0][1] = a.x*a.y*C - a.z*s;  R.m[0][2] = a.x*a.z*C + a.y*s;
        R.m[1][0] = a.y*a.x*C + a.z*s;  R.m[1][1] = c + a.y*a.y*C;      R.m[1][2] = a.y*a.z*C - a.x*s;
        R.m[2][0] = a.z*a.x*C - a.y*s;  R.m[2][1] = a.z*a.y*C + a.x*s;  R.m[2][2] = c + a.z*a.z*C;
        return R;
    }
};


// --- Quat
// Unit quaternion (w, x, y, z) of an active rotation (q and -q equal); q1*q2 = apply q2 then q1
struct Quat
{
    double w, x, y, z;

    Quat() : w(1), x(0), y(0), z(0) {}
    Quat(double W, double X, double Y, double Z) : w(W), x(X), y(Y), z(Z) {}

    static Quat axisAngle(const Vec3& axisIn, double angle)
    {
        Vec3 a = axisIn.normalised();
        double s = std::sin(0.5*angle);
        return Quat(std::cos(0.5*angle), a.x*s, a.y*s, a.z*s);
    }

    Quat conj() const { return Quat(w, -x, -y, -z); }

    double norm() const { return std::sqrt(w*w + x*x + y*y + z*z); }

    Quat normalised() const
    {
        double n = norm();
        return Quat(w/n, x/n, y/n, z/n);
    }

    double dot(const Quat& b) const { return w*b.w + x*b.x + y*b.y + z*b.z; }

    Quat operator*(const Quat& b) const
    {
        return Quat
        (
            w*b.w - x*b.x - y*b.y - z*b.z,
            w*b.x + x*b.w + y*b.z - z*b.y,
            w*b.y - x*b.z + y*b.w + z*b.x,
            w*b.z + x*b.y - y*b.x + z*b.w
        );
    }

    Quat operator-() const { return Quat(-w, -x, -y, -z); }

    // rotation angle in [0, pi]
    double angle() const
    {
        double a = std::fabs(w)/norm();
        if (a > 1.0) a = 1.0;
        return 2.0*std::acos(a);
    }

    Mat3 toMatrix() const
    {
        Quat q = normalised();
        Mat3 R;
        R.m[0][0] = 1 - 2*(q.y*q.y + q.z*q.z);
        R.m[0][1] = 2*(q.x*q.y - q.z*q.w);
        R.m[0][2] = 2*(q.x*q.z + q.y*q.w);
        R.m[1][0] = 2*(q.x*q.y + q.z*q.w);
        R.m[1][1] = 1 - 2*(q.x*q.x + q.z*q.z);
        R.m[1][2] = 2*(q.y*q.z - q.x*q.w);
        R.m[2][0] = 2*(q.x*q.z - q.y*q.w);
        R.m[2][1] = 2*(q.y*q.z + q.x*q.w);
        R.m[2][2] = 1 - 2*(q.x*q.x + q.y*q.y);
        return R;
    }

    // Shepperd's method (robust for all rotations)
    static Quat fromMatrix(const Mat3& R)
    {
        double t = R.trace();
        double q[4];
        if (t > R.m[0][0] && t > R.m[1][1] && t > R.m[2][2])
        {
            double s = std::sqrt(1.0 + t)*2.0;
            q[0] = 0.25*s;
            q[1] = (R.m[2][1] - R.m[1][2])/s;
            q[2] = (R.m[0][2] - R.m[2][0])/s;
            q[3] = (R.m[1][0] - R.m[0][1])/s;
        }
        else if (R.m[0][0] > R.m[1][1] && R.m[0][0] > R.m[2][2])
        {
            double s = std::sqrt(1.0 + R.m[0][0] - R.m[1][1] - R.m[2][2])*2.0;
            q[0] = (R.m[2][1] - R.m[1][2])/s;
            q[1] = 0.25*s;
            q[2] = (R.m[0][1] + R.m[1][0])/s;
            q[3] = (R.m[0][2] + R.m[2][0])/s;
        }
        else if (R.m[1][1] > R.m[2][2])
        {
            double s = std::sqrt(1.0 + R.m[1][1] - R.m[0][0] - R.m[2][2])*2.0;
            q[0] = (R.m[0][2] - R.m[2][0])/s;
            q[1] = (R.m[0][1] + R.m[1][0])/s;
            q[2] = 0.25*s;
            q[3] = (R.m[1][2] + R.m[2][1])/s;
        }
        else
        {
            double s = std::sqrt(1.0 + R.m[2][2] - R.m[0][0] - R.m[1][1])*2.0;
            q[0] = (R.m[1][0] - R.m[0][1])/s;
            q[1] = (R.m[0][2] + R.m[2][0])/s;
            q[2] = (R.m[1][2] + R.m[2][1])/s;
            q[3] = 0.25*s;
        }
        Quat out(q[0], q[1], q[2], q[3]);
        return out.normalised();
    }
};


// --- Euler conventions

//- Bunge Euler angles (rad) -> passive orientation matrix g (sample -> crystal)
inline Mat3 bungeToMatrix(double phi1, double Phi, double phi2)
{
    const double c1 = std::cos(phi1), s1 = std::sin(phi1);
    const double cP = std::cos(Phi),  sP = std::sin(Phi);
    const double c2 = std::cos(phi2), s2 = std::sin(phi2);
    Mat3 g;
    g.m[0][0] =  c1*c2 - s1*s2*cP;
    g.m[0][1] =  s1*c2 + c1*s2*cP;
    g.m[0][2] =  s2*sP;
    g.m[1][0] = -c1*s2 - s1*c2*cP;
    g.m[1][1] = -s1*s2 + c1*c2*cP;
    g.m[1][2] =  c2*sP;
    g.m[2][0] =  s1*sP;
    g.m[2][1] = -c1*sP;
    g.m[2][2] =  cP;
    return g;
}

//- Passive matrix g -> Bunge Euler angles (rad), phi1,phi2 in [0,2pi), Phi in [0,pi]
inline void matrixToBunge(const Mat3& g, double& phi1, double& Phi, double& phi2)
{
    double g33 = g.m[2][2];
    if (g33 > 1.0) g33 = 1.0;
    if (g33 < -1.0) g33 = -1.0;
    Phi = std::acos(g33);
    const double sP = std::sin(Phi);
    if (sP > 1e-9)
    {
        phi1 = std::atan2(g.m[2][0], -g.m[2][1]);
        phi2 = std::atan2(g.m[0][2],  g.m[1][2]);
    }
    else
    {
        // Phi = 0 or pi: only phi1 + phi2 (or phi1 - phi2) is defined; put it all in phi1
        phi1 = std::atan2(g.m[0][1], g.m[0][0]);
        phi2 = 0.0;
    }
    if (phi1 < 0) phi1 += 2*PI;
    if (phi2 < 0) phi2 += 2*PI;
}

//- Code angles (ROT.1, ROT.2, ROT.3) -> active matrix R = Rz(p3) Ry(p2) Rx(p1), as StressUpdate.H
inline Mat3 codeAnglesToMatrix(double p1, double p2, double p3)
{
    const double c1 = std::cos(p1), s1 = std::sin(p1);
    const double c2 = std::cos(p2), s2 = std::sin(p2);
    const double c3 = std::cos(p3), s3 = std::sin(p3);
    Mat3 R;
    R.m[0][0] =  c2*c3;
    R.m[0][1] =  s1*s2*c3 - c1*s3;
    R.m[0][2] =  c1*s2*c3 + s1*s3;
    R.m[1][0] =  c2*s3;
    R.m[1][1] =  s1*s2*s3 + c1*c3;
    R.m[1][2] =  c1*s2*s3 - s1*c3;
    R.m[2][0] = -s2;
    R.m[2][1] =  s1*c2;
    R.m[2][2] =  c1*c2;
    return R;
}

//- Active matrix R -> crystaldynamicsFoam angles (identical to rotationKernel in StressUpdate.H)
inline void matrixToCodeAngles(const Mat3& R, double& p1, double& p2, double& p3)
{
    if (std::fabs(R.m[2][0]) < 1.0 - 1e-12)
    {
        p2 = -std::asin(R.m[2][0]);
        const double ct = std::cos(p2);
        p1 = std::atan2(R.m[2][1]/ct, R.m[2][2]/ct);
        p3 = std::atan2(R.m[1][0]/ct, R.m[0][0]/ct);
    }
    else
    {
        p3 = 0.0;
        if (R.m[2][0] < 0.0)
        {
            p2 = 0.5*PI;
            p1 = std::atan2(R.m[0][1], R.m[0][2]);
        }
        else
        {
            p2 = -0.5*PI;
            p1 = std::atan2(-R.m[0][1], -R.m[0][2]);
        }
    }
}

//- Bunge (rad) -> code angles (rad):  R = g^T
inline void bungeToCodeAngles
(
    double phi1, double Phi, double phi2,
    double& p1, double& p2, double& p3
)
{
    matrixToCodeAngles(bungeToMatrix(phi1, Phi, phi2).T(), p1, p2, p3);
}

//- code angles (rad) -> Bunge (rad):  g = R^T
inline void codeAnglesToBunge
(
    double p1, double p2, double p3,
    double& phi1, double& Phi, double& phi2
)
{
    matrixToBunge(codeAnglesToMatrix(p1, p2, p3).T(), phi1, Phi, phi2);
}


// --- crystal symmetry

enum class Symmetry { none, cubic, hexagonal };

inline Symmetry symmetryFromName(const std::string& n)
{
    if (n == "cubic" || n == "m-3m" || n == "432") return Symmetry::cubic;
    if (n == "hexagonal" || n == "6/mmm" || n == "622") return Symmetry::hexagonal;
    if (n == "none" || n == "triclinic" || n == "1") return Symmetry::none;
    throw std::runtime_error("unknown symmetry '" + n + "' (cubic | hexagonal | none)");
}

inline std::string symmetryName(Symmetry s)
{
    switch (s)
    {
        case Symmetry::cubic: return "cubic";
        case Symmetry::hexagonal: return "hexagonal";
        default: return "none";
    }
}

//- Proper rotations of the Laue group (24 for cubic O, 12 for hexagonal D6, 1 for none)
inline const std::vector<Quat>& symmetryOps(Symmetry sym)
{
    static std::vector<Quat> cubic, hex, none;
    if (none.empty()) none.push_back(Quat());

    if (sym == Symmetry::cubic)
    {
        if (cubic.empty())
        {
            cubic.push_back(Quat());
            const Vec3 ax4[3] = {Vec3(1,0,0), Vec3(0,1,0), Vec3(0,0,1)};
            for (int a = 0; a < 3; a++)
                for (int k = 1; k <= 3; k++)
                    cubic.push_back(Quat::axisAngle(ax4[a], k*0.5*PI));
            const Vec3 ax3[4] = {Vec3(1,1,1), Vec3(1,-1,1), Vec3(-1,1,1), Vec3(1,1,-1)};
            for (int a = 0; a < 4; a++)
                for (int k = 1; k <= 2; k++)
                    cubic.push_back(Quat::axisAngle(ax3[a], k*2.0*PI/3.0));
            const Vec3 ax2[6] =
            {
                Vec3(1,1,0), Vec3(1,-1,0), Vec3(1,0,1),
                Vec3(1,0,-1), Vec3(0,1,1), Vec3(0,1,-1)
            };
            for (int a = 0; a < 6; a++)
                cubic.push_back(Quat::axisAngle(ax2[a], PI));
        }
        return cubic;
    }
    else if (sym == Symmetry::hexagonal)
    {
        if (hex.empty())
        {
            // Cartesian crystal frame: z || c, x || a1 = [2-1-10]
            for (int k = 0; k < 6; k++)
                hex.push_back(Quat::axisAngle(Vec3(0,0,1), k*PI/3.0));
            for (int k = 0; k < 6; k++)
            {
                const double a = k*PI/6.0;   // 2-fold axes in the basal plane every 30 deg
                hex.push_back(Quat::axisAngle(Vec3(std::cos(a), std::sin(a), 0), PI));
            }
        }
        return hex;
    }
    return none;
}

//- Minimum misorientation angle (rad) between two orientations R_a, R_b
//  (crystal->sample quaternions) accounting for crystal symmetry.
inline double misorientationAngle(const Quat& qa, const Quat& qb, Symmetry sym)
{
    const std::vector<Quat>& S = symmetryOps(sym);
    const Quat d = qb.conj()*qa;           // crystal-frame misorientation
    double best = 0.0;
    for (size_t i = 0; i < S.size(); i++)
    {
        const Quat c = d*S[i];             // one-sided suffices (angle is conjugation invariant)
        const double a = std::fabs(c.w);
        if (a > best) best = a;
    }
    if (best > 1.0) best = 1.0;
    return 2.0*std::acos(best);
}

//- The symmetry-equivalent of q (q*S, sign-adjusted) closest to ref
inline Quat nearestEquivalent(const Quat& q, const Quat& ref, Symmetry sym)
{
    const std::vector<Quat>& S = symmetryOps(sym);
    Quat best = q;
    double bestDot = -1.0;
    for (size_t i = 0; i < S.size(); i++)
    {
        Quat c = q*S[i];
        double d = c.dot(ref);
        if (std::fabs(d) > bestDot)
        {
            bestDot = std::fabs(d);
            best = (d < 0) ? -c : c;
        }
    }
    return best;
}

//- Symmetry-aware mean orientation of a set of quaternions (two-pass)
inline Quat meanOrientation(const std::vector<Quat>& qs, Symmetry sym)
{
    if (qs.empty()) throw std::runtime_error("meanOrientation: empty set");
    Quat ref = qs[0];
    for (int pass = 0; pass < 2; pass++)
    {
        double w = 0, x = 0, y = 0, z = 0;
        for (size_t i = 0; i < qs.size(); i++)
        {
            Quat c = nearestEquivalent(qs[i], ref, sym);
            w += c.w; x += c.x; y += c.y; z += c.z;
        }
        Quat m(w, x, y, z);
        if (m.norm() < 1e-12) return ref;
        ref = m.normalised();
    }
    return ref;
}

} // namespace ebsd

#endif
