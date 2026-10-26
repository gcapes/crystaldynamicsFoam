/*---------------------------------------------------------------------------*\
    ebsdTools :: ipfColour.H

    Inverse-pole-figure colouring (plain C++11): d_c = R^T d_s is reduced to the fundamental
    sector of the Laue group and mapped to RGB with the hsv (white-centre) or tsl key.
\*---------------------------------------------------------------------------*/

#ifndef ebsdTools_ipfColour_H
#define ebsdTools_ipfColour_H

#include "orientation.H"

namespace ebsd
{

struct RGB
{
    double r, g, b;
    RGB() : r(0), g(0), b(0) {}
    RGB(double R, double G, double B) : r(R), g(G), b(B) {}
};

enum class IpfKey { hsv, tsl };

inline IpfKey ipfKeyFromName(const std::string& n)
{
    if (n == "hsv" || n == "HSV" || n == "mtex") return IpfKey::hsv;
    if (n == "tsl" || n == "TSL" || n == "edax") return IpfKey::tsl;
    throw std::runtime_error("unknown IPF colour key '" + n + "' (hsv | tsl)");
}

//- Reduce a crystal direction to the fundamental sector of the Laue group
inline Vec3 toFundamentalSector(const Vec3& vIn, Symmetry sym)
{
    Vec3 v = vIn.normalised();
    if (sym == Symmetry::cubic)
    {
        double a[3] = {std::fabs(v.x), std::fabs(v.y), std::fabs(v.z)};
        std::sort(a, a + 3);                 // a0 <= a1 <= a2
        return Vec3(a[1], a[0], a[2]);       // z >= x >= y >= 0 : [001]-[101]-[111] triangle
    }
    else if (sym == Symmetry::hexagonal)
    {
        const double z = std::fabs(v.z);
        const double rho = std::sqrt(v.x*v.x + v.y*v.y);
        double th = std::atan2(v.y, v.x);
        th = std::fmod(th, PI/3.0);
        if (th < 0) th += PI/3.0;
        if (th > PI/6.0) th = PI/3.0 - th;   // fold onto 0..30 deg
        return Vec3(rho*std::cos(th), rho*std::sin(th), z);
    }
    return v;
}

//- Description of the fundamental sector as a spherical triangle
struct Sector
{
    Vec3 v[3];       // vertices: red, green, blue
    Vec3 n[3];       // plane normals of the great-circle edges v0-v1, v1-v2, v2-v0
    Vec3 centre;

    static const Sector& get(Symmetry sym)
    {
        static Sector cubic, hex, none;
        static bool init = false;
        if (!init)
        {
            cubic.v[0] = Vec3(0,0,1);
            cubic.v[1] = Vec3(1,0,1).normalised();
            cubic.v[2] = Vec3(1,1,1).normalised();
            cubic.finish();

            hex.v[0] = Vec3(0,0,1);
            hex.v[1] = Vec3(1,0,0);                                        // [2-1-10]
            hex.v[2] = Vec3(std::cos(PI/6.0), std::sin(PI/6.0), 0);        // [10-10]
            hex.finish();

            // 'none': whole sphere has no sector; use an octant as a stand-in
            none.v[0] = Vec3(0,0,1); none.v[1] = Vec3(1,0,0); none.v[2] = Vec3(0,1,0);
            none.finish();
            init = true;
        }
        if (sym == Symmetry::cubic) return cubic;
        if (sym == Symmetry::hexagonal) return hex;
        return none;
    }

    void finish()
    {
        for (int i = 0; i < 3; i++)
        {
            n[i] = v[i].cross(v[(i+1)%3]).normalised();
        }
        centre = (v[0] + v[1] + v[2]).normalised();
    }
};

//- TSL/EDAX style key
inline RGB ipfColourTSL(const Vec3& vSector, Symmetry sym)
{
    const Vec3& v = vSector;
    if (sym == Symmetry::hexagonal)
    {
        const double chi = std::acos(std::max(-1.0, std::min(1.0, v.z)));   // from [0001]
        const double eta = std::atan2(v.y, v.x);                             // 0..30 deg
        const double etaMax = PI/6.0, chiMax = 0.5*PI;
        double r = 1.0 - chi/chiMax;
        double b = (eta/etaMax)*(chi/chiMax);
        double g = (1.0 - eta/etaMax)*(chi/chiMax);
        double mx = std::max(r, std::max(g, b)); if (mx < 1e-12) mx = 1.0;
        return RGB(r/mx, g/mx, b/mx);
    }
    // cubic (and fallback)
    const double chi = std::acos(std::max(-1.0, std::min(1.0, v.z)));       // from [001]
    const double eta = std::atan2(v.y, v.x);                                 // 0..45 deg
    const double etaMax = PI/4.0;
    const double chiMax = std::atan2(1.0, std::cos(eta));                    // to the [101]-[111] edge
    double r = 1.0 - chi/chiMax;
    double b = (eta/etaMax)*(chi/chiMax);
    double g = (1.0 - eta/etaMax)*(chi/chiMax);
    double mx = std::max(r, std::max(g, b)); if (mx < 1e-12) mx = 1.0;
    r = std::max(0.0, std::min(1.0, r/mx));
    g = std::max(0.0, std::min(1.0, g/mx));
    b = std::max(0.0, std::min(1.0, b/mx));
    return RGB(r, g, b);
}

inline RGB hsvToRgb(double h, double s, double v)
{
    h = h - std::floor(h);
    const double h6 = h*6.0;
    const int i = static_cast<int>(std::floor(h6)) % 6;
    const double f = h6 - std::floor(h6);
    const double p = v*(1 - s), q = v*(1 - s*f), t = v*(1 - s*(1 - f));
    switch (i)
    {
        case 0: return RGB(v, t, p);
        case 1: return RGB(q, v, p);
        case 2: return RGB(p, v, t);
        case 3: return RGB(p, q, v);
        case 4: return RGB(t, p, v);
        default: return RGB(v, p, q);
    }
}

//- White-centre HSV key; saturationPower warps the radial saturation (1 = linear)
inline RGB ipfColourHSV(const Vec3& vSector, Symmetry sym, double saturationPower = 1.0)
{
    const Sector& S = Sector::get(sym);
    const Vec3& c = S.centre;

    // tangent basis at the centre, e1 towards the red vertex
    Vec3 e1 = (S.v[0] - c*(S.v[0].dot(c))).normalised();
    Vec3 e2 = c.cross(e1);

    auto azimuth = [&](const Vec3& u) -> double
    {
        double a = std::atan2(u.dot(e2), u.dot(e1));
        if (a < 0) a += 2*PI;
        if (a > 2*PI - 1e-9) a = 0.0;
        return a;
    };

    // orient e2 so the vertices are met in order red -> green -> blue
    if (azimuth(S.v[1]) > azimuth(S.v[2])) e2 = -e2;

    const double rhoV[4] = {0.0, azimuth(S.v[1]), azimuth(S.v[2]), 2*PI};
    const double hueV[4] = {0.0, 1.0/3.0, 2.0/3.0, 1.0};

    const Vec3 v = vSector.normalised();
    const double rho = azimuth(v);

    // piecewise-linear hue
    double hue = 0.0;
    for (int k = 0; k < 3; k++)
    {
        if (rho >= rhoV[k] && rho <= rhoV[k+1])
        {
            const double f = (rhoV[k+1] > rhoV[k]) ? (rho - rhoV[k])/(rhoV[k+1] - rhoV[k]) : 0.0;
            hue = hueV[k] + f*(hueV[k+1] - hueV[k]);
            break;
        }
    }

    // angular distance to centre, normalised by the distance to the sector boundary
    const double thetaV = std::acos(std::max(-1.0, std::min(1.0, v.dot(c))));
    const Vec3 t = e1*std::cos(rho) + e2*std::sin(rho);
    double thetaB = 1e300;
    for (int k = 0; k < 3; k++)
    {
        const double nc = S.n[k].dot(c), nt = S.n[k].dot(t);
        double th = std::atan2(-nc, nt);
        if (th < 0) th += PI;
        if (th > 1e-12 && th < thetaB) thetaB = th;
    }
    double r = (thetaB < 1e299) ? thetaV/thetaB : 0.0;
    r = std::max(0.0, std::min(1.0, r));
    r = std::pow(r, saturationPower);

    return hsvToRgb(hue, r, 1.0);
}

//- Colour for a crystal direction (any frame-reduced or not)
inline RGB ipfColour(const Vec3& crystalDir, Symmetry sym, IpfKey key, double satPower = 1.0)
{
    const Vec3 v = toFundamentalSector(crystalDir, sym);
    return (key == IpfKey::tsl) ? ipfColourTSL(v, sym) : ipfColourHSV(v, sym, satPower);
}

//- Colour for an orientation R (crystal->sample) viewed along sample direction d_s
inline RGB ipfColour
(
    const Mat3& R,
    const Vec3& sampleDir,
    Symmetry sym,
    IpfKey key,
    double satPower = 1.0
)
{
    const Vec3 dc = R.T()*sampleDir;      // = g . d_s
    return ipfColour(dc, sym, key, satPower);
}

} // namespace ebsd

#endif
