/*---------------------------------------------------------------------------*\
    ebsdTools :: ebsdMap.H

    EBSD map container and grid algorithms for ebsdToFoam (plain C++11): text reading,
    distance transform, grain reconstruction, graph colouring, boundary distance.
    Pixel (i,j) has index p = i + j*nx and map-frame coordinates (x0 + i*dx, y0 + j*dy).
\*---------------------------------------------------------------------------*/

#ifndef ebsdTools_ebsdMap_H
#define ebsdTools_ebsdMap_H

#include "orientation.H"

#include <fstream>
#include <sstream>
#include <iostream>
#include <map>
#include <set>
#include <deque>
#include <cctype>
#include <cstdlib>
#include <limits>

namespace ebsd
{

// --- helpers
inline std::string lower(std::string s)
{
    for (size_t i = 0; i < s.size(); i++) s[i] = static_cast<char>(std::tolower(s[i]));
    return s;
}

inline std::vector<std::string> splitWs(const std::string& line)
{
    std::vector<std::string> out;
    std::istringstream iss(line);
    std::string tok;
    while (iss >> tok) out.push_back(tok);
    return out;
}

inline bool isNumber(const std::string& s)
{
    char* end = nullptr;
    std::strtod(s.c_str(), &end);
    return end != s.c_str() && *end == '\0';
}


// --- distance transform

//- 1D squared distance transform (Felzenszwalb & Huttenlocher 2004), spacing h:
//  f in (0/INF or squared distances), d out, arg = index of the winning site
inline void dt1d
(
    const double* f, int n, double h,
    double* d, int* arg,
    std::vector<int>& v, std::vector<double>& z
)
{
    const double INF = 1e300;
    v.resize(n); z.resize(n + 1);
    int k = 0;
    v[0] = 0; z[0] = -INF; z[1] = INF;
    const double h2 = h*h;
    for (int q = 1; q < n; q++)
    {
        double s = ((f[q] + h2*q*q) - (f[v[k]] + h2*v[k]*v[k]))/(2.0*h2*(q - v[k]));
        while (s <= z[k])
        {
            k--;
            s = ((f[q] + h2*q*q) - (f[v[k]] + h2*v[k]*v[k]))/(2.0*h2*(q - v[k]));
        }
        k++;
        v[k] = q; z[k] = s; z[k+1] = INF;
    }
    k = 0;
    for (int q = 0; q < n; q++)
    {
        while (z[k+1] < q) k++;
        d[q] = h2*(q - v[k])*(q - v[k]) + f[v[k]];
        arg[q] = v[k];
    }
}

//- 2D exact Euclidean distance transform on a W x H grid (x fastest): nearest[p] = index
//  of the nearest feature pixel (-1 if none), dist[p] = distance in physical units
inline void edt2D
(
    int W, int H, double hx, double hy,
    const std::vector<char>& feature,
    std::vector<int>& nearest,
    std::vector<double>& dist
)
{
    const double BIG = 1e200;      // finite 'infinity' keeps the parabola arithmetic well defined
    const size_t N = static_cast<size_t>(W)*H;
    std::vector<double> d1(N), f(std::max(W, H)), d(std::max(W, H));
    std::vector<int> a1(N), arg(std::max(W, H)), v;
    std::vector<double> z;

    // pass 1: columns (along y)
    for (int i = 0; i < W; i++)
    {
        for (int j = 0; j < H; j++) f[j] = feature[i + static_cast<size_t>(j)*W] ? 0.0 : BIG;
        dt1d(f.data(), H, hy, d.data(), arg.data(), v, z);
        for (int j = 0; j < H; j++)
        {
            d1[i + static_cast<size_t>(j)*W] = d[j];
            a1[i + static_cast<size_t>(j)*W] = arg[j];    // nearest feature row in this column
        }
    }
    nearest.assign(N, -1);
    dist.assign(N, 0.0);
    // pass 2: rows (along x)
    for (int j = 0; j < H; j++)
    {
        for (int i = 0; i < W; i++) f[i] = d1[i + static_cast<size_t>(j)*W];
        dt1d(f.data(), W, hx, d.data(), arg.data(), v, z);
        for (int i = 0; i < W; i++)
        {
            const size_t p = i + static_cast<size_t>(j)*W;
            const int col = arg[i];
            const int row = a1[col + static_cast<size_t>(j)*W];
            if (d[i] < 0.5*BIG && feature[col + static_cast<size_t>(row)*W])
            {
                nearest[p] = col + row*W;
                dist[p] = std::sqrt(d[i]);
            }
            else
            {
                nearest[p] = -1;
                dist[p] = std::numeric_limits<double>::infinity();
            }
        }
    }
}


// --- EbsdMap
struct EbsdMap
{
    int nx = 0, ny = 0;
    double x0 = 0, y0 = 0, dx = 1, dy = 1;        // file units
    std::vector<double> phi1, Phi, phi2;         // Bunge angles, radians, per pixel
    std::vector<char> indexed;                   // measurement present
    std::vector<Quat> q;                         // crystal->sample (Euler frame), per pixel
    std::vector<int> grain;                      // grain label per pixel
    int nGrains = 0;
    std::size_t nPointsRead = 0;
    std::size_t nUnindexed = 0;

    std::size_t size() const { return static_cast<std::size_t>(nx)*ny; }
    std::size_t index(int i, int j) const { return i + static_cast<std::size_t>(j)*nx; }

    void buildQuaternions()
    {
        q.resize(size());
        for (std::size_t p = 0; p < size(); p++)
        {
            q[p] = Quat::fromMatrix(bungeToMatrix(phi1[p], Phi[p], phi2[p]).T());
        }
    }

    void setBungeFromQuat(std::size_t p)
    {
        matrixToBunge(q[p].toMatrix().T(), phi1[p], Phi[p], phi2[p]);
    }
};


struct ReadOptions
{
    std::string xName = "x";
    std::string yName = "y";
    std::string e1Name = "euler1";
    std::string e2Name = "euler2";
    std::string e3Name = "euler3";
    int xCol = -1, yCol = -1, e1Col = -1, e2Col = -1, e3Col = -1;   // >=0 overrides names
    bool degrees = true;
    bool zeroEulerIsUnindexed = true;
};


//- Read a whitespace/tab separated text export with an optional header line of column names.
inline EbsdMap readEbsdText(const std::string& file, const ReadOptions& opt, std::ostream& log)
{
    std::ifstream in(file.c_str());
    if (!in) throw std::runtime_error("cannot open EBSD file '" + file + "'");

    std::vector<std::string> header;
    std::vector<double> X, Y, E1, E2, E3;
    std::string line;
    int xc = opt.xCol, yc = opt.yCol, c1 = opt.e1Col, c2 = opt.e2Col, c3 = opt.e3Col;
    bool headerDone = false;
    std::size_t lineNo = 0;

    while (std::getline(in, line))
    {
        lineNo++;
        // strip CR and comments
        if (!line.empty() && line[line.size()-1] == '\r') line.erase(line.size()-1);
        if (line.empty()) continue;
        if (line[0] == '#' || line[0] == '%') continue;
        std::vector<std::string> tok = splitWs(line);
        if (tok.empty()) continue;

        if (!headerDone)
        {
            headerDone = true;
            if (!isNumber(tok[0]))
            {
                header = tok;
                if (xc < 0 || yc < 0 || c1 < 0 || c2 < 0 || c3 < 0)
                {
                    std::map<std::string, int> col;
                    for (size_t k = 0; k < header.size(); k++) col[lower(header[k])] = static_cast<int>(k);
                    auto find = [&](const std::string& nm, const char* alt1, const char* alt2) -> int
                    {
                        if (col.count(lower(nm))) return col[lower(nm)];
                        if (alt1 && col.count(alt1)) return col[alt1];
                        if (alt2 && col.count(alt2)) return col[alt2];
                        return -1;
                    };
                    if (xc < 0) xc = find(opt.xName, "x", nullptr);
                    if (yc < 0) yc = find(opt.yName, "y", nullptr);
                    if (c1 < 0) c1 = find(opt.e1Name, "euler1", "phi1");
                    if (c2 < 0) c2 = find(opt.e2Name, "euler2", "phi");
                    if (c3 < 0) c3 = find(opt.e3Name, "euler3", "phi2");
                    if (xc < 0 || yc < 0 || c1 < 0 || c2 < 0 || c3 < 0)
                    {
                        throw std::runtime_error
                        (
                            "could not identify the X, Y, Euler1, Euler2, Euler3 columns in the header '"
                          + line + "' of " + file
                        );
                    }
                }
                continue;
            }
            else
            {
                // no header: assume X Y E1 E2 E3 unless given
                if (xc < 0) xc = 0;
                if (yc < 0) yc = 1;
                if (c1 < 0) c1 = 2;
                if (c2 < 0) c2 = 3;
                if (c3 < 0) c3 = 4;
            }
        }

        const int need = std::max(std::max(xc, yc), std::max(c1, std::max(c2, c3)));
        if (static_cast<int>(tok.size()) <= need)
        {
            throw std::runtime_error("too few columns on line " + std::to_string(lineNo) + " of " + file);
        }
        X.push_back(std::strtod(tok[xc].c_str(), nullptr));
        Y.push_back(std::strtod(tok[yc].c_str(), nullptr));
        E1.push_back(std::strtod(tok[c1].c_str(), nullptr));
        E2.push_back(std::strtod(tok[c2].c_str(), nullptr));
        E3.push_back(std::strtod(tok[c3].c_str(), nullptr));
    }

    if (X.empty()) throw std::runtime_error("no data points read from " + file);

    // --- grid detection from the sorted unique coordinates
    auto uniqueSorted = [](std::vector<double> v, double& step)
    {
        std::sort(v.begin(), v.end());
        const double range = v.back() - v.front();
        const double tol = std::max(1e-9*std::max(range, 1.0), 1e-12);
        std::vector<double> u;
        for (size_t k = 0; k < v.size(); k++)
            if (u.empty() || v[k] - u.back() > tol) u.push_back(v[k]);
        step = 0.0;
        for (size_t k = 1; k < u.size(); k++)
        {
            const double s = u[k] - u[k-1];
            if (step == 0.0 || s < step) step = s;
        }
        if (step == 0.0) step = 1.0;
        return u;
    };

    EbsdMap m;
    std::vector<double> ux = uniqueSorted(X, m.dx);
    std::vector<double> uy = uniqueSorted(Y, m.dy);
    m.x0 = ux.front(); m.y0 = uy.front();
    m.nx = static_cast<int>(std::floor((ux.back() - ux.front())/m.dx + 0.5)) + 1;
    m.ny = static_cast<int>(std::floor((uy.back() - uy.front())/m.dy + 0.5)) + 1;
    if (static_cast<double>(m.nx)*m.ny > 4.0*X.size() + 16)
    {
        log << "WARNING: the detected grid (" << m.nx << " x " << m.ny
            << ") is much larger than the number of points (" << X.size()
            << "); the coordinates may not be on a regular grid." << std::endl;
    }

    const std::size_t N = m.size();
    m.phi1.assign(N, 0.0); m.Phi.assign(N, 0.0); m.phi2.assign(N, 0.0);
    m.indexed.assign(N, 0);
    m.nPointsRead = X.size();
    const double aScale = opt.degrees ? PI/180.0 : 1.0;

    std::size_t duplicates = 0;
    for (size_t k = 0; k < X.size(); k++)
    {
        const int i = static_cast<int>(std::floor((X[k] - m.x0)/m.dx + 0.5));
        const int j = static_cast<int>(std::floor((Y[k] - m.y0)/m.dy + 0.5));
        if (i < 0 || i >= m.nx || j < 0 || j >= m.ny) continue;
        const std::size_t p = m.index(i, j);
        if (opt.zeroEulerIsUnindexed && E1[k] == 0.0 && E2[k] == 0.0 && E3[k] == 0.0) continue;
        if (m.indexed[p]) duplicates++;
        m.phi1[p] = E1[k]*aScale;
        m.Phi[p]  = E2[k]*aScale;
        m.phi2[p] = E3[k]*aScale;
        m.indexed[p] = 1;
    }
    if (duplicates)
    {
        log << "WARNING: " << duplicates << " points mapped onto already occupied pixels" << std::endl;
    }
    m.nUnindexed = 0;
    for (std::size_t p = 0; p < N; p++) if (!m.indexed[p]) m.nUnindexed++;

    log << "EBSD file " << file << ": " << X.size() << " points, grid " << m.nx << " x " << m.ny
        << ", step (" << m.dx << ", " << m.dy << "), origin (" << m.x0 << ", " << m.y0 << "), "
        << m.nUnindexed << " empty/non-indexed pixels" << std::endl;
    if (!header.empty())
    {
        log << "  columns: X='" << header[xc] << "' Y='" << header[yc] << "' Euler=('"
            << header[c1] << "','" << header[c2] << "','" << header[c3] << "')" << std::endl;
    }
    m.buildQuaternions();
    return m;
}


//- Fill non-indexed pixels with the orientation of the nearest indexed pixel
inline void fillUnindexed(EbsdMap& m, std::ostream& log)
{
    if (m.nUnindexed == 0) return;
    std::vector<int> nearest; std::vector<double> dist;
    edt2D(m.nx, m.ny, m.dx, m.dy, m.indexed, nearest, dist);
    std::size_t filled = 0;
    for (std::size_t p = 0; p < m.size(); p++)
    {
        if (!m.indexed[p])
        {
            const int s = nearest[p];
            if (s < 0) throw std::runtime_error("fillUnindexed: no indexed pixels at all");
            m.phi1[p] = m.phi1[s]; m.Phi[p] = m.Phi[s]; m.phi2[p] = m.phi2[s];
            m.q[p] = m.q[s];
            filled++;
        }
    }
    log << "  filled " << filled << " non-indexed pixels by nearest neighbour" << std::endl;
}


//- Grain reconstruction: 4-connected flood fill, neighbours belong to the same grain when their
//  misorientation is below tol (rad).  Labels are then renumbered largest-first.
inline void segmentGrains(EbsdMap& m, double tol, Symmetry sym)
{
    const int nx = m.nx, ny = m.ny;
    const std::size_t N = m.size();
    std::vector<char> okx(N, 0), oky(N, 0);
    for (int j = 0; j < ny; j++)
        for (int i = 0; i < nx; i++)
        {
            const std::size_t p = m.index(i, j);
            if (i + 1 < nx) okx[p] = misorientationAngle(m.q[p], m.q[p+1], sym) < tol;
            if (j + 1 < ny) oky[p] = misorientationAngle(m.q[p], m.q[p+nx], sym) < tol;
        }

    m.grain.assign(N, -1);
    int ng = 0;
    std::deque<std::size_t> dq;
    for (std::size_t start = 0; start < N; start++)
    {
        if (m.grain[start] >= 0) continue;
        m.grain[start] = ng;
        dq.push_back(start);
        while (!dq.empty())
        {
            const std::size_t p = dq.front(); dq.pop_front();
            const int i = static_cast<int>(p % nx), j = static_cast<int>(p / nx);
            if (i + 1 < nx && okx[p] && m.grain[p+1] < 0)      { m.grain[p+1] = ng; dq.push_back(p+1); }
            if (i > 0 && okx[p-1] && m.grain[p-1] < 0)         { m.grain[p-1] = ng; dq.push_back(p-1); }
            if (j + 1 < ny && oky[p] && m.grain[p+nx] < 0)     { m.grain[p+nx] = ng; dq.push_back(p+nx); }
            if (j > 0 && oky[p-nx] && m.grain[p-nx] < 0)       { m.grain[p-nx] = ng; dq.push_back(p-nx); }
        }
        ng++;
    }
    m.nGrains = ng;
}

inline std::vector<int> grainSizes(const EbsdMap& m)
{
    std::vector<int> s(m.nGrains, 0);
    for (std::size_t p = 0; p < m.size(); p++) s[m.grain[p]]++;
    return s;
}

//- Renumber grains contiguously, largest first
inline void renumberGrains(EbsdMap& m)
{
    std::vector<int> sizes(m.nGrains, 0);
    for (std::size_t p = 0; p < m.size(); p++) if (m.grain[p] >= 0) sizes[m.grain[p]]++;
    std::vector<int> order;
    for (int g = 0; g < m.nGrains; g++) if (sizes[g] > 0) order.push_back(g);
    std::stable_sort(order.begin(), order.end(), [&](int a, int b) { return sizes[a] > sizes[b]; });
    std::vector<int> remap(m.nGrains, -1);
    for (size_t k = 0; k < order.size(); k++) remap[order[k]] = static_cast<int>(k);
    for (std::size_t p = 0; p < m.size(); p++) m.grain[p] = remap[m.grain[p]];
    m.nGrains = static_cast<int>(order.size());
}

//- Merge grains smaller than minPixels into the neighbouring grain sharing the longest boundary
inline void mergeSmallGrains(EbsdMap& m, int minPixels, std::ostream& log)
{
    if (minPixels <= 1) { renumberGrains(m); return; }
    const int nx = m.nx, ny = m.ny;
    int merged = 0;
    for (int iter = 0; iter < 1000; iter++)
    {
        std::vector<int> sizes(m.nGrains, 0);
        for (std::size_t p = 0; p < m.size(); p++) sizes[m.grain[p]]++;
        // boundary-length counts between small grains and their neighbours
        std::vector<std::map<int,int>> nb(m.nGrains);
        bool anySmall = false;
        for (int g = 0; g < m.nGrains; g++) if (sizes[g] > 0 && sizes[g] < minPixels) anySmall = true;
        if (!anySmall) break;
        for (int j = 0; j < ny; j++)
            for (int i = 0; i < nx; i++)
            {
                const std::size_t p = m.index(i, j);
                const int g = m.grain[p];
                const int di[4] = {1, -1, 0, 0}, dj[4] = {0, 0, 1, -1};
                for (int k = 0; k < 4; k++)
                {
                    const int i2 = i + di[k], j2 = j + dj[k];
                    if (i2 < 0 || i2 >= nx || j2 < 0 || j2 >= ny) continue;
                    const int g2 = m.grain[m.index(i2, j2)];
                    if (g2 != g && sizes[g] < minPixels) nb[g][g2]++;
                }
            }
        std::vector<int> target(m.nGrains, -1);
        bool changed = false;
        for (int g = 0; g < m.nGrains; g++)
        {
            if (sizes[g] == 0 || sizes[g] >= minPixels || nb[g].empty()) continue;
            int best = -1, bestCount = -1;
            for (std::map<int,int>::const_iterator it = nb[g].begin(); it != nb[g].end(); ++it)
            {
                // prefer merging into a grain that is itself not small, otherwise the longest boundary
                const int score = it->second + ((sizes[it->first] >= minPixels) ? 1000000 : 0);
                if (score > bestCount) { bestCount = score; best = it->first; }
            }
            target[g] = best;
            if (best >= 0) { changed = true; merged++; }
        }
        if (!changed) break;
        // resolve chains (a->b->c)
        for (int g = 0; g < m.nGrains; g++)
        {
            int t = target[g], guard = 0;
            while (t >= 0 && target[t] >= 0 && target[t] != t && guard++ < 100) t = target[t];
            if (t >= 0) target[g] = t;
        }
        for (std::size_t p = 0; p < m.size(); p++)
        {
            const int t = target[m.grain[p]];
            if (t >= 0) m.grain[p] = t;
        }
    }
    renumberGrains(m);
    log << "  merged " << merged << " grains smaller than " << minPixels << " pixels" << std::endl;
}

//- Symmetry-aware mean orientation of every grain
inline std::vector<Quat> grainMeanOrientations(const EbsdMap& m, Symmetry sym)
{
    std::vector<std::vector<Quat>> lists(m.nGrains);
    for (std::size_t p = 0; p < m.size(); p++) lists[m.grain[p]].push_back(m.q[p]);
    std::vector<Quat> means(m.nGrains);
    for (int g = 0; g < m.nGrains; g++) means[g] = meanOrientation(lists[g], sym);
    return means;
}

//- Replace pixel orientations that deviate from their grain mean by more than tol (rad)
inline std::size_t cleanOutliers(EbsdMap& m, const std::vector<Quat>& means, double tol, Symmetry sym)
{
    std::size_t n = 0;
    for (std::size_t p = 0; p < m.size(); p++)
    {
        if (misorientationAngle(m.q[p], means[m.grain[p]], sym) > tol)
        {
            m.q[p] = means[m.grain[p]];
            m.setBungeFromQuat(p);
            n++;
        }
    }
    return n;
}

//- Grain adjacency: grains closer than 'separation' pixels are adjacent; one neighbour set per grain
inline std::vector<std::set<int>> grainAdjacency(const EbsdMap& m, double separation)
{
    const int nx = m.nx, ny = m.ny;
    std::vector<std::set<int>> adj(m.nGrains);
    const int r = std::max(1, static_cast<int>(std::ceil(separation)));
    const double r2 = std::max(1.0, separation*separation);
    for (int j = 0; j < ny; j++)
        for (int i = 0; i < nx; i++)
        {
            const int g = m.grain[m.index(i, j)];
            // only boundary pixels need scanning
            bool boundary = false;
            if (i + 1 < nx && m.grain[m.index(i+1, j)] != g) boundary = true;
            if (i > 0 && m.grain[m.index(i-1, j)] != g) boundary = true;
            if (j + 1 < ny && m.grain[m.index(i, j+1)] != g) boundary = true;
            if (j > 0 && m.grain[m.index(i, j-1)] != g) boundary = true;
            if (!boundary) continue;
            for (int dj = -r; dj <= r; dj++)
                for (int di = -r; di <= r; di++)
                {
                    if (di*di + dj*dj > r2) continue;
                    const int i2 = i + di, j2 = j + dj;
                    if (i2 < 0 || i2 >= nx || j2 < 0 || j2 >= ny) continue;
                    const int g2 = m.grain[m.index(i2, j2)];
                    if (g2 != g) { adj[g].insert(g2); adj[g2].insert(g); }
                }
        }
    return adj;
}

//- Greedy (Welsh-Powell) colouring: grains ordered by size, each takes the smallest colour not
//  used by an adjacent grain.  Returns colour per grain and the number of colours.
inline std::vector<int> graphColouring(const EbsdMap& m, double separation, int& nColours)
{
    std::vector<std::set<int>> adj = grainAdjacency(m, separation);
    std::vector<int> sizes = grainSizes(m);
    std::vector<int> order(m.nGrains);
    for (int g = 0; g < m.nGrains; g++) order[g] = g;
    std::stable_sort(order.begin(), order.end(), [&](int a, int b) { return sizes[a] > sizes[b]; });
    std::vector<int> colour(m.nGrains, -1);
    nColours = 0;
    for (size_t k = 0; k < order.size(); k++)
    {
        const int g = order[k];
        std::set<int> used;
        for (std::set<int>::const_iterator it = adj[g].begin(); it != adj[g].end(); ++it)
            if (colour[*it] >= 0) used.insert(colour[*it]);
        int c = 0;
        while (used.count(c)) c++;
        colour[g] = c;
        if (c + 1 > nColours) nColours = c + 1;
    }
    return colour;
}


//- Grain-boundary network on a 2x refined grid (feature nodes halfway between 4-adjacent
//  pixels of different grains): distance to the nearest boundary and the grain across it
struct BoundaryDistance
{
    int W = 0, H = 0;                      // fine grid (2nx-1) x (2ny-1)
    double hx = 1, hy = 1;                 // fine spacings (file units)
    std::vector<int> nearest;              // nearest feature node per fine node
    std::vector<std::pair<int,int>> pairs; // grain pair per fine node (only meaningful on features)
    std::vector<double> fx, fy;            // fine-node coordinates in pixel units
    bool empty = true;

    void build(const EbsdMap& m)
    {
        W = 2*m.nx - 1; H = 2*m.ny - 1;
        hx = 0.5*m.dx; hy = 0.5*m.dy;
        const std::size_t N = static_cast<std::size_t>(W)*H;
        std::vector<char> feature(N, 0);
        pairs.assign(N, std::make_pair(-1, -1));
        std::size_t nf = 0;
        for (int j = 0; j < m.ny; j++)
            for (int i = 0; i < m.nx; i++)
            {
                const int g = m.grain[m.index(i, j)];
                if (i + 1 < m.nx)
                {
                    const int g2 = m.grain[m.index(i+1, j)];
                    if (g2 != g)
                    {
                        const std::size_t f = (2*i + 1) + static_cast<std::size_t>(2*j)*W;
                        feature[f] = 1; pairs[f] = std::make_pair(g, g2); nf++;
                    }
                }
                if (j + 1 < m.ny)
                {
                    const int g2 = m.grain[m.index(i, j+1)];
                    if (g2 != g)
                    {
                        const std::size_t f = (2*i) + static_cast<std::size_t>(2*j + 1)*W;
                        feature[f] = 1; pairs[f] = std::make_pair(g, g2); nf++;
                    }
                }
            }
        empty = (nf == 0);
        if (empty) return;
        std::vector<double> dist;
        edt2D(W, H, hx, hy, feature, nearest, dist);
    }

    //- Signed distance (file units, >0 inside the grain 'ownGrain') of a point at continuous
    //  pixel coordinates (u, v) and the grain across the nearest boundary.
    double query(const EbsdMap& m, double u, double v, int ownGrain, int& otherGrain) const
    {
        if (empty) { otherGrain = ownGrain; return 1e300; }
        int fi = static_cast<int>(std::floor(2.0*u + 0.5));
        int fj = static_cast<int>(std::floor(2.0*v + 0.5));
        fi = std::max(0, std::min(W - 1, fi));
        fj = std::max(0, std::min(H - 1, fj));
        const int f = nearest[fi + static_cast<std::size_t>(fj)*W];
        if (f < 0) { otherGrain = ownGrain; return 1e300; }
        const double px = (f % W)*0.5, py = (f / W)*0.5;      // feature position in pixel units
        const double ddx = (u - px)*m.dx, ddy = (v - py)*m.dy;
        const double d = std::sqrt(ddx*ddx + ddy*ddy);
        const std::pair<int,int>& pr = pairs[f];
        otherGrain = (pr.first == ownGrain) ? pr.second : pr.first;
        return d;
    }
};

} // namespace ebsd

#endif
