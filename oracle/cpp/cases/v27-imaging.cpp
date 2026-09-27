// Group 27 (imaging): AdaptiveSkeletonClimbing3.
//
// Everything here is integer bookkeeping or IEEE-exact arithmetic (+ - * /,
// sqrt in ComputeNormals, floor in the face subdivision, truncating
// conversions in GetGradient), so every case is compared bit for bit.
//
// Upstream is called through __declspec(noinline) wrappers (ORACLE.md, v23:
// the sign of zero of the ComputeNormals sums can reach an output). Large
// outputs are emitted as 32-bit FNV-1a digests computed identically on both
// sides; lists up to kFullLimit elements are emitted in full. Upstream's
// output order is deterministic (the octree recursion of Merge, the box
// order of Tessellate, the ear-clipping order and MakeUnique's
// first-occurrence numbering), so nothing is re-sorted: the verbatim order
// is a stronger comparison than a canonical sort.
#define ORACLE_FAMILY "v27-imaging"
#include "Oracle.h"

#include <Mathematics/AdaptiveSkeletonClimbing3.h>
#include <Mathematics/BSNumber.h>
#include <Mathematics/UIntegerAP32.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

using namespace gte;

namespace
{
    struct Digest
    {
        uint32_t h = 2166136261u;
        void word(uint32_t w)
        {
            for (int i = 0; i < 4; ++i)
            {
                h ^= (w >> (8 * i)) & 0xFFu;
                h *= 16777619u;
            }
        }
        void integer(int64_t v) { word(static_cast<uint32_t>(v)); }
        void real(double x)
        {
            uint64_t bits;
            std::memcpy(&bits, &x, sizeof(bits));
            word(static_cast<uint32_t>(bits));
            word(static_cast<uint32_t>(bits >> 32));
        }
    };

    // Lists of at most this many elements are emitted in full, longer ones
    // as a digest (the count is always emitted first, so both sides decide
    // identically).
    size_t const kFullLimit = 24;

    using ASC3Vertex = std::array<double, 3>;
    using ASC3Triangle = TriangleKey<true>;

    void OutVertices(oracle::Ctx& io, std::vector<ASC3Vertex> const& v)
    {
        io.outInt(v.size());
        if (v.size() <= kFullLimit)
        {
            for (auto const& p : v) { io.outReal(p[0]); io.outReal(p[1]); io.outReal(p[2]); }
        }
        else
        {
            Digest d;
            for (auto const& p : v) { d.real(p[0]); d.real(p[1]); d.real(p[2]); }
            io.outInt(d.h);
        }
    }

    void OutTriangles(oracle::Ctx& io, std::vector<ASC3Triangle> const& t)
    {
        io.outInt(t.size());
        if (t.size() <= kFullLimit)
        {
            for (auto const& k : t) { io.outInt(k.V[0]); io.outInt(k.V[1]); io.outInt(k.V[2]); }
        }
        else
        {
            Digest d;
            for (auto const& k : t) { d.integer(k.V[0]); d.integer(k.V[1]); d.integer(k.V[2]); }
            io.outInt(d.h);
        }
    }

    // The monoboxes, read back through the public PrintBoxes (the port keeps
    // them in the private mBoxes, which the replay reads). Each box is
    // x0, y0, z0, dx, dy, dz.
    template <typename T>
    std::vector<std::array<int32_t, 6>> GetBoxes(AdaptiveSkeletonClimbing3<T, double>& asc)
    {
        std::ostringstream stream;
        asc.PrintBoxes(stream);
        std::string const text = stream.str();
        std::vector<std::array<int32_t, 6>> boxes;
        std::array<int32_t, 6> box{};
        int field = 0;
        size_t pos = text.find('\n');
        while ((pos = text.find('=', pos)) != std::string::npos)
        {
            box[field] = static_cast<int32_t>(std::strtol(text.c_str() + pos + 1, nullptr, 10));
            ++pos;
            if (++field == 6)
            {
                boxes.push_back(box);
                field = 0;
            }
        }
        return boxes;
    }

    void OutBoxes(oracle::Ctx& io, std::vector<std::array<int32_t, 6>> const& boxes)
    {
        io.outInt(boxes.size());
        if (boxes.size() <= kFullLimit)
        {
            for (auto const& b : boxes) { for (auto c : b) { io.outInt(c); } }
        }
        else
        {
            Digest d;
            for (auto const& b : boxes) { for (auto c : b) { d.integer(c); } }
            io.outInt(d.h);
        }
    }

    void OutNormals(oracle::Ctx& io, std::vector<std::array<double, 3>> const& n)
    {
        io.outInt(n.size());
        if (n.size() <= kFullLimit)
        {
            for (auto const& p : n) { io.outReal(p[0]); io.outReal(p[1]); io.outReal(p[2]); }
        }
        else
        {
            Digest d;
            for (auto const& p : n) { d.real(p[0]); d.real(p[1]); d.real(p[2]); }
            io.outInt(d.h);
        }
    }

    // Everything one extraction produces.
    struct ASC3Result
    {
        std::vector<std::array<int32_t, 6>> boxes;
        std::vector<ASC3Vertex> v, vu;
        std::vector<ASC3Triangle> t, tu, to;
        std::vector<std::array<double, 3>> normals;
        // SaddleCheck on (vu, tu): faces checked, faces contradicting the
        // interpolant, ambiguous faces.
        std::array<int32_t, 3> saddle{ 0, 0, 0 };
    };

    void OutResult(oracle::Ctx& io, ASC3Result const& r)
    {
        OutBoxes(io, r.boxes);
        OutVertices(io, r.v);
        OutTriangles(io, r.t);
        OutVertices(io, r.vu);
        OutTriangles(io, r.tu);
        OutTriangles(io, r.to);
        OutNormals(io, r.normals);
        io.outInt(r.saddle[0]);
        io.outInt(r.saddle[1]);
        io.outInt(r.saddle[2]);
    }

    template <typename T>
    __declspec(noinline) void CallExtract(AdaptiveSkeletonClimbing3<T, double>& asc, double level,
        int32_t depth, std::vector<ASC3Vertex>& v, std::vector<ASC3Triangle>& t)
    {
        asc.Extract(level, depth, v, t);
    }

    template <typename T>
    __declspec(noinline) void CallMakeUnique(AdaptiveSkeletonClimbing3<T, double>& asc,
        std::vector<ASC3Vertex>& v, std::vector<ASC3Triangle>& t)
    {
        asc.MakeUnique(v, t);
    }

    template <typename T>
    __declspec(noinline) void CallOrient(AdaptiveSkeletonClimbing3<T, double>& asc,
        std::vector<ASC3Vertex>& v, std::vector<ASC3Triangle>& t, bool sameDir)
    {
        asc.OrientTriangles(v, t, sameDir);
    }

    template <typename T>
    __declspec(noinline) void CallNormals(AdaptiveSkeletonClimbing3<T, double>& asc,
        std::vector<ASC3Vertex> const& v, std::vector<ASC3Triangle> const& t,
        std::vector<std::array<double, 3>>& n)
    {
        asc.ComputeNormals(v, t, n);
    }
}

// ---- the face-saddle predicate (#544 in 3-D) ----

namespace
{
    using ASCNumber = BSNumber<UIntegerAP32>;

    // A unit face whose four edges are all crossed. Upstream's six
    // Get{X,Y,Z}{Min,Max}EdgesS read the corners as f00, f10 (first face
    // coordinate + 1), f01 (second + 1), f11 and pair the four crossings by
    // det = f00*f11 - f01*f10: det > 0 joins the two crossings next to
    // corner 00 and the two next to corner 11 (cutting those corners off),
    // det < 0 cuts off 10 and 01, det = 0 adds the plus-sign branch point.
    // The trilinear interpolant restricted to the face is bilinear, and its
    // topology is decided by dg = (f00-L)(f11-L) - (f01-L)(f10-L)
    // = det - L*(f00+f11-f01-f10): dg > 0 connects 00 and 11 through the
    // face (cut off 10 and 01), dg < 0 cuts off 00 and 11, dg = 0 is the
    // plus sign. Upstream is right iff sign(det) = -sign(dg). A corner is
    // "above" when static_cast<Real>(f) > level, as in SetLevel (a corner on
    // the level counts as below). Returns -1 when the face is not a
    // four-crossing face, 0 when upstream's pairing is right, 1 when not.
    // dg is evaluated exactly.
    int FaceDefect(int64_t f00, int64_t f10, int64_t f01, int64_t f11, double level)
    {
        bool a00 = static_cast<double>(f00) > level, a10 = static_cast<double>(f10) > level;
        bool a01 = static_cast<double>(f01) > level, a11 = static_cast<double>(f11) > level;
        if (!(a00 == a11 && a10 == a01 && a00 != a10))
        {
            return -1;
        }
        int64_t det = f00 * f11 - f01 * f10;
        ASCNumber dg = ASCNumber(det) - ASCNumber(level) * ASCNumber(f00 + f11 - f01 - f10);
        int sdg = dg.GetSign();
        int sdet = (det > 0 ? 1 : (det < 0 ? -1 : 0));
        return (sdet == -sdg ? 0 : 1);
    }

    // [0] four-crossing faces where upstream is right, [1] where it is not,
    // over every unit face of the grid.
    std::array<int, 2> CountFaces(std::vector<int64_t> const& px, int size, double level)
    {
        auto P = [&px, size](int x, int y, int z) { return px[x + size * (y + size * z)]; };
        std::array<int, 2> count{ 0, 0 };
        auto tally = [&count](int d) { if (d >= 0) { ++count[d]; } };
        for (int a = 0; a < size; ++a)
        {
            for (int j = 0; j + 1 < size; ++j)
            {
                for (int i = 0; i + 1 < size; ++i)
                {
                    // x = a face, (u, v) = (y, z) = (i, j)
                    tally(FaceDefect(P(a, i, j), P(a, i + 1, j), P(a, i, j + 1), P(a, i + 1, j + 1), level));
                    // y = a face, (u, v) = (x, z)
                    tally(FaceDefect(P(i, a, j), P(i + 1, a, j), P(i, a, j + 1), P(i + 1, a, j + 1), level));
                    // z = a face, (u, v) = (x, y)
                    tally(FaceDefect(P(i, j, a), P(i + 1, j, a), P(i, j + 1, a), P(i + 1, j + 1, a), level));
                }
            }
        }
        return count;
    }

    // An unrecorded image in [lo, hi], by mode:
    // 0 small random integers (ties, saddles, isolated voxels),
    // 1 a linear ramp (strongly monotone: merges completely, #194),
    // 2 a paraboloid bowl or dome (a sphere level set; a separable sum, so
    //   no four-crossing face),
    // 3 a constant,
    // 4 a signed checkerboard (saddle-rich),
    // 5 full-range values,
    // 6 a constant background with one or two isolated voxels,
    // 7 a quadric with cross terms (hyperboloids: saddles and tunnels),
    // 8 two spheres (the smaller of two bowls: two components or a neck).
    std::vector<int64_t> RawImage(oracle::Ctx& io, int size, int64_t lo, int64_t hi, int mode)
    {
        int64_t mid = (lo < 0 ? 0 : (lo + hi) / 2);
        int64_t span = hi - mid;
        int R = (io.rawInteger(0, 2) == 0 ? 2 : (io.rawInteger(0, 1) == 0 ? 5 : 100));
        std::array<int64_t, 12> c{};
        for (auto& e : c) { e = io.rawInteger(-3, 3); }
        int64_t cx = io.rawInteger(0, size - 1), cy = io.rawInteger(0, size - 1), cz = io.rawInteger(0, size - 1);
        int64_t dx = io.rawInteger(0, size - 1), dy = io.rawInteger(0, size - 1), dz = io.rawInteger(0, size - 1);
        int64_t sgn = (io.rawInteger(0, 1) == 0 ? 1 : -1);
        int64_t k0 = io.rawInteger(-5, 5), k1 = io.rawInteger(1, 6);
        std::vector<int64_t> px(static_cast<size_t>(size) * size * size);
        for (int z = 0; z < size; ++z)
        {
            for (int y = 0; y < size; ++y)
            {
                for (int x = 0; x < size; ++x)
                {
                    int64_t u = x - cx, v = y - cy, w = z - cz;
                    int64_t value;
                    switch (mode)
                    {
                    case 0: value = io.rawInteger(-R, R); break;
                    case 1: value = c[0] * x + c[1] * y + c[2] * z + k0; break;
                    case 2: value = sgn * (u * u + v * v + w * w) + k0; break;
                    case 3: value = k0; break;
                    case 4: value = ((x + y + z) % 2 == 0 ? 1 : -1) * io.rawInteger(1, R) + io.rawInteger(-1, 1); break;
                    case 5: value = static_cast<int64_t>(io.raw(-1.0, 1.0) * static_cast<double>(span)); break;
                    case 6:
                        value = k0 + ((x == cx && y == cy && z == cz) ? sgn * k1 : 0)
                            + ((x == dx && y == dy && z == dz) ? -sgn * k1 : 0);
                        break;
                    case 7:
                        value = c[0] * u * u + c[1] * v * v + c[2] * w * w + c[3] * u * v
                            + c[4] * v * w + c[5] * w * u + c[6] * x + c[7] * y + c[8] * z + k0;
                        break;
                    default:
                    {
                        int64_t s0 = u * u + v * v + w * w - k1 * k1;
                        int64_t s1 = (x - dx) * (x - dx) + (y - dy) * (y - dy) + (z - dz) * (z - dz) - k1;
                        value = sgn * std::min(s0, s1);
                        break;
                    }
                    }
                    px[x + size * (y + size * z)] = std::min(std::max(mid + value, lo), hi);
                }
            }
        }
        return px;
    }

    // An unrecorded level: half-integers (the documented use), reals, the
    // integers of the voxel range (a voxel on the level: upstream does not
    // assert in 3-D, SetLevel counts the voxel as below and the Get*Interp
    // return an edge endpoint), 0 and -0, values outside the range (nothing
    // extracted), quarter-integers, NaN and infinities (every comparison
    // false, or one-sided).
    double RawLevel(oracle::Ctx& io, std::vector<int64_t> const& px)
    {
        auto mm = std::minmax_element(px.begin(), px.end());
        double lo = static_cast<double>(*mm.first), hi = static_cast<double>(*mm.second);
        switch (io.rawInteger(0, 11))
        {
        case 0: case 1: case 2:
            return std::floor(io.raw(lo - 1.0, hi + 1.0)) + 0.5;
        case 3: case 4:
            return io.raw(lo - 0.5, hi + 0.5);
        case 5: case 6:
            return std::floor(io.raw(lo, hi + 1.0));
        case 7:
            return (io.rawInteger(0, 1) == 0 ? 0.0 : -0.0);
        case 8:
            return (io.rawInteger(0, 1) == 0 ? hi + io.raw(0.25, 3.0) : lo - io.raw(0.25, 3.0));
        case 9:
            return std::floor(4.0 * io.raw(lo - 1.0, hi + 1.0)) / 4.0;
        case 10:
            return 0.5;
        default:
        {
            int k = io.rawInteger(0, 2);
            return (k == 0 ? std::nan("") : (k == 1 ? std::numeric_limits<double>::infinity()
                : -std::numeric_limits<double>::infinity()));
        }
        }
    }
}

// ---- Extract, MakeUnique, OrientTriangles, ComputeNormals ----

namespace
{
    // The independent face-saddle check (the #544 check in 3-D), run on
    // upstream's own MakeUnique output and emitted, so that the deviation
    // case shows on the real build that upstream's mesh contradicts the
    // interpolant (the replay runs the same check on the port's mesh). For
    // every four-crossing unit face (finite level, no corner on the level)
    // the crossings E_v0 (edge 00-10), E_v1 (01-11), E_u0 (00-01), E_u1
    // (10-11) are located in the mesh by their coordinates, computed as
    // upstream's Get*Interp does. dg > 0 requires the mesh edges
    // E_v0-E_u1 and E_v1-E_u0 (corners 10 and 01 cut off), dg < 0 the edges
    // E_v0-E_u0 and E_v1-E_u1, and both that no vertex inside the face is
    // adjacent to all four (a plus-sign branch point); dg = 0 requires such
    // a vertex. A face with all four cuts is ambiguous: the ear clipping of
    // a box can put a diagonal into a face, and the diagonals complete the
    // other pairing (upstream does this on sound faces too). Returns
    // {faces checked, faces contradicting the interpolant, ambiguous faces}.
    std::array<int32_t, 3> SaddleCheck(std::vector<int64_t> const& px, int size, double level,
        std::vector<ASC3Vertex> const& vu, std::vector<ASC3Triangle> const& tu)
    {
        std::array<int32_t, 3> result{ 0, 0, 0 };
        if (!std::isfinite(level)) { return result; }
        std::map<ASC3Vertex, int32_t> index;
        for (size_t i = 0; i < vu.size(); ++i) { index.emplace(vu[i], static_cast<int32_t>(i)); }
        std::set<std::pair<int32_t, int32_t>> edges;
        std::vector<std::set<int32_t>> adjacent(vu.size());
        for (auto const& t : tu)
        {
            for (int k = 0; k < 3; ++k)
            {
                int32_t a = t.V[k], b = t.V[(k + 1) % 3];
                edges.insert({ std::min(a, b), std::max(a, b) });
                adjacent[a].insert(b);
                adjacent[b].insert(a);
            }
        }
        auto find = [&index](ASC3Vertex const& p) { auto it = index.find(p); return it == index.end() ? -1 : it->second; };
        auto hasEdge = [&edges](int32_t a, int32_t b)
        {
            return a >= 0 && b >= 0 && edges.count({ std::min(a, b), std::max(a, b) }) > 0;
        };
        auto P = [&px, size](int x, int y, int z) { return px[x + size * (y + size * z)]; };
        for (int axis = 0; axis < 3; ++axis)
        {
            for (int a = 0; a < size; ++a)
            {
                for (int j = 0; j + 1 < size; ++j)
                {
                    for (int i = 0; i + 1 < size; ++i)
                    {
                        // Face-local (u, v) = (i, j) along the two other axes,
                        // in upstream's order: (y, z), (x, z), (x, y).
                        int ku = (axis == 0 ? 1 : 0), kv = (axis == 2 ? 1 : 2);
                        auto corner = [&](int du, int dv)
                        {
                            std::array<int, 3> q{};
                            q[axis] = a; q[ku] = i + du; q[kv] = j + dv;
                            return P(q[0], q[1], q[2]);
                        };
                        int64_t f00 = corner(0, 0), f10 = corner(1, 0), f01 = corner(0, 1), f11 = corner(1, 1);
                        if (FaceDefect(f00, f10, f01, f11, level) < 0) { continue; }
                        double const g00 = static_cast<double>(f00), g10 = static_cast<double>(f10);
                        double const g01 = static_cast<double>(f01), g11 = static_cast<double>(f11);
                        if (g00 == level || g10 == level || g01 == level || g11 == level) { continue; }
                        auto crossing = [&](int k, double base, double h0, double h1, int du, int dv)
                        {
                            ASC3Vertex p{};
                            p[axis] = a; p[ku] = i + du; p[kv] = j + dv;
                            p[k] = base + (level - h0) / (h1 - h0);
                            return find(p);
                        };
                        int32_t ev0 = crossing(ku, i, g00, g10, 0, 0);
                        int32_t ev1 = crossing(ku, i, g01, g11, 0, 1);
                        int32_t eu0 = crossing(kv, j, g00, g01, 0, 0);
                        int32_t eu1 = crossing(kv, j, g10, g11, 1, 0);
                        int64_t det = f00 * f11 - f01 * f10;
                        int sdg = (ASCNumber(det) - ASCNumber(level) * ASCNumber(f00 + f11 - f01 - f10)).GetSign();
                        bool branch = false;
                        if (ev0 >= 0 && ev1 >= 0 && eu0 >= 0 && eu1 >= 0)
                        {
                            for (int32_t c : adjacent[ev0])
                            {
                                ASC3Vertex const& q = vu[c];
                                branch = branch || (q[axis] == a && i < q[ku] && q[ku] < i + 1
                                    && j < q[kv] && q[kv] < j + 1 && adjacent[ev1].count(c) > 0
                                    && adjacent[eu0].count(c) > 0 && adjacent[eu1].count(c) > 0);
                            }
                        }
                        bool cut00 = hasEdge(ev0, eu0), cut11 = hasEdge(ev1, eu1);
                        bool cut10 = hasEdge(ev0, eu1), cut01 = hasEdge(ev1, eu0);
                        ++result[0];
                        if (sdg != 0 && !branch && cut00 && cut11 && cut10 && cut01)
                        {
                            ++result[2];
                        }
                        else if (sdg > 0 ? branch || !(cut10 && cut01) :
                            (sdg < 0 ? branch || !(cut00 && cut11) : !branch))
                        {
                            ++result[1];
                        }
                    }
                }
            }
        }
        return result;
    }

    template <typename T>
    __declspec(noinline) std::vector<ASC3Result> RunASC(int32_t N, std::vector<int64_t> const& px,
        bool fixBoundary, std::vector<double> const& levels, std::vector<int32_t> const& depths,
        std::vector<bool> const& sameDirs)
    {
        std::vector<T> voxels(px.size());
        for (size_t i = 0; i < px.size(); ++i) { voxels[i] = static_cast<T>(px[i]); }
        AdaptiveSkeletonClimbing3<T, double> asc(N, voxels.data(), fixBoundary);
        std::vector<ASC3Result> results(levels.size());
        for (size_t k = 0; k < levels.size(); ++k)
        {
            ASC3Result& r = results[k];
            CallExtract(asc, levels[k], depths[k], r.v, r.t);
            r.boxes = GetBoxes(asc);
            r.vu = r.v;
            r.tu = r.t;
            CallMakeUnique(asc, r.vu, r.tu);
            r.to = r.tu;
            CallOrient(asc, r.vu, r.to, sameDirs[k]);
            CallNormals(asc, r.vu, r.to, r.normals);
            r.saddle = SaddleCheck(px, (1 << N) + 1, levels[k], r.vu, r.tu);
        }
        return results;
    }

    // N = 1 (a third of the records) or 2.
    int32_t DrawN(oracle::Ctx& io)
    {
        return static_cast<int32_t>(io.given(io.rawInteger(0, 2) == 0 ? 1 : 2));
    }

    // Voxels in [lo, hi] and 1 or 2 levels on which upstream's face pairing
    // is sound (FaceDefect; the port pairs by the interpolant and the
    // saddlePairing case covers the rest). Everything the acceptance test
    // depends on (image mode, image, levels) is redrawn per attempt; after
    // 64 attempts the image falls back to a paraboloid bowl at level 0.5,
    // which has no four-crossing face.
    void DrawSoundImage(oracle::Ctx& io, int32_t N, int64_t lo, int64_t hi,
        std::vector<int64_t>& px, std::vector<double>& levels)
    {
        int const size = (1 << N) + 1;
        bool sound = false;
        for (int attempt = 0; attempt < 64 && !sound; ++attempt)
        {
            int mode = io.rawInteger(0, 8);
            px = RawImage(io, size, lo, hi, mode);
            sound = true;
            for (auto& level : levels)
            {
                level = RawLevel(io, px);
                sound = sound && CountFaces(px, size, level)[1] == 0;
            }
        }
        if (!sound)
        {
            int64_t mid = (lo < 0 ? 0 : (lo + hi) / 2);
            for (int z = 0; z < size; ++z)
            {
                for (int y = 0; y < size; ++y)
                {
                    for (int x = 0; x < size; ++x)
                    {
                        px[x + size * (y + size * z)] = std::min(mid + (x - 1) * (x - 1)
                            + (y - 1) * (y - 1) + (z - 1) * (z - 1) - 2, hi);
                    }
                }
            }
            for (auto& level : levels) { level = static_cast<double>(mid) + 0.5; }
        }
        for (auto v : px) { io.given(static_cast<double>(v)); }
    }

    // The record layout shared by the extraction cases: N, fixBoundary, the
    // number of extractions, the voxels, then (level, depth, sameDir) per
    // extraction. Depth -2..N+1 (merging is blocked for the top 'depth'
    // octree levels; depth > N drops every mergeable leaf, see the report).
    template <typename T>
    void ASCCase(oracle::Ctx& io, int32_t N, int64_t lo, int64_t hi)
    {
        bool fixBoundary = (io.integer(0, 3) == 0);
        int32_t numExtracts = io.integer(1, 2);
        std::vector<int64_t> px;
        std::vector<double> levels(numExtracts);
        DrawSoundImage(io, N, lo, hi, px, levels);
        std::vector<int32_t> depths(levels.size());
        std::vector<bool> sameDirs(levels.size());
        for (size_t k = 0; k < levels.size(); ++k)
        {
            io.given(levels[k]);
            depths[k] = io.integer(-2, N + 1);
            sameDirs[k] = io.boolean();
        }
        std::vector<ASC3Result> results = RunASC<T>(N, px, fixBoundary, levels, depths, sameDirs);
        for (auto const& r : results) { OutResult(io, r); }
    }
}

// int32_t voxels, N = 1 or 2, one or two extractions on the same object (the
// merge trees are rebuilt per call), every depth, fixBoundary one record in
// four, levels on and between the sample values. Upstream-sound face
// pairings only (see DrawSoundImage).
ORACLE_CASE("AdaptiveSkeletonClimbing3.extract")
{
    int32_t N = DrawN(io);
    ASCCase<int32_t>(io, N, std::numeric_limits<int32_t>::min(), std::numeric_limits<int32_t>::max());
}

// The other voxel types upstream allows. uint32_t values stay below 2^31:
// above that, det = f00*f11 - f01*f10 overflows int64_t (undefined
// behaviour, as in AdaptiveSkeletonClimbing2).
ORACLE_CASE("AdaptiveSkeletonClimbing3.extract.types")
{
    int32_t type = io.integer(0, 4);
    int32_t N = DrawN(io);
    switch (type)
    {
    case 0: ASCCase<int8_t>(io, N, -128, 127); break;
    case 1: ASCCase<uint8_t>(io, N, 0, 255); break;
    case 2: ASCCase<int16_t>(io, N, -32768, 32767); break;
    case 3: ASCCase<uint16_t>(io, N, 0, 65535); break;
    default: ASCCase<uint32_t>(io, N, 0, int64_t(1) << 31); break;
    }
}

namespace
{
    // A saddle-rich image (random small values or a signed checkerboard)
    // with one level, optionally with a constructed face: the plus sign
    // (-p, p; q, -q, det = 0) at level 0 or at a small nonzero level, or a
    // face whose level is det/S (the interpolant's saddle on the level, up
    // to the rounding of the quotient).
    void RawSaddleImage(oracle::Ctx& io, int size, std::vector<int64_t>& px, double& level)
    {
        int mode = io.rawInteger(0, 3);
        int R = io.rawInteger(2, 9);
        px.assign(static_cast<size_t>(size) * size * size, 0);
        for (int z = 0; z < size; ++z)
        {
            for (int y = 0; y < size; ++y)
            {
                for (int x = 0; x < size; ++x)
                {
                    px[x + size * (y + size * z)] = (mode == 0 ? io.rawInteger(-R, R)
                        : ((x + y + z) % 2 == 0 ? 1 : -1) * io.rawInteger(1, R) + io.rawInteger(-2, 2));
                }
            }
        }
        level = (io.rawInteger(0, 4) == 0 ? std::floor(io.raw(-2.0, 3.0))
            : std::floor(io.raw(-2.0, 2.0)) + 0.5);
        // A random unit face: normal axis, plane and cell.
        int axis = io.rawInteger(0, 2);
        int a = io.rawInteger(0, size - 1), i = io.rawInteger(0, size - 2), j = io.rawInteger(0, size - 2);
        auto at = [&px, size, axis, a](int u, int v) -> int64_t&
        {
            int x = (axis == 0 ? a : u), y = (axis == 0 ? u : (axis == 1 ? a : v));
            int z = (axis == 2 ? a : v);
            return px[x + size * (y + size * z)];
        };
        int64_t& f00 = at(i, j);
        int64_t& f10 = at(i + 1, j);
        int64_t& f01 = at(i, j + 1);
        int64_t& f11 = at(i + 1, j + 1);
        if (mode == 2)
        {
            int64_t p = io.rawInteger(1, R), q = io.rawInteger(1, R);
            f00 = -p; f10 = p; f01 = q; f11 = -q;
            level = (io.rawInteger(0, 1) == 0 ? 0.0 : io.raw(-0.9, 0.9));
        }
        else if (mode == 3)
        {
            int64_t det = f00 * f11 - f01 * f10;
            int64_t S = f00 + f11 - f01 - f10;
            if (S != 0) { level = static_cast<double>(det) / static_cast<double>(S); }
        }
    }

    // wantDefect = false: at least one four-crossing face and upstream sound
    // on all of them. wantDefect = true: at least one face where upstream's
    // pairing contradicts the interpolant. Capped at 256 attempts; the
    // fallbacks, and record 0 of both cases (a pinned reproduction), are a
    // fixed image (every voxel -1 except 4 on the diagonal x = y of the
    // z = 0 layer, faces 4, -1; -1, 4 with det = 15, S = 10, saddle value
    // 1.5), sound at level 2.5 and defective at level 0.5.
    void DrawSaddleImage(oracle::Ctx& io, bool wantDefect, int32_t& N,
        std::vector<int64_t>& px, double& level)
    {
        N = static_cast<int32_t>(io.given(io.index() % 4 == 3 ? 2 : 1));
        int const size = (1 << N) + 1;
        bool ok = false;
        for (int attempt = 0; attempt < 256 && !ok && io.index() != 0; ++attempt)
        {
            RawSaddleImage(io, size, px, level);
            auto count = CountFaces(px, size, level);
            ok = (wantDefect ? count[1] > 0 : count[0] > 0 && count[1] == 0);
        }
        if (!ok)
        {
            px.assign(static_cast<size_t>(size) * size * size, 0);
            for (int z = 0; z < size; ++z)
            {
                for (int y = 0; y < size; ++y)
                {
                    for (int x = 0; x < size; ++x)
                    {
                        px[x + size * (y + size * z)] = (z == 0 && x == y ? 4 : -1);
                    }
                }
            }
            level = (wantDefect ? 0.5 : 2.5);
        }
    }

    void SaddleCase(oracle::Ctx& io, bool wantDefect)
    {
        int32_t N = 0;
        std::vector<int64_t> px;
        double level = 0.0;
        DrawSaddleImage(io, wantDefect, N, px, level);
        // Record 0 (the pinned reproduction): no fixBoundary, depth N.
        bool fixBoundary = (io.index() == 0 ? io.integer(1, 1) : io.integer(0, 3)) == 0;
        io.integer(1, 1);  // the number of extractions, the layout of ASCCase
        for (auto v : px) { io.given(static_cast<double>(v)); }
        io.given(level);
        std::vector<int32_t> depths{ io.index() == 0 ? io.integer(N, N) : io.integer(-2, N + 1) };
        std::vector<bool> sameDirs{ io.boolean() };
        std::vector<ASC3Result> results = RunASC<int32_t>(N, px, fixBoundary, { level }, depths, sameDirs);
        OutResult(io, results[0]);
    }
}

// Four-crossing faces on which upstream's det pairing agrees with the
// trilinear interpolant: disjoint pairs of both kinds and the plus-sign
// branch point (det = 0 at level 0), with levels on and between samples.
ORACLE_CASE("AdaptiveSkeletonClimbing3.extract.saddle") { SaddleCase(io, false); }

// Every record has a four-crossing face on which upstream's
// det = f00*f11 - f01*f10 pairing contradicts the interpolant (it ignores
// the level and its two disjoint pairings are swapped, #544 in 3-D). The
// port pairs by dg = det - L*(f00+f11-f01-f10), exactly.
ORACLE_CASE("AdaptiveSkeletonClimbing3.extract.saddlePairing") { SaddleCase(io, true); }

namespace
{
    // The voxel formula of the large case (mirrored in the replay), with
    // u = x - p[9], v = y - p[10], w = z - p[11]:
    // kind 0: a quadric p0 u^2 + p1 v^2 + p2 w^2 + p3 uv + p4 vw + p5 wu
    //         + p6 x + p7 y + p8 z (spheres, ellipsoids, hyperboloids,
    //         cylinders, planes);
    // kind 1: p0|u| + p1|v| + p2|w| + p3 min(|u|,|v|) + p4 max(|v|,|w|)
    //         + p6 x + p7 y + p8 z (octahedra, boxes: large flat faces);
    // kind 2: a torus (s^2 - 4 R^2 (a^2 + b^2)) * p5 with
    //         s = u^2 + v^2 + w^2 + R^2 - r^2, R = p0, r = p1 and (a, b) two
    //         of (u, v, w) chosen by p2 (a tunnel);
    // kind 3: two spheres p5 * min(|q - c|^2 - p0, |q - d|^2 - p1) with
    //         d = (p6, p7, p8);
    // kind 4: a cubic p0 uvw + p1 u^2 + p2 v^2 + p3 w^2 + p6 x + p7 y + p8 z
    //         (many saddles).
    int64_t LargeVoxel(int64_t kind, std::array<int64_t, 12> const& p, int64_t x, int64_t y, int64_t z)
    {
        int64_t u = x - p[9], v = y - p[10], w = z - p[11];
        int64_t au = (u < 0 ? -u : u), av = (v < 0 ? -v : v), aw = (w < 0 ? -w : w);
        switch (kind)
        {
        case 0:
            return p[0] * u * u + p[1] * v * v + p[2] * w * w + p[3] * u * v + p[4] * v * w
                + p[5] * w * u + p[6] * x + p[7] * y + p[8] * z;
        case 1:
            return p[0] * au + p[1] * av + p[2] * aw + p[3] * std::min(au, av)
                + p[4] * std::max(av, aw) + p[6] * x + p[7] * y + p[8] * z;
        case 2:
        {
            int64_t s = u * u + v * v + w * w + p[0] * p[0] - p[1] * p[1];
            int64_t a = (p[2] == 2 ? v : u), b = (p[2] == 0 ? v : w);
            return (s * s - 4 * p[0] * p[0] * (a * a + b * b)) * p[5];
        }
        case 3:
        {
            int64_t s0 = u * u + v * v + w * w - p[0];
            int64_t s1 = (x - p[6]) * (x - p[6]) + (y - p[7]) * (y - p[7]) + (z - p[8]) * (z - p[8]) - p[1];
            return p[5] * std::min(s0, s1);
        }
        default:
            return p[0] * u * v * w + p[1] * u * u + p[2] * v * v + p[3] * w * w
                + p[6] * x + p[7] * y + p[8] * z;
        }
    }

    void RawLargeParameters(oracle::Ctx& io, int size, int64_t& kind, std::array<int64_t, 12>& p)
    {
        kind = io.rawInteger(0, 4);
        for (int k = 0; k < 6; ++k) { p[k] = io.rawInteger(-3, 3); }
        for (int k = 6; k < 9; ++k) { p[k] = io.rawInteger(-10, 10); }
        for (int k = 9; k < 12; ++k) { p[k] = io.rawInteger(0, size - 1); }
        if (kind == 2)
        {
            p[0] = io.rawInteger(2, std::max(2, size / 3));
            p[1] = io.rawInteger(1, static_cast<int>(p[0]) - 1);
            p[2] = io.rawInteger(0, 2);
            p[5] = (io.rawInteger(0, 1) == 0 ? 1 : -1);
        }
        else if (kind == 3)
        {
            p[0] = io.rawInteger(1, size * size / 8);
            p[1] = io.rawInteger(1, size * size / 8);
            for (int k = 6; k < 9; ++k) { p[k] = io.rawInteger(0, size - 1); }
            p[5] = (io.rawInteger(0, 1) == 0 ? 1 : -1);
        }
    }
}

// 9^3, 17^3 and 33^3 images from a recorded closed-form integer formula
// (deep merge trees, large merged boxes, the face subdivision of
// Get*EdgesM and the fan triangulation RemoveTrianglesSE, tunnels, two
// components). The outputs are counts and digests (hundreds to thousands of
// triangles per record). Upstream-sound face pairings only (capped
// rejection over formula and level; fallback: a paraboloid bowl, which has
// no four-crossing face).
namespace
{
    // All border voxels on one side of the level (the level set does not
    // reach the image border) and at least one voxel on the other side.
    bool Interior(std::vector<int64_t> const& px, int size, double level)
    {
        int sides = 0, all = 0;
        for (int z = 0; z < size; ++z)
        {
            for (int y = 0; y < size; ++y)
            {
                for (int x = 0; x < size; ++x)
                {
                    int side = (static_cast<double>(px[x + size * (y + size * z)]) > level ? 1 : 2);
                    all |= side;
                    if (x == 0 || y == 0 || z == 0 || x == size - 1 || y == size - 1 || z == size - 1)
                    {
                        sides |= side;
                    }
                }
            }
        }
        return sides != 3 && all == 3;
    }

    void LargeCase(oracle::Ctx& io, bool interior)
    {
        int32_t N = static_cast<int32_t>(io.given(interior ? io.rawInteger(3, 4) : io.rawInteger(3, 5)));
        int const size = (1 << N) + 1;
        int64_t kind = 0;
        std::array<int64_t, 12> p{};
        std::vector<int64_t> px(static_cast<size_t>(size) * size * size);
        double level = 0.0;
        bool sound = false;
        for (int attempt = 0; attempt < (interior ? 64 : 32) && !sound; ++attempt)
        {
            RawLargeParameters(io, size, kind, p);
            if (interior)
            {
                for (int k = 9; k < 12; ++k) { p[k] = io.rawInteger(size / 4, 3 * size / 4); }
            }
            for (int z = 0; z < size; ++z)
            {
                for (int y = 0; y < size; ++y)
                {
                    for (int x = 0; x < size; ++x) { px[x + size * (y + size * z)] = LargeVoxel(kind, p, x, y, z); }
                }
            }
            auto mm = std::minmax_element(px.begin(), px.end());
            double lo = static_cast<double>(*mm.first), hi = static_cast<double>(*mm.second);
            int levelMode = io.rawInteger(0, 9);
            level = (levelMode < 7 ? std::floor(io.raw(lo, hi)) + 0.5
                : (levelMode < 9 ? std::floor(io.raw(lo, hi + 1.0)) : io.raw(lo, hi)));
            sound = CountFaces(px, size, level)[1] == 0 && (!interior || Interior(px, size, level));
        }
        if (!sound)
        {
            kind = 0;
            p = { 1, 1, 1, 0, 0, 0, 0, 0, 0, size / 2, size / 2, size / 2 };
            level = static_cast<double>((size / 4) * (size / 4)) + 0.5;
        }
        io.given(static_cast<double>(kind));
        for (auto v : p) { io.given(static_cast<double>(v)); }
        io.given(level);
        std::vector<int32_t> depths{ io.integer(-2, N + 1) };
        bool fixBoundary = (io.integer(0, 5) == 0);
        std::vector<bool> sameDirs{ io.boolean() };
        for (int z = 0; z < size; ++z)
        {
            for (int y = 0; y < size; ++y)
            {
                for (int x = 0; x < size; ++x) { px[x + size * (y + size * z)] = LargeVoxel(kind, p, x, y, z); }
            }
        }
        std::vector<ASC3Result> results = RunASC<int32_t>(N, px, fixBoundary, { level }, depths, sameDirs);
        OutResult(io, results[0]);
    }
}

ORACLE_CASE("AdaptiveSkeletonClimbing3.extract.large") { LargeCase(io, false); }

// The large case restricted to level sets that do not reach the image
// border (centers in the middle half, rejection on the border voxels): the
// mesh is expected to be closed, which the replay checks.
ORACLE_CASE("AdaptiveSkeletonClimbing3.extract.closedLarge") { LargeCase(io, true); }

// Level sets that do not reach the image border on 3^3 and 5^3 images: a
// constant border of value -sgn*(R+1) around random or checkerboard interior
// values in [-R, R] (isolated voxels, saddles, tunnels), levels between and
// (one in five) on the interior values. Upstream-sound face pairings only
// (capped rejection; fallback: a single interior voxel 5 in a border of -1,
// level 0.5). The layout is that of ASCCase.
ORACLE_CASE("AdaptiveSkeletonClimbing3.extract.closed")
{
    int32_t N = DrawN(io);
    int const size = (1 << N) + 1;
    bool fixBoundary = (io.integer(0, 3) == 0);
    io.integer(1, 1);  // the number of extractions
    std::vector<int64_t> px(static_cast<size_t>(size) * size * size);
    double level = 0.0;
    bool ok = false;
    for (int attempt = 0; attempt < 64 && !ok; ++attempt)
    {
        int64_t R = io.rawInteger(1, 9);
        int64_t sgn = (io.rawInteger(0, 1) == 0 ? 1 : -1);
        int mode = io.rawInteger(0, 1);
        for (int z = 0; z < size; ++z)
        {
            for (int y = 0; y < size; ++y)
            {
                for (int x = 0; x < size; ++x)
                {
                    bool border = (x == 0 || y == 0 || z == 0 || x == size - 1 || y == size - 1 || z == size - 1);
                    int64_t v = (mode == 0 ? io.rawInteger(-R, R)
                        : ((x + y + z) % 2 == 0 ? 1 : -1) * io.rawInteger(1, R));
                    px[x + size * (y + size * z)] = (border ? -sgn * (R + 1) : v);
                }
            }
        }
        level = (io.rawInteger(0, 4) == 0 ? std::floor(io.raw(-static_cast<double>(R), R + 1.0))
            : std::floor(io.raw(-R - 0.5, R + 0.5)) + 0.5);
        ok = CountFaces(px, size, level)[1] == 0 && Interior(px, size, level);
    }
    if (!ok)
    {
        std::fill(px.begin(), px.end(), -1);
        px[1 + size * (1 + size * 1)] = 5;
        level = 0.5;
    }
    for (auto v : px) { io.given(static_cast<double>(v)); }
    io.given(level);
    std::vector<int32_t> depths{ io.integer(-2, N + 1) };
    std::vector<bool> sameDirs{ io.boolean() };
    std::vector<ASC3Result> results = RunASC<int32_t>(N, px, fixBoundary, { level }, depths, sameDirs);
    OutResult(io, results[0]);
}

// #194 (preserved): an image that is strongly monotone along every axis (a
// ramp a*x + b*y + c*z + d with a nonzero coefficient, or a monotone cubic
// along one axis) merges into one box at the root, whose return value the
// root call discards, so at depth <= 0 the mesh is empty on both sides. The
// replay checks independently that the level set is not empty. N = 1..3,
// levels strictly inside the voxel range (half-integers, reals).
ORACLE_CASE("AdaptiveSkeletonClimbing3.extract.rootMonobox")
{
    int32_t N = static_cast<int32_t>(io.given(io.rawInteger(1, 2)));
    int const size = (1 << N) + 1;
    int64_t a = 0, b = 0, c = 0, cubic = 0;
    for (int attempt = 0; attempt < 64 && a == 0 && b == 0 && c == 0; ++attempt)
    {
        a = io.rawInteger(-3, 3);
        b = io.rawInteger(-3, 3);
        c = io.rawInteger(-3, 3);
    }
    if (a == 0 && b == 0 && c == 0) { a = 1; }
    cubic = io.rawInteger(0, 1);
    int64_t d = io.rawInteger(-5, 5);
    std::vector<int64_t> px(static_cast<size_t>(size) * size * size);
    for (int z = 0; z < size; ++z)
    {
        for (int y = 0; y < size; ++y)
        {
            for (int x = 0; x < size; ++x)
            {
                int64_t t = a * x + b * y + c * z;
                px[x + size * (y + size * z)] = (cubic != 0 && b == 0 && c == 0 ? t * t * t : t) + d;
            }
        }
    }
    bool fixBoundary = false;
    io.integer(1, 1);  // fixBoundary = false (the integer is 0 for true), the layout of ASCCase
    io.integer(1, 1);  // the number of extractions
    for (auto v : px) { io.given(static_cast<double>(v)); }
    auto mm = std::minmax_element(px.begin(), px.end());
    double lo = static_cast<double>(*mm.first), hi = static_cast<double>(*mm.second);
    double level = (io.rawInteger(0, 1) == 0 ? std::floor(io.raw(lo, hi)) + 0.5 : io.raw(lo, hi));
    io.given(level);
    std::vector<int32_t> depths{ io.integer(-2, 0) };
    std::vector<bool> sameDirs{ io.boolean() };
    std::vector<ASC3Result> results = RunASC<int32_t>(N, px, fixBoundary, { level }, depths, sameDirs);
    OutResult(io, results[0]);
}

// MakeUnique, OrientTriangles and ComputeNormals on arbitrary meshes over an
// arbitrary image: vertices on a quarter lattice in [-1.5, 2^N + 1] (so
// GetGradient's truncating conversion takes the -1 < x < 0 branch, and the
// zero gradient outside [0, 2^N) is reached), signed zeros, repeated
// vertices (a third of them re-use an earlier one with its zeros' signs
// flipped: std::map equivalence merges +0 and -0, the first one survives),
// triangles with repeated indices (degenerate: zero normal), rotations and
// duplicates of earlier triangles (MakeUnique keeps rotations apart),
// unused vertices (zero normals). Lists are emitted in full.
ORACLE_CASE("AdaptiveSkeletonClimbing3.meshOps")
{
    int32_t N = static_cast<int32_t>(io.given(io.rawInteger(1, 2)));
    int const size = (1 << N) + 1;
    std::vector<int32_t> voxels(static_cast<size_t>(size) * size * size);
    for (auto& v : voxels) { v = io.integer(-6, 6); }
    int32_t numVertices = io.integer(0, 12);
    std::vector<ASC3Vertex> vertices(numVertices);
    double const hi = static_cast<double>((1 << N) + 1);
    for (int32_t i = 0; i < numVertices; ++i)
    {
        int32_t kind = io.integer(0, 2);
        if (kind == 0 && i > 0)
        {
            int32_t j = io.integer(0, i - 1);
            vertices[i] = vertices[j];
            for (auto& c : vertices[i]) { if (c == 0.0) { c = -c; } }
        }
        else
        {
            for (auto& c : vertices[i])
            {
                double q = std::floor(io.raw(-1.5, hi) * 4.0) / 4.0;
                c = io.given(q == 0.0 && io.rawInteger(0, 1) == 0 ? -0.0 : q);
            }
        }
    }
    int32_t numTriangles = (numVertices == 0 ? io.integer(0, 0) : io.integer(0, 12));
    std::vector<ASC3Triangle> triangles(numTriangles);
    for (int32_t t = 0; t < numTriangles; ++t)
    {
        int32_t kind = io.integer(0, 3);
        if (kind == 0 && t > 0)
        {
            int32_t j = io.integer(0, t - 1);
            int32_t rotation = io.integer(0, 2);
            for (int k = 0; k < 3; ++k) { triangles[t].V[k] = triangles[j].V[(k + rotation) % 3]; }
        }
        else
        {
            for (int k = 0; k < 3; ++k) { triangles[t].V[k] = io.integer(0, numVertices - 1); }
        }
    }
    bool sameDir = io.boolean();
    AdaptiveSkeletonClimbing3<int32_t, double> asc(N, voxels.data());
    std::vector<ASC3Vertex> vu = vertices;
    std::vector<ASC3Triangle> tu = triangles;
    CallMakeUnique(asc, vu, tu);
    std::vector<ASC3Triangle> to = tu;
    CallOrient(asc, vu, to, sameDir);
    std::vector<std::array<double, 3>> normals;
    CallNormals(asc, vu, to, normals);
    std::vector<ASC3Triangle> ta = triangles;
    CallOrient(asc, vertices, ta, !sameDir);
    std::vector<std::array<double, 3>> na;
    CallNormals(asc, vertices, ta, na);
    OutVertices(io, vu);
    OutTriangles(io, tu);
    OutTriangles(io, to);
    OutNormals(io, normals);
    OutTriangles(io, ta);
    OutNormals(io, na);
}

// The constructor's LogError: N = 0 (a negative N makes upstream's 1 << N
// undefined behaviour) or a null voxel pointer. The other records (N = 1,
// valid pointer) extract a ramp normally.
ORACLE_CASE("AdaptiveSkeletonClimbing3.invalid")
{
    int32_t N = io.integer(0, 1);
    bool useNull = io.boolean();
    int const size = (1 << N) + 1;
    // a (x + y + z) + b y z, with (a, b) redrawn (at most 16 times, then
    // b = 0) until upstream's face pairing is sound at level 0.5: the image
    // 0, 1; 1, 0 of a = 1, b = -2 has a face saddle exactly on the level,
    // where the port's pairing deviates (saddlePairing covers that).
    std::vector<int64_t> px(static_cast<size_t>(size) * size * size);
    for (int attempt = 0; attempt < 17; ++attempt)
    {
        int64_t a = io.rawInteger(-3, 3), b = (attempt < 16 ? io.rawInteger(-3, 3) : 0);
        for (int z = 0; z < size; ++z)
        {
            for (int y = 0; y < size; ++y)
            {
                for (int x = 0; x < size; ++x) { px[x + size * (y + size * z)] = a * (x + y + z) + b * y * z; }
            }
        }
        if (CountFaces(px, size, 0.5)[1] == 0) { break; }
    }
    std::vector<int32_t> voxels(px.size());
    for (size_t i = 0; i < px.size(); ++i)
    {
        voxels[i] = static_cast<int32_t>(io.given(static_cast<double>(px[i])));
    }
    AdaptiveSkeletonClimbing3<int32_t, double> asc(N, useNull ? nullptr : voxels.data());
    std::vector<ASC3Vertex> v;
    std::vector<ASC3Triangle> t;
    CallExtract(asc, 0.5, -1, v, t);
    OutBoxes(io, GetBoxes(asc));
    OutVertices(io, v);
    OutTriangles(io, t);
}
