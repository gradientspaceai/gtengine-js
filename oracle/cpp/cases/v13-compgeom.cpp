// Verify group 13 (computational geometry): differential cases for
// SeparatePoints3.h and TriangulateCDT.h.
//
// INSTANTIATIONS. The C++ side instantiates exactly what the port implements:
//   SeparatePoints3<double, double>   ComputeType is vestigial upstream (the
//                                     current ConvexHull3 is templated on the
//                                     input type only); the port drops it.
//   TriangulateCDT<double>            the replacement class; the deprecated
//                                     TriangulateCDT<InputType, ComputeType>
//                                     specialization is not ported.
//
// COMPARISON POLICY. Everything compared here is combinatorial (booleans,
// indices, counts) or computed with + - * / sqrt on doubles (the separating
// plane is UnitCross and Dot), so every case is declared { exact: true }.
// Nothing on a compared path calls the C math library.
//
// CANONICALIZATION. TriangulateCDT fills PolygonTreeEx::allTriangles from
// ConstrainedDelaunay2::GetIndices() and outsideTriangles from
// ETManifoldMesh::GetTriangles(), both std::unordered_map iteration orders;
// the port enumerates them in sorted-key order. Those two lists are sorted
// by their (restored) vertex tuple on both sides. Every other list is
// emitted verbatim: Node::triangulation is filled from a
// std::set<TriangleKey<true>>, whose order is defined, and the interior,
// exterior and inside lists are built from the per-node triangulations in
// node order.
//
// SeparatePoints3 iterates a std::set<std::pair<size_t,size_t>> (defined
// order) and ConvexHull3::GetHull(). v11 showed that the hull's triangle SET
// is identical to the port's for at most 32 points (the sets here have at
// most 8), but its ORDER is std::unordered_map order and, for more than four
// points, even differs between two calls in one process (heap-address
// dependent, see TraceSeparate3). The query returns the first separating face
// in that order, so the SeparatePoints3 cases are built on an
// order-independent soundness predicate and emit the plane only where it is
// order independent; see SeparatePoints3.compute.
//
// One io draw per C++ statement: MSVC evaluates function arguments right to
// left, so every generated value goes into a named local before it is used.
#define ORACLE_FAMILY "v13-compgeom"
#include "Oracle.h"

#include <Mathematics/ArbitraryPrecision.h>
#include <Mathematics/ConvexHull3.h>
#include <Mathematics/SeparatePoints3.h>
#include <Mathematics/TriangulateCDT.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

using namespace gte;

namespace
{
    // BSNumber<UIntegerAP32> grows as needed and the expressions below use
    // only + - *, so these signs are exact for any double input. They are
    // used by the generator probes and by the independent reference checks.
    using Exact = BSNumber<UIntegerAP32>;
    using E3 = std::array<Exact, 3>;

    // Record an integer that was not produced by io.integer().
    void GivenInt(oracle::Ctx& io, int32_t v)
    {
        io.given(static_cast<double>(v));
    }

    int32_t ExactOrient2(Vector2<double> const& A, Vector2<double> const& B,
        Vector2<double> const& C)
    {
        Exact ax(A[0]), ay(A[1]), bx(B[0]), by(B[1]), cx(C[0]), cy(C[1]);
        Exact x0 = bx - ax, y0 = by - ay, x1 = cx - ax, y1 = cy - ay;
        Exact det = x0 * y1 - x1 * y0;
        return det.GetSign();
    }

    E3 ToE3(Vector3<double> const& p)
    {
        return E3{ Exact(p[0]), Exact(p[1]), Exact(p[2]) };
    }

    E3 SubE3(E3 const& a, E3 const& b)
    {
        return E3{ a[0] - b[0], a[1] - b[1], a[2] - b[2] };
    }

    E3 CrossE3(E3 const& u, E3 const& v)
    {
        return E3{
            u[1] * v[2] - u[2] * v[1],
            u[2] * v[0] - u[0] * v[2],
            u[0] * v[1] - u[1] * v[0] };
    }

    Exact DotE3(E3 const& u, E3 const& v)
    {
        return u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
    }

    // Sign of Dot(Cross(B - A, C - A), D - A): positive when D is on the side
    // the normal of the counterclockwise triangle <A,B,C> points to.
    int32_t ExactOrient3(E3 const& A, E3 const& B, E3 const& C, E3 const& D)
    {
        return DotE3(CrossE3(SubE3(B, A), SubE3(C, A)), SubE3(D, A)).GetSign();
    }

    // Sign of the 2D orientation of the projections of A, B, C that drop
    // coordinate 'drop'.
    int32_t ExactOrient2Drop(E3 const& A, E3 const& B, E3 const& C, int32_t drop)
    {
        int32_t i0 = (drop + 1) % 3, i1 = (drop + 2) % 3;
        Exact det = (B[i0] - A[i0]) * (C[i1] - A[i1]) - (B[i1] - A[i1]) * (C[i0] - A[i0]);
        return det.GetSign();
    }

    // Closed intersection of two coplanar segments in the projection that
    // drops coordinate 'drop' (the projected plane is nondegenerate).
    bool OnSegmentDrop(E3 const& P, E3 const& Q, E3 const& X, int32_t drop)
    {
        // X is collinear with PQ; test that it lies in the bounding box.
        for (int32_t k = 0; k < 3; ++k)
        {
            if (k == drop)
            {
                continue;
            }
            Exact lo = (P[k] < Q[k] ? P[k] : Q[k]);
            Exact hi = (P[k] < Q[k] ? Q[k] : P[k]);
            if (X[k] < lo || hi < X[k])
            {
                return false;
            }
        }
        return true;
    }

    bool SegmentsMeetDrop(E3 const& P, E3 const& Q, E3 const& A, E3 const& B, int32_t drop)
    {
        int32_t d0 = ExactOrient2Drop(P, Q, A, drop);
        int32_t d1 = ExactOrient2Drop(P, Q, B, drop);
        int32_t d2 = ExactOrient2Drop(A, B, P, drop);
        int32_t d3 = ExactOrient2Drop(A, B, Q, drop);
        if (d0 * d1 < 0 && d2 * d3 < 0)
        {
            return true;
        }
        return (d0 == 0 && OnSegmentDrop(P, Q, A, drop))
            || (d1 == 0 && OnSegmentDrop(P, Q, B, drop))
            || (d2 == 0 && OnSegmentDrop(A, B, P, drop))
            || (d3 == 0 && OnSegmentDrop(A, B, Q, drop));
    }

    // Closed segment <P,Q> against the closed nondegenerate triangle <A,B,C>.
    bool SegmentMeetsTriangle(E3 const& P, E3 const& Q, E3 const& A, E3 const& B,
        E3 const& C)
    {
        int32_t sP = ExactOrient3(A, B, C, P);
        int32_t sQ = ExactOrient3(A, B, C, Q);
        if (sP * sQ > 0)
        {
            return false;
        }
        if (sP == 0 && sQ == 0)
        {
            // Coplanar: project along a coordinate in which the triangle
            // normal is nonzero.
            E3 normal = CrossE3(SubE3(B, A), SubE3(C, A));
            int32_t drop = (normal[0].GetSign() != 0 ? 0 : (normal[1].GetSign() != 0 ? 1 : 2));
            std::array<E3 const*, 2> ends{ &P, &Q };
            for (auto end : ends)
            {
                int32_t o0 = ExactOrient2Drop(A, B, *end, drop);
                int32_t o1 = ExactOrient2Drop(B, C, *end, drop);
                int32_t o2 = ExactOrient2Drop(C, A, *end, drop);
                if ((o0 >= 0 && o1 >= 0 && o2 >= 0) || (o0 <= 0 && o1 <= 0 && o2 <= 0))
                {
                    return true;
                }
            }
            return SegmentsMeetDrop(P, Q, A, B, drop) || SegmentsMeetDrop(P, Q, B, C, drop)
                || SegmentsMeetDrop(P, Q, C, A, drop);
        }
        // The segment meets the plane in exactly one point; it is in the
        // closed triangle when the line PQ passes the three edges on one side.
        int32_t v0 = ExactOrient3(P, Q, A, B);
        int32_t v1 = ExactOrient3(P, Q, B, C);
        int32_t v2 = ExactOrient3(P, Q, C, A);
        return (v0 >= 0 && v1 >= 0 && v2 >= 0) || (v0 <= 0 && v1 <= 0 && v2 <= 0);
    }
}

// ---- SeparatePoints3 -----------------------------------------------------

namespace
{
    // Verbatim copies of upstream's private OnSameSide and WhichSide.
    int32_t FloatOnSameSide(Plane3<double> const& plane, size_t numTriangles,
        size_t const* indices, Vector3<double> const* points)
    {
        size_t posSide = 0, negSide = 0;
        for (size_t t = 0; t < numTriangles; ++t)
        {
            for (size_t i = 0; i < 3; ++i)
            {
                size_t v = indices[3 * t + i];
                double c0 = Dot(plane.normal, points[v]);
                if (c0 > plane.constant)
                {
                    ++posSide;
                }
                else if (c0 < plane.constant)
                {
                    ++negSide;
                }
                if (posSide && negSide)
                {
                    return 0;
                }
            }
        }
        return (posSide ? +1 : -1);
    }

    int32_t FloatWhichSide(Plane3<double> const& plane, size_t numTriangles,
        size_t const* indices, Vector3<double> const* points)
    {
        for (size_t t = 0; t < numTriangles; ++t)
        {
            for (size_t i = 0; i < 3; ++i)
            {
                size_t v = indices[3 * t + i];
                double c0 = Dot(plane.normal, points[v]);
                if (c0 > plane.constant)
                {
                    return +1;
                }
                if (c0 < plane.constant)
                {
                    return -1;
                }
            }
        }
        return 0;
    }

    // The same two functions with the exact side of a point: the sign of
    // Dot(N, Q - O) for the unnormalized normal N of the candidate plane
    // through O. This is the port's predicate.
    struct ExactPlane3
    {
        E3 normal;
        E3 origin;
    };

    int32_t ExactSideOf(ExactPlane3 const& plane, E3 const& q)
    {
        return DotE3(plane.normal, SubE3(q, plane.origin)).GetSign();
    }

    int32_t ExactOnSameSide(ExactPlane3 const& plane, std::vector<size_t> const& hull,
        std::vector<E3> const& points)
    {
        size_t posSide = 0, negSide = 0;
        for (size_t v : hull)
        {
            int32_t s = ExactSideOf(plane, points[v]);
            if (s > 0)
            {
                ++posSide;
            }
            else if (s < 0)
            {
                ++negSide;
            }
            if (posSide && negSide)
            {
                return 0;
            }
        }
        return (posSide ? +1 : -1);
    }

    int32_t ExactWhichSide(ExactPlane3 const& plane, std::vector<size_t> const& hull,
        std::vector<E3> const& points)
    {
        for (size_t v : hull)
        {
            int32_t s = ExactSideOf(plane, points[v]);
            if (s != 0)
            {
                return s;
            }
        }
        return 0;
    }

    std::vector<E3> ToE3s(std::vector<Vector3<double>> const& pts)
    {
        std::vector<E3> out(pts.size());
        for (size_t i = 0; i < pts.size(); ++i)
        {
            out[i] = ToE3(pts[i]);
        }
        return out;
    }

    ExactPlane3 ExactFacePlane(std::vector<E3> const& pts, size_t i0, size_t i1, size_t i2)
    {
        return ExactPlane3{ CrossE3(SubE3(pts[i1], pts[i0]), SubE3(pts[i2], pts[i0])), pts[i0] };
    }

    // Stages of SeparatePoints3::operator(), in the order upstream reaches them.
    enum SepStage : int32_t
    {
        STAGE_DEGENERATE = 0,   // a hull has dimension below 3 (early 'false')
        STAGE_FACE0 = 1,        // separated by a face plane of hull 0
        STAGE_FACE1 = 2,        // separated by a face plane of hull 1
        STAGE_EDGES = 3,        // separated by an edge-edge cross-product plane
        STAGE_NONE = 4          // every candidate failed ('false')
    };

    struct SepTrace
    {
        // True when upstream's float control flow provably equals the
        // port's exact one in every face order (TraceSeparate3).
        bool sound;
        int32_t stage;
        bool separated;
        Plane3<double> plane;
        // Whether upstream's returned plane is independent of its face
        // order and equal to the port's (TraceSeparate3); 'plane' is then
        // that plane. Always true for the non-face stages.
        bool determined;
    };

    // The float sign pattern of the owner hull's vertices against a
    // candidate face plane of that hull: whether any vertex classifies as
    // strictly positive and whether any classifies as strictly negative.
    // Upstream's WhichSide returns the first nonzero sign in hull order.
    void FloatOwnerSigns(Plane3<double> const& plane, std::vector<size_t> const& hull,
        std::vector<Vector3<double>> const& points, bool& hasPos, bool& hasNeg)
    {
        hasPos = false;
        hasNeg = false;
        for (size_t v : hull)
        {
            double c0 = Dot(plane.normal, points[v]);
            hasPos = hasPos || (c0 > plane.constant);
            hasNeg = hasNeg || (c0 < plane.constant);
        }
    }

    // One candidate face of one face loop, classified both ways.
    struct FaceClass
    {
        bool exactSeparates;    // the port's test: exact OnSameSide == +1
        bool floatPossible;     // upstream's float test passes in SOME order
        bool floatRobust;       // upstream's float test passes in EVERY order
        Plane3<double> plane;   // upstream's float plane of the face
    };

    // Upstream's algorithm and the port's evaluated side by side, as a
    // predicate on the input that does not depend on the order in which
    // upstream visits the hull faces and vertices.
    //
    // Why order independence is needed: upstream visits the faces, and in
    // WhichSide the vertices, in ConvexHull3::GetHull() order, which is
    // std::unordered_map order, and ConvexHull3::Hull3 seeds each
    // visible-region search from Vertex::TAdjacent, a std::unordered_set of
    // Triangle pointers. For more than four points the same input therefore
    // yields different hull orders in two calls within one process (measured
    // on 2000 uniform 8-point sets: 884 differ, and the returned plane of
    // SeparatePoints3 differs between two calls on 97), so no replica can
    // follow upstream's order.
    //
    // A face test is 'side_other != 0 && side_owner * side_other <= 0'.
    // OnSameSide's value is order independent (0 when both signs occur,
    // otherwise the sign that occurs, -1 when none). WhichSide's is the first
    // nonzero sign in visiting order: exactly it is always -1 for the owning
    // hull (outward face normals), in floating point it can be any nonzero
    // sign that occurs among the owner's vertices (0 when none occurs). So a
    // face's float test passes in some order when side_other != 0 and the
    // owner's float signs include -side_other or are all zero, and in every
    // order when they exclude side_other.
    //
    // The input is 'sound' when upstream's boolean result equals the port's
    // and upstream's plane is a genuine separator, in every order:
    //  - if some face separates exactly (the port returns true from a face
    //    loop), some face's float test must pass in every order (upstream
    //    then returns true, from that face or an earlier one), and every
    //    face whose float test can pass in a loop upstream can reach must
    //    separate exactly (otherwise upstream can return the right boolean
    //    with a non-separating plane, the #348 defect again; the first deep
    //    run without this condition showed it on 240 of 2000 records);
    //  - otherwise no face's float test may pass in any order, and the
    //    edge-edge loop, which runs over a std::set in a defined order, is
    //    replicated with the float and exact OnSameSide values compared
    //    along the loop until it returns; its plane is then bit-identical to
    //    upstream's (normal and constant).
    // A face stage's plane is 'determined' when the returning loop has
    // exactly one exactly-separating face, it is also the loop's only face
    // whose float test can pass, it passes in every order, and (for the loop
    // of hull 1) no face of hull 0 can pass in any order. Upstream then
    // returns that face's float plane in every order, and the port returns
    // the same face.
    //
    // With withFloat == false the float tests are skipped and the trace is
    // the port's algorithm alone.
    SepTrace TraceSeparate3(std::vector<Vector3<double>> const& pts0,
        std::vector<Vector3<double>> const& pts1, bool withFloat = true)
    {
        SepTrace unsound{ false, -1, false, Plane3<double>{}, false };
        ConvexHull3<double> ch0{};
        ch0(pts0.size(), pts0.data(), 0);
        if (ch0.GetDimension() != 3)
        {
            return SepTrace{ true, STAGE_DEGENERATE, false, Plane3<double>{}, true };
        }
        ConvexHull3<double> ch1{};
        ch1(pts1.size(), pts1.data(), 0);
        if (ch1.GetDimension() != 3)
        {
            return SepTrace{ true, STAGE_DEGENERATE, false, Plane3<double>{}, true };
        }

        auto const& hull0 = ch0.GetHull();
        auto const& hull1 = ch1.GetHull();
        size_t numTriangles0 = hull0.size() / 3;
        size_t numTriangles1 = hull1.size() / 3;
        std::vector<E3> ex0 = ToE3s(pts0);
        std::vector<E3> ex1 = ToE3s(pts1);

        std::array<std::vector<Vector3<double>> const*, 2> pts{ &pts0, &pts1 };
        std::array<std::vector<E3> const*, 2> ex{ &ex0, &ex1 };
        std::array<std::vector<size_t> const*, 2> hull{ &hull0, &hull1 };
        std::array<size_t, 2> numTriangles{ numTriangles0, numTriangles1 };
        std::array<std::vector<FaceClass>, 2> faces{};
        for (size_t owner = 0; owner < 2; ++owner)
        {
            size_t other = 1 - owner;
            auto const& ownerPts = *pts[owner];
            auto const& ownerHull = *hull[owner];
            for (size_t i = 0; i < numTriangles[owner]; ++i)
            {
                size_t i0 = ownerHull[3 * i];
                size_t i1 = ownerHull[3 * i + 1];
                size_t i2 = ownerHull[3 * i + 2];
                FaceClass fc{ false, false, false, Plane3<double>{} };
                ExactPlane3 eplane = ExactFacePlane(*ex[owner], i0, i1, i2);
                int32_t eOther = ExactOnSameSide(eplane, *hull[other], *ex[other]);
                if (eOther != 0)
                {
                    int32_t eOwner = ExactWhichSide(eplane, ownerHull, *ex[owner]);
                    fc.exactSeparates = (eOwner * eOther <= 0);
                }
                if (withFloat)
                {
                    fc.plane = Plane3<double>({ ownerPts[i0], ownerPts[i1], ownerPts[i2] });
                    int32_t fOther = FloatOnSameSide(fc.plane, numTriangles[other],
                        hull[other]->data(), pts[other]->data());
                    if (fOther != 0)
                    {
                        bool hasPos = false, hasNeg = false;
                        FloatOwnerSigns(fc.plane, ownerHull, ownerPts, hasPos, hasNeg);
                        bool allZero = !hasPos && !hasNeg;
                        fc.floatPossible = allZero || (fOther > 0 ? hasNeg : hasPos);
                        fc.floatRobust = allZero || (fOther > 0 ? !hasPos : !hasNeg);
                    }
                }
                faces[owner].push_back(fc);
            }
        }

        std::array<int32_t, 2> numExact{ 0, 0 }, numPossible{ 0, 0 }, numRobust{ 0, 0 };
        std::array<int32_t, 2> firstExact{ -1, -1 };
        for (size_t h = 0; h < 2; ++h)
        {
            for (size_t i = 0; i < faces[h].size(); ++i)
            {
                FaceClass const& fc = faces[h][i];
                if (fc.exactSeparates)
                {
                    if (numExact[h] == 0)
                    {
                        firstExact[h] = static_cast<int32_t>(i);
                    }
                    ++numExact[h];
                }
                numPossible[h] += (fc.floatPossible ? 1 : 0);
                numRobust[h] += (fc.floatRobust ? 1 : 0);
            }
        }

        if (numExact[0] + numExact[1] > 0)
        {
            size_t loop = (numExact[0] > 0 ? 0 : 1);
            int32_t stage = (loop == 0 ? STAGE_FACE0 : STAGE_FACE1);
            FaceClass const& first = faces[loop][static_cast<size_t>(firstExact[loop])];
            if (!withFloat)
            {
                return SepTrace{ true, stage, true, Plane3<double>{}, numExact[loop] == 1 };
            }
            if (numRobust[0] + numRobust[1] == 0)
            {
                return unsound;
            }
            // Every face upstream may return from must separate exactly:
            // a face that passes the float test in some order without
            // separating exactly is the #348 defect returning a wrong plane
            // with the right boolean. Loop 1 is always visited; loop 2 only
            // in orders in which no face of loop 1 passes, which exist only
            // when no face of loop 1 passes in every order.
            for (size_t h = 0; h < 2; ++h)
            {
                if (h == 1 && numRobust[0] > 0)
                {
                    break;
                }
                for (auto const& fc : faces[h])
                {
                    if (fc.floatPossible && !fc.exactSeparates)
                    {
                        return unsound;
                    }
                }
            }
            bool determined = numExact[loop] == 1 && numPossible[loop] == 1
                && first.floatRobust && (loop == 0 || numPossible[0] == 0);
            return SepTrace{ true, stage, true, first.plane, determined };
        }
        if (withFloat && numPossible[0] + numPossible[1] > 0)
        {
            return unsound;
        }

        std::array<std::set<std::pair<size_t, size_t>>, 2> edgeSet{};
        for (size_t h = 0; h < 2; ++h)
        {
            auto const& hh = *hull[h];
            for (size_t i = 0; i < numTriangles[h]; ++i)
            {
                size_t i0 = hh[3 * i];
                size_t i1 = hh[3 * i + 1];
                size_t i2 = hh[3 * i + 2];
                edgeSet[h].insert(std::make_pair(i0, i1));
                edgeSet[h].insert(std::make_pair(i0, i2));
                edgeSet[h].insert(std::make_pair(i1, i2));
            }
        }

        for (auto const& e0 : edgeSet[0])
        {
            Vector3<double> diff0 = pts0[e0.second] - pts0[e0.first];
            E3 ediff0 = SubE3(ex0[e0.second], ex0[e0.first]);
            for (auto const& e1 : edgeSet[1])
            {
                Vector3<double> diff1 = pts1[e1.second] - pts1[e1.first];
                ExactPlane3 eplane{ CrossE3(ediff0, SubE3(ex1[e1.second], ex1[e1.first])),
                    ex0[e0.first] };
                int32_t eSide0 = ExactOnSameSide(eplane, hull0, ex0);
                int32_t eSide1 = ExactOnSameSide(eplane, hull1, ex1);
                Plane3<double> fplane{};
                if (withFloat)
                {
                    fplane.normal = UnitCross(diff0, diff1);
                    fplane.constant = Dot(fplane.normal, pts0[e0.first]);
                    int32_t fSide0 = FloatOnSameSide(fplane, numTriangles0, hull0.data(),
                        pts0.data());
                    int32_t fSide1 = FloatOnSameSide(fplane, numTriangles1, hull1.data(),
                        pts1.data());
                    // The loop continues or returns on the sign of the
                    // product alone.
                    if ((fSide0 * fSide1 < 0) != (eSide0 * eSide1 < 0))
                    {
                        return unsound;
                    }
                }
                if (eSide0 * eSide1 < 0)
                {
                    return SepTrace{ true, STAGE_EDGES, true, fplane, true };
                }
            }
        }
        return SepTrace{ true, STAGE_NONE, false, Plane3<double>{}, true };
    }

    // ---- independent reference (shares no code with the query) ----------

    // Closed intersection of two convex polyhedra given by their point sets
    // and hull triangles, exactly: a point of one set in the other closed
    // hull, or a hull edge of one meeting a closed hull triangle of the
    // other. The vertices of the intersection of two convex polytopes are
    // vertices of either or edge-face intersections, so this is complete.
    bool PointInClosedHull(E3 const& q, std::vector<E3> const& pts,
        std::vector<size_t> const& hull)
    {
        for (size_t t = 0; t + 2 < hull.size(); t += 3)
        {
            if (ExactOrient3(pts[hull[t]], pts[hull[t + 1]], pts[hull[t + 2]], q) > 0)
            {
                return false;
            }
        }
        return true;
    }

    bool ClosedHullsMeet(std::vector<E3> const& pts0, std::vector<size_t> const& hull0,
        std::vector<E3> const& pts1, std::vector<size_t> const& hull1)
    {
        for (auto const& q : pts0)
        {
            if (PointInClosedHull(q, pts1, hull1))
            {
                return true;
            }
        }
        for (auto const& q : pts1)
        {
            if (PointInClosedHull(q, pts0, hull0))
            {
                return true;
            }
        }
        std::array<std::vector<E3> const*, 2> pts{ &pts0, &pts1 };
        std::array<std::vector<size_t> const*, 2> hull{ &hull0, &hull1 };
        for (size_t a = 0; a < 2; ++a)
        {
            auto const& pa = *pts[a];
            auto const& ha = *hull[a];
            auto const& pb = *pts[1 - a];
            auto const& hb = *hull[1 - a];
            for (size_t t = 0; t + 2 < ha.size(); t += 3)
            {
                for (size_t j = 0; j < 3; ++j)
                {
                    E3 const& P = pa[ha[t + j]];
                    E3 const& Q = pa[ha[t + (j + 1) % 3]];
                    for (size_t u = 0; u + 2 < hb.size(); u += 3)
                    {
                        if (SegmentMeetsTriangle(P, Q, pb[hb[u]], pb[hb[u + 1]], pb[hb[u + 2]]))
                        {
                            return true;
                        }
                    }
                }
            }
        }
        return false;
    }

    // The largest residual by which the rounded plane fails to put one set
    // on its nonpositive side and the other on its nonnegative side.
    double SeparationResidual(std::vector<Vector3<double>> const& pts0,
        std::vector<Vector3<double>> const& pts1, Plane3<double> const& plane)
    {
        double maxNeg0 = 0.0, maxPos0 = 0.0, maxNeg1 = 0.0, maxPos1 = 0.0;
        for (auto const& p : pts0)
        {
            double s = Dot(plane.normal, p) - plane.constant;
            maxPos0 = std::max(maxPos0, s);
            maxNeg0 = std::max(maxNeg0, -s);
        }
        for (auto const& p : pts1)
        {
            double s = Dot(plane.normal, p) - plane.constant;
            maxPos1 = std::max(maxPos1, s);
            maxNeg1 = std::max(maxNeg1, -s);
        }
        // Set 0 below and set 1 above, or the other way round.
        return std::min(std::max(maxPos0, maxNeg1), std::max(maxNeg0, maxPos1));
    }

    double MaxAbsCoordinate(std::vector<Vector3<double>> const& pts0,
        std::vector<Vector3<double>> const& pts1)
    {
        double scale = 1.0;
        for (auto const* set : { &pts0, &pts1 })
        {
            for (auto const& p : *set)
            {
                for (int32_t j = 0; j < 3; ++j)
                {
                    scale = std::max(scale, std::fabs(p[j]));
                }
            }
        }
        return scale;
    }

    // true: the returned plane separates to 1e-12 relative and is unit
    // length. false with two 3D hulls: the closed hulls meet (exact).
    bool SeparationReference(std::vector<Vector3<double>> const& pts0,
        std::vector<Vector3<double>> const& pts1, bool separated,
        Plane3<double> const& plane)
    {
        if (separated)
        {
            double tol = 1e-12 * MaxAbsCoordinate(pts0, pts1);
            double unit = Dot(plane.normal, plane.normal) - 1.0;
            return std::fabs(unit) <= 1e-12 && SeparationResidual(pts0, pts1, plane) <= tol;
        }
        ConvexHull3<double> ch0{};
        ch0(pts0.size(), pts0.data(), 0);
        ConvexHull3<double> ch1{};
        ch1(pts1.size(), pts1.data(), 0);
        if (ch0.GetDimension() != 3 || ch1.GetDimension() != 3)
        {
            return true;
        }
        return ClosedHullsMeet(ToE3s(pts0), ch0.GetHull(), ToE3s(pts1), ch1.GetHull());
    }

    // ---- generators (raw, unrecorded) ------------------------------------

    void RawLattice3(oracle::Ctx& io, size_t n, int32_t lo, int32_t hi,
        std::vector<Vector3<double>>& pts)
    {
        pts.resize(n);
        for (size_t i = 0; i < n; ++i)
        {
            for (int32_t j = 0; j < 3; ++j)
            {
                pts[i][j] = static_cast<double>(io.rawInteger(lo, hi));
            }
        }
    }

    void RawUniform3(oracle::Ctx& io, size_t n, double lo, double hi,
        std::vector<Vector3<double>>& pts)
    {
        pts.resize(n);
        for (size_t i = 0; i < n; ++i)
        {
            for (int32_t j = 0; j < 3; ++j)
            {
                pts[i][j] = io.raw(lo, hi);
            }
        }
    }

    void Translate(std::vector<Vector3<double>>& pts, Vector3<double> const& shift)
    {
        for (auto& p : pts)
        {
            p = p + shift;
        }
    }

    // Apply the same random signed axis permutation to both sets. Lattice
    // values stay exact.
    void RawSignedPermutation(oracle::Ctx& io, std::vector<Vector3<double>>& pts0,
        std::vector<Vector3<double>>& pts1)
    {
        std::array<int32_t, 3> perm{ 0, 1, 2 };
        for (int32_t i = 2; i > 0; --i)
        {
            int32_t j = io.rawInteger(0, i);
            std::swap(perm[static_cast<size_t>(i)], perm[static_cast<size_t>(j)]);
        }
        std::array<double, 3> sign{};
        for (int32_t j = 0; j < 3; ++j)
        {
            sign[static_cast<size_t>(j)] = (io.rawInteger(0, 1) != 0 ? -1.0 : 1.0);
        }
        for (auto* set : { &pts0, &pts1 })
        {
            for (auto& p : *set)
            {
                Vector3<double> q = p;
                for (int32_t j = 0; j < 3; ++j)
                {
                    p[j] = sign[static_cast<size_t>(j)] * q[perm[static_cast<size_t>(j)]];
                }
            }
        }
    }

    // Two tetrahedra with skew edges: A = {(-s,0,0), (s,0,0), (0,-t,-h),
    // (0,t,-h)} has its top edge on the x-axis, B = {(0,-u,d), (0,u,d),
    // (-v,0,d+k), (v,0,d+k)} its bottom edge parallel to the y-axis at height
    // d. d < 0 overlaps, d = 0 touches at the crossing point, d > 0 is
    // separated by the cross product of the two edges (or by a face plane
    // when the tetrahedra are flat enough). Extra points are midpoints of
    // two base vertices (exact halves). 'jitter' perturbs every coordinate.
    void RawCrossedTetrahedra(oracle::Ctx& io, size_t n0, size_t n1, bool jitter,
        std::vector<Vector3<double>>& pts0, std::vector<Vector3<double>>& pts1)
    {
        double s = static_cast<double>(io.rawInteger(1, 3));
        double t = static_cast<double>(io.rawInteger(1, 3));
        double h = static_cast<double>(io.rawInteger(1, 3));
        double u = static_cast<double>(io.rawInteger(1, 3));
        double v = static_cast<double>(io.rawInteger(1, 3));
        double k = static_cast<double>(io.rawInteger(1, 3));
        double d = static_cast<double>(io.rawInteger(-1, 2));
        std::array<Vector3<double>, 4> a{ {
            { -s, 0.0, 0.0 }, { s, 0.0, 0.0 }, { 0.0, -t, -h }, { 0.0, t, -h } } };
        std::array<Vector3<double>, 4> b{ {
            { 0.0, -u, d }, { 0.0, u, d }, { -v, 0.0, d + k }, { v, 0.0, d + k } } };
        pts0.assign(a.begin(), a.end());
        pts1.assign(b.begin(), b.end());
        while (pts0.size() < n0)
        {
            size_t i = static_cast<size_t>(io.rawInteger(0, 3));
            size_t j = static_cast<size_t>(io.rawInteger(0, 3));
            pts0.push_back((a[i] + a[j]) * 0.5);
        }
        while (pts1.size() < n1)
        {
            size_t i = static_cast<size_t>(io.rawInteger(0, 3));
            size_t j = static_cast<size_t>(io.rawInteger(0, 3));
            pts1.push_back((b[i] + b[j]) * 0.5);
        }
        if (jitter)
        {
            for (auto* set : { &pts0, &pts1 })
            {
                for (auto& p : *set)
                {
                    for (int32_t j = 0; j < 3; ++j)
                    {
                        p[j] += io.raw(-0.2, 0.2);
                    }
                }
            }
        }
        RawSignedPermutation(io, pts0, pts1);
        if (io.rawInteger(0, 1) != 0)
        {
            std::swap(pts0, pts1);
        }
    }

    // Points on two spheres whose gap is small (negative gaps overlap).
    void RawNearTangentSpheres(oracle::Ctx& io, size_t n0, size_t n1,
        std::vector<Vector3<double>>& pts0, std::vector<Vector3<double>>& pts1)
    {
        auto rawUnit = [&io]()
        {
            Vector3<double> w{ 1.0, 0.0, 0.0 };
            for (int32_t attempt = 0; attempt < 32; ++attempt)
            {
                Vector3<double> c{ io.raw(-1.0, 1.0), 0.0, 0.0 };
                c[1] = io.raw(-1.0, 1.0);
                c[2] = io.raw(-1.0, 1.0);
                double len = Length(c);
                if (len >= 0.1 && len <= 1.0)
                {
                    w = c / len;
                    break;
                }
            }
            return w;
        };
        Vector3<double> c0{ io.raw(-2.0, 2.0), 0.0, 0.0 };
        c0[1] = io.raw(-2.0, 2.0);
        c0[2] = io.raw(-2.0, 2.0);
        double r0 = io.raw(1.0, 2.0);
        double r1 = io.raw(1.0, 2.0);
        double gap = io.raw(-0.5, 0.5);
        Vector3<double> dir = rawUnit();
        Vector3<double> c1 = c0 + (r0 + r1 + gap) * dir;
        pts0.resize(n0);
        for (size_t i = 0; i < n0; ++i)
        {
            pts0[i] = c0 + r0 * rawUnit();
        }
        pts1.resize(n1);
        for (size_t i = 0; i < n1; ++i)
        {
            pts1[i] = c1 + r1 * rawUnit();
        }
    }

    // Axis-aligned lattice boxes that share a face, an edge or a vertex, are
    // a unit apart, or overlap by a unit; each set is a random subset of its
    // box's corners (a set of four corners of one face is coplanar, which is
    // the degenerate early return). With 'shear' both sets are mapped by the
    // same integer shear, so faces and shared features are no longer axis
    // aligned (upstream's rounded face normals then decide the touching
    // classifications, and the soundness probe rejects many of them).
    void RawTouchingBoxes(oracle::Ctx& io, size_t n0, size_t n1,
        std::vector<Vector3<double>>& pts0, std::vector<Vector3<double>>& pts1)
    {
        std::array<double, 3> size0{}, size1{}, shift{};
        for (size_t j = 0; j < 3; ++j)
        {
            size0[j] = static_cast<double>(io.rawInteger(1, 3));
            size1[j] = static_cast<double>(io.rawInteger(1, 3));
        }
        // contact: 0 face, 1 edge, 2 vertex, 3 gap, 4 overlap
        int32_t contact = io.rawInteger(0, 4);
        int32_t numAxes = (contact == 0 || contact >= 3 ? 1 : (contact == 1 ? 2 : 3));
        for (size_t j = 0; j < 3; ++j)
        {
            if (static_cast<int32_t>(j) < numAxes)
            {
                shift[j] = size0[j] + (contact == 3 ? 1.0 : (contact == 4 ? -1.0 : 0.0));
            }
            else
            {
                shift[j] = static_cast<double>(io.rawInteger(-1, 1));
            }
        }
        auto corners = [&io](std::array<double, 3> const& size, std::array<double, 3> const& at,
            size_t n, std::vector<Vector3<double>>& pts)
        {
            std::array<int32_t, 8> order{ 0, 1, 2, 3, 4, 5, 6, 7 };
            for (int32_t i = 7; i > 0; --i)
            {
                int32_t j = io.rawInteger(0, i);
                std::swap(order[static_cast<size_t>(i)], order[static_cast<size_t>(j)]);
            }
            pts.resize(n);
            for (size_t i = 0; i < n; ++i)
            {
                int32_t c = order[i % 8];
                for (size_t j = 0; j < 3; ++j)
                {
                    pts[i][static_cast<int32_t>(j)] = at[j] + ((c >> j) & 1 ? size[j] : 0.0);
                }
            }
        };
        std::array<double, 3> origin{ 0.0, 0.0, 0.0 };
        corners(size0, origin, n0, pts0);
        corners(size1, shift, n1, pts1);
        if (io.rawInteger(0, 2) == 0)
        {
            int32_t a01 = io.rawInteger(-1, 1);
            int32_t a02 = io.rawInteger(-1, 1);
            int32_t a12 = io.rawInteger(-1, 1);
            for (auto* set : { &pts0, &pts1 })
            {
                for (auto& p : *set)
                {
                    p[0] = p[0] + a01 * p[1] + a02 * p[2];
                    p[1] = p[1] + a12 * p[2];
                }
            }
        }
        RawSignedPermutation(io, pts0, pts1);
    }

    // One or both sets have a hull of dimension below 3: all equal, a
    // lattice line, a lattice plane with a non-axis normal, or at most three
    // points; kind 4 puts both sets in the same plane.
    void RawDegenerate(oracle::Ctx& io, size_t n0, size_t n1,
        std::vector<Vector3<double>>& pts0, std::vector<Vector3<double>>& pts1)
    {
        int32_t kind = io.rawInteger(0, 4);
        RawLattice3(io, n0, -3, 3, pts0);
        RawLattice3(io, n1, -3, 3, pts1);
        double bx = static_cast<double>(io.rawInteger(-3, 3));
        double by = static_cast<double>(io.rawInteger(-3, 3));
        double bz = static_cast<double>(io.rawInteger(-3, 3));
        double dx = static_cast<double>(io.rawInteger(-2, 2));
        double dy = static_cast<double>(io.rawInteger(1, 2));
        double dz = static_cast<double>(io.rawInteger(-2, 2));
        double a = static_cast<double>(io.rawInteger(1, 3));
        double b = static_cast<double>(io.rawInteger(-3, 3));
        double c = static_cast<double>(io.rawInteger(-3, 3));
        auto flatten = [&](std::vector<Vector3<double>>& pts)
        {
            for (auto& p : pts)
            {
                if (kind == 0)
                {
                    p = { bx, by, bz };
                }
                else if (kind == 1)
                {
                    double k = static_cast<double>(io.rawInteger(-3, 3));
                    p = { bx + k * dx, by + k * dy, bz + k * dz };
                }
                else if (kind >= 2)
                {
                    // The lattice plane through (bx,by,bz) spanned by
                    // (b,-a,0) and (c,0,-a).
                    double s = static_cast<double>(io.rawInteger(-3, 3));
                    double t = static_cast<double>(io.rawInteger(-3, 3));
                    p = { bx + s * b + t * c, by - s * a, bz - t * a };
                }
            }
        };
        bool first = (io.rawInteger(0, 1) != 0);
        if (kind == 4)
        {
            flatten(pts0);
            flatten(pts1);
        }
        else if (kind < 3)
        {
            flatten(first ? pts0 : pts1);
        }
        // kind 3 keeps the lattice sets; they are degenerate only when a set
        // has at most three points or happens to be flat.
    }

    // mode 0: two lattice sets in [-3,3]^3, the second shifted by a lattice
    //         vector in [-6,6]^3
    // mode 1: two uniform sets in [-4,4]^3, integer shift in [-8,8]^3
    // mode 2: crossed tetrahedra (edge-edge separation, touching, overlap),
    //         lattice or jittered
    // mode 3: near-tangent spheres (uniform)
    // mode 4: touching lattice boxes (face, edge, vertex contact), sheared
    //         on a third of the records
    // mode 5: degenerate sets (dimension 0, 1, 2)
    void RawSeparatePair3(oracle::Ctx& io, int32_t mode, size_t n0, size_t n1,
        std::vector<Vector3<double>>& pts0, std::vector<Vector3<double>>& pts1)
    {
        if (mode == 0 || mode == 1)
        {
            if (mode == 0)
            {
                RawLattice3(io, n0, -3, 3, pts0);
                RawLattice3(io, n1, -3, 3, pts1);
            }
            else
            {
                RawUniform3(io, n0, -4.0, 4.0, pts0);
                RawUniform3(io, n1, -4.0, 4.0, pts1);
            }
            int32_t range = (mode == 0 ? 6 : 8);
            Vector3<double> shift{ 0.0, 0.0, 0.0 };
            for (int32_t j = 0; j < 3; ++j)
            {
                shift[j] = static_cast<double>(io.rawInteger(-range, range));
            }
            Translate(pts1, shift);
        }
        else if (mode == 2)
        {
            bool jitter = (io.rawInteger(0, 1) != 0);
            RawCrossedTetrahedra(io, n0, n1, jitter, pts0, pts1);
        }
        else if (mode == 3)
        {
            RawNearTangentSpheres(io, n0, n1, pts0, pts1);
        }
        else if (mode == 4)
        {
            RawTouchingBoxes(io, n0, n1, pts0, pts1);
        }
        else
        {
            RawDegenerate(io, n0, n1, pts0, pts1);
        }
    }

    void FallbackCubes(size_t n0, size_t n1, std::vector<Vector3<double>>& pts0,
        std::vector<Vector3<double>>& pts1)
    {
        // Two axis-aligned unit cubes far apart: every side test is on small
        // integers against an axis-aligned unit normal.
        pts0.resize(n0);
        pts1.resize(n1);
        for (size_t i = 0; i < n0; ++i)
        {
            pts0[i] = { static_cast<double>(i & 1), static_cast<double>((i >> 1) & 1),
                static_cast<double>((i >> 2) & 1) };
        }
        for (size_t i = 0; i < n1; ++i)
        {
            pts1[i] = { 10.0 + static_cast<double>(i & 1), static_cast<double>((i >> 1) & 1),
                static_cast<double>((i >> 2) & 1) };
        }
        // Four corners (i < 4) of a cube are coplanar; lift the fourth.
        if (n0 == 4)
        {
            pts0[3][2] = 1.0;
        }
        if (n1 == 4)
        {
            pts1[3][2] = 1.0;
        }
    }

    void GivePoints3(oracle::Ctx& io, std::vector<Vector3<double>> const& pts)
    {
        for (auto const& p : pts)
        {
            io.givenVec<3>(p);
        }
    }

    // The returned plane is determined (independent of the face order) on
    // the degenerate, edge-edge and no-separation stages and on a face stage
    // with exactly one separating face.
    bool PlaneDetermined(SepTrace const& trace)
    {
        return trace.determined;
    }

    // Run upstream and check it against the trace: a mismatch would mean the
    // probe does not describe upstream, and becomes a throw record that the
    // replay reports.
    bool RunSeparatePoints3(std::vector<Vector3<double>> const& pts0,
        std::vector<Vector3<double>> const& pts1, SepTrace const& trace,
        Plane3<double>& plane)
    {
        SeparatePoints3<double, double> query{};
        bool separated = query(pts0.size(), pts0.data(), pts1.size(), pts1.data(), plane);
        if (trace.sound)
        {
            if (separated != trace.separated || (separated && PlaneDetermined(trace) &&
                (plane.normal != trace.plane.normal || plane.constant != trace.plane.constant)))
            {
                throw std::logic_error("SeparatePoints3 trace does not replicate upstream");
            }
        }
        return separated;
    }
}

// SeparatePoints3<double, double>::operator() on inputs where upstream's
// floating-point answer provably equals the port's exact one in every face
// order (TraceSeparate3: the boolean always, and every plane upstream can
// return is a genuine separator). The capped rejection loop redraws both
// sets on every attempt; measured acceptance per draw over 2000 records:
// lattice 35%, uniform 13%, crossed tetrahedra 52%, near-tangent spheres
// 11%, touching boxes 94%, degenerate 95%, with 2 fallbacks (the far-apart
// cubes) in 2000 records. The rejected draws are exactly the inputs on
// which upstream's result depends on its hull order: a face of the owner
// hull whose own vertex rounds to the wrong side (#348 finding 1).
//
// ORDER DEPENDENCE. The face loops visit ConvexHull3::GetHull() triangles
// in std::unordered_map order (heap-address dependent; the port enumerates
// them in sorted-key order) and return the FIRST face that separates. When
// a loop has more than one separating face, which plane upstream returns is
// not reproducible (group report, "Upstream bug suspects"); 'separated' is
// order independent. The plane is therefore emitted only when it is
// determined: a face stage whose loop has exactly one face that can pass,
// or the edge-edge stage (reached only when no face separates; its
// std::set order is defined). The reference check below covers the
// undetermined records on each side.
//
// The stage the trace reached (degenerate, face of hull 0, face of hull 1,
// edge-edge plane, none) and the 'determined' flag are recorded as inputs
// after the points. They only steer which outputs are emitted: 'separated'
// always; the plane's normal and constant on a determined true return; the
// plane's origin only when a face plane separated. On the edge-edge path
// upstream assigns normal and constant directly and leaves 'origin' holding
// the last face plane it built, the inconsistency the port fixes (#348,
// finding 2); that field is compared in
// SeparatePoints3.deviation.edgeAxisOrigin instead. On a false return
// upstream leaves the caller's plane holding the last candidate, which its
// documentation does not promise; the port returns a default plane, and the
// plane is not emitted.
//
// The last output is an independent reference check computed on each side
// from its own result: on a true return the plane is unit length and puts
// one set on each closed side (residual <= 1e-12 * max|coordinate|); on a
// false return with two 3D hulls the closed hulls meet, exhibited exactly (a
// point of one set in the other hull, or a hull edge meeting a hull
// triangle).
ORACLE_CASE("SeparatePoints3.compute")
{
    int32_t mode = io.index() % 6;
    int32_t minCount = (mode == 5 ? 1 : 4);
    int32_t n0 = io.rawInteger(minCount, 8);
    int32_t n1 = io.rawInteger(minCount, 8);
    std::vector<Vector3<double>> pts0{}, pts1{};
    SepTrace trace{ false, -1, false, Plane3<double>{}, false };
    for (int32_t attempt = 0; attempt < 64 && !trace.sound; ++attempt)
    {
        pts0.clear();
        pts1.clear();
        RawSeparatePair3(io, mode, static_cast<size_t>(n0), static_cast<size_t>(n1),
            pts0, pts1);
        trace = TraceSeparate3(pts0, pts1);
    }
    if (!trace.sound)
    {
        FallbackCubes(static_cast<size_t>(n0), static_cast<size_t>(n1), pts0, pts1);
        trace = TraceSeparate3(pts0, pts1);
    }
    // The crossed-tetrahedra generator may swap the sets, so the sizes are
    // recorded after the draw.
    GivenInt(io, static_cast<int32_t>(pts0.size()));
    GivenInt(io, static_cast<int32_t>(pts1.size()));
    GivePoints3(io, pts0);
    GivePoints3(io, pts1);
    GivenInt(io, trace.stage);
    bool determined = PlaneDetermined(trace);
    GivenInt(io, determined ? 1 : 0);

    Plane3<double> plane{};
    bool separated = RunSeparatePoints3(pts0, pts1, trace, plane);
    io.outBool(separated);
    if (separated && determined)
    {
        io.outVec(plane.normal);
        io.outReal(plane.constant);
        if (trace.stage == STAGE_FACE0 || trace.stage == STAGE_FACE1)
        {
            io.outVec(plane.origin);
        }
    }
    io.outBool(SeparationReference(pts0, pts1, separated, plane));
}

namespace
{
    // A verbatim replica of upstream's first face loop (the faces of hull 0)
    // in ConvexHull3::GetHull() order: true when some face of hull 0 passes
    // the float test 'side1 != 0 && side0 * side1 <= 0'. Upstream returns
    // true at the first such face.
    bool FloatFirstLoopSeparates(std::vector<Vector3<double>> const& pts0,
        std::vector<Vector3<double>> const& pts1)
    {
        ConvexHull3<double> ch0{};
        ch0(pts0.size(), pts0.data(), 0);
        ConvexHull3<double> ch1{};
        ch1(pts1.size(), pts1.data(), 0);
        if (ch0.GetDimension() != 3 || ch1.GetDimension() != 3)
        {
            return false;
        }
        auto const& hull0 = ch0.GetHull();
        auto const& hull1 = ch1.GetHull();
        for (size_t i = 0; i + 2 < hull0.size(); i += 3)
        {
            Plane3<double> fplane({ pts0[hull0[i]], pts0[hull0[i + 1]], pts0[hull0[i + 2]] });
            int32_t side1 = FloatOnSameSide(fplane, hull1.size() / 3, hull1.data(), pts1.data());
            if (side1)
            {
                int32_t side0 = FloatWhichSide(fplane, hull0.size() / 3, hull0.data(), pts0.data());
                if (side0 * side1 <= 0)
                {
                    return true;
                }
            }
        }
        return false;
    }
}

// DELIBERATE DEVIATION (#348, finding 1). OnSameSide and WhichSide compare
// Dot(normal, Q) with the rounded constant Dot(normal, P0) of the rounded
// plane Plane3({P0,P1,P2}), so a vertex of the hull's own candidate face can
// classify as strictly positive; when WhichSide meets that vertex before any
// strictly negative one it returns +1 instead of -1, and 'side0 * side1 <= 0'
// reports two overlapping clouds as separated. The port evaluates the side
// with the exact sign of Dot(Cross(P1-P0, P2-P0), Q-P0).
//
// DETERMINISM. Whether WhichSide meets the rounded vertex first depends on
// the order of ConvexHull3::GetHull(), which for more than four points
// varies between calls in one process (ConvexHull3::Hull3 starts each
// visible-region search from Vertex::TAdjacent, a std::unordered_set of
// Triangle pointers, so the hull mesh's std::unordered_map is filled in a
// heap-address-dependent order; see the group report). For exactly four
// points Hull3 inserts the four faces of the initial tetrahedron in a fixed
// sequence and never searches, so the order is reproducible. Set 0 is
// therefore a tetrahedron, and a draw is accepted only when the defect fires
// in upstream's FIRST loop (the faces of hull 0; the order of hull 1 matters
// there only through OnSameSide, whose value is order independent), as
// replicated by FloatFirstLoopSeparates on the same deterministic order.
//
// The generator draws a tetrahedron and a cloud of 4 to 8 points in
// [-1,1]^3 (the cloud shifted by at most 0.5 per coordinate) and accepts a
// draw only when (a) the exact algorithm (the port's) reports no separation,
// (b) the closed hulls meet, shown exactly by the independent reference,
// (c) the float first loop separates and upstream reports a separation, and
// (d) upstream's plane fails to separate by more than 1e-12 (it puts points
// strictly on both sides). Every record is therefore a genuinely wrong
// upstream answer, and every record deviates.
ORACLE_CASE("SeparatePoints3.deviation.roundoff")
{
    int32_t n1 = io.rawInteger(4, 8);
    std::vector<Vector3<double>> pts0{}, pts1{};
    bool accepted = false;
    for (int32_t attempt = 0; attempt < 4096 && !accepted; ++attempt)
    {
        RawUniform3(io, 4, -1.0, 1.0, pts0);
        RawUniform3(io, static_cast<size_t>(n1), -1.0, 1.0, pts1);
        Vector3<double> shift{ io.raw(-0.5, 0.5), 0.0, 0.0 };
        shift[1] = io.raw(-0.5, 0.5);
        shift[2] = io.raw(-0.5, 0.5);
        Translate(pts1, shift);
        SepTrace exact = TraceSeparate3(pts0, pts1, false);
        if (exact.stage != STAGE_NONE || !FloatFirstLoopSeparates(pts0, pts1))
        {
            continue;
        }
        Plane3<double> plane{};
        SeparatePoints3<double, double> probe{};
        bool separated = probe(pts0.size(), pts0.data(), pts1.size(), pts1.data(), plane);
        if (!separated)
        {
            continue;
        }
        double tol = 1e-12 * MaxAbsCoordinate(pts0, pts1);
        accepted = SeparationResidual(pts0, pts1, plane) > tol
            && SeparationReference(pts0, pts1, false, plane);
    }
    if (!accepted)
    {
        // A draw this generator accepted (record 0 of a 200-record run);
        // 17 significant digits round-trip exactly.
        pts0 = {
            { -0.030824049258111996, 0.31266741903663320, 0.18246756566068445 },
            { -0.54321918448773987, 0.060158433097447617, 0.052621387469387049 },
            { 0.98731848468561356, 0.43104447012556624, 0.25548570324684983 },
            { -0.96693191883509089, 0.056435527527842977, -0.24337962436569960 } };
        pts1 = {
            { 0.37316851385309235, 0.51908704770480940, 0.10903053815676567 },
            { 0.75694897211238088, -0.16034819722610727, 0.58566420501042160 },
            { 0.94991695522726094, 0.86534417964425214, -0.70869795673856595 },
            { -0.24294187389809563, 0.53794457510720350, -0.29747934732221337 },
            { 0.12714361127570861, -0.20217162565416691, -1.0380191039479851 },
            { -0.17244120207973623, 0.084590629352267244, -0.98061615406497193 },
            { 0.14465206170111722, 0.044064600399401299, -0.14247425194123087 } };
    }
    GivenInt(io, static_cast<int32_t>(pts0.size()));
    GivenInt(io, static_cast<int32_t>(pts1.size()));
    GivePoints3(io, pts0);
    GivePoints3(io, pts1);

    Plane3<double> plane{};
    SeparatePoints3<double, double> query{};
    bool separated = query(pts0.size(), pts0.data(), pts1.size(), pts1.data(), plane);
    io.outBool(separated);
}

// DELIBERATE DEVIATION (#348, finding 2). On the edge-edge path upstream
// writes separatingPlane.normal and .constant and leaves .origin holding
// the origin of the last face plane of hull 1 it constructed; the port
// builds the plane from normal and constant, so origin = constant * normal.
// The records are crossed tetrahedra separated (or touching, d = 0) by the
// cross product of their skew edges (RawCrossedTetrahedra, lattice or
// jittered), accepted only when the trace is sound and reaches the edge-edge
// stage, so 'separated', the normal and the constant agree bit for bit and
// only the origin differs. Residue: a record agrees when the plane passes
// through the origin (constant 0, so the port's origin is 0 * normal) and
// upstream's stale plane does too, with the same zero signs.
//
// DETERMINISM. The stale origin is that of the LAST face of hull 1 in
// ConvexHull3::GetHull() order, which for more than four points varies
// between runs (heap-address dependent; a first version of this case with
// up to 8 points in set 1 produced 7 different records of 2000 in two
// identical deep runs). Set 1 is therefore always a tetrahedron, whose hull
// order is reproducible.
ORACLE_CASE("SeparatePoints3.deviation.edgeAxisOrigin")
{
    int32_t n0 = io.rawInteger(4, 8);
    int32_t n1 = 4;
    std::vector<Vector3<double>> pts0{}, pts1{};
    SepTrace trace{ false, -1, false, Plane3<double>{}, false };
    bool accepted = false;
    for (int32_t attempt = 0; attempt < 256 && !accepted; ++attempt)
    {
        bool jitter = (io.rawInteger(0, 1) != 0);
        RawCrossedTetrahedra(io, static_cast<size_t>(n0), static_cast<size_t>(n1), jitter,
            pts0, pts1);
        if (pts1.size() != 4)
        {
            std::swap(pts0, pts1);
        }
        trace = TraceSeparate3(pts0, pts1);
        accepted = trace.sound && trace.stage == STAGE_EDGES;
    }
    if (!accepted)
    {
        // Long skew slabs (test/SeparatePoints3.test.ts): no face separates.
        pts0.clear();
        pts1.clear();
        for (size_t i = 0; i < 8; ++i)
        {
            pts0.push_back({ (i & 4) ? 4.0 : -4.0, (i & 2) ? 0.25 : -0.25, (i & 1) ? 0.25 : -0.25 });
            pts1.push_back({ (i & 1) ? 0.25 : -0.25, (i & 4) ? 4.0 : -4.0, (i & 2) ? 2.0 : 1.5 });
        }
        trace = TraceSeparate3(pts0, pts1);
    }
    GivenInt(io, static_cast<int32_t>(pts0.size()));
    GivenInt(io, static_cast<int32_t>(pts1.size()));
    GivePoints3(io, pts0);
    GivePoints3(io, pts1);

    Plane3<double> plane{};
    bool separated = RunSeparatePoints3(pts0, pts1, trace, plane);
    io.outBool(separated);
    io.outVec(plane.normal);
    io.outReal(plane.constant);
    io.outVec(plane.origin);
}

// ---- TriangulateCDT ------------------------------------------------------

namespace
{
    struct RawTree
    {
        std::vector<int32_t> polygon;
        std::vector<RawTree> child;
    };

    struct CdtInput
    {
        std::vector<Vector2<double>> pool;
        RawTree root;
    };

    int32_t const gCircle5[12][2] =
    {
        { 5, 0 }, { 4, 3 }, { 3, 4 }, { 0, 5 }, { -3, 4 }, { -4, 3 },
        { -5, 0 }, { -4, -3 }, { -3, -4 }, { 0, -5 }, { 3, -4 }, { 4, -3 }
    };

    int32_t const gSmallDirs[8][2] =
    {
        { 1, 0 }, { 1, 1 }, { 0, 1 }, { -1, 1 },
        { -1, 0 }, { -1, -1 }, { 0, -1 }, { 1, -1 }
    };

    std::vector<int32_t> AddPolygon(std::vector<Vector2<double>>& pool,
        std::vector<Vector2<double>> const& coords)
    {
        std::vector<int32_t> polygon{};
        for (auto const& p : coords)
        {
            polygon.push_back(static_cast<int32_t>(pool.size()));
            pool.push_back(p);
        }
        return polygon;
    }

    std::vector<Vector2<double>> Reversed(std::vector<Vector2<double>> coords)
    {
        std::reverse(coords.begin(), coords.end());
        return coords;
    }

    // The 12 lattice directions of the circle of radius 5 scaled by 'scale',
    // counterclockwise, about (cx, cy): convex.
    std::vector<Vector2<double>> Circle12(double cx, double cy, double scale)
    {
        std::vector<Vector2<double>> P{};
        for (int32_t k = 0; k < 12; ++k)
        {
            P.push_back({ cx + scale * gCircle5[k][0], cy + scale * gCircle5[k][1] });
        }
        return P;
    }

    // A counterclockwise star of 3 to 8 vertices about (cx, cy) in the 8
    // lattice directions scaled by 1 or 'maxScale'; its vertices are within
    // maxScale*sqrt(2) of the centre.
    std::vector<Vector2<double>> SmallStar(oracle::Ctx& io, double cx, double cy,
        int32_t maxScale)
    {
        int32_t numVertices = io.rawInteger(3, 8);
        int32_t start = io.rawInteger(0, 7);
        std::vector<Vector2<double>> ccw{};
        for (int32_t k = 0; k < numVertices; ++k)
        {
            int32_t d = (start + (k * 8) / numVertices) % 8;
            double scale = static_cast<double>(io.rawInteger(1, maxScale));
            ccw.push_back({ cx + scale * gSmallDirs[d][0], cy + scale * gSmallDirs[d][1] });
        }
        return ccw;
    }

    // A simple counterclockwise lattice polygon: an axis-aligned rectangle
    // with every lattice point of its boundary as a vertex (long collinear
    // runs), or a star in the 12 circle directions with radii 5, 10 or 15.
    std::vector<Vector2<double>> SimpleLatticePolygon(oracle::Ctx& io)
    {
        std::vector<Vector2<double>> polygon{};
        if (io.rawInteger(0, 1) == 0)
        {
            int32_t w = io.rawInteger(1, 4);
            int32_t h = io.rawInteger(1, 4);
            double ox = static_cast<double>(io.rawInteger(-3, 3));
            double oy = static_cast<double>(io.rawInteger(-3, 3));
            for (int32_t x = 0; x < w; ++x)
            {
                polygon.push_back({ ox + x, oy });
            }
            for (int32_t y = 0; y < h; ++y)
            {
                polygon.push_back({ ox + w, oy + y });
            }
            for (int32_t x = w; x > 0; --x)
            {
                polygon.push_back({ ox + x, oy + h });
            }
            for (int32_t y = h; y > 0; --y)
            {
                polygon.push_back({ ox, oy + y });
            }
            return polygon;
        }
        int32_t n = io.rawInteger(3, 12);
        int32_t start = io.rawInteger(0, 11);
        std::vector<int32_t> pick{};
        for (int32_t k = 0; k < 12 && static_cast<int32_t>(pick.size()) < n; ++k)
        {
            if (io.rawInteger(0, 11) < n)
            {
                pick.push_back((start + k) % 12);
            }
        }
        if (pick.size() < 3)
        {
            pick = { start % 12, (start + 4) % 12, (start + 8) % 12 };
        }
        // Close every angular gap wider than 150 degrees (5 directions) with
        // the direction in its middle, so that the polygon is star-shaped
        // about the origin and hence simple. (A gap wider than 180 degrees
        // produced a self-intersecting polygon, which upstream documents as
        // unsupported; it threw "Unexpected condition." on one deep-run
        // record, and the port threw too.)
        for (size_t i = 0; i < pick.size(); ++i)
        {
            int32_t next = pick[(i + 1) % pick.size()];
            int32_t gap = (next - pick[i] + 12) % 12;
            if (gap > 5)
            {
                pick.insert(pick.begin() + static_cast<std::ptrdiff_t>(i + 1),
                    (pick[i] + gap / 2) % 12);
                i = static_cast<size_t>(-1);  // rescan from the start
            }
        }
        for (int32_t k : pick)
        {
            double scale = static_cast<double>(io.rawInteger(1, 3));
            polygon.push_back({ scale * gCircle5[k][0], scale * gCircle5[k][1] });
        }
        return polygon;
    }

    // mode 0: one simple polygon, no children.
    void BuildSimple(oracle::Ctx& io, CdtInput& in)
    {
        in.root.polygon = AddPolygon(in.pool, SimpleLatticePolygon(io));
    }

    // mode 1: the 12-direction outer polygon with per-vertex radius 10 or 15
    // (every outer edge is at least 9.66 from the origin) and one hole about
    // a lattice centre of length <= 4 or two holes about (-4,0) and (4,0); a
    // hole star has radius <= 2*sqrt(2), so the holes are inside and
    // disjoint.
    void BuildHoles(oracle::Ctx& io, CdtInput& in)
    {
        std::vector<Vector2<double>> outer{};
        for (int32_t k = 0; k < 12; ++k)
        {
            double scale = static_cast<double>(io.rawInteger(2, 3));
            outer.push_back({ scale * gCircle5[k][0], scale * gCircle5[k][1] });
        }
        in.root.polygon = AddPolygon(in.pool, outer);
        if (io.rawInteger(0, 1) == 0)
        {
            double cx = static_cast<double>(io.rawInteger(-2, 2));
            double cy = static_cast<double>(io.rawInteger(-2, 2));
            RawTree hole{};
            hole.polygon = AddPolygon(in.pool, Reversed(SmallStar(io, cx, cy, 2)));
            in.root.child.push_back(hole);
        }
        else
        {
            for (double cx : { -4.0, 4.0 })
            {
                RawTree hole{};
                hole.polygon = AddPolygon(in.pool, Reversed(SmallStar(io, cx, 0.0, 2)));
                in.root.child.push_back(hole);
            }
        }
    }

    // mode 2: concentric nesting two to four levels deep: outer radius 15,
    // hole radius 10, island radius 5, a hole in the island (a star of
    // radius <= 2*sqrt(2) < 4.83), plus optionally a small hole (radius
    // <= sqrt(2)) in the annulus between 10 and 14.49 at (12,0) and a small
    // island in the annulus between 5 and 9.66 at (0,7).
    void BuildNested(oracle::Ctx& io, CdtInput& in)
    {
        int32_t depth = io.rawInteger(2, 4);
        bool extraHole = (io.rawInteger(0, 1) != 0);
        bool extraIsland = (depth >= 2 && io.rawInteger(0, 1) != 0);
        in.root.polygon = AddPolygon(in.pool, Circle12(0.0, 0.0, 3.0));
        RawTree hole{};
        hole.polygon = AddPolygon(in.pool, Reversed(Circle12(0.0, 0.0, 2.0)));
        if (depth >= 3)
        {
            RawTree island{};
            island.polygon = AddPolygon(in.pool, Circle12(0.0, 0.0, 1.0));
            if (depth >= 4)
            {
                RawTree inner{};
                inner.polygon = AddPolygon(in.pool, Reversed(SmallStar(io, 0.0, 0.0, 2)));
                island.child.push_back(inner);
            }
            hole.child.push_back(island);
        }
        if (extraIsland)
        {
            RawTree island{};
            island.polygon = AddPolygon(in.pool, SmallStar(io, 0.0, 7.0, 1));
            hole.child.push_back(island);
        }
        in.root.child.push_back(hole);
        if (extraHole)
        {
            RawTree small{};
            small.polygon = AddPolygon(in.pool, Reversed(SmallStar(io, 12.0, 0.0, 1)));
            in.root.child.push_back(small);
        }
    }

    // mode 3: coincident configurations on a small lattice, offset by a
    // lattice vector. Indices are shared explicitly where polygons share a
    // vertex.
    //   kind 0: a hole sharing a vertex (the same index) with the outer square
    //   kind 1: a hole vertex in the interior of an outer edge (the outer
    //           edge is split: vertex-edge)
    //   kind 2: a U-shaped outer polygon and a hole sharing the notch's bottom
    //           edge (edge-edge, the edge interior to the CDT hull)
    //   kind 3: two holes sharing an edge, or only a vertex
    //   kind 4: an island touching its hole at a vertex
    //   kind 5: a pinched outer polygon that visits one vertex twice
    //   kind 6: collinear vertices on both the outer polygon and the hole
    //   kind 7: the U-shaped outer polygon of kind 2 and a hole sharing part
    //           of the notch's bottom edge (the outer edge is split at the
    //           hole's vertices: vertex-edge plus edge-edge)
    void BuildCoincident(oracle::Ctx& io, CdtInput& in)
    {
        int32_t kind = io.rawInteger(0, 7);
        double ox = static_cast<double>(io.rawInteger(-3, 3));
        double oy = static_cast<double>(io.rawInteger(-3, 3));
        auto P = [&in, ox, oy](double x, double y)
        {
            in.pool.push_back({ ox + x, oy + y });
            return static_cast<int32_t>(in.pool.size() - 1);
        };
        if (kind == 0)
        {
            // Outer square (0,0)-(6,6); hole triangle with the corner (0,0).
            int32_t a = P(0, 0), b = P(6, 0), c = P(6, 6), d = P(0, 6);
            in.root.polygon = { a, b, c, d };
            double hx = static_cast<double>(io.rawInteger(2, 4));
            double hy = static_cast<double>(io.rawInteger(2, 4));
            RawTree hole{};
            int32_t e = P(hx, hy - 1.0), f = P(hx - 1.0, hy);
            hole.polygon = { a, f, e };
            in.root.child.push_back(hole);
        }
        else if (kind == 1)
        {
            // The hole's first vertex lies on the outer edge (6,0)-(6,6).
            int32_t a = P(0, 0), b = P(6, 0), c = P(6, 6), d = P(0, 6);
            in.root.polygon = { a, b, c, d };
            double y = static_cast<double>(io.rawInteger(1, 5));
            RawTree hole{};
            int32_t e = P(6, y), f = P(4, y - 1.0), g = P(4, y + 1.0);
            hole.polygon = { e, f, g };
            in.root.child.push_back(hole);
        }
        else if (kind == 2)
        {
            int32_t p0 = P(0, 0), p1 = P(9, 0), p2 = P(9, 9), p3 = P(6, 9);
            int32_t p4 = P(6, 3), p5 = P(3, 3), p6 = P(3, 9), p7 = P(0, 9);
            in.root.polygon = { p0, p1, p2, p3, p4, p5, p6, p7 };
            double apex = static_cast<double>(io.rawInteger(1, 2));
            RawTree hole{};
            int32_t q = P(static_cast<double>(io.rawInteger(3, 6)), apex);
            hole.polygon = { p5, p4, q };
            in.root.child.push_back(hole);
        }
        else if (kind == 3)
        {
            // Outer square (0,0)-(8,8); holes left and right of x = 4.
            int32_t a = P(0, 0), b = P(8, 0), c = P(8, 8), d = P(0, 8);
            in.root.polygon = { a, b, c, d };
            bool shareEdge = (io.rawInteger(0, 1) != 0);
            int32_t m0 = P(4, 3), m1 = (shareEdge ? P(4, 5) : -1);
            RawTree left{}, right{};
            int32_t l0 = P(2, 4);
            int32_t r0 = P(6, 4);
            if (shareEdge)
            {
                left.polygon = { m0, l0, m1 };      // clockwise
                right.polygon = { m1, r0, m0 };     // clockwise
            }
            else
            {
                int32_t l1 = P(2, 2), r1 = P(6, 2);
                left.polygon = { m0, l1, l0 };      // clockwise
                right.polygon = { m0, r0, r1 };     // clockwise
            }
            in.root.child.push_back(left);
            in.root.child.push_back(right);
        }
        else if (kind == 4)
        {
            // Outer square (0,0)-(10,10), square hole (2,2)-(8,8), and an
            // island triangle touching the hole at (2,2).
            int32_t a = P(0, 0), b = P(10, 0), c = P(10, 10), d = P(0, 10);
            in.root.polygon = { a, b, c, d };
            int32_t h0 = P(2, 2), h1 = P(2, 8), h2 = P(8, 8), h3 = P(8, 2);
            RawTree hole{};
            hole.polygon = { h0, h1, h2, h3 };
            double s = static_cast<double>(io.rawInteger(2, 4));
            RawTree island{};
            int32_t i1 = P(2 + s, 3), i2 = P(3, 2 + s);
            island.polygon = { h0, i1, i2 };
            hole.child.push_back(island);
            in.root.child.push_back(hole);
        }
        else if (kind == 5)
        {
            // Two squares sharing the corner m, traversed as one polygon.
            double s = static_cast<double>(io.rawInteger(1, 3));
            int32_t a = P(0, 0), b = P(s, 0), m = P(s, s), c = P(2 * s, s);
            int32_t d = P(2 * s, 2 * s), e = P(s, 2 * s), f = P(0, s);
            in.root.polygon = { a, b, m, c, d, e, m, f };
            if (io.rawInteger(0, 1) != 0 && s >= 2)
            {
                RawTree hole{};
                int32_t g = P(1, 0.5 * s), h = P(1, 0.5 * s + 0.5), j = P(1.5, 0.5 * s);
                hole.polygon = { g, h, j };     // clockwise
                in.root.child.push_back(hole);
            }
        }
        else if (kind == 7)
        {
            int32_t p0 = P(0, 0), p1 = P(9, 0), p2 = P(9, 9), p3 = P(6, 9);
            int32_t p4 = P(6, 3), p5 = P(3, 3), p6 = P(3, 9), p7 = P(0, 9);
            in.root.polygon = { p0, p1, p2, p3, p4, p5, p6, p7 };
            double a = static_cast<double>(io.rawInteger(3, 5));
            double b = a + static_cast<double>(io.rawInteger(1, static_cast<int32_t>(6.0 - a)));
            double c = static_cast<double>(io.rawInteger(3, 6));
            double apex = static_cast<double>(io.rawInteger(1, 2));
            // Reuse the outer vertex when the hole's vertex coincides with it.
            int32_t ha = (a == 3.0 ? p5 : P(a, 3));
            int32_t hb = (b == 6.0 ? p4 : P(b, 3));
            int32_t hc = P(c, apex);
            RawTree hole{};
            hole.polygon = { ha, hb, hc };      // clockwise
            in.root.child.push_back(hole);
        }
        else
        {
            std::vector<Vector2<double>> outer{}, inner{};
            int32_t w = io.rawInteger(4, 6);
            int32_t h = io.rawInteger(4, 6);
            for (int32_t x = 0; x < w; ++x) { outer.push_back({ ox + x, oy }); }
            for (int32_t y = 0; y < h; ++y) { outer.push_back({ ox + w, oy + y }); }
            for (int32_t x = w; x > 0; --x) { outer.push_back({ ox + x, oy + h }); }
            for (int32_t y = h; y > 0; --y) { outer.push_back({ ox, oy + y }); }
            // Clockwise hole (1,1)-(w-1,h-1) with its boundary lattice points.
            for (int32_t y = 1; y < h - 1; ++y) { inner.push_back({ ox + 1, oy + y }); }
            for (int32_t x = 1; x < w - 1; ++x) { inner.push_back({ ox + x, oy + h - 1 }); }
            for (int32_t y = h - 1; y > 1; --y) { inner.push_back({ ox + w - 1, oy + y }); }
            for (int32_t x = w - 1; x > 1; --x) { inner.push_back({ ox + x, oy + 1 }); }
            in.root.polygon = AddPolygon(in.pool, outer);
            RawTree hole{};
            hole.polygon = AddPolygon(in.pool, inner);
            in.root.child.push_back(hole);
        }
    }

    // mode 5: non-lattice. A star-shaped outer polygon with 6 to 9 vertices
    // at radii in [8,10] and jittered angles (consecutive angles differ by at
    // most 90 degrees, so every outer edge is at least 8*cos(45) = 5.66 from
    // the origin); optionally a hole of 6 to 8 vertices at radii in [4,5]
    // (edges at least 2.83 from the origin) and an island at radii in [1,2.5].
    std::vector<Vector2<double>> UniformStar(oracle::Ctx& io, int32_t k, double rmin,
        double rmax)
    {
        std::vector<Vector2<double>> P{};
        for (int32_t i = 0; i < k; ++i)
        {
            double angle = 6.28318530717958647692 * (i + io.raw(0.0, 0.5)) / k;
            double radius = io.raw(rmin, rmax);
            P.push_back({ radius * std::cos(angle), radius * std::sin(angle) });
        }
        return P;
    }

    void BuildUniform(oracle::Ctx& io, CdtInput& in)
    {
        int32_t depth = io.rawInteger(1, 3);
        in.root.polygon = AddPolygon(in.pool, UniformStar(io, io.rawInteger(6, 9), 8.0, 10.0));
        if (depth >= 2)
        {
            RawTree hole{};
            hole.polygon = AddPolygon(in.pool,
                Reversed(UniformStar(io, io.rawInteger(6, 8), 4.0, 5.0)));
            if (depth >= 3)
            {
                RawTree island{};
                island.polygon = AddPolygon(in.pool, UniformStar(io, io.rawInteger(3, 6), 1.0, 2.5));
                hole.child.push_back(island);
            }
            in.root.child.push_back(hole);
        }
    }

    void VisitTree(RawTree& node, std::vector<RawTree*>& out)
    {
        out.push_back(&node);
        for (auto& c : node.child)
        {
            VisitTree(c, out);
        }
    }

    // Replace about half of the polygon references by fresh pool entries
    // with equal coordinates (a zero coordinate becomes -0 half the time,
    // which std::map<Vector2> and the port's string key both treat as +0).
    void DuplicateSome(oracle::Ctx& io, CdtInput& in)
    {
        std::vector<RawTree*> nodes{};
        VisitTree(in.root, nodes);
        for (auto* node : nodes)
        {
            for (auto& index : node->polygon)
            {
                if (io.rawInteger(0, 1) == 0)
                {
                    Vector2<double> p = in.pool[static_cast<size_t>(index)];
                    for (int32_t j = 0; j < 2; ++j)
                    {
                        if (p[j] == 0.0 && io.rawInteger(0, 1) != 0)
                        {
                            p[j] = -0.0;
                        }
                    }
                    index = static_cast<int32_t>(in.pool.size());
                    in.pool.push_back(p);
                }
            }
        }
    }

    // Unused pool points and a random permutation of the pool, so that the
    // CDT's point order (first encounter in breadth-first order) differs
    // from the pool order and the restore step has real work to do.
    void ShuffleAndPad(oracle::Ctx& io, CdtInput& in)
    {
        int32_t numExtra = io.rawInteger(0, 3);
        for (int32_t k = 0; k < numExtra; ++k)
        {
            double x = static_cast<double>(io.rawInteger(-12, 12));
            double y = static_cast<double>(io.rawInteger(-12, 12));
            in.pool.push_back({ x, y });
        }
        size_t n = in.pool.size();
        std::vector<int32_t> perm(n);
        for (size_t i = 0; i < n; ++i)
        {
            perm[i] = static_cast<int32_t>(i);
        }
        for (size_t i = n - 1; i > 0; --i)
        {
            size_t j = static_cast<size_t>(io.rawInteger(0, static_cast<int32_t>(i)));
            std::swap(perm[i], perm[j]);
        }
        // Old index perm[i] moves to position i.
        std::vector<int32_t> newIndex(n);
        std::vector<Vector2<double>> pool(n);
        for (size_t i = 0; i < n; ++i)
        {
            newIndex[static_cast<size_t>(perm[i])] = static_cast<int32_t>(i);
            pool[i] = in.pool[static_cast<size_t>(perm[i])];
        }
        in.pool = pool;
        std::vector<RawTree*> nodes{};
        VisitTree(in.root, nodes);
        for (auto* node : nodes)
        {
            for (auto& index : node->polygon)
            {
                index = newIndex[static_cast<size_t>(index)];
            }
        }
    }

    // The unique points TriangulateCDT hands to ConstrainedDelaunay2: the
    // first occurrence of each coordinate pair in breadth-first node order,
    // with std::map<Vector2<double>> equality (as RemapPolygonTree).
    std::vector<Vector2<double>> CdtPoints(CdtInput const& in)
    {
        std::map<Vector2<double>, int32_t> seen{};
        std::vector<Vector2<double>> points{};
        std::vector<RawTree const*> queue{ &in.root };
        for (size_t head = 0; head < queue.size(); ++head)
        {
            for (int32_t index : queue[head]->polygon)
            {
                Vector2<double> const& p = in.pool[static_cast<size_t>(index)];
                if (seen.find(p) == seen.end())
                {
                    seen.insert(std::make_pair(p, static_cast<int32_t>(points.size())));
                    points.push_back(p);
                }
            }
            for (auto const& c : queue[head]->child)
            {
                queue.push_back(&c);
            }
        }
        return points;
    }

    // Delaunay2<T>::operator() trusts the floating-point IntrinsicsVector2
    // classification (finding #391, fixed in the port with an exact
    // predicate); ConstrainedDelaunay2 inherits it. This is v09/v11's exact
    // separator of the point sets on which that classification is sound.
    bool Sound2(std::vector<Vector2<double>> const& pts)
    {
        if (pts.size() < 3)
        {
            return false;
        }
        IntrinsicsVector2<double> info(static_cast<int32_t>(pts.size()), pts.data(), 0.0);
        if (info.dimension != 2)
        {
            return false;
        }
        int32_t sign = ExactOrient2(pts[info.extreme[0]], pts[info.extreme[1]],
            pts[info.extreme[2]]);
        if (sign == 0)
        {
            return false;
        }
        return info.extremeCCW == (sign > 0);
    }

    // A lattice square with a square hole, sound by construction.
    void FallbackTree(CdtInput& in)
    {
        in.pool.clear();
        in.root = RawTree{};
        in.root.polygon = AddPolygon(in.pool, { { 0.0, 0.0 }, { 4.0, 0.0 }, { 4.0, 4.0 },
            { 0.0, 4.0 } });
        RawTree hole{};
        hole.polygon = AddPolygon(in.pool, { { 1.0, 1.0 }, { 1.0, 3.0 }, { 3.0, 3.0 },
            { 3.0, 1.0 } });
        in.root.child.push_back(hole);
    }

    void GiveTree(oracle::Ctx& io, RawTree const& node)
    {
        GivenInt(io, static_cast<int32_t>(node.polygon.size()));
        for (int32_t index : node.polygon)
        {
            GivenInt(io, index);
        }
        GivenInt(io, static_cast<int32_t>(node.child.size()));
        for (auto const& c : node.child)
        {
            GiveTree(io, c);
        }
    }

    void GiveCdtInput(oracle::Ctx& io, CdtInput const& in)
    {
        GivenInt(io, static_cast<int32_t>(in.pool.size()));
        for (auto const& p : in.pool)
        {
            io.givenVec<2>(p);
        }
        GiveTree(io, in.root);
    }

    std::shared_ptr<PolygonTree> ToPolygonTree(RawTree const& node)
    {
        auto tree = std::make_shared<PolygonTree>();
        tree->polygon = node.polygon;
        for (auto const& c : node.child)
        {
            tree->child.push_back(ToPolygonTree(c));
        }
        return tree;
    }

    // ---- emission ----------------------------------------------------------

    void EmitTriangles(oracle::Ctx& io, std::vector<std::array<int32_t, 3>> const& tris)
    {
        io.outInt(tris.size());
        for (auto const& t : tris)
        {
            io.outInt(t[0]);
            io.outInt(t[1]);
            io.outInt(t[2]);
        }
    }

    void EmitNodeIndices(oracle::Ctx& io, std::vector<size_t> const& indices)
    {
        io.outInt(indices.size());
        for (size_t i : indices)
        {
            io.outInt(i);
        }
    }

    // allTriangles and outsideTriangles follow std::unordered_map order
    // upstream; they are sorted by their stored (restored) vertex tuple.
    void EmitSortedTriangles(oracle::Ctx& io, std::vector<std::array<int32_t, 3>> tris)
    {
        std::sort(tris.begin(), tris.end());
        EmitTriangles(io, tris);
    }

    void EmitTreeEx(oracle::Ctx& io, PolygonTreeEx const& tree)
    {
        io.outInt(tree.nodes.size());
        for (auto const& node : tree.nodes)
        {
            io.outInt(node.self);
            io.outInt(node.chirality);
            io.outInt(node.parent == std::numeric_limits<size_t>::max()
                ? static_cast<int64_t>(-1) : static_cast<int64_t>(node.parent));
            io.outInt(node.minChild);
            io.outInt(node.supChild);
            io.outInt(node.polygon.size());
            for (int32_t v : node.polygon)
            {
                io.outInt(v);
            }
            EmitTriangles(io, node.triangulation);
        }
        EmitTriangles(io, tree.interiorTriangles);
        EmitNodeIndices(io, tree.interiorNodeIndices);
        EmitTriangles(io, tree.exteriorTriangles);
        EmitNodeIndices(io, tree.exteriorNodeIndices);
        EmitTriangles(io, tree.insideTriangles);
        EmitNodeIndices(io, tree.insideNodeIndices);
        EmitSortedTriangles(io, tree.outsideTriangles);
        EmitSortedTriangles(io, tree.allTriangles);
    }

    // ---- independent reference (exact) -------------------------------------

    using E2 = std::array<Exact, 2>;

    // Twice the signed area of the polygon (exact shoelace).
    Exact TwiceArea(std::vector<E2> const& pool, std::vector<int32_t> const& polygon)
    {
        Exact sum(0);
        size_t n = polygon.size();
        for (size_t i = 0, j = n - 1; i < n; j = i++)
        {
            E2 const& a = pool[static_cast<size_t>(polygon[j])];
            E2 const& b = pool[static_cast<size_t>(polygon[i])];
            sum = sum + (a[0] * b[1] - b[0] * a[1]);
        }
        return sum;
    }

    int32_t SignE(Exact const& x)
    {
        return x.GetSign();
    }

    int32_t Orient2E(E2 const& a, E2 const& b, E2 const& c)
    {
        return SignE((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]));
    }

    // Winding number of the polygon about q, where the polygon vertices are
    // scaled by 3 and q is the sum of a triangle's three vertices (three
    // times its centroid), so everything stays exact.
    int32_t Winding3(std::vector<E2> const& pool, std::vector<int32_t> const& polygon,
        E2 const& q)
    {
        Exact three(3);
        int32_t winding = 0;
        size_t n = polygon.size();
        for (size_t i = 0, j = n - 1; i < n; j = i++)
        {
            E2 const& pa = pool[static_cast<size_t>(polygon[j])];
            E2 const& pb = pool[static_cast<size_t>(polygon[i])];
            E2 a{ three * pa[0], three * pa[1] };
            E2 b{ three * pb[0], three * pb[1] };
            if (!(q[1] < a[1]))   // a.y <= q.y
            {
                if (q[1] < b[1] && Orient2E(a, b, q) > 0)
                {
                    ++winding;
                }
            }
            else if (!(q[1] < b[1]) && Orient2E(a, b, q) < 0)
            {
                --winding;
            }
        }
        return winding;
    }

    std::array<int32_t, 3> Canonical(std::array<int32_t, 3> const& t)
    {
        // Rotate the minimum first; keeps the orientation.
        size_t k = (t[0] <= t[1] && t[0] <= t[2] ? 0 : (t[1] <= t[2] ? 1 : 2));
        return { t[k], t[(k + 1) % 3], t[(k + 2) % 3] };
    }

    // Checks on the output tree, exact:
    //  (1) every node's triangles have the node's chirality strictly, and
    //      their exact signed areas sum to the node polygon's shoelace area
    //      plus its children's (which have the opposite orientation);
    //  (2) the centroid of every node triangle is inside the node polygon
    //      (winding number == chirality) and outside every child polygon
    //      (winding number 0); every outside triangle's centroid has winding
    //      number 0 with respect to the root polygon;
    //  (3) every edge of every output polygon is an edge of allTriangles;
    //  (4) interior + exterior = inside (in order), inside and outside
    //      partition allTriangles (as sets of oriented triangles).
    bool TreeReference(std::vector<Vector2<double>> const& poolD, PolygonTreeEx const& tree)
    {
        std::vector<E2> pool(poolD.size());
        for (size_t i = 0; i < poolD.size(); ++i)
        {
            pool[i] = E2{ Exact(poolD[i][0]), Exact(poolD[i][1]) };
        }
        auto sum3 = [&pool](std::array<int32_t, 3> const& t)
        {
            E2 const& a = pool[static_cast<size_t>(t[0])];
            E2 const& b = pool[static_cast<size_t>(t[1])];
            E2 const& c = pool[static_cast<size_t>(t[2])];
            return E2{ a[0] + b[0] + c[0], a[1] + b[1] + c[1] };
        };

        std::set<std::pair<int32_t, int32_t>> allEdges{};
        std::multiset<std::array<int32_t, 3>> all{}, parts{};
        for (auto const& t : tree.allTriangles)
        {
            all.insert(Canonical(t));
            for (size_t j = 0; j < 3; ++j)
            {
                int32_t u = t[j], v = t[(j + 1) % 3];
                allEdges.insert(std::make_pair(std::min(u, v), std::max(u, v)));
            }
        }

        for (auto const& node : tree.nodes)
        {
            Exact expected = TwiceArea(pool, node.polygon);
            for (size_t c = node.minChild; c < node.supChild; ++c)
            {
                expected = expected + TwiceArea(pool, tree.nodes[c].polygon);
            }
            Exact total(0);
            for (auto const& t : node.triangulation)
            {
                E2 const& a = pool[static_cast<size_t>(t[0])];
                E2 const& b = pool[static_cast<size_t>(t[1])];
                E2 const& c = pool[static_cast<size_t>(t[2])];
                if (Orient2E(a, b, c) != node.chirality)
                {
                    return false;
                }
                total = total + ((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]));
                E2 q = sum3(t);
                if (Winding3(pool, node.polygon, q) != node.chirality)
                {
                    return false;
                }
                for (size_t ci = node.minChild; ci < node.supChild; ++ci)
                {
                    if (Winding3(pool, tree.nodes[ci].polygon, q) != 0)
                    {
                        return false;
                    }
                }
            }
            if (total != expected)
            {
                return false;
            }
            size_t n = node.polygon.size();
            for (size_t i = 0, j = n - 1; i < n; j = i++)
            {
                int32_t u = node.polygon[j], v = node.polygon[i];
                if (allEdges.find(std::make_pair(std::min(u, v), std::max(u, v))) == allEdges.end())
                {
                    return false;
                }
            }
        }

        for (auto const& t : tree.outsideTriangles)
        {
            if (Winding3(pool, tree.nodes[0].polygon, sum3(t)) != 0)
            {
                return false;
            }
            parts.insert(Canonical(t));
        }

        std::vector<std::array<int32_t, 3>> joined{};
        size_t ie = 0, ee = 0;
        for (size_t k = 0; k < tree.insideTriangles.size(); ++k)
        {
            size_t nIndex = tree.insideNodeIndices[k];
            if (tree.nodes[nIndex].chirality > 0)
            {
                if (ie >= tree.interiorTriangles.size()
                    || tree.interiorTriangles[ie] != tree.insideTriangles[k]
                    || tree.interiorNodeIndices[ie] != nIndex)
                {
                    return false;
                }
                ++ie;
            }
            else
            {
                if (ee >= tree.exteriorTriangles.size()
                    || tree.exteriorTriangles[ee] != tree.insideTriangles[k]
                    || tree.exteriorNodeIndices[ee] != nIndex)
                {
                    return false;
                }
                ++ee;
            }
            // Exterior triangles are stored clockwise; compare the
            // counterclockwise version with allTriangles.
            auto t = tree.insideTriangles[k];
            if (tree.nodes[nIndex].chirality < 0)
            {
                std::swap(t[1], t[2]);
            }
            parts.insert(Canonical(t));
        }
        if (ie != tree.interiorTriangles.size() || ee != tree.exteriorTriangles.size())
        {
            return false;
        }
        return parts == all;
    }

    // Draw a polygon tree of the given mode on which upstream's CDT seed
    // classification is sound, recording nothing. Capped, everything redrawn.
    CdtInput DrawCdtInput(oracle::Ctx& io, int32_t mode)
    {
        CdtInput in{};
        for (int32_t attempt = 0; attempt < 64; ++attempt)
        {
            in = CdtInput{};
            // Mode 4 draws the coincident shapes on three records of four:
            // only a vertex referenced twice (shared by two polygons, or
            // visited twice by one) can be referenced through two indices
            // with equal coordinates, which is what RemapPolygonTree's
            // duplicate branch needs.
            int32_t shape = mode;
            if (mode == 4)
            {
                shape = (io.rawInteger(0, 3) == 0 ? io.rawInteger(0, 2) : 3);
            }
            if (shape == 0) { BuildSimple(io, in); }
            else if (shape == 1) { BuildHoles(io, in); }
            else if (shape == 2) { BuildNested(io, in); }
            else if (shape == 3) { BuildCoincident(io, in); }
            else { BuildUniform(io, in); }
            if (mode == 4)
            {
                DuplicateSome(io, in);
            }
            ShuffleAndPad(io, in);
            if (Sound2(CdtPoints(in)))
            {
                return in;
            }
        }
        FallbackTree(in);
        return in;
    }
}

// TriangulateCDT<double>: both operator() overloads (the vector overload on
// even records, the pointer overload on odd ones; the port has one
// compute(points, tree)), and every member of the output PolygonTreeEx:
// per node self, chirality, parent (-1 for the root), minChild, supChild,
// the output polygon (with edge splits) and the triangulation; interior,
// exterior and inside triangles with their node indices; outside triangles
// and allTriangles (sorted, see the header comment).
//
//   mode 0: one simple lattice polygon (rectangles with collinear runs,
//           12-direction stars)
//   mode 1: an outer polygon with one or two holes
//   mode 2: concentric nesting 2 to 4 levels deep, extra hole and island
//   mode 3: coincident configurations (shared vertices and edges, a vertex
//           on another polygon's edge, partly shared edges, pinched
//           polygon, collinear vertices)
//   mode 4: mostly mode 3, else modes 0-2, with about half of the polygon
//           references replaced by duplicated pool points, so that a shared
//           vertex is referenced through two indices with equal coordinates
//           (RemapPolygonTree's duplicate branch and its overwrite quirk,
//           #348 finding 3, preserved), including -0 for +0
// A shared edge of two polygons must not lie on the convex hull of the
// referenced points: upstream then throws, which has its own case below.
//   mode 5: non-lattice star polygons (outer, hole, island)
// Every record has 0-3 unused pool points and a shuffled pool.
//
// The unique points given to ConstrainedDelaunay2 must pass Sound2 (the
// #391 restriction inherited from Delaunay2); CdtPoints reproduces
// RemapPolygonTree's order exactly.
//
// The last output is an independent exact reference check computed by each
// side on its own output (TreeReference).
ORACLE_CASE("TriangulateCDT.compute")
{
    int32_t mode = io.index() % 6;
    CdtInput in = DrawCdtInput(io, mode);
    GiveCdtInput(io, in);

    auto inputTree = ToPolygonTree(in.root);
    PolygonTreeEx outputTree{};
    TriangulateCDT<double> triangulator{};
    if (io.index() % 2 == 0)
    {
        triangulator(in.pool, inputTree, outputTree);
    }
    else
    {
        triangulator(in.pool.size(), in.pool.data(), inputTree, outputTree);
    }
    EmitTreeEx(io, outputTree);
    io.outBool(TreeReference(in.pool, outputTree));
}

// Throw parity for TriangulateCDT's input asserts:
//   kind 0: fewer than 3 input points ("Invalid argument.")
//   kind 1: at least 3 input points but fewer than 3 distinct points
//           referenced by the tree ("Invalid polygon tree."): the polygon
//           repeats indices, or uses distinct indices of equal points
// Every record is a throw record.
ORACLE_CASE("TriangulateCDT.compute.invalidThrows")
{
    int32_t kind = io.index() % 2;
    CdtInput in{};
    if (kind == 0)
    {
        int32_t n = io.rawInteger(1, 2);
        for (int32_t i = 0; i < n; ++i)
        {
            double x = static_cast<double>(io.rawInteger(-3, 3));
            double y = static_cast<double>(io.rawInteger(-3, 3));
            in.pool.push_back({ x, y });
        }
        int32_t m = io.rawInteger(3, 4);
        for (int32_t i = 0; i < m; ++i)
        {
            in.root.polygon.push_back(io.rawInteger(0, n - 1));
        }
    }
    else
    {
        int32_t n = io.rawInteger(3, 6);
        double x0 = static_cast<double>(io.rawInteger(-3, 3));
        double y0 = static_cast<double>(io.rawInteger(-3, 3));
        double x1 = x0 + static_cast<double>(io.rawInteger(1, 3));
        double y1 = y0 + static_cast<double>(io.rawInteger(-3, 3));
        for (int32_t i = 0; i < n; ++i)
        {
            if (io.rawInteger(0, 1) != 0)
            {
                in.pool.push_back({ x0, y0 });
            }
            else
            {
                in.pool.push_back({ x1, y1 });
            }
        }
        // At most two distinct coordinates can be referenced.
        int32_t m = io.rawInteger(3, 5);
        for (int32_t i = 0; i < m; ++i)
        {
            in.root.polygon.push_back(io.rawInteger(0, n - 1));
        }
    }
    GiveCdtInput(io, in);

    auto inputTree = ToPolygonTree(in.root);
    PolygonTreeEx outputTree{};
    TriangulateCDT<double> triangulator{};
    triangulator(in.pool, inputTree, outputTree);
    EmitTreeEx(io, outputTree);
}

// Throw parity for a documented-supported input that upstream rejects (new
// suspect, group report): a hole that shares an edge, or part of an edge,
// with the outer polygon where that edge lies on the convex hull of the
// referenced points. ClassifyDFS processes the hole first and extracts its
// triangles; ETManifoldMesh::Remove then deletes the shared edge, which had
// no triangle on its other side, and the outer polygon's pass asserts
// "Unexpected condition." at 'emap.find(ekey) != emap.end()'. The port
// preserves the behavior (graph.getEdge returns null, same assert). The
// header says "The algorithm supports coincident vertex-edge and coincident
// edge-edge configurations"; the same configuration with the shared edge
// inside the hull (TriangulateCDT.compute, mode 3 kinds 2 and 7) succeeds.
//
// The outer polygon is a lattice rectangle w x h (its boundary lattice
// points as vertices on half the records); the hole is a clockwise
// triangle or quadrilateral with an edge on the rectangle's left side, from
// y = a to y = b, 0 <= a < b <= h. The whole configuration is rotated by a
// multiple of 90 degrees (exactly). Every record is a throw record.
ORACLE_CASE("TriangulateCDT.compute.hullSharedEdgeThrows")
{
    CdtInput in{};
    int32_t w = io.rawInteger(2, 5);
    int32_t h = io.rawInteger(2, 5);
    bool runs = (io.rawInteger(0, 1) != 0);
    int32_t a = io.rawInteger(0, h - 1);
    int32_t b = io.rawInteger(a + 1, h);
    int32_t depth = io.rawInteger(1, w);
    bool quad = (io.rawInteger(0, 1) != 0);
    int32_t rotation = io.rawInteger(0, 3);
    double ox = static_cast<double>(io.rawInteger(-3, 3));
    double oy = static_cast<double>(io.rawInteger(-3, 3));
    auto rotate = [rotation, ox, oy](double x, double y)
    {
        for (int32_t r = 0; r < rotation; ++r)
        {
            double t = x;
            x = -y;
            y = t;
        }
        return Vector2<double>{ ox + x, oy + y };
    };
    std::vector<Vector2<double>> outer{};
    if (runs)
    {
        for (int32_t x = 0; x < w; ++x) { outer.push_back(rotate(x, 0)); }
        for (int32_t y = 0; y < h; ++y) { outer.push_back(rotate(w, y)); }
        for (int32_t x = w; x > 0; --x) { outer.push_back(rotate(x, h)); }
        for (int32_t y = h; y > 0; --y) { outer.push_back(rotate(0, y)); }
    }
    else
    {
        outer = { rotate(0, 0), rotate(w, 0), rotate(w, h), rotate(0, h) };
    }
    in.root.polygon = AddPolygon(in.pool, outer);
    std::vector<Vector2<double>> hole{ rotate(0, a), rotate(0, b) };
    double mid = 0.5 * (a + b);
    if (quad)
    {
        hole.push_back(rotate(depth, b));
        hole.push_back(rotate(depth, a));
    }
    else
    {
        hole.push_back(rotate(depth, mid));
    }
    RawTree child{};
    child.polygon = AddPolygon(in.pool, hole);   // clockwise
    in.root.child.push_back(child);
    ShuffleAndPad(io, in);
    if (!Sound2(CdtPoints(in)))
    {
        // Not expected for these lattice sets (the #391 restriction).
        in = CdtInput{};
        in.root.polygon = AddPolygon(in.pool, { { 0.0, 0.0 }, { 4.0, 0.0 }, { 4.0, 3.0 },
            { 0.0, 3.0 } });
        RawTree fallback{};
        fallback.polygon = AddPolygon(in.pool, { { 0.0, 1.0 }, { 0.0, 2.0 }, { 1.0, 1.5 } });
        in.root.child.push_back(fallback);
    }
    GiveCdtInput(io, in);

    auto inputTree = ToPolygonTree(in.root);
    PolygonTreeEx outputTree{};
    TriangulateCDT<double> triangulator{};
    triangulator(in.pool, inputTree, outputTree);
    EmitTreeEx(io, outputTree);
}
