// Verify group 40 (numerical): differential cases for LinearSystem.h,
// GaussNewtonMinimizer.h and LevenbergMarquardtMinimizer.h.
//
// Everything in this family is + - * / and sqrt: closed-form 2x2/3x3/4x4
// inverses, Gaussian elimination with full pivoting, the tridiagonal
// recurrences, conjugate gradient (sqrt only, in the residual test) and the
// Cholesky-based normal-equation solves of the two minimizers. The residual
// functions F and the Jacobians J that drive the minimizers are written here
// and in the TypeScript replay as the same polynomial expressions in the same
// association, so no case in this file needs a tolerance: every case is
// declared { exact: true }.
//
// Generator notes:
//  * Matrices are drawn in three modes - uniform, small integer lattice and
//    an aimed singular construction (one row a linear combination of the
//    others) - so the "not invertible" branch of every solver is reached.
//    Every mode records the same number of doubles in the same (row-major)
//    order, so the replay does not need to know the mode.
//  * The tridiagonal solvers get aimed modes for diagonal[0] == 0 and for an
//    exactly zero pivot in the recurrence.
//  * The conjugate gradient cases vary the tolerance over 0, a realistic
//    1e-8 and a huge value that stops the loop at the first test, and vary
//    maxIterations over 0..6, so the returned iteration count covers the
//    "never converged" (maxIterations + 1) and "converged" outcomes.
//
// LevenbergMarquardtMinimizer::DoIteration builds -J^T*F from the member mF,
// which holds F at the previously *rejected* candidate whenever the inner
// lambda-adjustment loop, or the Gauss-Newton fallback after it, runs
// DoIteration more than once for the same pCurrent. The port re-evaluates F
// at pCurrent (issue #261, "fixed"). Upstream and the port therefore agree
// exactly on the runs in which no outer iteration repeats a DoIteration, and
// that is exactly the runs whose result.numAdjustments is 0 for every prefix
// of the iteration. SoundLMIterations below runs upstream itself with
// maxIterations = 1, 2, ... and keeps the longest sound prefix (main cases)
// or the first unsound one (the deviation case).
#define ORACLE_FAMILY "v40-numerical"
#include "Oracle.h"

#include <Mathematics/GaussNewtonMinimizer.h>
#include <Mathematics/LevenbergMarquardtMinimizer.h>
#include <Mathematics/LinearSystem.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <map>
#include <vector>

using namespace gte;

namespace
{
    // ---- shared matrix generators -----------------------------------------
    // Mode 0 uniform, mode 1 integer lattice, mode 2 an exactly singular
    // matrix whose last row is a linear combination of the earlier rows.
    // Every mode records N*N doubles in row-major order.

    void FillMatrix(oracle::Ctx& io, int32_t n, int mode, std::vector<double>& a)
    {
        a.resize(static_cast<size_t>(n) * static_cast<size_t>(n));
        if (mode == 0)
        {
            for (size_t i = 0; i < a.size(); ++i) { a[i] = io.real(-3.0, 3.0); }
        }
        else if (mode == 1)
        {
            for (size_t i = 0; i < a.size(); ++i) { a[i] = io.lattice(-2, 2); }
        }
        else
        {
            std::vector<double> raw(a.size());
            for (size_t i = 0; i < raw.size(); ++i) { raw[i] = io.raw(-3.0, 3.0); }
            std::vector<double> coeff(static_cast<size_t>(n) - 1);
            for (size_t i = 0; i < coeff.size(); ++i) { coeff[i] = io.raw(-2.0, 2.0); }
            int32_t last = n - 1;
            for (int32_t c = 0; c < n; ++c)
            {
                double sum = 0.0;
                for (int32_t r = 0; r < last; ++r)
                {
                    sum += coeff[static_cast<size_t>(r)]
                        * raw[static_cast<size_t>(r) * static_cast<size_t>(n)
                            + static_cast<size_t>(c)];
                }
                raw[static_cast<size_t>(last) * static_cast<size_t>(n)
                    + static_cast<size_t>(c)] = sum;
            }
            for (size_t i = 0; i < raw.size(); ++i) { a[i] = io.given(raw[i]); }
        }
    }

    // A symmetric positive-definite NxN matrix on a small lattice: the
    // off-diagonal entries are integers and the diagonal is made dominant.
    // Recorded row-major, N*N doubles.
    void FillSymmetric(oracle::Ctx& io, int32_t n, std::vector<double>& a)
    {
        a.assign(static_cast<size_t>(n) * static_cast<size_t>(n), 0.0);
        for (int32_t r = 0; r < n; ++r)
        {
            for (int32_t c = r; c < n; ++c)
            {
                double value = (r == c
                    ? static_cast<double>(io.rawInteger(1, 3)) + 2.0 * n
                    : static_cast<double>(io.rawInteger(-2, 2)));
                a[static_cast<size_t>(r) * static_cast<size_t>(n)
                    + static_cast<size_t>(c)] = value;
                a[static_cast<size_t>(c) * static_cast<size_t>(n)
                    + static_cast<size_t>(r)] = value;
            }
        }
        for (size_t i = 0; i < a.size(); ++i) { a[i] = io.given(a[i]); }
    }

    // The conjugate gradient tolerance, one of four regimes.
    double CGTolerance(oracle::Ctx& io, int mode)
    {
        if (mode == 0) { return io.given(0.0); }
        if (mode == 1) { return io.given(1e-8); }
        if (mode == 2) { return io.real(1e-3, 1e-1); }
        return io.given(1e12);
    }
}

// ---- LinearSystem: the closed-form fixed-size solvers -----------------------
// Inverse(A, &invertible) resolves to the closed-form overload of
// Matrix2x2.h / Matrix3x3.h / Matrix4x4.h, which LinearSystem.h includes, not
// to the Gaussian-elimination template of Matrix.h.

ORACLE_CASE("LinearSystem.solve.2x2")
{
    int mode = io.index() % 3;
    std::vector<double> a;
    FillMatrix(io, 2, mode, a);
    Vector2<double> b = io.vec<2>(-3.0, 3.0);
    Matrix2x2<double> A{};
    for (int32_t r = 0; r < 2; ++r)
    {
        for (int32_t c = 0; c < 2; ++c)
        {
            A(r, c) = a[static_cast<size_t>(2 * r + c)];
        }
    }
    Vector2<double> x{};
    bool invertible = LinearSystem<double>::Solve(A, b, x);
    io.outBool(invertible);
    io.outVec(x);
}

ORACLE_CASE("LinearSystem.solve.3x3")
{
    int mode = io.index() % 3;
    std::vector<double> a;
    FillMatrix(io, 3, mode, a);
    Vector3<double> b = io.vec<3>(-3.0, 3.0);
    Matrix3x3<double> A{};
    for (int32_t r = 0; r < 3; ++r)
    {
        for (int32_t c = 0; c < 3; ++c)
        {
            A(r, c) = a[static_cast<size_t>(3 * r + c)];
        }
    }
    Vector3<double> x{};
    bool invertible = LinearSystem<double>::Solve(A, b, x);
    io.outBool(invertible);
    io.outVec(x);
}

ORACLE_CASE("LinearSystem.solve.4x4")
{
    int mode = io.index() % 3;
    std::vector<double> a;
    FillMatrix(io, 4, mode, a);
    Vector4<double> b = io.vec<4>(-3.0, 3.0);
    Matrix4x4<double> A{};
    for (int32_t r = 0; r < 4; ++r)
    {
        for (int32_t c = 0; c < 4; ++c)
        {
            A(r, c) = a[static_cast<size_t>(4 * r + c)];
        }
    }
    Vector4<double> x{};
    bool invertible = LinearSystem<double>::Solve(A, b, x);
    io.outBool(invertible);
    io.outVec(x);
}

// ---- LinearSystem: Gaussian elimination ------------------------------------
// X is emitted on the singular path too: GaussianElimination zero-fills it
// before returning false, so the value is defined and is compared bit for bit.

ORACLE_CASE("LinearSystem.solve.nxn")
{
    int mode = io.index() % 3;
    int32_t n = static_cast<int32_t>(io.integer(2, 5));
    std::vector<double> a;
    FillMatrix(io, n, mode, a);
    std::vector<double> b(static_cast<size_t>(n));
    for (int32_t i = 0; i < n; ++i) { b[static_cast<size_t>(i)] = io.real(-3.0, 3.0); }
    std::vector<double> x(static_cast<size_t>(n), 0.0);
    bool invertible = LinearSystem<double>::Solve(n, a.data(), b.data(), x.data());
    io.outBool(invertible);
    for (int32_t i = 0; i < n; ++i) { io.outReal(x[static_cast<size_t>(i)]); }
}

ORACLE_CASE("LinearSystem.solve.nxm")
{
    int mode = io.index() % 3;
    int32_t n = static_cast<int32_t>(io.integer(2, 5));
    int32_t m = static_cast<int32_t>(io.integer(1, 3));
    std::vector<double> a;
    FillMatrix(io, n, mode, a);
    std::vector<double> b(static_cast<size_t>(n) * static_cast<size_t>(m));
    for (size_t i = 0; i < b.size(); ++i) { b[i] = io.real(-3.0, 3.0); }
    std::vector<double> y(b.size(), 0.0);
    bool invertible = LinearSystem<double>::Solve(n, m, a.data(), b.data(), y.data());
    io.outBool(invertible);
    for (size_t i = 0; i < y.size(); ++i) { io.outReal(y[i]); }
}

// ---- LinearSystem: the tridiagonal solvers ---------------------------------
// X is emitted on the failure path too. Upstream's X is caller-allocated and
// the function writes into it in a fixed order before returning false, so the
// partially written array is deterministic; both sides start from a
// zero-filled X of length N. Mode 2 forces diagonal[0] == 0 (the early
// return) and mode 3 forces an exactly zero pivot at the first step of the
// recurrence.

ORACLE_CASE("LinearSystem.solveTridiagonal")
{
    int mode = io.index() % 4;
    int32_t n = static_cast<int32_t>(io.integer(2, 6));
    std::vector<double> sub(static_cast<size_t>(n) - 1);
    std::vector<double> diag(static_cast<size_t>(n));
    std::vector<double> super(static_cast<size_t>(n) - 1);
    if (mode == 0)
    {
        for (size_t i = 0; i < sub.size(); ++i) { sub[i] = io.real(-2.0, 2.0); }
        for (size_t i = 0; i < diag.size(); ++i) { diag[i] = io.real(2.0, 5.0); }
        for (size_t i = 0; i < super.size(); ++i) { super[i] = io.real(-2.0, 2.0); }
    }
    else if (mode == 1)
    {
        for (size_t i = 0; i < sub.size(); ++i) { sub[i] = io.lattice(-2, 2); }
        for (size_t i = 0; i < diag.size(); ++i) { diag[i] = io.lattice(1, 4); }
        for (size_t i = 0; i < super.size(); ++i) { super[i] = io.lattice(-2, 2); }
    }
    else if (mode == 2)
    {
        for (size_t i = 0; i < sub.size(); ++i) { sub[i] = io.lattice(-2, 2); }
        diag[0] = io.given(0.0);
        for (size_t i = 1; i < diag.size(); ++i) { diag[i] = io.lattice(-3, 3); }
        for (size_t i = 0; i < super.size(); ++i) { super[i] = io.lattice(-2, 2); }
    }
    else
    {
        // diagonal[0] = 1 makes invExpr exactly 1, so tmp[0] = superdiagonal[0]
        // and the next pivot is diagonal[1] - subdiagonal[0]*superdiagonal[0],
        // which is exactly zero for the recorded diagonal[1].
        for (size_t i = 0; i < sub.size(); ++i) { sub[i] = io.lattice(-2, 2); }
        std::vector<double> rest(diag.size() - 2);
        for (size_t i = 0; i < rest.size(); ++i) { rest[i] = io.raw(-3.0, 3.0); }
        std::vector<double> sup(super.size());
        for (size_t i = 0; i < sup.size(); ++i)
        {
            sup[i] = static_cast<double>(io.rawInteger(-2, 2));
        }
        diag[0] = io.given(1.0);
        diag[1] = io.given(sub[0] * sup[0]);
        for (size_t i = 2; i < diag.size(); ++i) { diag[i] = io.given(rest[i - 2]); }
        for (size_t i = 0; i < super.size(); ++i) { super[i] = io.given(sup[i]); }
    }
    std::vector<double> b(static_cast<size_t>(n));
    for (int32_t i = 0; i < n; ++i) { b[static_cast<size_t>(i)] = io.real(-3.0, 3.0); }
    std::vector<double> x(static_cast<size_t>(n), 0.0);
    bool solved = LinearSystem<double>::SolveTridiagonal(n, sub.data(), diag.data(),
        super.data(), b.data(), x.data());
    io.outBool(solved);
    for (int32_t i = 0; i < n; ++i) { io.outReal(x[static_cast<size_t>(i)]); }
}

ORACLE_CASE("LinearSystem.solveConstantTridiagonal")
{
    int mode = io.index() % 4;
    int32_t n = static_cast<int32_t>(io.integer(2, 6));
    double sub = 0.0, diag = 0.0, super = 0.0;
    if (mode == 0)
    {
        sub = io.real(-2.0, 2.0);
        diag = io.real(2.0, 5.0);
        super = io.real(-2.0, 2.0);
    }
    else if (mode == 1)
    {
        sub = io.lattice(-2, 2);
        diag = io.lattice(1, 4);
        super = io.lattice(-2, 2);
    }
    else if (mode == 2)
    {
        sub = io.lattice(-2, 2);
        diag = io.given(0.0);
        super = io.lattice(-2, 2);
    }
    else
    {
        // subdiagonal = diagonal = superdiagonal = c makes the first pivot
        // c - c*(c/c) exactly zero.
        double c = static_cast<double>(io.rawInteger(1, 4));
        sub = io.given(c);
        diag = io.given(c);
        super = io.given(c);
    }
    std::vector<double> b(static_cast<size_t>(n));
    for (int32_t i = 0; i < n; ++i) { b[static_cast<size_t>(i)] = io.real(-3.0, 3.0); }
    std::vector<double> x(static_cast<size_t>(n), 0.0);
    bool solved = LinearSystem<double>::SolveConstantTridiagonal(n, sub, diag,
        super, b.data(), x.data());
    io.outBool(solved);
    for (int32_t i = 0; i < n; ++i) { io.outReal(x[static_cast<size_t>(i)]); }
}

// ---- LinearSystem: conjugate gradient --------------------------------------

ORACLE_CASE("LinearSystem.solveSymmetricCG.dense")
{
    int mode = io.index() % 4;
    int32_t n = static_cast<int32_t>(io.integer(2, 5));
    std::vector<double> a;
    FillSymmetric(io, n, a);
    std::vector<double> b(static_cast<size_t>(n));
    for (int32_t i = 0; i < n; ++i) { b[static_cast<size_t>(i)] = io.real(-3.0, 3.0); }
    int maxIterations = io.integer(0, 6);
    double tolerance = CGTolerance(io, mode);
    std::vector<double> x(static_cast<size_t>(n), 0.0);
    uint32_t iterations = LinearSystem<double>::SolveSymmetricCG(n, a.data(),
        b.data(), x.data(), static_cast<uint32_t>(maxIterations), tolerance);
    io.outInt(iterations);
    for (int32_t i = 0; i < n; ++i) { io.outReal(x[static_cast<size_t>(i)]); }
}

// The sparse form stores only one of (i,j) and (j,i). The presence of each
// off-diagonal entry is recorded as a 0/1 integer, so the replay reads the
// same fixed layout without knowing the generator mode.
ORACLE_CASE("LinearSystem.solveSymmetricCG.sparse")
{
    int mode = io.index() % 4;
    int32_t n = static_cast<int32_t>(io.integer(2, 5));
    LinearSystem<double>::SparseMatrix A{};
    for (int32_t r = 0; r < n; ++r)
    {
        double value = io.given(static_cast<double>(io.rawInteger(1, 3)) + 2.0 * n);
        A[{ r, r }] = value;
    }
    for (int32_t r = 0; r < n; ++r)
    {
        for (int32_t c = r + 1; c < n; ++c)
        {
            bool present = io.boolean();
            double value = io.lattice(-2, 2);
            if (present) { A[{ r, c }] = value; }
        }
    }
    std::vector<double> b(static_cast<size_t>(n));
    for (int32_t i = 0; i < n; ++i) { b[static_cast<size_t>(i)] = io.real(-3.0, 3.0); }
    int maxIterations = io.integer(0, 6);
    double tolerance = CGTolerance(io, mode);
    std::vector<double> x(static_cast<size_t>(n), 0.0);
    uint32_t iterations = LinearSystem<double>::SolveSymmetricCG(n, A, b.data(),
        x.data(), static_cast<uint32_t>(maxIterations), tolerance);
    io.outInt(iterations);
    for (int32_t i = 0; i < n; ++i) { io.outReal(x[static_cast<size_t>(i)]); }
}

// A preserved upstream defect (issue #261): with B = 0 the first iteration
// computes alpha = 0/0 and overwrites the exact solution X = 0 with NaN. The
// iteration count stays meaningful, and the NaN pattern is emitted as
// booleans because the harness treats any NaN as equal to any NaN.
ORACLE_CASE("LinearSystem.solveSymmetricCG.zeroRHS")
{
    int32_t n = static_cast<int32_t>(io.integer(2, 5));
    std::vector<double> a;
    FillSymmetric(io, n, a);
    std::vector<double> b(static_cast<size_t>(n));
    for (int32_t i = 0; i < n; ++i) { b[static_cast<size_t>(i)] = io.given(0.0); }
    int maxIterations = io.integer(0, 4);
    std::vector<double> x(static_cast<size_t>(n), 0.0);
    uint32_t iterations = LinearSystem<double>::SolveSymmetricCG(n, a.data(),
        b.data(), x.data(), static_cast<uint32_t>(maxIterations), 1e-8);
    io.outInt(iterations);
    for (int32_t i = 0; i < n; ++i)
    {
        io.outBool(x[static_cast<size_t>(i)] != x[static_cast<size_t>(i)]);
    }
}

// ---- the minimizers: shared residual problems ------------------------------
// Every F and J below is a polynomial in the parameters, written here and in
// the TypeScript replay as the same expression in the same association, so
// the whole minimizer path (Cholesky factor, the two triangular solves, the
// normal equations) is + - * / and sqrt.

namespace
{
    // Rosenbrock: F = (a*(p1 - p0^2), b - p0), minimized at (b, b^2) with
    // F = 0. The Gauss-Newton step from a random start frequently increases
    // the error, which is what drives the Levenberg-Marquardt adjustments.
    struct Rosenbrock
    {
        double a, b;

        void Eval(GVector<double> const& p, GVector<double>& f) const
        {
            f[0] = a * (p[1] - p[0] * p[0]);
            f[1] = b - p[0];
        }

        void Jacobian(GVector<double> const& p, GMatrix<double>& j) const
        {
            j(0, 0) = -2.0 * a * p[0];
            j(0, 1) = a;
            j(1, 0) = -1.0;
            j(1, 1) = 0.0;
        }
    };

    // Draws a, b and the initial guess; 6 recorded doubles.
    Rosenbrock MakeRosenbrock(oracle::Ctx& io, GVector<double>& p0)
    {
        Rosenbrock problem{};
        problem.a = io.real(5.0, 15.0);
        problem.b = io.real(0.5, 1.5);
        p0.SetSize(2);
        p0[0] = io.real(-3.0, 3.0);
        p0[1] = io.real(-3.0, 3.0);
        return problem;
    }

    // Circle fit by algebraic residuals: F_i = |Q_i - C|^2 - r^2 with the
    // parameters p = (C.x, C.y, r).
    struct CircleFit
    {
        std::vector<Vector2<double>> points;

        void Eval(GVector<double> const& p, GVector<double>& f) const
        {
            for (size_t i = 0; i < points.size(); ++i)
            {
                double dx = points[i][0] - p[0];
                double dy = points[i][1] - p[1];
                f[static_cast<int32_t>(i)] = dx * dx + dy * dy - p[2] * p[2];
            }
        }

        void Jacobian(GVector<double> const& p, GMatrix<double>& j) const
        {
            for (size_t i = 0; i < points.size(); ++i)
            {
                double dx = points[i][0] - p[0];
                double dy = points[i][1] - p[1];
                int32_t r = static_cast<int32_t>(i);
                j(r, 0) = -2.0 * dx;
                j(r, 1) = -2.0 * dy;
                j(r, 2) = -2.0 * p[2];
            }
        }

        // The normal equations, accumulated point by point. This is the
        // "J plus" formulation and is deliberately a different accumulation
        // order from MultiplyATB, which sweeps (r, c) outside and the points
        // inside.
        void Normal(GVector<double> const& p, GMatrix<double>& jtj,
            GVector<double>& negJTF) const
        {
            for (int32_t r = 0; r < 3; ++r)
            {
                negJTF[r] = 0.0;
                for (int32_t c = 0; c < 3; ++c) { jtj(r, c) = 0.0; }
            }
            for (size_t i = 0; i < points.size(); ++i)
            {
                double dx = points[i][0] - p[0];
                double dy = points[i][1] - p[1];
                double fi = dx * dx + dy * dy - p[2] * p[2];
                double d[3] = { -2.0 * dx, -2.0 * dy, -2.0 * p[2] };
                for (int32_t r = 0; r < 3; ++r)
                {
                    for (int32_t c = 0; c < 3; ++c)
                    {
                        jtj(r, c) += d[r] * d[c];
                    }
                    negJTF[r] -= fi * d[r];
                }
            }
        }
    };

    // Mode 0 lattice points, mode 1 uniform, mode 2 points on a circle whose
    // parameters are dyadic (so the sample points are exact) with a small
    // dyadic perturbation. Records 1 + 2*n doubles in every mode, then the
    // 3 components of the initial guess.
    CircleFit MakeCircleFit(oracle::Ctx& io, int mode, GVector<double>& p0)
    {
        CircleFit problem{};
        int n = io.integer(5, 8);
        problem.points.resize(static_cast<size_t>(n));
        if (mode == 0)
        {
            for (int i = 0; i < n; ++i) { problem.points[static_cast<size_t>(i)] = io.latticeVec<2>(-4, 4); }
        }
        else if (mode == 1)
        {
            for (int i = 0; i < n; ++i) { problem.points[static_cast<size_t>(i)] = io.vec<2>(-4.0, 4.0); }
        }
        else
        {
            double cx = static_cast<double>(io.rawInteger(-4, 4)) * 0.5;
            double cy = static_cast<double>(io.rawInteger(-4, 4)) * 0.5;
            for (int i = 0; i < n; ++i)
            {
                // A Pythagorean-like offset on a dyadic grid keeps the points
                // exactly representable without calling sin or cos.
                double ox = static_cast<double>(io.rawInteger(-8, 8)) * 0.25;
                double oy = static_cast<double>(io.rawInteger(-8, 8)) * 0.25;
                Vector2<double> q{ cx + ox, cy + oy };
                problem.points[static_cast<size_t>(i)] = io.givenVec(q);
            }
        }
        p0.SetSize(3);
        p0[0] = io.real(-2.0, 2.0);
        p0[1] = io.real(-2.0, 2.0);
        p0[2] = io.real(0.5, 4.0);
        return problem;
    }

    // A linear model y = p0*t + p1 whose Jacobian columns are chosen to be
    // full rank (mode 2), duplicated (mode 1) or one of them identically zero
    // (mode 0). Modes 0 and 1 make J^T*J singular, so CholeskyDecomposition::
    // Factor fails and the minimizer returns on its first iteration.
    struct LinearModel
    {
        int mode;
        std::vector<double> t, y;

        void Eval(GVector<double> const& p, GVector<double>& f) const
        {
            for (size_t i = 0; i < t.size(); ++i)
            {
                double model = (mode == 0 ? p[0] * t[i]
                    : mode == 1 ? (p[0] + p[1]) * t[i]
                    : p[0] * t[i] + p[1]);
                f[static_cast<int32_t>(i)] = model - y[i];
            }
        }

        void Jacobian(GVector<double> const&, GMatrix<double>& j) const
        {
            for (size_t i = 0; i < t.size(); ++i)
            {
                int32_t r = static_cast<int32_t>(i);
                j(r, 0) = t[i];
                j(r, 1) = (mode == 0 ? 0.0 : mode == 1 ? t[i] : 1.0);
            }
        }
    };

    // Records 1 + 2*n + 2 doubles. The samples are on a lattice, so in mode 2
    // the residual is exactly zero at the generating parameters; atSolution
    // starts the minimizer there.
    LinearModel MakeLinearModel(oracle::Ctx& io, int mode, GVector<double>& p0,
        bool atSolution = false)
    {
        LinearModel problem{};
        problem.mode = mode;
        int n = io.integer(3, 6);
        problem.t.resize(static_cast<size_t>(n));
        problem.y.resize(static_cast<size_t>(n));
        double slope = static_cast<double>(io.rawInteger(-3, 3));
        double intercept = static_cast<double>(io.rawInteger(-3, 3));
        for (int i = 0; i < n; ++i)
        {
            problem.t[static_cast<size_t>(i)] = io.lattice(-4, 4);
        }
        for (int i = 0; i < n; ++i)
        {
            double ti = problem.t[static_cast<size_t>(i)];
            problem.y[static_cast<size_t>(i)] = io.given(slope * ti + intercept);
        }
        p0.SetSize(2);
        if (atSolution)
        {
            p0[0] = io.given(slope);
            p0[1] = io.given(intercept);
        }
        else
        {
            p0[0] = io.lattice(-3, 3);
            p0[1] = io.lattice(-3, 3);
        }
        return problem;
    }

    // The number of Levenberg-Marquardt iterations on these inputs for which
    // upstream never repeats a DoIteration for the same pCurrent, capped at
    // maxCap (wantSound), or the first prefix length at which it does
    // (!wantSound). See the issue #261 note at the top of the file.
    size_t SoundLMIterations(LevenbergMarquardtMinimizer<double>& minimizer,
        GVector<double> const& p0, double updateLengthTolerance,
        double errorDifferenceTolerance, double lambdaFactor,
        double lambdaAdjust, size_t maxAdjustments, size_t maxCap,
        bool wantSound)
    {
        size_t chosen = (wantSound ? 0 : maxCap);
        for (size_t m = 1; m <= maxCap; ++m)
        {
            auto probe = minimizer(p0, m, updateLengthTolerance,
                errorDifferenceTolerance, lambdaFactor, lambdaAdjust,
                maxAdjustments);
            if (wantSound)
            {
                if (probe.numAdjustments != 0) { break; }
                chosen = m;
            }
            else if (probe.numAdjustments != 0)
            {
                chosen = m;
                break;
            }
        }
        return chosen;
    }

    void EmitGN(oracle::Ctx& io, GaussNewtonMinimizer<double>::Result const& r)
    {
        io.outInt(r.numIterations);
        io.outBool(r.converged);
        io.outReal(r.minError);
        io.outReal(r.minErrorDifference);
        io.outReal(r.minUpdateLength);
        for (int32_t i = 0; i < r.minLocation.GetSize(); ++i)
        {
            io.outReal(r.minLocation[i]);
        }
    }

    void EmitLM(oracle::Ctx& io, LevenbergMarquardtMinimizer<double>::Result const& r)
    {
        io.outInt(r.numIterations);
        io.outInt(r.numAdjustments);
        io.outBool(r.converged);
        io.outReal(r.minError);
        io.outReal(r.minErrorDifference);
        io.outReal(r.minUpdateLength);
        for (int32_t i = 0; i < r.minLocation.GetSize(); ++i)
        {
            io.outReal(r.minLocation[i]);
        }
    }
}

// ---- GaussNewtonMinimizer --------------------------------------------------

ORACLE_CASE("GaussNewtonMinimizer.minimize.rosenbrock")
{
    GVector<double> p0{};
    Rosenbrock problem = MakeRosenbrock(io, p0);
    int maxIterations = io.integer(0, 6);
    double updateLengthTolerance = (io.index() % 3 == 0 ? io.real(1e-8, 1e-4)
        : io.given(0.0));
    double errorDifferenceTolerance = (io.index() % 3 == 1 ? io.real(1e-8, 1e-4)
        : io.given(0.0));
    GaussNewtonMinimizer<double> minimizer(2, 2,
        [&problem](GVector<double> const& p, GVector<double>& f) { problem.Eval(p, f); },
        [&problem](GVector<double> const& p, GMatrix<double>& j) { problem.Jacobian(p, j); });
    auto result = minimizer(p0, static_cast<size_t>(maxIterations),
        updateLengthTolerance, errorDifferenceTolerance);
    EmitGN(io, result);
}

ORACLE_CASE("GaussNewtonMinimizer.minimize.circleFit")
{
    int mode = io.index() % 3;
    GVector<double> p0{};
    CircleFit problem = MakeCircleFit(io, mode, p0);
    int numF = static_cast<int>(problem.points.size());
    int maxIterations = io.integer(0, 5);
    double updateLengthTolerance = (io.index() % 4 == 0 ? io.real(1e-8, 1e-4)
        : io.given(0.0));
    double errorDifferenceTolerance = (io.index() % 4 == 1 ? io.real(1e-8, 1e-4)
        : io.given(0.0));
    GaussNewtonMinimizer<double> minimizer(3, numF,
        [&problem](GVector<double> const& p, GVector<double>& f) { problem.Eval(p, f); },
        [&problem](GVector<double> const& p, GMatrix<double>& j) { problem.Jacobian(p, j); });
    auto result = minimizer(p0, static_cast<size_t>(maxIterations),
        updateLengthTolerance, errorDifferenceTolerance);
    EmitGN(io, result);
}

// The "J plus" constructor: the caller supplies J^T*J and -J^T*F directly.
ORACLE_CASE("GaussNewtonMinimizer.minimize.jPlus")
{
    int mode = io.index() % 3;
    GVector<double> p0{};
    CircleFit problem = MakeCircleFit(io, mode, p0);
    int numF = static_cast<int>(problem.points.size());
    int maxIterations = io.integer(0, 5);
    GaussNewtonMinimizer<double> minimizer(3, numF,
        [&problem](GVector<double> const& p, GVector<double>& f) { problem.Eval(p, f); },
        [&problem](GVector<double> const& p, GMatrix<double>& jtj, GVector<double>& negJTF)
        { problem.Normal(p, jtj, negJTF); });
    auto result = minimizer(p0, static_cast<size_t>(maxIterations), 0.0, 0.0);
    EmitGN(io, result);
}

// Modes 0 and 1 give J^T*J a zero eigenvalue, so CholeskyDecomposition::
// Factor fails on the first iteration and the minimizer returns with
// numIterations = 1, converged = false and minErrorDifference = DBL_MAX.
ORACLE_CASE("GaussNewtonMinimizer.minimize.rankDeficient")
{
    int mode = io.index() % 3;
    GVector<double> p0{};
    LinearModel problem = MakeLinearModel(io, mode, p0);
    int numF = static_cast<int>(problem.t.size());
    int maxIterations = io.integer(1, 4);
    GaussNewtonMinimizer<double> minimizer(2, numF,
        [&problem](GVector<double> const& p, GVector<double>& f) { problem.Eval(p, f); },
        [&problem](GVector<double> const& p, GMatrix<double>& j) { problem.Jacobian(p, j); });
    auto result = minimizer(p0, static_cast<size_t>(maxIterations), 0.0, 0.0);
    EmitGN(io, result);
}

// F(p0) = 0 exactly: minError is 0, no later error can be smaller, so the
// loop runs to maxIterations with converged = false, minUpdateLength = 0 and
// minErrorDifference = DBL_MAX.
ORACLE_CASE("GaussNewtonMinimizer.minimize.zeroResidual")
{
    GVector<double> p0{};
    LinearModel problem = MakeLinearModel(io, 2, p0, true);
    int numF = static_cast<int>(problem.t.size());
    int maxIterations = io.integer(1, 4);
    GaussNewtonMinimizer<double> minimizer(2, numF,
        [&problem](GVector<double> const& p, GVector<double>& f) { problem.Eval(p, f); },
        [&problem](GVector<double> const& p, GMatrix<double>& j) { problem.Jacobian(p, j); });
    auto result = minimizer(p0, static_cast<size_t>(maxIterations), 0.0, 0.0);
    EmitGN(io, result);
}

// Throw parity for LogAssert(mNumPDimensions > 0 && mNumFDimensions > 0).
ORACLE_CASE("GaussNewtonMinimizer.construct.invalidDimensions")
{
    int numP = io.integer(0, 2);
    int numF = (numP > 0 ? io.given(0.0) : io.integer(0, 2));
    GaussNewtonMinimizer<double> minimizer(numP, static_cast<int32_t>(numF),
        [](GVector<double> const&, GVector<double>& f) { f[0] = 0.0; },
        [](GVector<double> const&, GMatrix<double>& j) { j(0, 0) = 0.0; });
    io.outInt(minimizer.GetNumPDimensions());
    io.outInt(minimizer.GetNumFDimensions());
}

// ---- LevenbergMarquardtMinimizer -------------------------------------------
// The main cases restrict maxIterations to the longest prefix on which
// upstream never rebuilds -J^T*F from a stale residual (issue #261); see the
// note at the top of the file. Every fourth record passes a nonpositive
// lambdaFactor, which upstream turns into a single Gauss-Newton adjustment.

ORACLE_CASE("LevenbergMarquardtMinimizer.minimize.rosenbrock")
{
    GVector<double> p0{};
    Rosenbrock problem = MakeRosenbrock(io, p0);
    int maxAdjustments = io.integer(1, 3);
    double lambdaFactor = (io.index() % 4 == 0 ? io.real(-1.0, 0.0)
        : io.real(1e-4, 1e-1));
    double lambdaAdjust = io.real(2.0, 10.0);
    LevenbergMarquardtMinimizer<double> minimizer(2, 2,
        [&problem](GVector<double> const& p, GVector<double>& f) { problem.Eval(p, f); },
        [&problem](GVector<double> const& p, GMatrix<double>& j) { problem.Jacobian(p, j); });
    size_t maxIterations = SoundLMIterations(minimizer, p0, 0.0, 0.0,
        lambdaFactor, lambdaAdjust, static_cast<size_t>(maxAdjustments), 5, true);
    io.given(static_cast<double>(maxIterations));
    auto result = minimizer(p0, maxIterations, 0.0, 0.0, lambdaFactor,
        lambdaAdjust, static_cast<size_t>(maxAdjustments));
    EmitLM(io, result);
}

ORACLE_CASE("LevenbergMarquardtMinimizer.minimize.circleFit")
{
    int mode = io.index() % 3;
    GVector<double> p0{};
    CircleFit problem = MakeCircleFit(io, mode, p0);
    int numF = static_cast<int>(problem.points.size());
    int maxAdjustments = io.integer(1, 3);
    double lambdaFactor = (io.index() % 4 == 0 ? io.real(-1.0, 0.0)
        : io.real(1e-4, 1e-1));
    double lambdaAdjust = io.real(2.0, 10.0);
    double updateLengthTolerance = (io.index() % 5 == 0 ? io.real(1e-8, 1e-4)
        : io.given(0.0));
    LevenbergMarquardtMinimizer<double> minimizer(3, numF,
        [&problem](GVector<double> const& p, GVector<double>& f) { problem.Eval(p, f); },
        [&problem](GVector<double> const& p, GMatrix<double>& j) { problem.Jacobian(p, j); });
    size_t maxIterations = SoundLMIterations(minimizer, p0,
        updateLengthTolerance, 0.0, lambdaFactor, lambdaAdjust,
        static_cast<size_t>(maxAdjustments), 5, true);
    io.given(static_cast<double>(maxIterations));
    auto result = minimizer(p0, maxIterations, updateLengthTolerance, 0.0,
        lambdaFactor, lambdaAdjust, static_cast<size_t>(maxAdjustments));
    EmitLM(io, result);
}

// The deliberate port fix of issue #261: maxIterations is the first prefix
// length at which upstream repeats a DoIteration for the same pCurrent and
// therefore builds the step from the residual at the previously rejected
// candidate.
ORACLE_CASE("LevenbergMarquardtMinimizer.minimize.staleResidual.deviation")
{
    GVector<double> p0{};
    Rosenbrock problem = MakeRosenbrock(io, p0);
    int maxAdjustments = io.integer(1, 3);
    double lambdaFactor = io.real(1e-4, 1e-1);
    double lambdaAdjust = io.real(2.0, 10.0);
    LevenbergMarquardtMinimizer<double> minimizer(2, 2,
        [&problem](GVector<double> const& p, GVector<double>& f) { problem.Eval(p, f); },
        [&problem](GVector<double> const& p, GMatrix<double>& j) { problem.Jacobian(p, j); });
    size_t maxIterations = SoundLMIterations(minimizer, p0, 0.0, 0.0,
        lambdaFactor, lambdaAdjust, static_cast<size_t>(maxAdjustments), 12, false);
    io.given(static_cast<double>(maxIterations));
    auto result = minimizer(p0, maxIterations, 0.0, 0.0, lambdaFactor,
        lambdaAdjust, static_cast<size_t>(maxAdjustments));
    EmitLM(io, result);
}

// The "J plus" path is unaffected by issue #261: the callback computes both
// J^T*J and -J^T*F from pCurrent itself, so no stale residual can leak in and
// the generator is unrestricted.
ORACLE_CASE("LevenbergMarquardtMinimizer.minimize.jPlus")
{
    int mode = io.index() % 3;
    GVector<double> p0{};
    CircleFit problem = MakeCircleFit(io, mode, p0);
    int numF = static_cast<int>(problem.points.size());
    int maxIterations = io.integer(0, 5);
    int maxAdjustments = io.integer(1, 3);
    double lambdaFactor = (io.index() % 4 == 0 ? io.real(-1.0, 0.0)
        : io.real(1e-4, 1e-1));
    double lambdaAdjust = io.real(2.0, 10.0);
    LevenbergMarquardtMinimizer<double> minimizer(3, numF,
        [&problem](GVector<double> const& p, GVector<double>& f) { problem.Eval(p, f); },
        [&problem](GVector<double> const& p, GMatrix<double>& jtj, GVector<double>& negJTF)
        { problem.Normal(p, jtj, negJTF); });
    auto result = minimizer(p0, static_cast<size_t>(maxIterations), 0.0, 0.0,
        lambdaFactor, lambdaAdjust, static_cast<size_t>(maxAdjustments));
    EmitLM(io, result);
}

// CholeskyDecomposition::Factor fails on the damped normal equations too when
// a Jacobian column is identically zero and lambda is zero; with a positive
// lambda the damping makes J^T*J invertible, so both branches are reached.
ORACLE_CASE("LevenbergMarquardtMinimizer.minimize.rankDeficient")
{
    int mode = io.index() % 3;
    GVector<double> p0{};
    LinearModel problem = MakeLinearModel(io, mode, p0);
    int numF = static_cast<int>(problem.t.size());
    int maxAdjustments = io.integer(1, 3);
    double lambdaFactor = (io.index() % 2 == 0 ? io.given(0.0)
        : io.real(1e-3, 1e-1));
    double lambdaAdjust = io.real(2.0, 10.0);
    LevenbergMarquardtMinimizer<double> minimizer(2, numF,
        [&problem](GVector<double> const& p, GVector<double>& f) { problem.Eval(p, f); },
        [&problem](GVector<double> const& p, GMatrix<double>& j) { problem.Jacobian(p, j); });
    size_t maxIterations = SoundLMIterations(minimizer, p0, 0.0, 0.0,
        lambdaFactor, lambdaAdjust, static_cast<size_t>(maxAdjustments), 4, true);
    io.given(static_cast<double>(maxIterations));
    auto result = minimizer(p0, maxIterations, 0.0, 0.0, lambdaFactor,
        lambdaAdjust, static_cast<size_t>(maxAdjustments));
    EmitLM(io, result);
}

// Throw parity for LogAssert(mNumPDimensions > 0 && mNumFDimensions > 0).
ORACLE_CASE("LevenbergMarquardtMinimizer.construct.invalidDimensions")
{
    int numP = io.integer(0, 2);
    int numF = (numP > 0 ? io.given(0.0) : io.integer(0, 2));
    LevenbergMarquardtMinimizer<double> minimizer(numP, static_cast<int32_t>(numF),
        [](GVector<double> const&, GVector<double>& f) { f[0] = 0.0; },
        [](GVector<double> const&, GMatrix<double>& j) { j(0, 0) = 0.0; });
    io.outInt(minimizer.GetNumPDimensions());
    io.outInt(minimizer.GetNumFDimensions());
}
