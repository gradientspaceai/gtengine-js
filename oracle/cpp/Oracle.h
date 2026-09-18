// gtengine-js C++ oracle harness.
//
// Each ORACLE_CASE generates its own inputs through a Ctx (every generated
// value is recorded bit-for-bit as an input), runs the upstream GTE code and
// records the outputs. The TypeScript side replays the recorded inputs through
// the port and compares outputs. See ORACLE.md.
#pragma once

#include <Mathematics/Vector2.h>
#include <Mathematics/Vector3.h>
#include <Mathematics/Vector4.h>
#include <Mathematics/Matrix.h>
#include <Mathematics/GVector.h>
#include <Mathematics/GMatrix.h>
#include <cstdint>
#include <cstring>
#include <map>
#include <string>
#include <type_traits>
#include <vector>

namespace oracle
{
    inline uint64_t SplitMix64(uint64_t& state)
    {
        uint64_t z = (state += 0x9E3779B97F4A7C15ull);
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
        return z ^ (z >> 31);
    }

    inline uint64_t Fnv1a(std::string const& s)
    {
        uint64_t h = 0xCBF29CE484222325ull;
        for (unsigned char c : s) { h ^= c; h *= 0x100000001B3ull; }
        return h;
    }

    class Ctx
    {
    public:
        Ctx(uint64_t seed, int index) : mState(seed), mIndex(index) {}

        // The record number, 0 <= index < n. Use it to switch generator
        // modes deterministically (e.g. every 4th record on a lattice).
        int index() const { return mIndex; }

        // ---- inputs (everything returned here is recorded) ----
        double real(double lo, double hi) { return given(raw(lo, hi)); }

        // Integer in [lo, hi], inclusive.
        int integer(int lo, int hi)
        {
            int v = rawInteger(lo, hi);
            mIn.push_back(static_cast<double>(v));
            return v;
        }

        bool boolean() { return integer(0, 1) != 0; }

        // Integer-valued double in [lo, hi]. Lattice inputs reach the
        // degenerate branches (parallel, touching, coincident) with exactly
        // representable arithmetic.
        double lattice(int lo, int hi) { return static_cast<double>(integer(lo, hi)); }

        // Record a derived value as an input (e.g. a normalized direction).
        double given(double x) { mIn.push_back(x); return x; }

        template <int N>
        gte::Vector<N, double> vec(double lo, double hi)
        {
            gte::Vector<N, double> v;
            for (int i = 0; i < N; ++i) { v[i] = real(lo, hi); }
            return v;
        }

        template <int N>
        gte::Vector<N, double> latticeVec(int lo, int hi)
        {
            gte::Vector<N, double> v;
            for (int i = 0; i < N; ++i) { v[i] = lattice(lo, hi); }
            return v;
        }

        // Unit-length vector. Only the normalized components are recorded.
        template <int N>
        gte::Vector<N, double> unit()
        {
            gte::Vector<N, double> v;
            double len = 0.0;
            do
            {
                for (int i = 0; i < N; ++i) { v[i] = raw(-1.0, 1.0); }
                len = gte::Length(v);
            } while (len < 0.1 || len > 1.0);
            gte::Normalize(v);
            return givenVec(v);
        }

        // Nonzero integer-valued vector (not normalized).
        template <int N>
        gte::Vector<N, double> latticeDir(int lo, int hi)
        {
            gte::Vector<N, double> v;
            bool zero = true;
            do
            {
                zero = true;
                for (int i = 0; i < N; ++i)
                {
                    v[i] = static_cast<double>(rawInteger(lo, hi));
                    zero = zero && v[i] == 0.0;
                }
            } while (zero);
            return givenVec(v);
        }

        template <int N>
        gte::Vector<N, double> givenVec(gte::Vector<N, double> const& v)
        {
            for (int i = 0; i < N; ++i) { mIn.push_back(v[i]); }
            return v;
        }

        // Unrecorded draws, for building derived inputs that are then
        // recorded with given()/givenVec().
        double raw(double lo, double hi)
        {
            double u = static_cast<double>(SplitMix64(mState) >> 11) * 0x1.0p-53;
            return lo + (hi - lo) * u;
        }

        int rawInteger(int lo, int hi)
        {
            uint64_t span = static_cast<uint64_t>(static_cast<int64_t>(hi) - lo + 1);
            return lo + static_cast<int>(SplitMix64(mState) % span);
        }

        // ---- outputs ----
        void outReal(double x) { mOut.push_back(x); }
        void outBool(bool b) { mOut.push_back(b ? 1.0 : 0.0); }

        template <typename I, typename = std::enable_if_t<std::is_integral<I>::value>>
        void outInt(I i) { mOut.push_back(static_cast<double>(i)); }

        template <int N>
        void outVec(gte::Vector<N, double> const& v)
        {
            for (int i = 0; i < N; ++i) { mOut.push_back(v[i]); }
        }

        // Row-major.
        template <int R, int C>
        void outMat(gte::Matrix<R, C, double> const& m)
        {
            for (int r = 0; r < R; ++r)
            {
                for (int c = 0; c < C; ++c) { mOut.push_back(m(r, c)); }
            }
        }

        // Size first, then the elements.
        void outGVec(gte::GVector<double> const& v)
        {
            outInt(v.GetSize());
            for (int i = 0; i < v.GetSize(); ++i) { mOut.push_back(v[i]); }
        }

        // Rows, columns, then the elements row-major.
        void outGMat(gte::GMatrix<double> const& m)
        {
            outInt(m.GetNumRows());
            outInt(m.GetNumCols());
            for (int r = 0; r < m.GetNumRows(); ++r)
            {
                for (int c = 0; c < m.GetNumCols(); ++c) { mOut.push_back(m(r, c)); }
            }
        }

        std::vector<double> const& inputs() const { return mIn; }
        std::vector<double> const& outputs() const { return mOut; }

    private:
        uint64_t mState;
        int mIndex;
        std::vector<double> mIn, mOut;
    };

    using CaseFn = void (*)(Ctx&);

    struct CaseEntry
    {
        std::string family, name;
        CaseFn fn;
    };

    inline std::vector<CaseEntry>& Registry()
    {
        static std::vector<CaseEntry> registry;
        return registry;
    }

    struct Registrar
    {
        Registrar(char const* family, char const* name, CaseFn fn)
        {
            Registry().push_back({ family, name, fn });
        }
    };
}

// Define ORACLE_FAMILY (a string literal, the golden file basename) before
// the first ORACLE_CASE of a translation unit. 'name' is the case name the
// TypeScript side claims, conventionally "<Header>.<method>[.<variant>]".
// Inside the body the context is available as 'io'.
#define ORACLE_CONCAT2(a, b) a##b
#define ORACLE_CONCAT(a, b) ORACLE_CONCAT2(a, b)
#define ORACLE_CASE(name) \
    static void ORACLE_CONCAT(OracleCaseFn, __LINE__)(oracle::Ctx&); \
    static oracle::Registrar ORACLE_CONCAT(OracleCaseReg, __LINE__)( \
        ORACLE_FAMILY, name, ORACLE_CONCAT(OracleCaseFn, __LINE__)); \
    static void ORACLE_CONCAT(OracleCaseFn, __LINE__)(oracle::Ctx& io)
