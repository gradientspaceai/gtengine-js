// gtengine-js C++ oracle driver. See ORACLE.md.
//
//   oracle.exe --out <dir> [--n <records per case>] [--filter <substring>] [--list]
//
// Writes one golden file per family: <dir>/<family>.txt
#include "Oracle.h"
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <iostream>

#ifndef ORACLE_UPSTREAM_COMMIT
#define ORACLE_UPSTREAM_COMMIT "unknown"
#endif

static void WriteHex(std::ofstream& f, char tag, std::vector<double> const& v)
{
    f << tag;
    char buf[20];
    for (double d : v)
    {
        uint64_t bits;
        std::memcpy(&bits, &d, sizeof(bits));
        std::snprintf(buf, sizeof(buf), " %016llx", static_cast<unsigned long long>(bits));
        f << buf;
    }
    f << '\n';
}

int main(int argc, char** argv)
{
    std::string outDir = ".", filter;
    int n = 20;
    bool list = false;
    for (int i = 1; i < argc; ++i)
    {
        std::string a = argv[i];
        if (a == "--out" && i + 1 < argc) { outDir = argv[++i]; }
        else if (a == "--n" && i + 1 < argc) { n = std::atoi(argv[++i]); }
        else if (a == "--filter" && i + 1 < argc) { filter = argv[++i]; }
        else if (a == "--list") { list = true; }
        else { std::cerr << "unknown argument: " << a << "\n"; return 2; }
    }

    auto cases = oracle::Registry();
    std::stable_sort(cases.begin(), cases.end(), [](auto const& a, auto const& b)
        { return a.family != b.family ? a.family < b.family : a.name < b.name; });
    for (size_t i = 1; i < cases.size(); ++i)
    {
        if (cases[i].family == cases[i - 1].family && cases[i].name == cases[i - 1].name)
        {
            std::cerr << "duplicate case: " << cases[i].family << "/" << cases[i].name << "\n";
            return 2;
        }
    }

    // A filter selects whole families so a golden file is never partial.
    std::map<std::string, std::vector<oracle::CaseEntry>> families;
    for (auto const& c : cases)
    {
        if (filter.empty() || c.family.find(filter) != std::string::npos)
        {
            families[c.family].push_back(c);
        }
    }

    size_t numRecords = 0, numThrew = 0;
    for (auto const& fam : families)
    {
        if (list)
        {
            for (auto const& c : fam.second) { std::cout << fam.first << "/" << c.name << "\n"; }
            continue;
        }
        std::string path = outDir + "/" + fam.first + ".txt";
        std::ofstream f(path, std::ios::binary);
        if (!f) { std::cerr << "cannot write " << path << "\n"; return 1; }
        f << "# gtengine-js oracle golden v1\n";
        f << "# upstream davideberly/GeometricTools " << ORACLE_UPSTREAM_COMMIT << "\n";
        f << "# compiler MSVC " << _MSC_FULL_VER << " x64 /O2 /fp:precise\n";
        f << "# records-per-case " << n << "\n";
        for (auto const& c : fam.second)
        {
            f << "case " << c.name << "\n";
            uint64_t seed = oracle::Fnv1a(c.family + "/" + c.name);
            for (int i = 0; i < n; ++i)
            {
                uint64_t s = seed + 0xD1B54A32D192ED03ull * static_cast<uint64_t>(i + 1);
                oracle::Ctx io(oracle::SplitMix64(s), i);
                bool threw = false;
                try { c.fn(io); }
                catch (std::exception const&) { threw = true; }
                WriteHex(f, 'i', io.inputs());
                if (threw) { f << "o !\n"; ++numThrew; }
                else { WriteHex(f, 'o', io.outputs()); }
                ++numRecords;
            }
        }
        std::cout << fam.first << ": " << fam.second.size() << " cases\n";
    }
    if (!list)
    {
        std::cout << numRecords << " records (" << numThrew << " threw)\n";
    }
    return 0;
}
