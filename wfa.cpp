#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <unistd.h>   // for getopt
#include <iomanip>
#include <cmath>
#include <cctype>
#include <algorithm>
#include "bindings/cpp/WFAligner.hpp"

// g++ -O3 -std=c++11 wavefront_twoseq_stats_time.cpp \
//     WFA2-lib/bindings/cpp/WFAligner.cpp \
//     -I./WFA2-lib -I./WFA2-lib/bindings/cpp \
//     ./WFA2-lib/build/libwfa2.a -o wavefront_twoseq_stats_time

using namespace std;
using namespace wfa;

// ---------- Utils for stats ----------
static inline double mean_of(const vector<double>& v) {
    if (v.empty()) return NAN;
    double s = 0.0; for (double x : v) s += x;
    return s / v.size();
}
static inline double var_of(const vector<double>& v, bool unbiased=true) {
    size_t n = v.size(); if (n < 2) return 0.0;
    double m = mean_of(v), s2 = 0.0;
    for (double x : v) { double d = x - m; s2 += d*d; }
    return unbiased ? s2 / double(n-1) : s2 / double(n);
}
static inline double sd_of(const vector<double>& v) {
    return std::sqrt(var_of(v));
}

// ---------- Original helpers ----------
void print_usage(const char* prog) {
    cerr << "Usage: " << prog << " [options] <two-seq-fasta>\n"
         << "Options:\n"
         << "  -x <int>        mismatch penalty (default 6)\n"
         << "  -O <int>        gap opening1 penalty (default 4)\n"
         << "  -E <int>        gap extension1 penalty (default 2)\n"
         << "  -o <int>        gap opening2 penalty (default 100)\n"
         << "  -e <int>        gap extension2 penalty (default 1)\n"
         << "  -c [sam|full]   CIGAR format (default full)\n"
         << "  -u <float>      mutation rate \xCE\xBC (default 3e-8)\n"
         << "  -W <int>        window size for p-distance metrics (0=off)\n"
         << "  -S <int>        slide size (default W/2)\n"
         << "  -h              show this help message\n\n"
         << "Output columns:\n"
         << "  header1, header2, CIGAR, [identity if full],\n"
         << "  total_length, total_substitutions, transitions, transversions,\n"
         << "  raw_d, raw_T, JC69_d, JC69_T, K2P_d, K2P_T,\n"
         << "  [if -W>0: win_n, win_mean, win_overdisp]\n";
}

bool read_two_fasta(const string& path, string& h1, string& s1, string& h2, string& s2) {
    ifstream in(path);
    if (!in) return false;
    string line, header, seq;
    vector<pair<string,string>> recs;
    while (getline(in, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            if (!header.empty()) {
                recs.emplace_back(header, seq);
                seq.clear();
                if (recs.size() >= 2) break;
            }
            header = line.substr(1);
            size_t tpos = header.find('\t');
            if (tpos != string::npos) header.erase(tpos);
        } else {
            seq += line;
        }
    }
    if (!header.empty() && recs.size() < 2)
        recs.emplace_back(header, seq);
    if (recs.size() != 2) return false;
    h1 = recs[0].first;  s1 = recs[0].second;
    h2 = recs[1].first;  s2 = recs[1].second;
    return true;
}

double compute_gc_identity(const string& cigar) {
    long matches = 0, mismatches = 0, gaps = 0;
    int  num = 0;
    for (char c : cigar) {
        if (isdigit(c)) {
            num = num * 10 + (c - '0');
        } else {
            switch (c) {
                case '=': matches    += num; break;
                case 'X': mismatches += num; break;
                case 'D': case 'I': gaps += 1;   break;
                default:               break;
            }
            num = 0;
        }
    }
    return double(matches) / (matches + mismatches + gaps);
}

void count_substitutions(
    const string& cigar,
    const string& s1,
    const string& s2,
    long &total_subs,
    long &transitions,
    long &transversions
) {
    size_t i = 0, j = 0, pos = 0;
    total_subs = transitions = transversions = 0;
    while (pos < cigar.size()) {
        int num = 0;
        while (pos < cigar.size() && isdigit(cigar[pos])) {
            num = num * 10 + (cigar[pos++] - '0');
        }
        if (pos >= cigar.size()) break;
        char op = cigar[pos++];
        if (op == '=') {
            i += num;  j += num;
        } else if (op == 'X') {
            for (int k = 0; k < num; ++k) {
                char a = toupper(s1[i++]);
                char b = toupper(s2[j++]);
                ++total_subs;
                // transitions: A<->G or C<->T
                if ((a=='A' && b=='G') || (a=='G' && b=='A') ||
                    (a=='C' && b=='T') || (a=='T' && b=='C')) {
                        ++transitions;
                } else {
                        ++transversions;
                }
            }
        } else if (op == 'D') {
            i += num;
        } else if (op == 'I') {
            j += num;
        }
    }
}

// ---------- Window p-distances from full CIGAR ----------
static vector<int> mismatch_series_from_cigar(const string& cigar) {
    // 1 element per aligned PAIR (= or X). 0 for match, 1 for mismatch. Gaps are skipped.
    vector<int> mm;
    size_t pos = 0;
    while (pos < cigar.size()) {
        int num = 0;
        while (pos < cigar.size() && isdigit(cigar[pos])) {
            num = num * 10 + (cigar[pos++] - '0');
        }
        if (pos >= cigar.size()) break;
        char op = cigar[pos++];
        if (op == '=') {
            mm.insert(mm.end(), (size_t)num, 0);
        } else if (op == 'X') {
            mm.insert(mm.end(), (size_t)num, 1);
        } else {
            // I or D -> skip; not aligned pairs
        }
    }
    return mm;
}
static vector<double> sliding_pdist(const vector<int>& mm, int W, int S) {
    vector<double> out;
    if (W <= 0 || mm.empty()) return out;
    if (S <= 0) S = max(1, W/2);
    // prefix sums
    vector<int> pref(mm.size()+1, 0);
    for (size_t i=0;i<mm.size();++i) pref[i+1] = pref[i] + mm[i];
    for (int start = 0; start + W <= (int)mm.size(); start += S) {
        int end = start + W;
        int mism = pref[end] - pref[start];
        out.push_back( double(mism) / double(W) );
    }
    if (out.empty() && (int)mm.size() > 0) {
        int W2 = (int)mm.size();
        int mism = pref[W2] - pref[0];
        out.push_back(double(mism) / double(W2));
    }
    return out;
}

// ---------- Main ----------
int main(int argc, char** argv) {
    // default penalties and format
    int mismatch  = 6;
    int go1       = 4;
    int ge1       = 2;
    int go2       = 100;
    int ge2       = 1;
    string format = "full";
    double miu    = 3e-8;

    // window options
    int W = 0;   // window size (0=disabled)
    int S = 0;   // slide (0 => W/2)

    // parse options
    int opt;
    while ((opt = getopt(argc, argv, "x:O:E:o:e:c:u:W:S:h")) != -1) {
        switch (opt) {
            case 'x': mismatch = stoi(optarg);   break;
            case 'O': go1      = stoi(optarg);   break;
            case 'E': ge1      = stoi(optarg);   break;
            case 'o': go2      = stoi(optarg);   break;
            case 'e': ge2      = stoi(optarg);   break;
            case 'c': format   = optarg;         break;
            case 'u': miu      = stod(optarg);   break;
            case 'W': W        = stoi(optarg);   break;
            case 'S': S        = stoi(optarg);   break;
            case 'h':
            default:
                print_usage(argv[0]);
                return (opt=='h'? 0 : 1);
        }
    }
    if (format!="full" && format!="sam") {
        cerr << "Error: unknown CIGAR format '" << format << "'\n";
        print_usage(argv[0]);
        return 1;
    }
    if (optind + 1 != argc) {
        print_usage(argv[0]);
        return 1;
    }
    const string fasta_path = argv[optind];

    // read two FASTA records
    string h1, s1, h2, s2;
    if (!read_two_fasta(fasta_path, h1, s1, h2, s2)) {
        cerr << "Error: need exactly two sequences in " << fasta_path << "\n";
        return 1;
    }

    // run the 2-piece gap-affine aligner
    WFAlignerGapAffine2Pieces aligner(
        mismatch, go1, ge1, go2, ge2,
        WFAligner::Alignment, WFAligner::MemoryHigh
    );
    int status = aligner.alignEnd2End(s1, s2);
    if (status < 0) {
        cerr << "Alignment failed (status=" << status << ")\n";
        return 1;
    }

    // get CIGAR
    string cigar = aligner.getCIGAR(format=="full");

    // compute substitution statistics
    long total_subs, trans_count, transv_count;
    count_substitutions(cigar, s1, s2, total_subs, trans_count, transv_count);

    // compute distances
    double tot_len = double(s1.size());
    double raw_d    = (trans_count + transv_count) / tot_len;
    double P        = trans_count / tot_len;
    double Q        = transv_count / tot_len;
    double JC69_d   = -0.75 * log(1.0 - (4.0 * raw_d / 3.0));
    double K2P_d    = -0.5 * log((1.0 - 2.0*P - Q) * sqrt(1.0 - 2.0*Q));

    long raw_T  = lround(raw_d  / (2.0 * miu));
    long JC69_T = lround(JC69_d / (2.0 * miu));
    long K2P_T  = lround(K2P_d  / (2.0 * miu));

    // ---------- Sliding-window metrics (reduced) ----------
    vector<double> win_p;
    int used_W = W, used_S = S;
    if (format=="full" && W > 0) {
        vector<int> mm = mismatch_series_from_cigar(cigar);
        if (used_S <= 0) used_S = max(1, used_W/2);
        win_p = sliding_pdist(mm, used_W, used_S);
    }

    size_t win_n = win_p.size();
    double win_mean = NAN;
    double win_overdisp = NAN;

    if (win_n > 0) {
        win_mean = mean_of(win_p);
        // Overdispersion vs binomial sampling at window size used_W:
        // expected variance of p-hat in each window = p(1-p)/W
        if (used_W > 0) {
            double expected = win_mean * (1.0 - win_mean) / double(used_W);
            if (expected <= 0.0) {
                win_overdisp = 0.0;
            } else {
                win_overdisp = var_of(win_p) / expected; // >1 => extra-binomial heterogeneity
            }
        }
    }

    // output results
    cout << h1 << "\t" << h2 << "\t" << cigar;
    if (format=="full") {
        double id = compute_gc_identity(cigar);
        cout << "\t" << fixed << setprecision(6) << id;
    }
    cout << "\t"
         << long(tot_len) << "\t"
         << total_subs << "\t"
         << trans_count << "\t"
         << transv_count << "\t"
         << fixed << setprecision(6)
         << raw_d    << "\t"
         << raw_T    << "\t"
         << JC69_d   << "\t"
         << JC69_T   << "\t"
         << K2P_d    << "\t"
         << K2P_T;

    if (format=="full" && W > 0) {
        cout << "\t" << win_n
             << "\t" << setprecision(6) << win_mean
             << "\t" << setprecision(6) << win_overdisp;
    }
    cout << "\n";
    return 0;
}
