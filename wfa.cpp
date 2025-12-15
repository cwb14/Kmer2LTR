#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <unistd.h>   // getopt
#include <iomanip>
#include <cmath>
#include <cctype>
#include <algorithm>
#include <limits>
#include "bindings/cpp/WFAligner.hpp"

// g++ -O3 -std=c++11 wavefront_twoseq_stats_time.cpp \
//     WFA2-lib/bindings/cpp/WFAligner.cpp \
//     -I./WFA2-lib -I./WFA2-lib/bindings/cpp \
//     ./WFA2-lib/build/libwfa2.a -o wavefront_twoseq_stats_time

// -K is added to test if its better to include k bases in alignment length or to exclude them... benchmarking suggests its better to exclude (use K). I guess because k bases are, by definition, always matching... including those bases in aln length appears to artificially inflate seq identity.
// In sort, its likely best to use '-k 5 -K' although the k [int] is not extensively tested. Including '-k [int]' appears to be important with empirical data where the LTR-RT may not have been precisely provided from start to end. Any extra bp tacked onto the LTR-RT appear to inflate divergence since theyre likely not perfect identity. With '-k [int]', we can exclude counting them.  

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
         << "  -k <int>        trim unreliable ends: require k consecutive matches to enter/exit reliable region (default 0=off)\n"
         << "  -K              exclude the k boundary-match columns from the reliable region metrics (requires -k>0 and successful trimming)\n"
         << "  -d              debug: print a 3-line alignment view to stderr (also prints 5' and 3' k-mers around boundaries)\n"
         << "  -h              show this help message\n\n"
         << "Output columns:\n"
         << "  header1, header2, CIGAR, [identity if full],\n"
         << "  reliable_pairs, total_substitutions, transitions, transversions,\n"
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

// ---------- Debug alignment pretty-printer ----------
struct AlignView {
    string a1, mid, a2;
};

static AlignView build_alignment_view(const string& cigar,
                                     const string& s1,
                                     const string& s2) {
    // Supports WFA "full" CIGAR (= X I D) and SAM-style (M = X I D).
    AlignView v;
    v.a1.reserve(cigar.size());
    v.a2.reserve(cigar.size());
    v.mid.reserve(cigar.size());

    size_t i = 0, j = 0, pos = 0;
    while (pos < cigar.size()) {
        int num = 0;
        while (pos < cigar.size() && isdigit((unsigned char)cigar[pos])) {
            num = num * 10 + (cigar[pos++] - '0');
        }
        if (pos >= cigar.size()) break;
        char op = cigar[pos++];

        if (op == '=' || op == 'X' || op == 'M') {
            for (int k = 0; k < num; ++k) {
                char a = (i < s1.size()) ? s1[i++] : '-';
                char b = (j < s2.size()) ? s2[j++] : '-';
                char ua = toupper((unsigned char)a);
                char ub = toupper((unsigned char)b);

                v.a1.push_back(ua);
                v.a2.push_back(ub);

                if (op == '=') {
                    v.mid.push_back('|');
                } else if (op == 'X') {
                    v.mid.push_back('*');
                } else { // 'M' (SAM): decide by base comparison
                    v.mid.push_back(ua == ub ? '|' : '*');
                }
            }
        } else if (op == 'D') {
            for (int k = 0; k < num; ++k) {
                char a = (i < s1.size()) ? toupper((unsigned char)s1[i++]) : '-';
                v.a1.push_back(a);
                v.mid.push_back(' ');
                v.a2.push_back('-');
            }
        } else if (op == 'I') {
            for (int k = 0; k < num; ++k) {
                char b = (j < s2.size()) ? toupper((unsigned char)s2[j++]) : '-';
                v.a1.push_back('-');
                v.mid.push_back(' ');
                v.a2.push_back(b);
            }
        } else {
            // Unknown op: ignore
        }
    }
    return v;
}

static void print_alignment_view(const AlignView& v, ostream& os, size_t width=80) {
    for (size_t p = 0; p < v.a1.size(); p += width) {
        os << v.a1.substr(p, width) << "\n";
        os << v.mid.substr(p, width) << "\n";
        os << v.a2.substr(p, width) << "\n\n";
    }
}

// ---------- Reliable region detection (k consecutive matches) ----------
struct ReliableBounds {
    // Alignment-column coordinates (including gaps) in [start_col, end_col_excl)
    size_t start_col = 0;
    size_t end_col_excl = 0;
    bool   enabled = false;
    bool   found = false;
};

static inline bool is_match_op(char op,
                               bool format_full,
                               const string& s1, const string& s2,
                               size_t i, size_t j) {
    if (format_full) {
        return op == '=';
    } else {
        // SAM: 'M' can be match or mismatch; '=' may appear too depending on source
        if (op == '=') return true;
        if (op == 'X') return false;
        if (op == 'M') {
            if (i >= s1.size() || j >= s2.size()) return false;
            char a = toupper((unsigned char)s1[i]);
            char b = toupper((unsigned char)s2[j]);
            return a == b;
        }
        return false;
    }
}

static ReliableBounds find_reliable_bounds_k(const string& cigar,
                                            const string& s1,
                                            const string& s2,
                                            bool format_full,
                                            int k) {
    ReliableBounds rb;
    rb.end_col_excl = 0;
    rb.start_col = 0;
    rb.enabled = (k > 0);

    // If k disabled, whole alignment is reliable (we still need total alignment length in cols)
    // We'll compute alignment length in cols regardless.
    auto compute_alignment_cols = [&]() -> size_t {
        size_t cols = 0;
        size_t pos = 0;
        while (pos < cigar.size()) {
            int num = 0;
            while (pos < cigar.size() && isdigit((unsigned char)cigar[pos])) {
                num = num * 10 + (cigar[pos++] - '0');
            }
            if (pos >= cigar.size()) break;
            char op = cigar[pos++];
            if (op=='=' || op=='X' || op=='M' || op=='I' || op=='D') cols += (size_t)num;
        }
        return cols;
    };

    const size_t aln_cols = compute_alignment_cols();
    rb.end_col_excl = aln_cols;

    if (!rb.enabled) {
        rb.found = true;
        rb.start_col = 0;
        rb.end_col_excl = aln_cols;
        return rb;
    }

    // Forward scan: find first run of k consecutive matches on aligned PAIRS only (not gaps).
    size_t i = 0, j = 0, pos = 0, col = 0;
    int run = 0;
    bool have_start = false;
    size_t start_col = 0;

    while (pos < cigar.size() && !have_start) {
        int num = 0;
        while (pos < cigar.size() && isdigit((unsigned char)cigar[pos])) {
            num = num * 10 + (cigar[pos++] - '0');
        }
        if (pos >= cigar.size()) break;
        char op = cigar[pos++];

        if (op=='=' || op=='X' || op=='M') {
            for (int t=0; t<num; ++t) {
                bool m = is_match_op(op, format_full, s1, s2, i, j);
                // advance indices for aligned pair
                if (i < s1.size()) ++i; else {}
                if (j < s2.size()) ++j; else {}
                // alignment column advances
                ++col;

                if (m) {
                    ++run;
                    if (run >= k) {
                        start_col = col - (size_t)k; // first col of this k-run (0-based)
                        have_start = true;
                        break;
                    }
                } else {
                    run = 0;
                }
            }
        } else if (op=='D') {
            // gap in query; NOT an aligned pair; reset run
            i = min(s1.size(), i + (size_t)num);
            col += (size_t)num;
            run = 0;
        } else if (op=='I') {
            // gap in ref; NOT an aligned pair; reset run
            j = min(s2.size(), j + (size_t)num);
            col += (size_t)num;
            run = 0;
        } else {
            // unknown op; ignore
        }
    }

    // Backward scan: find last run of k consecutive matches (aligned pairs only), from the end.
    // Easiest: build AlignView mid string once and scan it, treating '|' as match, '*' mismatch,
    // and ' ' (gap) as reset (not aligned pairs).
    AlignView av = build_alignment_view(cigar, s1, s2);
    // av.mid length == alignment columns
    size_t end_excl = av.mid.size();
    int rrun = 0;
    bool have_end = false;
    size_t end_col_excl = end_excl;

    for (size_t idx = end_excl; idx-- > 0; ) {
        char m = av.mid[idx];
        if (m == '|') {
            ++rrun;
            if (rrun >= k) {
                // idx is the first position of the reverse k-run (from the end),
                // and (idx + k - 1) is the last match position of that run.
                size_t last_match_pos = idx + (size_t)k - 1;
                end_col_excl = last_match_pos + 1; // exclusive
                have_end = true;
                break;
            }
        } else if (m == '*') {
            rrun = 0;
        } else {
            // gap column
            rrun = 0;
        }
    }

    rb.found = have_start && have_end && (start_col < end_col_excl);
    if (rb.found) {
        rb.start_col = start_col;
        rb.end_col_excl = end_col_excl;
    } else {
        // If we can't find k consecutive matches on one/both ends, fall back to whole alignment.
        rb.start_col = 0;
        rb.end_col_excl = aln_cols;
        rb.found = false; // indicates fallback
    }
    return rb;
}

// ---------- Stats within reliable region ----------
struct RegionStats {
    long pairs = 0;          // aligned pairs (match+mismatch)
    long matches = 0;
    long mismatches = 0;
    long gap_cols = 0;       // alignment columns that are gaps (I or D)
    long total_subs = 0;
    long transitions = 0;
    long transversions = 0;

    vector<int> mm_series;   // 0/1 for aligned pairs only (for window pdist), within region
};

static RegionStats compute_region_stats(const string& cigar,
                                        const string& s1,
                                        const string& s2,
                                        bool format_full,
                                        size_t start_col,
                                        size_t end_col_excl) {
    RegionStats st;
    st.mm_series.reserve(end_col_excl > start_col ? (end_col_excl - start_col) : 0);

    size_t i = 0, j = 0, pos = 0, col = 0;

    while (pos < cigar.size() && col < end_col_excl) {
        int num = 0;
        while (pos < cigar.size() && isdigit((unsigned char)cigar[pos])) {
            num = num * 10 + (cigar[pos++] - '0');
        }
        if (pos >= cigar.size()) break;
        char op = cigar[pos++];

        if (op=='=' || op=='X' || op=='M') {
            for (int t=0; t<num; ++t) {
                bool in_region = (col >= start_col && col < end_col_excl);
                bool match = is_match_op(op, format_full, s1, s2, i, j);

                char a = (i < s1.size()) ? toupper((unsigned char)s1[i]) : 'N';
                char b = (j < s2.size()) ? toupper((unsigned char)s2[j]) : 'N';

                // advance indices and column
                if (i < s1.size()) ++i;
                if (j < s2.size()) ++j;
                ++col;

                if (!in_region) continue;

                ++st.pairs;
                if (match) {
                    ++st.matches;
                    st.mm_series.push_back(0);
                } else {
                    ++st.mismatches;
                    st.mm_series.push_back(1);
                    ++st.total_subs;
                    // transitions: A<->G or C<->T
                    if ((a=='A' && b=='G') || (a=='G' && b=='A') ||
                        (a=='C' && b=='T') || (a=='T' && b=='C')) {
                        ++st.transitions;
                    } else {
                        ++st.transversions;
                    }
                }
            }
        } else if (op=='D') {
            for (int t=0; t<num; ++t) {
                bool in_region = (col >= start_col && col < end_col_excl);
                if (i < s1.size()) ++i;
                ++col;
                if (in_region) ++st.gap_cols;
            }
        } else if (op=='I') {
            for (int t=0; t<num; ++t) {
                bool in_region = (col >= start_col && col < end_col_excl);
                if (j < s2.size()) ++j;
                ++col;
                if (in_region) ++st.gap_cols;
            }
        } else {
            // ignore unknown
        }
    }
    return st;
}

static double compute_identity_region(long matches, long mismatches, long gap_cols) {
    double denom = double(matches) + double(mismatches) + double(gap_cols);
    if (denom <= 0.0) return NAN;
    return double(matches) / denom;
}

// ---------- Window p-distances on aligned pairs only ----------
static vector<double> sliding_pdist(const vector<int>& mm, int W, int S) {
    vector<double> out;
    if (W <= 0 || mm.empty()) return out;
    if (S <= 0) S = max(1, W/2);
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

static void print_boundary_kmers(const AlignView& av,
                                 size_t start_col, size_t end_col_excl,
                                 int k, ostream& os) {
    if (k <= 0 || av.mid.empty()) return;

    os << "---- Reliable region (alignment cols) ----\n";
    os << "start_col=" << start_col << ", end_col_excl=" << end_col_excl
       << ", length=" << (end_col_excl > start_col ? (end_col_excl - start_col) : 0) << "\n";

    // 5' k columns of reliable region
    size_t left_len = min<size_t>((size_t)k, (end_col_excl > start_col ? end_col_excl - start_col : 0));
    if (left_len > 0) {
        os << "[5' first " << left_len << " cols]\n";
        os << av.a1.substr(start_col, left_len) << "\n";
        os << av.mid.substr(start_col, left_len) << "\n";
        os << av.a2.substr(start_col, left_len) << "\n";
    }

    // 3' k columns of reliable region
    size_t right_len = left_len;
    if (right_len > 0) {
        size_t right_start = end_col_excl - right_len;
        os << "[3' last " << right_len << " cols]\n";
        os << av.a1.substr(right_start, right_len) << "\n";
        os << av.mid.substr(right_start, right_len) << "\n";
        os << av.a2.substr(right_start, right_len) << "\n";
    }
    os << "-----------------------------------------\n";
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

    // boundary trimming
    int ktrim = 0; // 0=off

    bool debug = false;
    bool exclude_k = false;

    // parse options
    int opt;

    while ((opt = getopt(argc, argv, "x:O:E:o:e:c:u:W:S:k:Kdh")) != -1) {

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
            case 'k': ktrim    = stoi(optarg);   break;
            case 'K': exclude_k = true;          break;
            case 'd': debug    = true;           break;
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
    if (ktrim < 0) ktrim = 0;

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
    const bool format_full = (format == "full");

    // determine reliable region bounds in alignment columns
    ReliableBounds rb = find_reliable_bounds_k(cigar, s1, s2, format_full, ktrim);

    // Effective bounds used for stats/identity/divergence
    size_t eff_start = rb.start_col;
    size_t eff_end_excl = rb.end_col_excl;

    // Optionally exclude the k boundary-match columns from metrics
    if (exclude_k && ktrim > 0 && rb.found) {
        // Remove k columns from each side (these are guaranteed matches by construction)
        eff_start = min(eff_start + (size_t)ktrim, eff_end_excl);
        if (eff_end_excl >= (size_t)ktrim) {
            eff_end_excl = max(eff_start, eff_end_excl - (size_t)ktrim);
        } else {
            eff_end_excl = eff_start;
        }
    }

    // build alignment view for debug (and boundary k-mers)
    AlignView av;
    if (debug || ktrim > 0) {
        av = build_alignment_view(cigar, s1, s2);
    }

    // debug: print alignment view (to stderr, so stdout stays parseable)
    if (debug) {
        cerr << "==== Alignment (" << h1 << " vs " << h2 << ") ====\n";
        print_alignment_view(av, cerr, 80);
        if (ktrim > 0) {
            if (!rb.found) {
                cerr << "Note: -k " << ktrim << " could not find k consecutive matches on one/both ends; using full alignment.\n";
            }
            print_boundary_kmers(av, rb.start_col, rb.end_col_excl, ktrim, cerr);
        }
        cerr << "==== End alignment ====\n";
    }

    // compute stats within reliable region (or full if ktrim=0 / fallback)
    RegionStats st = compute_region_stats(cigar, s1, s2, format_full, eff_start, eff_end_excl);

    // denominator for divergence: aligned PAIRS only within reliable region
    double denom = double(st.pairs);
    double raw_d = NAN, P = NAN, Q = NAN, JC69_d = NAN, K2P_d = NAN;
    long raw_T = 0, JC69_T = 0, K2P_T = 0;

    if (denom > 0) {
        raw_d = double(st.transitions + st.transversions) / denom;
        P     = double(st.transitions) / denom;
        Q     = double(st.transversions) / denom;

        // Guard against invalid logs due to saturation / rounding
        double jc_arg = 1.0 - (4.0 * raw_d / 3.0);
        if (jc_arg > 0.0) JC69_d = -0.75 * log(jc_arg);

        double k2_a = (1.0 - 2.0*P - Q);
        double k2_b = (1.0 - 2.0*Q);
        if (k2_a > 0.0 && k2_b > 0.0) K2P_d = -0.5 * log(k2_a * sqrt(k2_b));

        if (std::isfinite(raw_d))  raw_T  = lround(raw_d  / (2.0 * miu));
        if (std::isfinite(JC69_d)) JC69_T = lround(JC69_d / (2.0 * miu));
        if (std::isfinite(K2P_d))  K2P_T  = lround(K2P_d  / (2.0 * miu));
    }

    // identity (only printed for full CIGAR output; uses region)
    double id = NAN;
    if (format_full) {
        id = compute_identity_region(st.matches, st.mismatches, st.gap_cols);
    }

    // ---------- Sliding-window metrics (on aligned pairs only, within reliable region) ----------
    vector<double> win_p;
    int used_W = W, used_S = S;
    if (W > 0) {
        // windowing now applies to the aligned-pair mismatch series within region
        if (used_S <= 0) used_S = max(1, used_W/2);
        win_p = sliding_pdist(st.mm_series, used_W, used_S);
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
        cout << "\t" << fixed << setprecision(6) << id;
    }
    cout << "\t"
         << st.pairs << "\t"
         << st.total_subs << "\t"
         << st.transitions << "\t"
         << st.transversions << "\t"
         << fixed << setprecision(6)
         << (std::isfinite(raw_d)  ? raw_d  : NAN) << "\t"
         << raw_T    << "\t"
         << (std::isfinite(JC69_d) ? JC69_d : NAN) << "\t"
         << JC69_T   << "\t"
         << (std::isfinite(K2P_d)  ? K2P_d  : NAN) << "\t"
         << K2P_T;

    if (W > 0) {
        cout << "\t" << win_n
             << "\t" << setprecision(6) << win_mean
             << "\t" << setprecision(6) << win_overdisp;
    }
    cout << "\n";
    return 0;
}
