#include <iostream>
#include <vector>
#include <numeric>
#include <map>

#include <boost/optional.hpp>
#include <boost/regex.hpp>

#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TH1.h>

#include "catstr.hh"
#include "timed_counter.hh"
#include "program_options.hh"
#include "binner_root.hh"
#include "tc_msg.hh"
// #include "math.hh"

using std::cout;
using std::cerr;
using std::endl;
using std::get;
using namespace ivanp;

std::string get_data_str(const char* name) {
  static const boost::regex re("^(.*/)?hgam_(data_h\\d+\\l?(_.+)?).root$");
  boost::cmatch match;
  if (boost::regex_search(name,name+std::strlen(name),match,re)) {
    return match[2];
  } else throw std::runtime_error("data regex failed to match");
}
std::string get_mc_str(const char* name) {
  static const boost::regex re("^(.*/)?hgam_mc_(h\\d+\\l?_.+).root$");
  boost::cmatch match;
  if (boost::regex_search(name,name+std::strlen(name),match,re)) {
    return match[2];
  } else throw std::runtime_error("mc regex failed to match");
}

double weight, lumi = 0;

struct bin {
  static std::map<std::string,unsigned> m;
  std::vector<double> v;
  static unsigned i;
  void operator()() { v[i] += weight; }
  static std::vector<bin*> all;
  bin() { v.reserve(2); all.push_back(this); }
};
std::map<std::string,unsigned> bin::m;
unsigned bin::i;
std::vector<bin*> bin::all;

namespace ivanp { namespace root {
template <> struct bin_converter<bin> {
  inline auto weight(const bin& b) const noexcept { return b.v[bin::i]; }
};
}}
struct all_mc {
  inline auto weight(const bin& b) const noexcept {
    return std::accumulate(++b.v.begin(),b.v.end(),0.);
  };
};

template <typename... Axes>
using hist_t = ivanp::binner<bin,
  std::tuple<ivanp::axis_spec<Axes>...>>;
using hist = hist_t<ivanp::ref_axis<double>>;

std::map<std::string,std::unique_ptr<ivanp::abstract_axis<double>>> axes;

void binning_parser(char* str, decltype(axes)& x) {
  char* it = str;
  while (*it != ':') ++it;
  std::string name(str,it-str);
  std::vector<double> edges;
  bool uniform = false;
  str = ++it;
  for (char c; ; ) {
    c = *it;
    if (c==',' || c=='\0' || c==':') {
      if (c==':') {
        if (edges.size()==0) {
          uniform = true;
        } else throw std::runtime_error("to many \':\' in axis definition");
      }
      *it = '\0';
      edges.push_back(std::atof(str));
      if (c=='\0') break;
      str = ++it;
    } else ++it;
  }
  if (uniform) {
    if (edges.size()!=3 || edges[1] > edges[2])
      throw std::runtime_error("bad axis definition");
    axes.emplace(
      std::move(name),
      new ivanp::uniform_axis<double,true>(edges[0],edges[1],edges[2])
    );
  } else {
    axes.emplace(
      std::move(name),
      new ivanp::container_axis<std::vector<double>,true>(std::move(edges))
    );
  }
}

int main(int argc, char* argv[]) {
  std::vector<const char*> ifnames;
  const char* ofname;

  try {
    using namespace ivanp::po;
    if (program_options()
      (ifnames,'i',"input files",req(),pos())
      (ofname,'o',"output file name",req())
      (axes,{"-b","--binning"},"",multi(),binning_parser)
      .parse(argc,argv,true)) return 0;
  } catch (const std::exception& e) {
    cerr << e << endl;
    return 1;
  }

  // { auto it = axes.find("pT_yy");
  //   if (it == axes.end()) axes.emplace(
  //     "pT_yy", new ivanp::uniform_axis<double,true>(100,0.,500.)
  //   );
  // }

  hist h_pT_yy("pT_yy",axes.at("pT_yy").get());

  bool is_mc;
  TNamed *lumi_str;
  std::string type;
  for (const char* ifname : ifnames) {
    TFile f(ifname);
    cout << f.GetName() << endl;
    if (f.IsZombie()) return 1;

    lumi_str = nullptr;
    f.GetObject("Lumi",lumi_str);
    if (lumi_str) {
      cout << "Lumi: " << lumi_str->GetTitle() << endl;
      lumi += std::atof(lumi_str->GetTitle());
      type = get_data_str(ifname);
      is_mc = false;
      weight = 1;
    } else {
      type = get_mc_str(ifname);
      is_mc = true;
    }
    cout << type << endl;

    TTreeReader reader;
    reader.SetTree("hgam",&f);
    boost::optional<TTreeReaderValue<Float_t>> _weight;
    if (is_mc) _weight.emplace(reader,"weight");
    TTreeReaderValue<Float_t>
      _m_yy  (reader,"m_yy"),
      _pT_yy (reader,"pT_yy");

    auto it = bin::m.find(type);
    if (it == bin::m.end()) {
      const auto n = bin::m.size();
      it = bin::m.emplace(type,n).first;
      for (auto& b : bin::all) b->v.push_back(0);
    }
    bin::i = it->second;

    using tc = ivanp::timed_counter<Long64_t>;
    for (tc ent(reader.GetEntries(true)); reader.Next(); ++ent) {
      const auto m_yy = *_m_yy;
      const bool in_window = (121.<m_yy && m_yy<129.);
      if (is_mc) { // signal from MC
        if (!in_window) continue;
        weight = **_weight;
      } else { // background from data
        if (in_window) continue;
      }
      h_pT_yy(*_pT_yy);
    }
  }

  TFile fout(ofname,"recreate");
  cout << fout.GetName() << endl;
  if (fout.IsZombie()) return 1;

  // cout << bin::m.size() << endl;

  using ivanp::root::to_root;
  // using ivanp::root::slice_to_root;

  cout << "\nSaving histograms:\n";

  for (const auto& m : bin::m) {
    cout << m.first << endl;
    bin::i = m.second;
    fout.mkdir(m.first.c_str())->cd();

    for (auto& h : hist::all) to_root(*h,h.name);
  }

  fout.mkdir("all_MC")->cd();
  for (auto& h : hist::all) to_root(*h,h.name,all_mc{});

  fout.Write();
}
