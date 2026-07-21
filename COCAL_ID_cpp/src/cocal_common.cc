#include "cocal_id_cpp_internal.hh"
#include "cctk.h"

namespace cocal_id_cpp {

[[noreturn]] void abort_cocal(const string &msg) {
  CCTK_WARN(CCTK_WARN_ABORT, msg.c_str());
  throw std::runtime_error(msg);
}

vector<string> data_tokens(string line) {
  const auto hash = line.find('#');
  if (hash != string::npos)
    line = line.substr(0, hash);
  const auto colon = line.find(':');
  if (colon != string::npos)
    line = line.substr(0, colon);
  std::istringstream iss(line);
  vector<string> tokens;
  string tok;
  while (iss >> tok)
    tokens.push_back(tok);
  return tokens;
}

vector<string> next_tokens(std::istream &is) {
  string line;
  while (std::getline(is, line)) {
    auto tokens = data_tokens(line);
    if (!tokens.empty())
      return tokens;
  }
  return {};
}

namespace {
string context_suffix(const string &context) {
  return context.empty() ? string{} : " while reading " + context;
}
} // namespace

double parse_real(string tok, const string &context) {
  for (char &c : tok) {
    if (c == 'd' || c == 'D')
      c = 'e';
  }
  size_t parsed = 0;
  try {
    const double value = std::stod(tok, &parsed);
    if (parsed != tok.size())
      abort_cocal("Malformed COCAL real token '" + tok + "'" + context_suffix(context));
    return value;
  } catch (const std::exception &) {
    abort_cocal("Could not parse COCAL real token '" + tok + "'" + context_suffix(context));
  }
}

int parse_int(const string &tok, const string &context) {
  size_t parsed = 0;
  try {
    const int value = std::stoi(tok, &parsed);
    if (parsed != tok.size())
      abort_cocal("Malformed COCAL integer token '" + tok + "'" + context_suffix(context));
    return value;
  } catch (const std::exception &) {
    abort_cocal("Could not parse COCAL integer token '" + tok + "'" + context_suffix(context));
  }
}

RnsData rns;
BnsData bns;
BhtData bht;
std::mutex rns_mutex;
std::mutex bns_mutex;
std::mutex bht_mutex;

} // namespace cocal_id_cpp
