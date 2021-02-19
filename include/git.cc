#include "git.h"

bool GitMetadata::Populated() { return true; }
bool GitMetadata::AnyUncommittedChanges() { return true; }
std::string GitMetadata::AuthorName() { return "schd"; }
std::string GitMetadata::AuthorEmail() { return "sacha.duverger@inrae.fr"; }
std::string GitMetadata::CommitSHA1() { return "5ad8c545a3c20ee4d43d89e53013080b9512dce5"; }
std::string GitMetadata::CommitDate() { return "2021-02-19 19:05:51 +0100"; }
std::string GitMetadata::CommitSubject() { return ":hammer: Update git.cc"; }
std::string GitMetadata::CommitBody() { return ""; }
std::string GitMetadata::Describe() { return "5ad8c545"; }
