#include "git.h"

bool GitMetadata::Populated() { return true; }
bool GitMetadata::AnyUncommittedChanges() { return false; }
std::string GitMetadata::AuthorName() { return "schd"; }
std::string GitMetadata::AuthorEmail() { return "sacha.duverger@inrae.fr"; }
std::string GitMetadata::CommitSHA1() { return "ffbec37e9c680ba62ac9b0b5f94bac70de6ee0e6";}
std::string GitMetadata::CommitDate() { return "2021-02-19 18:59:00 +0100"; }
std::string GitMetadata::CommitSubject() { return ":hammer: Manually format git.cc.in"; }
std::string GitMetadata::CommitBody() { return ""; }
std::string GitMetadata::Describe() { return "ffbec37e"; }

