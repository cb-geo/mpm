#include "git.h"

bool GitMetadata::Populated() {
    return true;
}
bool GitMetadata::AnyUncommittedChanges() {
    return true;
}
std::string GitMetadata::AuthorName() {
    return "schd";
}
std::string GitMetadata::AuthorEmail() {
    return "sacha.duverger@inrae.fr";
}
std::string GitMetadata::CommitSHA1() {
    return "e4c44e702eca1d5ecb667573a60fee852a119357";
}
std::string GitMetadata::CommitDate() {
    return "2021-02-19 10:30:55 +0100";
}
std::string GitMetadata::CommitSubject() {
    return ":construction: Add Andrew Hardin's scripts for embed git metadata";
}
std::string GitMetadata::CommitBody() {
    return "";
}
std::string GitMetadata::Describe() {
    return "e4c44e70";
}

