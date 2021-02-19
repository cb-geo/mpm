#include "git.h"

bool GitMetadata::Populated() {
    return true;
}
bool GitMetadata::AnyUncommittedChanges() {
    return false;
}
std::string GitMetadata::AuthorName() {
    return "schd";
}
std::string GitMetadata::AuthorEmail() {
    return "sacha.duverger@inrae.fr";
}
std::string GitMetadata::CommitSHA1() {
    return "4a2bc8fc1d997ae0c906c4442d4a7210cc94082d";
}
std::string GitMetadata::CommitDate() {
    return "2021-02-19 13:59:25 +0100";
}
std::string GitMetadata::CommitSubject() {
    return ":wrench: clang format";
}
std::string GitMetadata::CommitBody() {
    return "";
}
std::string GitMetadata::Describe() {
    return "4a2bc8fc";
}

