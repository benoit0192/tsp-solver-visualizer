#include "utils.hpp"


bool
isValidFilePath(
    const std::string& filepath
)
{
    fs::path path(filepath);
    return fs::exists(path) && fs::is_regular_file(path);
}

bool
isValidDirectoryPath(
    const std::string& dirPath
)
{
    fs::path path(dirPath);
    return fs::exists(path) && fs::is_directory(path);
}

