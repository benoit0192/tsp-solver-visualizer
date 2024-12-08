#ifndef UTILS_H
#define UTILS_H

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <filesystem>

#include <GL/glew.h>
#include <GLFW/glfw3.h>


namespace fs = std::filesystem;

bool
isValidFilePath(
    const std::string& filepath
);

bool
isValidDirectoryPath(
    const std::string& dirPath
);
std::pair<GLint, std::string>
loadShaderSource(
    const std::string& filepath
);

std::pair<GLint, GLuint>
compileShaderFromFile(
    GLenum type,
    const std::string& filepath
);

GLint
addShaderToProgram(
    const GLuint shaderProg,
    const std::string&& shaderPath,
    const GLenum type
);

#endif // !UTILS_H
