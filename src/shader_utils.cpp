#include "utils.hpp"


std::pair<GLint, std::string>
loadShaderSource(
    const std::string& filepath
)
{
    std::ifstream file(filepath);
    if (!file)
    {
        std::cerr << "ERROR: Could not open shader file: "
                  << filepath
                  << '\n';
        return {0, ""};
    }

    std::stringstream buffer;
    buffer << file.rdbuf();
    return {1, buffer.str()};
}

std::pair<GLint, GLuint>
compileShaderFromFile(
    GLenum type,
    const std::string& filepath
)
{
    auto [success, source] = loadShaderSource(filepath);
    if (!success) {
        return {success, 0};
    }
    const char* sourceCStr = source.c_str();

    GLuint shader = glCreateShader(type);
    glShaderSource(shader, 1, &sourceCStr, NULL);
    glCompileShader(shader);

    glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
    if (!success)
    {
        char infoLog[512];
        glGetShaderInfoLog(shader, 512, NULL, infoLog);
        std::cerr << "ERROR::SHADER_COMPILATION_FAILED\n"
                  << infoLog
                  << '\n';
    }

    return {success, shader};
}

GLint
addShaderToProgram(
    const GLuint shaderProg,
    const std::string&& shaderPath,
    const GLenum type
)
{
    auto [success, shader] = compileShaderFromFile(type, shaderPath);
    if (!success)   return success;

    glAttachShader(shaderProg, shader);
    glLinkProgram(shaderProg);

    glGetProgramiv(shaderProg, GL_LINK_STATUS, &success);
    if (!success)
    {
        char infoLog[512];
        glGetProgramInfoLog(shaderProg, 512, NULL, infoLog);
        std::cerr << "ERROR::PROGRAM_LINKING_FAILED\n"
                  << infoLog
                  << '\n';
    }

    glDeleteShader(shader);

    return success;
}
