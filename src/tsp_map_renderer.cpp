#include <vector>
#include <climits>
#include <fstream>
#include <sstream>
#include <iostream>
#include <filesystem>

#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>

#include "geo.hpp"
#include "utils.hpp"


#define WINDOW_WIDTH  600
#define WINDOW_HEIGHT 600

#define MIN_CANVAS_WIDTH  400
#define MIN_CANVAS_HEIGHT 400


enum ShaderType
{
    SHADER_VERTICES,
    SHADER_EDGES,
    SHADER_COUNT
};

typedef struct AppState
{
    size_t windowWidth    =600;
    size_t windowHeight   =600;
    size_t minCanvasWidth =400;
    size_t minCanvasHeight=400;
    GLuint shaderProgs[SHADER_COUNT]={0};
} AppState;

void
keyCallback(
    GLFWwindow* window,
    int key,
    [[maybe_unused]]int scancode,
    int action,
    [[maybe_unused]]int mods
)
{
    if(action==GLFW_PRESS) {
        if(key==GLFW_KEY_Q ||  key==GLFW_KEY_ESCAPE) {
            glfwSetWindowShouldClose(window, GLFW_TRUE);
        }
    }
}

void
framebufferSizeCallback(
    GLFWwindow* window,
    int win_w,
    int win_h
)
{
    AppState* state = static_cast<AppState*>(glfwGetWindowUserPointer(window));

    int draw_w = std::max(win_w, (int)state->minCanvasWidth);
    int draw_h = std::max(win_h, (int)state->minCanvasHeight);
    int mn = std::min(draw_w, draw_h);
    /* Square drawing window */
    draw_w = mn;
    draw_h = mn;
    /* Drawing window centering */
    int xOff = (win_w - draw_w) * 0.5;
    int yOff = (win_h - draw_h) * 0.5;
    glViewport(xOff, yOff, draw_w, draw_h);

    // Update the uniform scale factor in the vertex shader
    float scaleFactor = (mn < (int)state->minCanvasWidth)
                                    ? 1.0
                                    : (float)mn / state->minCanvasWidth;
    GLuint shPg = state->shaderProgs[SHADER_VERTICES];
    glUseProgram(shPg);
    GLint scaleUniform = glGetUniformLocation(shPg, "f_ScaleFactor");
    glUniform1f(scaleUniform, scaleFactor);
}

typedef struct gl_Vertices
{
    GLuint VAO=0;
    GLuint VBO=0;
    size_t count=0;
} gl_Vertices;

GLint
initVertices(
    gl_Vertices& gl_vs,
    std::vector<GLfloat>& vs
)
{
    if (vs.size()%3 != 0) {
        std::cerr << "ERROR: Expected vertices size to be a multiple of 2.\n";
        return 0;
    }
    gl_vs.count=vs.size()/3;

    glGenVertexArrays(1, &gl_vs.VAO);
    glBindVertexArray(gl_vs.VAO);

    glGenBuffers(1, &gl_vs.VBO);
    glBindBuffer(GL_ARRAY_BUFFER, gl_vs.VBO);
    glBufferData(GL_ARRAY_BUFFER,
                 vs.size()*sizeof(GLfloat),
                 vs.data(),
                 GL_STATIC_DRAW);

    glVertexAttribPointer(0,
                          3,
                          GL_FLOAT,
                          GL_FALSE,
                          3*sizeof(GL_FLOAT),
                          (GLvoid*)0);
    glEnableVertexAttribArray(0);
    /* Reset */
    glBindVertexArray(0);

    return 1;
}

typedef struct gl_Edges
{
    GLuint VAO=0;
    GLuint VBO=0;
    size_t count=0;
} gl_Edges;

GLint
initEdges(
    gl_Edges& gl_es,
    std::vector<GLfloat>& es
)
{
    //WARNING: Would adjust accordingly if dealing with 3D data
    if (es.size()%4 != 0) {
        std::cerr << "ERROR: Expected number of edges to be a multiple of 4.\n";
        return 0;
    }

    std::vector<GLfloat> quadEdges;
    //WARNING: Would need to increase the looping offset if dealing with 3D data
    for (size_t i=0; i<es.size(); i+=4) {
        std::vector<GLfloat> quadVertices = generatePrismVertices(
                                           glm::vec3(es[i+0], es[i+1], 0.0), // No 3d data yet
                                           glm::vec3(es[i+2], es[i+3], 0.0),
                                           0.03f,
                                           0.03f);
        quadEdges.insert(quadEdges.end(),
                         quadVertices.begin(),
                         quadVertices.end());
    }
    // Nums of triangles (2 triangles per edge)
    gl_es.count=quadEdges.size()*0.5;

    glGenVertexArrays(1, &gl_es.VAO);
    glBindVertexArray(gl_es.VAO);

    glGenBuffers(1, &gl_es.VBO);
    glBindBuffer(GL_ARRAY_BUFFER, gl_es.VBO);
    glBufferData(GL_ARRAY_BUFFER,
                 quadEdges.size()*sizeof(GLfloat),
                 quadEdges.data(),
                 GL_STATIC_DRAW);

    glVertexAttribPointer(0,
                          3,
                          GL_FLOAT,
                          GL_FALSE,
                          3*sizeof(GL_FLOAT),
                          (GLvoid*)0);
    glEnableVertexAttribArray(0);
    /* Reset */
    glBindVertexArray(0);

    return 1;
}

void
renderVertices(
    AppState& appState,
    gl_Vertices& gl_vs
)
{
    glUseProgram(appState.shaderProgs[SHADER_VERTICES]);
    glBindVertexArray(gl_vs.VAO);
    glDrawArrays(GL_POINTS, 0, gl_vs.count);
}

void
renderEdges(
    AppState& appState,
    gl_Edges& gl_es
)
{
    glUseProgram(appState.shaderProgs[SHADER_EDGES]);
    glBindVertexArray(gl_es.VAO);
    glDrawArrays(GL_TRIANGLES, 0, gl_es.count);
}

bool
loadData(
    std::string filename,
    std::vector<std::vector<float>>& nodes_dst,
    std::vector<int>& path_dst
)
{
    std::ifstream inFile(filename);
    if (!inFile) {
        std::cerr << "ERROR: Failed to open the file" << '\n';
        return false;
    }
    constexpr size_t MAX_DIM=3; // Support up to 3D
    std::string line;
    size_t curLineNum=1;
    while (std::getline(inFile, line)) {
        if (line.starts_with('#')) {
            //Do nothing
        }
        /* Decode nodes */
        else if (line.contains(':')) {
            std::stringstream ss(line);
            std::vector<float> x_nd;
            std::string compStr;
            size_t dimCnt=0;
            while (std::getline(ss, compStr, ':') && dimCnt < MAX_DIM) {
                try {
                    float comp = std::stof(compStr);
                    x_nd.push_back(comp);
                } catch (const std::exception& e) {
                    std::cerr << "ERROR: Failed to parse '"
                              << compStr << "' at line "
                              << curLineNum << ": " << e.what();
                    return false;
                }
                dimCnt++;
            }
            if (!nodes_dst.empty() && nodes_dst.back().size() != x_nd.size()) {
                std::cerr << "ERROR: Inconsistent dimension size ("
                          << nodes_dst.back().size() << " != " << x_nd.size()
                          << ") at line " << curLineNum;
                return false;
            }
            if (x_nd.size() < 2) {
                std::cerr << "ERROR: Node dimension too small at line " 
                          << curLineNum << '\n';
                return false;
            }
            nodes_dst.emplace_back(std::move(x_nd));
        }
        /* Decode path */
        else if (line.contains('-')) {
            std::istringstream iss(line);
            std::string numStr;
            while (std::getline(iss, numStr, '-')) {
                try {
                    int num = std::stoi(numStr);
                    path_dst.push_back(num);
                } catch (const std::exception& e) {
                    std::cerr << "ERROR: Failed to parse '"
                              << numStr << "' at line "
                              << curLineNum << ": " << e.what();
                    return false;
                }
            }
        }
        curLineNum++;
    }
    if (nodes_dst.empty()) {
        std::cerr << "ERROR: Number of nodes is empty.\n";
        return false;
    }
    if (path_dst.empty()) {
        std::cerr << "ERROR: Path is empty.\n";
        return false;
    }
    return true;
}
/**
 * This function processes a set of nodes represented as vectors of floats.
 * It normalizes each node's components to a [-1, 1] range while applying
 * an optional padding to shrink the range.
 */
std::vector<float>
normNodes(
    std::vector<std::vector<float>>& nodes,
    float padding = 0.0
)
{
    if (nodes.empty())
        return {};
    const size_t dim = nodes.front().size();
    constexpr float MIN_PAD = 0.0, MAX_PAD = 0.4;
    if (padding < MIN_PAD || MAX_PAD < padding) {
        std::cerr << "Padding exceeds the allowed range [0.0, 0.4]."
                  << "It will be adjusted to fit within the allowed range.\n";
        padding = std::clamp(padding, MIN_PAD, MAX_PAD);
    }
    using Limits = std::numeric_limits<float>;
    std::vector<std::pair<float, float>> mn_mx_pairs(dim,
                                            {Limits::max(), Limits::lowest()});
    for (auto& x_nd : nodes) {
        for (size_t i=0; i<dim; ++i) {
            auto& [mn, mx] = mn_mx_pairs[i];
            mn = std::min(mn, x_nd[i]);
            mx = std::max(mx, x_nd[i]);
        }
    }
    const float scaleFactor = 2.0f * (1 - 2 * padding);
    const float offset      = 1.0f - 2 * padding;
    std::vector<float> procNodes;
    for (auto& x_nd : nodes) {
        for (size_t i=0; i<x_nd.size(); ++i) {
            auto& [mn, mx] = mn_mx_pairs[i];
            float nx  = x_nd[i] - mn;
            if (mx - mn > 0.0f) {
                nx = nx * scaleFactor / (mx - mn) - offset;
            }
            procNodes.push_back(nx);
        }
    }
    return procNodes;
}

std::vector<float>
extractEdgeNodes(
    std::vector<int>& path,
    std::vector<float>& nodes1d
)
{
    std::vector<float> edgeNodes;
    for (size_t i=0; i<path.size()-1; i++) {
        int e1Ix = path[i], e2Ix = path[i+1];
        /* edgeNodes.push_back(nodes1d[2*e1Ix]); */
        /* edgeNodes.push_back(nodes1d[2*e1Ix+1]); */
        /* edgeNodes.push_back(nodes1d[2*e2Ix]); */
        /* edgeNodes.push_back(nodes1d[2*e2Ix+1]); */
        edgeNodes.push_back(nodes1d[3*e1Ix]);
        edgeNodes.push_back(nodes1d[3*e1Ix+1]);
        edgeNodes.push_back(nodes1d[3*e2Ix]);
        edgeNodes.push_back(nodes1d[3*e2Ix+1]);
    }
    return edgeNodes;
}

GLFWwindow*
initOpenGLWindow(void)
{
    if(!glfwInit()) {
        return nullptr;
    }

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    GLFWwindow* window = glfwCreateWindow(WINDOW_WIDTH,
                                          WINDOW_HEIGHT,
                                          "TSL Path Map",
                                          NULL,
                                          NULL);
    if(!window) {
        glfwTerminate();
        return nullptr;
    }

    glfwMakeContextCurrent(window);
    glEnable(GL_PROGRAM_POINT_SIZE);

    if(glewInit() != GLEW_OK) {
        std::cerr << "Failed to initialize GLEW" << '\n';
        return nullptr;
    }

    return window;
}

void
cleanOpenGLContext(GLFWwindow *window)
{
    glfwDestroyWindow(window);
    glfwTerminate();
}

int main(int argc, char *argv[])
{
    std::cout << "== TSP Map Rendering ==" << '\n';

    std::string usage = std::string("Usage:\n")
         + "  ./map_renderer <DATA_FILEPATH>\n"
         + "  <DATA_FILEPATH> is the data filepath generated by the TSP solver program.\n";

    if (argc != 2) {
        std::cerr << "ERROR: Unexpected number of input arguments.\n";
        std::cout << usage << '\n';
        return -1;
    }
    std::string filepath = argv[1];
    if(!isValidFilePath(filepath)) {
        std::cerr << "ERROR: Provided data filepath does not exist:\n"
                  << "       " << filepath << '\n';
        return -1;
    }

    std::vector<std::vector<float>> nodes;
    std::vector<int> path;
    GLint success;
    success = loadData(filepath, nodes, path);
    if (!success) {
        return -1;
    }
    std::vector<float> vertices = normNodes(nodes, 0.1);
    std::vector<float> edges    = extractEdgeNodes(path, vertices);

    AppState appState = {
            .windowWidth=WINDOW_WIDTH,
            .windowHeight=WINDOW_HEIGHT,
            .minCanvasWidth=MIN_CANVAS_WIDTH,
            .minCanvasHeight=MIN_CANVAS_HEIGHT,
    };
    GLFWwindow* window = initOpenGLWindow();
    if (!window) {
        cleanOpenGLContext(window);
        return -1;
    }
    const GLubyte* version = glGetString(GL_VERSION);
    std::cout << "OpenGL Version: " << version << '\n';

    glfwSetWindowUserPointer(window, &appState);

    GLuint shaderProg = glCreateProgram();
    success = addShaderToProgram(shaderProg,
                                 "shader/vertex_shader_nodes.glsl",
                                 GL_VERTEX_SHADER);
    if (!success) {
        cleanOpenGLContext(window);
        return -1;
    }
    success = addShaderToProgram(shaderProg,
                                 "shader/fragment_shader_nodes.glsl",
                                 GL_FRAGMENT_SHADER);
    if (!success) {
        cleanOpenGLContext(window);
        return -1;
    }
    appState.shaderProgs[SHADER_VERTICES] = shaderProg;

    GLuint shaderProgLines = glCreateProgram();
    success = addShaderToProgram(shaderProgLines,
                                 "shader/vertex_shader_edges.glsl",
                                 GL_VERTEX_SHADER);
    if (!success) {
        cleanOpenGLContext(window);
        return -1;
    }
    success = addShaderToProgram(shaderProgLines,
                                 "shader/fragment_shader_edges.glsl",
                                 GL_FRAGMENT_SHADER);
    if (!success) {
        cleanOpenGLContext(window);
        return -1;
    }
    appState.shaderProgs[SHADER_EDGES] = shaderProgLines;
    /* Vertices data */
    gl_Vertices gl_vs;
    success = initVertices(gl_vs, vertices);
    if (!success) {
        cleanOpenGLContext(window);
        return -1;
    }
    /* Edges data */
    gl_Edges gl_es;
    success = initEdges(gl_es, edges);
    if (!success) {
        cleanOpenGLContext(window);
        return -1;
    }
    /* Callbacks */
    glfwSetFramebufferSizeCallback(window, framebufferSizeCallback);
    glfwSetKeyCallback(window, keyCallback);

    /* Call the resize window callback at least once */
    int win_w, win_h;
    glfwGetFramebufferSize(window, &win_w, &win_h);
    framebufferSizeCallback(window, win_w, win_h);

    /* Rendering loop */
    while(!glfwWindowShouldClose(window)) {
        glClear(GL_COLOR_BUFFER_BIT);

        renderVertices(appState, gl_vs);
        renderEdges(appState, gl_es);

        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    glDeleteVertexArrays(1, &gl_vs.VAO);
    glDeleteBuffers(1, &gl_vs.VBO);

    glDeleteVertexArrays(1, &gl_es.VAO);
    glDeleteBuffers(1, &gl_vs.VBO);

    cleanOpenGLContext(window);

    return 0;
}
