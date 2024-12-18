#include <vector>
#include <climits>
#include <sstream>
#include <iostream>

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
    size_t windowWidth    = 600;
    size_t windowHeight   = 600;
    size_t minCanvasWidth = 400;
    size_t minCanvasHeight= 400;
    GLuint shaderProgs[SHADER_COUNT]={0};
} AppState;

typedef struct gl_Data
{
    GLuint VAO=0;
    GLuint VBO=0;
    GLuint EBO=0;
    size_t count=0;     // Number of vertices
    size_t vSize=0;     // Total size of flattened vertices data
    size_t iSize=0;     // Total size of flattened indices data
    size_t compSize=0;  // Components per vertex (Vec2 or Vec3)
} gl_Data;


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
std::vector<std::vector<float>>
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
    std::vector<std::vector<float>> procNodes;
    procNodes.reserve(nodes.size() * dim);
    for (auto& x_nd : nodes) {
        std::vector<float> x_nd_n;
        x_nd.reserve(x_nd.size());
        for (size_t i=0; i<x_nd.size(); ++i) {
            auto& [mn, mx] = mn_mx_pairs[i];
            float nx  = x_nd[i] - mn;
            if (mx - mn > 0.0f) {
                nx = nx * scaleFactor / (mx - mn) - offset;
            }
            x_nd_n.push_back(nx);
        }
        procNodes.push_back(std::move(x_nd_n));
    }
    return procNodes;
}

std::vector<float>
extractEdgeNodes(
    const std::vector<int>& path,
    const std::vector<std::vector<float>>& centers
)
{
    size_t dim = 0;
    if (!centers.empty()) {
        dim = centers.front().size();
    }
    std::vector<float> edgeNodes;
    edgeNodes.reserve((path.size()-1) * 2 * dim);
    for (size_t i=0; i<path.size()-1; i++) {
        int e1Ix = path[i], e2Ix = path[i+1];

        edgeNodes.insert(edgeNodes.end(),
                         centers[e1Ix].begin(),
                         centers[e1Ix].end());

        edgeNodes.insert(edgeNodes.end(),
                         centers[e2Ix].begin(),
                         centers[e2Ix].end());
    }
    return edgeNodes;
}

void
preprocessEdges(
    const std::vector<int>& path,
    const std::vector<std::vector<float>>& node_centers,
    std::vector<float>& edge_vertices,
    std::vector<unsigned int>& edge_indices
)
{
    std::vector<float> edges = extractEdgeNodes(path, node_centers);

    size_t dim = 0;
    if (!node_centers.empty()) {
        dim = node_centers.front().size();
    }
    const float width = 0.05;
    generatePrisms(edges, width, width, dim, edge_vertices, edge_indices);
}

gl_Data
initGLData(
    const std::vector<GLfloat>& vertices,
    const std::vector<GLuint>& indices,
    size_t compSize
)
{

    gl_Data gl_data  = {};
    gl_data.compSize = compSize;
    gl_data.vSize    = vertices.size();
    gl_data.iSize    = indices.size();
    gl_data.count    = vertices.size() / compSize;

    glGenVertexArrays(1, &gl_data.VAO);
    glBindVertexArray(gl_data.VAO);

    glGenBuffers(1, &gl_data.VBO);
    glBindBuffer(GL_ARRAY_BUFFER, gl_data.VBO);
    glBufferData(GL_ARRAY_BUFFER,
                 vertices.size() * sizeof(GLfloat),
                 vertices.data(),
                 GL_STATIC_DRAW);

    glGenBuffers(1, &gl_data.EBO);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, gl_data.EBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER,
                 indices.size() * sizeof(GLuint),
                 indices.data(),
                 GL_STATIC_DRAW);

    glVertexAttribPointer(0,
                          compSize,
                          GL_FLOAT,
                          GL_FALSE,
                          compSize * sizeof(GL_FLOAT),
                          (GLvoid*)0);
    glEnableVertexAttribArray(0);

    glBindVertexArray(0);

    return gl_data;
}

void
clearGLData(
    gl_Data& gl_data
)
{
    glDeleteVertexArrays(1, &gl_data.VAO);
    glDeleteBuffers(1, &gl_data.VBO);
    glDeleteBuffers(1, &gl_data.EBO);
}

void
renderVertices(
    AppState& appState,
    gl_Data& gl_vs
)
{
    glUseProgram(appState.shaderProgs[SHADER_VERTICES]);
    glBindVertexArray(gl_vs.VAO);
    glDrawElements(GL_TRIANGLES, gl_vs.iSize, GL_UNSIGNED_INT, 0);
}

void
renderEdges(
    AppState& appState,
    gl_Data& gl_es
)
{
    glUseProgram(appState.shaderProgs[SHADER_EDGES]);
    glBindVertexArray(gl_es.VAO);
    glDrawElements(GL_TRIANGLES, gl_es.iSize, GL_UNSIGNED_INT, 0);
}

bool
initShaders(
    AppState& appState
)
{
    /* Nodes shaders */
    GLuint shaderProg = glCreateProgram();
    if (!addShaderToProgram(shaderProg,
                            "shader/vertex_shader_nodes.glsl",
                            GL_VERTEX_SHADER))
    {
        return false;
    }
    if (!addShaderToProgram(shaderProg,
                            "shader/fragment_shader_nodes.glsl",
                            GL_FRAGMENT_SHADER))
    {
        return false;
    }
    appState.shaderProgs[SHADER_VERTICES] = shaderProg;

    /* Edges shaders */
    GLuint shaderProgLines = glCreateProgram();
    if (!addShaderToProgram(shaderProgLines,
                            "shader/vertex_shader_edges.glsl",
                            GL_VERTEX_SHADER))
    {
        return false;
    }
    if (!addShaderToProgram(shaderProgLines,
                            "shader/fragment_shader_edges.glsl",
                            GL_FRAGMENT_SHADER))
    {
        return false;
    }
    appState.shaderProgs[SHADER_EDGES] = shaderProgLines;

    return true;
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

bool parseArgs(
    int argc,
    char *argv[],
    std::string& filepath
)
{
    std::string usage = std::string("Usage:\n")
         + "  ./map_renderer <DATA_FILEPATH>\n"
         + "  <DATA_FILEPATH> is the data filepath generated by the TSP solver program.\n";

    if (argc != 2) {
        std::cerr << "ERROR: Unexpected number of input arguments.\n";
        std::cout << usage << '\n';
        return false;
    }
    filepath = argv[1];
    if(!isValidFilePath(filepath)) {
        std::cerr << "ERROR: Provided data filepath does not exist:\n"
                  << "       " << filepath << '\n';
        return false;
    }
    return true;
}

int main(int argc, char *argv[])
{
    std::cout << "== TSP Map Rendering ==" << '\n';

    std::string filepath;
    if (!parseArgs(argc, argv, filepath)) {
        return -1;
    }

    std::vector<std::vector<float>> node_centers;
    std::vector<int> path;
    if (!loadData(filepath, node_centers, path)) {
        return -1;
    }

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

    if (!initShaders(appState)) {
        cleanOpenGLContext(window);
        return -1;
    }

    /* Init nodes data */
    node_centers = normNodes(node_centers, 0.1);
    std::vector<float> node_vertices;
    std::vector<unsigned int> node_indices;
    float radius=0.07;
    generateSpheres(radius, 20, 20, node_centers, node_vertices, node_indices);
    gl_Data gl_vs = initGLData(node_vertices, node_indices, node_centers.front().size());

    /* Init edges data */
    std::vector<float> edge_vertices;
    std::vector<unsigned int> edge_indices;
    preprocessEdges(path, node_centers, edge_vertices, edge_indices);
    gl_Data gl_es = initGLData(edge_vertices, edge_indices, node_centers.front().size());

    /* Callbacks */
    glfwSetFramebufferSizeCallback(window, framebufferSizeCallback);
    glfwSetKeyCallback(window, keyCallback);

    /* Call the resize window callback at least once */
    int win_w, win_h;
    glfwGetFramebufferSize(window, &win_w, &win_h);
    framebufferSizeCallback(window, win_w, win_h);

    /* Rendering loop */
    while (!glfwWindowShouldClose(window)) {
        glClear(GL_COLOR_BUFFER_BIT);

        renderVertices(appState, gl_vs);
        renderEdges(appState, gl_es);

        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    clearGLData(gl_vs);
    clearGLData(gl_es);

    cleanOpenGLContext(window);

    return 0;
}
