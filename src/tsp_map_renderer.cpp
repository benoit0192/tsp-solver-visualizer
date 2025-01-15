#include <vector>
#include <thread>
#include <climits>
#include <sstream>
#include <iostream>

#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/quaternion.hpp>
#include <glm/gtc/matrix_transform.hpp>

#include "geo.hpp"
#include "utils.hpp"


enum ShaderType
{
    NODES_SHADER,
    EDGES_SHADER,
    SHADER_COUNT
};

typedef struct WindowConfig {
    size_t windowWidth;
    size_t windowHeight;
} WindowConfig;

typedef struct RenderState {
    GLuint shaderProgs[SHADER_COUNT];
    float aspectRatio;
} RenderState;

typedef struct gl_Data
{
    GLuint VAO;
    GLuint VBO;         // Vertices buffer
    GLuint NBO;         // Normals buffer
    GLuint EBO;         // Indices buffer
    size_t count;       // Number of vertices
    size_t vSize;       // Total size of flattened vertices data
    size_t iSize;       // Total size of flattened indices data
    size_t compSize;    // Components per vertex (Vec2 or Vec3)
} gl_Data;

typedef struct CamParams
{
    float fov;
    float near;
    float far;
    glm::vec3 pos;     // Camera position
    glm::vec3 target;  // Point camera looks at
    glm::vec3 upAxis;  // World up vector
} CamParams;

typedef struct InteractionState {
    glm::mat4 rotationMatrix;
    glm::vec3 lastPosOnSphere;
    bool isDragging;
    float zoomSpeed;
    CamParams* camParams;
    float maxCamDist;
    float minCamDist;
    bool isRotating;
} InteractionState;

typedef struct AppContext {
    WindowConfig* windowConfig;
    RenderState* renderState;
    InteractionState* interactionState;
} AppContext;

typedef struct SphereParams {
    const float radius;
    const int stacks;
    const int slices;
} SphereParams;

typedef struct PrismParams {
    const float edgeWidth;
    const int numEdgeSamples;
} PrismParams;


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
    AppContext* ctx = static_cast<AppContext*>(glfwGetWindowUserPointer(window));
    glViewport(0, 0, win_w, win_h);
    ctx->renderState->aspectRatio = float(win_w) / float(win_h);
}

/* Mouse button callback */
void
mouseButtonCallback(
    GLFWwindow* window,
    int button,
    int action,
    [[maybe_unused]]int mods
)
{
    AppContext* ctx = static_cast<AppContext*>(glfwGetWindowUserPointer(window));
    InteractionState* interaction = ctx->interactionState;

    /* Stop the initial rotation animation */
    interaction->isRotating = false;

    if (button == GLFW_MOUSE_BUTTON_LEFT) {
        if (action == GLFW_PRESS) {
            interaction->isDragging = true;

            double xpos, ypos;
            glfwGetCursorPos(window, &xpos, &ypos);

            /* Map to sphere */
            int width, height;
            glfwGetWindowSize(window, &width, &height);
            interaction->lastPosOnSphere = mapToSphere(xpos, ypos, width, height);
        }
        else if (action == GLFW_RELEASE) {
            interaction->isDragging = false;
        }
    }
}

/* Mouse drag callback */
void
cursorPositionCallback(
    GLFWwindow* window,
    double xpos,
    double ypos
)
{
    AppContext* ctx = static_cast<AppContext*>(glfwGetWindowUserPointer(window));
    InteractionState* interaction = ctx->interactionState;

    if (!interaction->isDragging) return;

    /* Map current mouse position to sphere */
    int width, height;
    glfwGetWindowSize(window, &width, &height);
    glm::vec3 curPos = mapToSphere(xpos, ypos, width, height);

    /* Compute rotation axis */
    glm::vec3 axis = glm::cross(interaction->lastPosOnSphere, curPos);
    if (glm::length(axis) > 1e-6) { // Avoid numerical instability
        axis = glm::normalize(axis);

        float angle = std::acos(glm::clamp(glm::dot(interaction->lastPosOnSphere,
                                                    curPos),
                                           -1.0f, 1.0f));

        glm::quat rotation = glm::angleAxis(angle, axis);
        interaction->rotationMatrix = glm::toMat4(rotation) * interaction->rotationMatrix;

        interaction->lastPosOnSphere = curPos;
    }
}

/* Mouse scroll callback */
void
scrollCallback(
    GLFWwindow* window,
    [[maybe_unused]] double xOffset,
    double yOffset
)
{
    AppContext* ctx = static_cast<AppContext*>(glfwGetWindowUserPointer(window));
    InteractionState* interaction = ctx->interactionState;
    CamParams* cam = interaction->camParams;

    glm::vec3 zoomDirection = glm::normalize(cam->target - cam->pos);

    /* Update camera position */
    if (yOffset > 0) {
        cam->pos += zoomDirection * interaction->zoomSpeed;
    }
    else {
        cam->pos -= zoomDirection * interaction->zoomSpeed;
    }

    float distance = glm::length(cam->pos - cam->target);
    float maxCamDist = interaction->maxCamDist;
    float minCamDist = interaction->minCamDist;

    /* Clip camera position */
    if (distance > maxCamDist) {
        cam->pos = cam->target - zoomDirection * maxCamDist;
    }
    else if (distance < minCamDist) {
        cam->pos = cam->target - zoomDirection * minCamDist;
    }
}

void
updateRotationAnimation(
    InteractionState& interaction,
    float angularVelocity)
{
    if (!interaction.isRotating)
        return;
    glm::mat4& rotationMatrix = interaction.rotationMatrix;
    /* Rotation matrix around the camera y-axis */
    glm::mat4 rotation = glm::rotate(glm::mat4(1.0f),
                                     angularVelocity,
                                     glm::vec3(0.0f, 1.0f, 0.0f));
    /* Accumulate the rotation */
    rotationMatrix *= rotation;
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
            /* Do nothing */
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
    return true;
}
/**
 * This function processes a set of nodes represented as vectors of floats.
 * It normalizes each node's components to a [-1, 1] range.
 */
std::vector<std::vector<float>>
normNodes(
    std::vector<std::vector<float>>& nodes
)
{
    if (nodes.empty())
        return {};

    size_t dim = nodes.front().size();

    std::vector<float> centroid(dim, 0.0f);
    for (size_t i = 0; i < nodes.size(); ++i) {
        for (size_t d = 0; d < dim; ++d) {
            centroid[d] += (nodes[i][d] - centroid[d]) / (i + 1);
        }
    }

    std::vector<std::vector<float>> procNodes;
    procNodes.reserve(nodes.size() * dim);
    float maxDist = 0.0;
    for (const auto& x_nd : nodes) {
        float distSqr = 0.0;
        std::vector<float> x_nd_c;
        x_nd_c.reserve(x_nd.size());
        for (size_t d = 0; d < dim; ++d) {
            x_nd_c.push_back(x_nd[d] - centroid[d]);
            distSqr += x_nd_c[d] * x_nd_c[d];
        }
        procNodes.push_back(x_nd_c); // Centered point
        maxDist = std::max(maxDist, std::sqrt(distSqr));
    }

    if (maxDist > 0.0f) {
        for (auto& x_nd : procNodes) {
            for (size_t d = 0; d < dim; ++d) {
                x_nd[d] /= maxDist;
            }
        }
    }
    return procNodes;
}

std::vector<float>
extractEdgeNodes(
    const std::vector<int>& path,
    const std::vector<std::vector<float>>& centers
)
{
    if (path.empty())
        return {};

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
    PrismParams& prismParams,
    const float jointRadius,
    std::vector<float>& edge_vertices,
    std::vector<float>& edge_normals,
    std::vector<unsigned int>& edge_indices
)
{
    std::vector<float> edges = extractEdgeNodes(path, node_centers);

    size_t dim = 0;
    if (!node_centers.empty()) {
        dim = node_centers.front().size();
    }

    generatePrisms(edges, dim, prismParams.edgeWidth, prismParams.edgeWidth,
                   prismParams.numEdgeSamples, jointRadius,
                   edge_vertices, edge_normals, edge_indices);
}

gl_Data
initGLData(
    const std::vector<GLfloat>& vertices,
    const std::vector<GLfloat>& normals,
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

    /* Vertex Buffer Object (VBO) */
    glGenBuffers(1, &gl_data.VBO);
    glBindBuffer(GL_ARRAY_BUFFER, gl_data.VBO);
    glBufferData(GL_ARRAY_BUFFER,
                 vertices.size() * sizeof(GLfloat),
                 vertices.data(),
                 GL_STATIC_DRAW);

    glVertexAttribPointer(0,
                          compSize,
                          GL_FLOAT,
                          GL_FALSE,
                          compSize * sizeof(GL_FLOAT),
                          (GLvoid*)0);
    glEnableVertexAttribArray(0);

    /* Normal Buffer Object (NBO) */
    glGenBuffers(1, &gl_data.NBO);
    glBindBuffer(GL_ARRAY_BUFFER, gl_data.NBO);
    glBufferData(GL_ARRAY_BUFFER,
                 normals.size() * sizeof(GLfloat),
                 normals.data(),
                 GL_STATIC_DRAW);

    glVertexAttribPointer(1,
            compSize,
            GL_FLOAT,
            GL_FALSE,
            compSize * sizeof(GL_FLOAT),
            (GLvoid*)0);
    glEnableVertexAttribArray(1);

    /* Element Buffer Object (EBO) */
    glGenBuffers(1, &gl_data.EBO);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, gl_data.EBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER,
                 indices.size() * sizeof(GLuint),
                 indices.data(),
                 GL_STATIC_DRAW);

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
    glDeleteBuffers(1, &gl_data.NBO);
    glDeleteBuffers(1, &gl_data.EBO);
}

void
getCameraMatrices(
    RenderState& state,
    CamParams& cam,
    InteractionState& interaction,
    glm::mat4& view,
    glm::mat4& projection
)
{
    /* World to camera transform */
    view = glm::lookAt(cam.pos, cam.target, cam.upAxis);

    glm::mat4 translationToOrigin = glm::translate(glm::mat4(1.0f), -cam.pos);
    glm::mat4 translationBack     = glm::translate(glm::mat4(1.0f), cam.pos);
    /* Apply the rotation relative to the world origin */
    view = translationBack * interaction.rotationMatrix * translationToOrigin * view;

    projection = glm::perspective(glm::radians(cam.fov),
                                  state.aspectRatio,
                                  cam.near,
                                  cam.far);
}

void
renderNodes(
    RenderState& state,
    gl_Data& gl_vs,
    CamParams& cam,
    InteractionState& interaction
)
{
    glm::mat4 model = glm::mat4(1.0f);

    glm::mat4 view, projection;
    getCameraMatrices(state, cam, interaction, view, projection);

    GLuint shProg = state.shaderProgs[NODES_SHADER];
    glUseProgram(shProg);

    GLint modelLoc = glGetUniformLocation(shProg, "model");
    glUniformMatrix4fv(modelLoc, 1, GL_FALSE, glm::value_ptr(model));
    GLuint viewLoc = glGetUniformLocation(shProg, "view");
    glUniformMatrix4fv(viewLoc, 1, GL_FALSE, glm::value_ptr(view));
    GLint projectionLoc = glGetUniformLocation(shProg, "projection");
    glUniformMatrix4fv(projectionLoc, 1, GL_FALSE, glm::value_ptr(projection));
    GLint camPosLoc = glGetUniformLocation(shProg, "camPos");
    glUniform3fv(camPosLoc, 1, glm::value_ptr(cam.pos));

    glBindVertexArray(gl_vs.VAO);
    glDrawElements(GL_TRIANGLES, gl_vs.iSize, GL_UNSIGNED_INT, 0);
}

void
renderEdges(
    RenderState& state,
    gl_Data& gl_es,
    CamParams& cam,
    InteractionState& interaction
)
{
    /* Object to world transform */
    glm::mat4 model = glm::mat4(1.0f);

    glm::mat4 view, projection;
    getCameraMatrices(state, cam, interaction, view, projection);

    GLuint shProg = state.shaderProgs[EDGES_SHADER];
    glUseProgram(shProg);

    GLint modelLoc = glGetUniformLocation(shProg, "model");
    glUniformMatrix4fv(modelLoc, 1, GL_FALSE, glm::value_ptr(model));
    GLuint viewLoc = glGetUniformLocation(shProg, "view");
    glUniformMatrix4fv(viewLoc, 1, GL_FALSE, glm::value_ptr(view));
    GLint projectionLoc = glGetUniformLocation(shProg, "projection");
    glUniformMatrix4fv(projectionLoc, 1, GL_FALSE, glm::value_ptr(projection));

    glBindVertexArray(gl_es.VAO);

    glDrawElements(GL_TRIANGLES, gl_es.iSize, GL_UNSIGNED_INT, 0);
}

bool
initShaders(
    RenderState& state
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
    state.shaderProgs[NODES_SHADER] = shaderProg;

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
    state.shaderProgs[EDGES_SHADER] = shaderProgLines;

    return true;
}

GLFWwindow*
initOpenGLWindow(
    WindowConfig& windowConfig
)
{
    if(!glfwInit()) {
        return nullptr;
    }

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    /* glfwWindowHint(GLFW_DEPTH_BITS, 24); */

    GLFWwindow* window = glfwCreateWindow(windowConfig.windowWidth,
                                          windowConfig.windowHeight,
                                          "TSL Path Map",
                                          NULL,
                                          NULL);
    if(!window) {
        glfwTerminate();
        return nullptr;
    }

    glfwMakeContextCurrent(window);

    glEnable(GL_CULL_FACE);
    glFrontFace(GL_CCW); // Front face are counter-clock wise formed
    /* glDisable(GL_CULL_FACE); */

    glEnable(GL_DEPTH_TEST);
    /* glDepthFunc(GL_LESS); */

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

    WindowConfig windowConfig = {
        .windowWidth     = 600,
        .windowHeight    = 600,
    };

    RenderState renderState = {
        .shaderProgs = {0},
        .aspectRatio = 1.0f,
    };

    /* Init camera */
    CamParams cam = {
        .fov  = 90.0f,
        .near = 0.01f,
        .far  = 10.0f,
        .pos    = glm::vec3(0.0f, 0.0f, -2.0f), // Camera position
        .target = glm::vec3(0.0f, 0.0f, 0.0f),  // Point camera looks at
        .upAxis = glm::vec3(0.0f, 1.0f, 0.0f)   // World up vector
    };

    InteractionState interaction = {
        .rotationMatrix = glm::mat4(1.0f),
        .lastPosOnSphere = glm::vec3(0.0f),
        .isDragging = false,
        .zoomSpeed = 0.1f,
        .camParams = &cam,
        .maxCamDist = 2.0f,
        .minCamDist = 0.1f,
        .isRotating = true,
    };

    AppContext appContext = {
        .windowConfig = &windowConfig,
        .renderState = &renderState,
        .interactionState = &interaction,
    };

    GLFWwindow* window = initOpenGLWindow(windowConfig);
    if (!window) {
        cleanOpenGLContext(window);
        return -1;
    }

    const GLubyte* version = glGetString(GL_VERSION);
    std::cout << "OpenGL Version: " << version << '\n';

    glfwSetWindowUserPointer(window, &appContext);

    if (!initShaders(renderState)) {
        cleanOpenGLContext(window);
        return -1;
    }

    /* Init nodes data */
    node_centers = normNodes(node_centers);
    std::vector<float> node_vertices, node_normals;
    std::vector<unsigned int> node_indices;
    SphereParams sphereParams = {
        .radius = 0.07f,
        .stacks = 20,
        .slices = 20,
    };
    generateSpheres(sphereParams.radius, sphereParams.stacks, sphereParams.slices,
                    node_centers, node_vertices, node_normals, node_indices);
    gl_Data gl_vs = initGLData(node_vertices, node_normals, node_indices, 3);

    /* Init edges data */
    std::vector<float> edge_vertices, edge_normals;
    std::vector<unsigned int> edge_indices;

    PrismParams prismParams = {
        .edgeWidth = sphereParams.radius * 0.7f,
        .numEdgeSamples = 20,
    };
    preprocessEdges(path, node_centers, prismParams, sphereParams.radius,
                    edge_vertices, edge_normals, edge_indices);
    gl_Data gl_es = initGLData(edge_vertices, edge_normals, edge_indices, 3);

    /* Callbacks */
    glfwSetFramebufferSizeCallback(window, framebufferSizeCallback);
    glfwSetKeyCallback(window, keyCallback);
    glfwSetMouseButtonCallback(window, mouseButtonCallback);
    glfwSetCursorPosCallback(window, cursorPositionCallback);
    glfwSetScrollCallback(window, scrollCallback);

    /* Call the resize window callback at least once */
    int win_w, win_h;
    glfwGetFramebufferSize(window, &win_w, &win_h);
    framebufferSizeCallback(window, win_w, win_h);

    /* Background color */
    glClearColor(0.2f, 0.3f, 0.3f, 1.0f);

    const int TARGET_FPS = 60; // Desired FPS
    const float FRAME_DURATION = 1.0f / TARGET_FPS;
    auto lastTime = std::chrono::high_resolution_clock::now();
    float elapsedTime = 0.0f;
    float angle = -1.0f;

    /* Rendering loop */
    while (!glfwWindowShouldClose(window)) {
        auto frameStartTime = std::chrono::high_resolution_clock::now();
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        updateRotationAnimation(interaction, angle * elapsedTime);
        renderEdges(renderState, gl_es, cam, interaction);
        renderNodes(renderState, gl_vs, cam, interaction);

        glfwSwapBuffers(window);
        glfwPollEvents();

        /* Fix FPS */
        auto frameEndTime = std::chrono::high_resolution_clock::now();
        float frameTime = std::chrono::duration<float>(frameEndTime - frameStartTime).count();
        if (frameTime < FRAME_DURATION) {
            std::this_thread::sleep_for(std::chrono::duration<float>(FRAME_DURATION - frameTime));
        }

        auto currentTime = std::chrono::high_resolution_clock::now();
        elapsedTime = std::chrono::duration<float>(currentTime - lastTime).count();
        lastTime = currentTime;
    }

    clearGLData(gl_vs);
    clearGLData(gl_es);

    cleanOpenGLContext(window);

    return 0;
}
