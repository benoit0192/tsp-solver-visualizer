#include "geo.hpp"


glm::vec2
computeOrthogonal2D(
    const glm::vec2& dir
)
{
    return glm::vec2(-dir.y, dir.x);
}

glm::vec3
computeOrthogonal3D(
    const glm::vec3& dir
)
{
    glm::vec3 arbitrary = glm::abs(dir.x) > glm::abs(dir.y) ?
                          glm::vec3(0.0f, 1.0f, 0.0f) : glm::vec3(1.0f, 0.0f, 0.0f);
    return glm::normalize(glm::cross(dir, arbitrary));
}

std::vector<GLfloat>
generateQuadVertices(
    const glm::vec2& point1,
    const glm::vec2& point2,
    float lineWidth
)
{
    glm::vec2 dir = point2 - point1;
    dir = glm::normalize(dir);

    glm::vec2 perp = computeOrthogonal2D(dir);
    glm::vec2 offset = perp * lineWidth * 0.5f;

    std::vector<GLfloat> quadVertices = {
        // First triangle
        point1.x - offset.x, point1.y - offset.y,
        point2.x - offset.x, point2.y - offset.y,
        point1.x + offset.x, point1.y + offset.y,

        // Second triangle
        point1.x + offset.x, point1.y + offset.y,
        point2.x + offset.x, point2.y + offset.y,
        point2.x - offset.x, point2.y - offset.y
    };

    return quadVertices;
}

void
generatePrisms(
    const std::vector<float>& edges,
    const float width,
    const float height,
    const size_t dim,
    std::vector<float>& vertices,
    std::vector<unsigned int>& indices
)
{
    size_t offset = 2 * dim;
    if (edges.size() % offset != 0) {
        throw std::runtime_error("Invalid edge data: size mismatch.");
    }

    unsigned int vertexOffset = 0;
    for (size_t i = 0; i < edges.size(); i += offset) {
        glm::vec3 p1(0.0f), p2(0.0f);
        for (size_t d = 0; d < dim; ++d) {
            p1[d] = edges[i + d];
            p2[d] = edges[i + dim + d];
        }
        glm::vec3 dir = p2 - p1;

        dir = glm::normalize(dir);

        glm::vec3 orth1 = computeOrthogonal3D(dir);
        glm::vec3 orth2 = glm::normalize(glm::cross(dir, orth1));

        // Scale orthogonal vectors to form the width and height of the prism
        glm::vec3 offset1 = orth1 * width  * 0.5f;
        glm::vec3 offset2 = orth2 * height * 0.5f;

        // Front face (centered at point1)
        glm::vec3 frontTopLeft     = p1 + offset1 + offset2;
        glm::vec3 frontTopRight    = p1 - offset1 + offset2;
        glm::vec3 frontBottomLeft  = p1 + offset1 - offset2;
        glm::vec3 frontBottomRight = p1 - offset1 - offset2;

        // Back face (centered at point2)
        glm::vec3 backTopLeft     = p2 + offset1 + offset2;
        glm::vec3 backTopRight    = p2 - offset1 + offset2;
        glm::vec3 backBottomLeft  = p2 + offset1 - offset2;
        glm::vec3 backBottomRight = p2 - offset1 - offset2;

        vertices.insert(vertices.end(), {
            // Front face
            frontTopLeft.x   , frontTopLeft.y   , frontTopLeft.z,
            frontTopRight.x  , frontTopRight.y  , frontTopRight.z,
            frontBottomLeft.x, frontBottomLeft.y, frontBottomLeft.z,
            frontBottomRight.x, frontBottomRight.y, frontBottomRight.z,

            // Back face
            backTopLeft.x   , backTopLeft.y   , backTopLeft.z,
            backTopRight.x  , backTopRight.y  , backTopRight.z,
            backBottomLeft.x, backBottomLeft.y, backBottomLeft.z,
            backBottomRight.x, backBottomRight.y, backBottomRight.z,
        });

        // Front face (vertices 0-3)
        indices.insert(indices.end(), {vertexOffset + 0,
                                       vertexOffset + 2,
                                       vertexOffset + 1,
                                       vertexOffset + 1,
                                       vertexOffset + 2,
                                       vertexOffset + 3});
        // Back face (vertices 4-7)
        indices.insert(indices.end(), {vertexOffset + 4,
                                       vertexOffset + 5,
                                       vertexOffset + 6,
                                       vertexOffset + 5,
                                       vertexOffset + 7,
                                       vertexOffset + 6});

        // Side faces (connecting front and back vertices)
        for (int j = 0; j < 4; ++j) {
            unsigned int front     = vertexOffset + j;               // Current front vertex
            unsigned int back      = vertexOffset + j + 4;           // Corresponding back vertex
            unsigned int nextFront = vertexOffset + (j + 1) % 4;     // Next front vertex (wraps around)
            unsigned int nextBack  = vertexOffset + (j + 1) % 4 + 4; // Corresponding next back vertex

            indices.insert(indices.end(), {front, back, nextBack});
            indices.insert(indices.end(), {front, nextBack, nextFront});
        }
        vertexOffset += 8;
    }
}

void
generateSpheres(
    float radius,
    int stacks,
    int slices,
    const std::vector<std::vector<float>>& centers,
    std::vector<float>& vertices,
    std::vector<unsigned int>& indices
)
{
    size_t dim = 0;
    if (!centers.empty()) {
        dim = centers.front().size();
    }
    vertices.reserve(centers.size() * (stacks + 1) * (slices + 1) * 3);
    indices.reserve(centers.size() * stacks * slices * 2 * 3);
    unsigned int vertexOffset = 0;
    for (const auto& center : centers) {
        for (int i = 0; i <= stacks; ++i) {
            float phi = M_PI * i / stacks;
            for (int j = 0; j <= slices; ++j) {
                float theta = 2.0f * M_PI * j / slices;

                float x = radius * sin(phi) * cos(theta) + center[0];
                float y = radius * cos(phi) + center[1];
                float z = radius * sin(phi) * sin(theta) + (dim==3) ? center[2]
                                                                    : 0.0;
                vertices.push_back(x);
                vertices.push_back(y);
                vertices.push_back(z);
            }
        }

        for (int i = 0; i < stacks; ++i) {
            for (int j = 0; j < slices; ++j) {
                unsigned int first = vertexOffset + i * (slices + 1) + j;
                unsigned int second = first + slices + 1;

                // Two triangles for the current quad
                indices.push_back(first);
                indices.push_back(second);
                indices.push_back(first + 1);

                indices.push_back(second);
                indices.push_back(second + 1);
                indices.push_back(first + 1);
            }
        }
        vertexOffset += (stacks + 1) * (slices + 1);
    }
}
