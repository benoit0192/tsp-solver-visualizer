#include "geo.hpp"


glm::vec3
computeOrthogonal3D(
    const glm::vec3& dir
)
{
    glm::vec3 arbitrary = glm::abs(dir.x) > glm::abs(dir.y) ?
                          glm::vec3(0.0f, 1.0f, 0.0f) : glm::vec3(1.0f, 0.0f, 0.0f);
    return glm::normalize(glm::cross(dir, arbitrary));
}

void generatePrisms(
    const std::vector<float>& edges,
    const float width,
    const float height,
    const size_t dim,
    std::vector<float>& vertices,
    std::vector<unsigned int>& indices
)
{
    size_t offset = 2 * dim; // assuming 2 points per edge
    if (edges.size() % offset != 0) {
        throw std::runtime_error("invalid edge data: size mismatch.");
    }

    unsigned int vertexOffset = 0;
    for (size_t i = 0; i < edges.size(); i += offset) {
        glm::vec3 p1(0.0f), p2(0.0f);
        for (size_t d = 0; d < dim; ++d) {
            p1[d] = edges[i + d];
            p2[d] = edges[i + dim + d];
        }

        /* direction of the edge */
        glm::vec3 dir = p2 - p1;
        dir = glm::normalize(dir);

        /* compute orthogonal vectors for the face's width and height */
        glm::vec3 orth1 = computeOrthogonal3D(dir);
        glm::vec3 orth2 = glm::normalize(glm::cross(dir, orth1));

        /* scale orthogonal vectors for width and height */
        glm::vec3 offset1 = orth1 * width * 0.5f;
        glm::vec3 offset2 = orth2 * height * 0.5f;

        /* front face vertices (centered at p1) */
        glm::vec3 frontTopLeft     = p1 + offset1 + offset2;
        glm::vec3 frontTopRight    = p1 - offset1 + offset2;
        glm::vec3 frontBottomLeft  = p1 + offset1 - offset2;
        glm::vec3 frontBottomRight = p1 - offset1 - offset2;

        /* back face vertices (centered at p2) */
        glm::vec3 backTopLeft     = p2 + offset1 + offset2;
        glm::vec3 backTopRight    = p2 - offset1 + offset2;
        glm::vec3 backBottomLeft  = p2 + offset1 - offset2;
        glm::vec3 backBottomRight = p2 - offset1 - offset2;

        /* add vertices to the array */
        vertices.insert(vertices.end(), {
            /* front face */
            frontTopLeft.x, frontTopLeft.y, frontTopLeft.z,
            frontTopRight.x, frontTopRight.y, frontTopRight.z,
            frontBottomLeft.x, frontBottomLeft.y, frontBottomLeft.z,
            frontBottomRight.x, frontBottomRight.y, frontBottomRight.z,

            /* back face */
            backTopLeft.x, backTopLeft.y, backTopLeft.z,
            backTopRight.x, backTopRight.y, backTopRight.z,
            backBottomLeft.x, backBottomLeft.y, backBottomLeft.z,
            backBottomRight.x, backBottomRight.y, backBottomRight.z,
        });


        /* indices for the front face */
        indices.insert(indices.end(), {
            vertexOffset + 0, vertexOffset + 2, vertexOffset + 1,
            vertexOffset + 1, vertexOffset + 2, vertexOffset + 3
        });

        /* indices for the back face */
        indices.insert(indices.end(), {
            vertexOffset + 4, vertexOffset + 5, vertexOffset + 6,
            vertexOffset + 5, vertexOffset + 7, vertexOffset + 6
        });

        /* side faces (connecting front and back vertices) */
        indices.insert(indices.end(), {
            vertexOffset + 0, vertexOffset + 4, vertexOffset + 1,
            vertexOffset + 1, vertexOffset + 4, vertexOffset + 5
        });

        indices.insert(indices.end(), {
            vertexOffset + 1, vertexOffset + 5, vertexOffset + 7,
            vertexOffset + 1, vertexOffset + 3, vertexOffset + 7
        });

        indices.insert(indices.end(), {
            vertexOffset + 2, vertexOffset + 6, vertexOffset + 3,
            vertexOffset + 3, vertexOffset + 6, vertexOffset + 7
        });

        indices.insert(indices.end(), {
            vertexOffset + 0, vertexOffset + 2, vertexOffset + 6,
            vertexOffset + 6, vertexOffset + 4, vertexOffset + 0
        });

        vertexOffset += 8;  // 8 vertices per prism (4 for front, 4 for back)
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
                float z = radius * sin(phi) * sin(theta) + ((dim==3) ? center[2]
                                                                     : 0.0);
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

/* Map a 2D vector captured from the mouse position onto half of a unit sphere */
glm::vec3
mapToSphere(
    float x,
    float y,
    int width,
    int height
)
{
    /* Normalize coordinates (-1,1) */
    float nx = (2.0f * x - width) / width;
    /* Flip y-axis (-1: bottom, 1: top) */
    float ny = (height - 2.0f * y) / height;

    float length = nx * nx + ny * ny;

    /* Map to the sphere */
    if (length > 1.0f) {
        float norm = 1.0f / std::sqrt(length);
        return glm::vec3(nx * norm, ny * norm, 0.0f);
    }
    else {
        float z = std::sqrt(1.0f - length);
        return glm::vec3(nx, ny, z);
    }
}
