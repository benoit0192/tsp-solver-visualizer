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
    std::vector<float>& normals,
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

        std::vector<glm::vec3> prismVertices = {
            /* Front face vertices */
            frontTopLeft, frontTopRight, frontBottomLeft, frontBottomRight,     // 0,1,2,3
            /* Back face vertices */
            backTopLeft, backTopRight, backBottomLeft, backBottomRight,         // 4,5,6,7
            /* Left face vertices */
            backTopLeft, frontTopLeft, backBottomLeft, frontBottomLeft,         // 8,9,10,11
            /* Right face vertices */
            frontTopRight, backTopRight, frontBottomRight, backBottomRight,     // 12,13,14,15
            /* Top face vertices */
            backTopLeft, backTopRight, frontTopLeft, frontTopRight,             // 16,17,18,19
            /* Bottom face vertices */
            frontBottomLeft, frontBottomRight, backBottomLeft, backBottomRight, // 20,21,22,23
        };

        /* Add vertices to the array */
        for (const auto& v : prismVertices) {
            vertices.insert(vertices.end(), { v.x, v.y, v.z });
        }

        /* Face normals */
        glm::vec3 frontNormal = glm::normalize(glm::cross(frontTopRight - frontTopLeft,
                                                          frontBottomLeft - frontTopLeft));
        glm::vec3 backNormal = glm::normalize(glm::cross(backBottomLeft - backTopLeft,
                                                         backTopRight - backTopLeft));
        glm::vec3 leftNormal = glm::normalize(glm::cross(frontBottomLeft - frontTopLeft,
                                                         backTopLeft - frontTopLeft));
        glm::vec3 rightNormal = glm::normalize(glm::cross(backTopRight - frontTopRight,
                                                          frontBottomRight - frontTopRight));
        glm::vec3 topNormal = glm::normalize(glm::cross(frontTopLeft - backTopLeft,
                                                        backTopRight - backTopLeft));
        glm::vec3 bottomNormal = glm::normalize(glm::cross(backBottomLeft - backTopLeft,
                                                           frontBottomLeft - backTopLeft));

        std::vector<glm::vec3> faceNormals = {
            frontNormal, backNormal,
            leftNormal, rightNormal,
            topNormal, bottomNormal
        };

        for (size_t i = 0; i < 6; ++i) { // Loop over the 6 faces
            const glm::vec3& normal = faceNormals[i];
            for (size_t j = 0; j < 4; ++j) { // Each face has 4 vertices
                normals.insert(normals.end(), { normal.x, normal.y, normal.z });
            }
        }

        indices.insert(indices.end(), {
            /* indices for the front face */
            vertexOffset + 0, vertexOffset + 2, vertexOffset + 1,
            vertexOffset + 1, vertexOffset + 2, vertexOffset + 3,

            /* indices for the back face */
            vertexOffset + 4, vertexOffset + 5, vertexOffset + 6,
            vertexOffset + 5, vertexOffset + 7, vertexOffset + 6,

            /* indices for the left face */
            vertexOffset + 8, vertexOffset + 10, vertexOffset + 9,
            vertexOffset + 9, vertexOffset + 10, vertexOffset + 11,

            /* indices for the right face */
            vertexOffset + 12, vertexOffset + 14, vertexOffset + 13,
            vertexOffset + 13, vertexOffset + 14, vertexOffset + 15,

            /* indices for the top face */
            vertexOffset + 16, vertexOffset + 18, vertexOffset + 17,
            vertexOffset + 17, vertexOffset + 18, vertexOffset + 19,

            /* indices for the bottom face */
            vertexOffset + 20, vertexOffset + 22, vertexOffset + 21,
            vertexOffset + 21, vertexOffset + 22, vertexOffset + 23
        });

        vertexOffset += 24;  // 24 vertices per prism (4 for each face)
    }
}

void
generateSpheres(
    float radius,
    int stacks,
    int slices,
    const std::vector<std::vector<float>>& centers,
    std::vector<float>& vertices,
    std::vector<float>& normals,
    std::vector<unsigned int>& indices
)
{
    size_t dim = 0;
    if (!centers.empty()) {
        dim = centers.front().size();
    }
    vertices.reserve(centers.size() * (stacks + 1) * (slices + 1) * 3);
    normals.reserve(centers.size() * (stacks + 1) * (slices + 1) * 3);
    indices.reserve(centers.size() * stacks * slices * 2 * 3);
    unsigned int vertexOffset = 0;
    for (const auto& center : centers) {
        for (int i = 0; i <= stacks; ++i) {
            float phi = M_PI * i / stacks;
            for (int j = 0; j <= slices; ++j) {
                float theta = 2.0f * M_PI * j / slices;

                float nx = sin(phi) * cos(theta);
                float ny = cos(phi);
                float nz = sin(phi) * sin(theta);

                float x = radius * nx + center[0];
                float y = radius * ny + center[1];
                float z = radius * nz + ((dim==3) ? center[2] : 0.0);

                vertices.push_back(x);
                vertices.push_back(y);
                vertices.push_back(z);

                normals.push_back(nx);
                normals.push_back(ny);
                normals.push_back(nz);
            }
        }

        for (int i = 0; i < stacks; ++i) {
            for (int j = 0; j < slices; ++j) {
                unsigned int first = vertexOffset + i * (slices + 1) + j;
                unsigned int second = first + slices + 1;

                // Two triangles for the current quad
                indices.push_back(first);
                indices.push_back(first + 1);
                indices.push_back(second);

                indices.push_back(second);
                indices.push_back(first + 1);
                indices.push_back(second + 1);
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
