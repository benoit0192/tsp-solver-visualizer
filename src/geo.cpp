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
    const size_t dim,
    const float width,
    const float height,
    const int numEdgeSamples,
    const float jointRadius,
    std::vector<float>& vertices,
    std::vector<float>& normals,
    std::vector<unsigned int>& indices,
    std::vector<float>& weights
)
{
    size_t offset = 2 * dim; // An edge is defined by its 2 far end points
    if (edges.size() % offset != 0) {
        throw std::runtime_error("invalid edge data: size mismatch.");
    }

    /* The function finds the vertex of the front/back face of a projected
     * prism onto the surface of a connected sphere. The projection is
     * performed for points contained within the sphere. The point is
     * projected along the prism's main axis direction at a distance 'dist'.
     * We solve the quadratic equation: |point + dist*dir - center|**2 = R**2 */
    auto adjustToSphere = [](const glm::vec3& point,
                             const glm::vec3& center,
                             const float R,
                             const glm::vec3 dir)
    {
        float distToSphere = glm::length(point - center);
        if (distToSphere < R) {
            glm::vec3 m = point - center;
            glm::vec3 dir_unit = glm::normalize(dir);
            float b = glm::dot(m,dir_unit);
            float c = glm::dot(m,m) - R * R;
            float discriminant = b * b - c;
            if (discriminant < 0.0f)
                return point;
            float dist = -b + std::sqrt(discriminant);
            glm::vec3 p_sph = point + dir_unit * dist;
            return p_sph;
        }
        return point;
    };

    auto getVertex = [&](unsigned int index) -> glm::vec3 {
                            return glm::vec3(vertices[index * 3 + 0],
                                             vertices[index * 3 + 1],
                                             vertices[index * 3 + 2]);
    };

    float totalPathLen = 0.0f;
    for (size_t i = 0; i < edges.size(); i += offset) {
        glm::vec3 p1(0.0f), p2(0.0f);
        for (size_t d = 0; d < dim; ++d) {
            p1[d] = edges[i + d];
            p2[d] = edges[i + dim + d];
        }
        float edgeLen = glm::length(p1 - p2);
        /* If the edge's end joints are overlapping, do not create an edge */
        if (edgeLen <= 2 * jointRadius)
            continue;
        totalPathLen += edgeLen;
    }

    unsigned int vertexOffset = 0;
    float accumulatedLen = 0.0f;
    for (size_t i = 0; i < edges.size(); i += offset) {
        glm::vec3 p1(0.0f), p2(0.0f);
        for (size_t d = 0; d < dim; ++d) {
            p1[d] = edges[i + d];
            p2[d] = edges[i + dim + d];
        }
        float edgeLen = glm::length(p1 - p2);
        /* If the edge's end joints are overlapping, do not create an edge */
        if (edgeLen <= 2 * jointRadius)
            continue;

        float startGrad = accumulatedLen / totalPathLen;
        float endGrad   = (accumulatedLen + edgeLen) / totalPathLen;
        accumulatedLen += edgeLen;

        /* direction of the edge */
        glm::vec3 dir = glm::normalize(p2 - p1);

        /* compute orthogonal vectors for the face's width and height */
        glm::vec3 orth1 = computeOrthogonal3D(dir);
        glm::vec3 orth2 = glm::normalize(glm::cross(dir, orth1));

        /* scale orthogonal vectors for width and height */
        glm::vec3 offset1 = orth1 * width * 0.5f;
        glm::vec3 offset2 = orth2 * height * 0.5f;

        std::vector<glm::vec3> frontEdgePoints, backEdgePoints;
        std::vector<float> frontGradWeights, backGradWeights;
        /* Loop over each front/back face edges (4 edges per face) */
        for (int i = 0; i <= numEdgeSamples; ++i) {
            /* Ensures t varies from -1.0 to 1.0 */
            float t = -1.0f + 2.0f * (i * 1.0f / numEdgeSamples);

            /* Compute edge points for the front face */
            glm::vec3 frontTop    = p1 + t * offset1 + offset2;
            glm::vec3 frontBottom = p1 + t * offset1 - offset2;
            glm::vec3 frontLeft   = p1 - offset1 + t * offset2;
            glm::vec3 frontRight  = p1 + offset1 + t * offset2;

            /* Adjust points to fit the sphere */
            frontEdgePoints.push_back(adjustToSphere(frontTop, p1, jointRadius, dir));
            frontEdgePoints.push_back(adjustToSphere(frontBottom, p1, jointRadius, dir));
            frontEdgePoints.push_back(adjustToSphere(frontLeft, p1, jointRadius, dir));
            frontEdgePoints.push_back(adjustToSphere(frontRight, p1, jointRadius, dir));
            frontGradWeights.insert(frontGradWeights.end(),
                                    {startGrad, startGrad, startGrad, startGrad});
            /* Compute edge points for the back face */
            glm::vec3 backTop    = p2 + t * offset1 + offset2;
            glm::vec3 backBottom = p2 + t * offset1 - offset2;
            glm::vec3 backLeft   = p2 - offset1 + t * offset2;
            glm::vec3 backRight  = p2 + offset1 + t * offset2;

            /* Adjust points to fit the sphere */
            backEdgePoints.push_back(adjustToSphere(backTop, p2, jointRadius, -dir));
            backEdgePoints.push_back(adjustToSphere(backBottom, p2, jointRadius, -dir));
            backEdgePoints.push_back(adjustToSphere(backLeft, p2, jointRadius, -dir));
            backEdgePoints.push_back(adjustToSphere(backRight, p2, jointRadius, -dir));
            backGradWeights.insert(backGradWeights.end(),
                                   {endGrad, endGrad, endGrad, endGrad});
        }

        for (const auto& vert : frontEdgePoints) {
            vertices.insert(vertices.end(), {vert.x, vert.y, vert.z});
        }
        for (const auto& vert : backEdgePoints) {
            vertices.insert(vertices.end(), {vert.x, vert.y, vert.z});
        }
        weights.insert(weights.end(),
                       frontGradWeights.begin(), frontGradWeights.end());
        weights.insert(weights.end(),
                       backGradWeights.begin(), backGradWeights.end());
        /* 4 points per sample per face (top_i, bottom_i, left_i, right_i) */
        size_t verticesPerSample = 4;

        normals.resize(vertices.size());

        for (size_t j = 0; j < (size_t)numEdgeSamples; ++j) {

            /* Offset indices for the front and back faces */
            unsigned int frontOffset = vertexOffset;
            unsigned int backOffset  = vertexOffset +
                                       (numEdgeSamples + 1) * verticesPerSample;

            for (size_t k = 0; k < verticesPerSample; ++k) {
                unsigned int v0 = frontOffset + j * verticesPerSample + k;
                unsigned int v1 = frontOffset + (j + 1) * verticesPerSample + k;
                unsigned int v2 = backOffset + j * verticesPerSample + k;
                unsigned int v3 = backOffset + (j + 1) * verticesPerSample + k;

                /* Apply mirroring for opposite faces to correct CCW winding */
                glm::vec3 n;
                if (k % 2 == 0)
                {
                    indices.push_back(v0);
                    indices.push_back(v2);
                    indices.push_back(v1);

                    indices.push_back(v1);
                    indices.push_back(v2);
                    indices.push_back(v3);

                    n = glm::normalize(glm::cross(getVertex(v2) - getVertex(v0),
                                                  getVertex(v1) - getVertex(v0)));
                }
                else {
                    indices.push_back(v0);
                    indices.push_back(v1);
                    indices.push_back(v2);

                    indices.push_back(v2);
                    indices.push_back(v1);
                    indices.push_back(v3);

                    n = glm::normalize(glm::cross(getVertex(v1) - getVertex(v0),
                                                  getVertex(v2) - getVertex(v0)));
                }
                for (unsigned int v : {v0, v1, v2, v3}) {
                    normals[v * 3 + 0] = n.x;
                    normals[v * 3 + 1] = n.y;
                    normals[v * 3 + 2] = n.z;
                }
            }
        }
        vertexOffset += 2 * (numEdgeSamples + 1) * verticesPerSample; // Front and back faces
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
                float z = radius * nz + ((dim==3) ? center[2] : 0.0f);

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

                /* Two triangles for the current quad */
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
