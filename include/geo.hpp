#ifndef GEO_H
#define GEO_H

#include <vector>
#include <iostream>
#include <glm/glm.hpp>
#include <GLFW/glfw3.h>

glm::vec2
computeOrthogonal2D(
    const glm::vec2& dir
);

glm::vec3
computeOrthogonal3D(
    const glm::vec3& dir
);

std::vector<GLfloat>
generateQuadVertices(
    const glm::vec2& point1,
    const glm::vec2& point2,
    float lineWidth
);

void
generatePrisms(
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
);

void
generateSpheres(
    float radius,
    int stacks,
    int slices,
    const std::vector<std::vector<float>>& centers,
    std::vector<float>& vertices,
    std::vector<float>& normals,
    std::vector<unsigned int>& indices
);

glm::vec3
mapToSphere(
    float x,
    float y,
    int width,
    int height
);

#endif // !GEO_H
