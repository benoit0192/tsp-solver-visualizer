#ifndef GEO_H
#define GEO_H

#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <vector>

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

std::vector<GLfloat>
generatePrismVertices(
    const glm::vec3& point1,
    const glm::vec3& point2,
    float width,
    float height
);

#endif // !GEO_H
