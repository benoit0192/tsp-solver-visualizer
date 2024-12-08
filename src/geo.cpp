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

std::vector<GLfloat>
generatePrismVertices(
    const glm::vec3& point1,
    const glm::vec3& point2,
    float width,
    float height
)
{
    glm::vec3 dir = point2 - point1;
    dir = glm::normalize(dir);

    glm::vec3 orth1 = computeOrthogonal3D(dir);
    glm::vec3 orth2 = glm::normalize(glm::cross(dir, orth1));

    // Scale orthogonal vectors to form the width and height of the prism
    glm::vec3 offset1 = orth1 * width  * 0.5f;
    glm::vec3 offset2 = orth2 * height * 0.5f;

    // Front face (centered at point1)
    glm::vec3 frontTopLeft     = point1 + offset1 + offset2;
    glm::vec3 frontTopRight    = point1 - offset1 + offset2;
    glm::vec3 frontBottomLeft  = point1 + offset1 - offset2;
    glm::vec3 frontBottomRight = point1 - offset1 - offset2;

    // Back face (centered at point2)
    glm::vec3 backTopLeft     = point2 + offset1 + offset2;
    glm::vec3 backTopRight    = point2 - offset1 + offset2;
    glm::vec3 backBottomLeft  = point2 + offset1 - offset2;
    glm::vec3 backBottomRight = point2 - offset1 - offset2;

    std::vector<GLfloat> prismVertices = {
        // Front face
        frontTopLeft.x   , frontTopLeft.y   , frontTopLeft.z,
        frontBottomLeft.x, frontBottomLeft.y, frontBottomLeft.z,
        frontTopRight.x  , frontTopRight.y  , frontTopRight.z,

        frontTopRight.x   , frontTopRight.y   , frontTopRight.z,
        frontBottomLeft.x , frontBottomLeft.y , frontBottomLeft.z,
        frontBottomRight.x, frontBottomRight.y, frontBottomRight.z,

        // Back face
        backTopLeft.x   , backTopLeft.y   , backTopLeft.z,
        backTopRight.x  , backTopRight.y  , backTopRight.z,
        backBottomLeft.x, backBottomLeft.y, backBottomLeft.z,

        backTopRight.x   , backTopRight.y   , backTopRight.z,
        backBottomRight.x, backBottomRight.y, backBottomRight.z,
        backBottomLeft.x , backBottomLeft.y , backBottomLeft.z,

        // Left face
        frontTopLeft.x  , frontTopLeft.y  , frontTopLeft.z,
        backBottomLeft.x, backBottomLeft.y, backBottomLeft.z,
        backTopLeft.x   , backTopLeft.y   , backTopLeft.z,

        frontTopLeft.x   , frontTopLeft.y   , frontTopLeft.z,
        frontBottomLeft.x, frontBottomLeft.y, frontBottomLeft.z,
        backBottomLeft.x , backBottomLeft.y , backBottomLeft.z,

        // Right face
        frontTopRight.x  , frontTopRight.y  , frontTopRight.z,
        backTopRight.x   , backTopRight.y   , backTopRight.z,
        backBottomRight.x, backBottomRight.y, backBottomRight.z,

        frontTopRight.x   , frontTopRight.y   , frontTopRight.z,
        backBottomRight.x , backBottomRight.y , backBottomRight.z,
        frontBottomRight.x, frontBottomRight.y, frontBottomRight.z,

        // Top face
        frontTopLeft.x, frontTopLeft.y, frontTopLeft.z,
        backTopRight.x, backTopRight.y, backTopRight.z,
        backTopLeft.x , backTopLeft.y , backTopLeft.z,

        frontTopLeft.x , frontTopLeft.y , frontTopLeft.z,
        frontTopRight.x, frontTopRight.y, frontTopRight.z,
        backTopRight.x , backTopRight.y , backTopRight.z,

        // Bottom face
        frontBottomLeft.x, frontBottomLeft.y, frontBottomLeft.z,
        backBottomLeft.x , backBottomLeft.y , backBottomLeft.z,
        backBottomRight.x, backBottomRight.y, backBottomRight.z,

        frontBottomLeft.x , frontBottomLeft.y , frontBottomLeft.z,
        backBottomRight.x , backBottomRight.y , backBottomRight.z,
        frontBottomRight.x, frontBottomRight.y, frontBottomRight.z
    };

    return prismVertices;
}
