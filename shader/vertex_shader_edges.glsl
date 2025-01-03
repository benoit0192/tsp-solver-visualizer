#version 330 core

layout(location = 0) in vec3 aPos;
out vec3 vertexPos;

uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;

void main()
{
    vec4 viewPos = view * model * vec4(aPos, 1.0);
    gl_Position = projection * viewPos;
    vertexPos = viewPos.xyz;
}
