#version 330 core

layout(location = 0) in vec2 aPos;

uniform float f_ScaleFactor = 1.0;

void main()
{
    gl_Position = vec4(aPos, 0.0, 1.0);
    gl_PointSize = 25.0 * f_ScaleFactor;
}
