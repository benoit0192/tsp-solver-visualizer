#version 330 core

in vec3 vertexPos;
out vec4 FragColor;
uniform float fogNearDist = -1.0;
uniform float fogFarDist = 6.0;

void main()
{
    vec4 fogColor    = vec4(0.0, 0.0, 0.0, 1.0);
    vec4 objectColor = vec4(1.0, 1.0, 1.0, 1.0);
    float dist       = length(vertexPos);
    float fogFactor  = smoothstep(fogNearDist, fogFarDist, dist);
    vec4 color       = mix(objectColor, fogColor, fogFactor);
    FragColor        = color;
}
