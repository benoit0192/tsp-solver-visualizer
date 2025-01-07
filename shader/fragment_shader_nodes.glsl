#version 330 core

in vec3 fragPos;
in vec3 fragNormal;

out vec4 FragColor;

uniform vec3 lightDir = normalize(vec3(0.5, -0.5, 0.5));
uniform vec3 lightColor = vec3(1.0, 1.0, 1.0);
uniform vec3 objectColor = vec3(0.9, 0.9, 0.9);
uniform vec3 ambientColor = vec3(0.8, 0.8, 0.8);

void main()
{
    vec3 norm = normalize(fragNormal);

    float diff = max(dot(norm, -lightDir), 0.0);
    vec3 diffuse = diff * lightColor;

    vec3 ambient = ambientColor * lightColor;

    vec3 finalColor = (ambient + diffuse) * objectColor;

    FragColor = vec4(finalColor, 1.0);
}
