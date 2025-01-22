#version 330 core

in vec3 fragPos;
in vec3 fragNormal;

out vec4 FragColor;

uniform vec3 objectColor  = vec3(0.85, 0.85, 0.85);
uniform vec3 ambientColor = vec3(0.8, 0.8, 0.8);
uniform vec3 lightColor   = vec3(1.0, 1.0, 1.0);
uniform vec3 lightDir     = normalize(vec3(0.5, -0.5, 0.5));

uniform vec3 camPos;
uniform float shininess = 32.0;

void main()
{
    vec3 norm = normalize(fragNormal);

    // Diffuse component
    float diff   = max(dot(norm, -lightDir), 0.0);
    vec3 diffuse = diff * lightColor;

    // Ambient component
    vec3 ambient = ambientColor * lightColor;

    // Specular component
    vec3 viewDir    = normalize(camPos - fragPos);
    vec3 halfwayDir = normalize(viewDir - lightDir);
    float spec      = pow(max(dot(norm, halfwayDir), 0.0), shininess);
    vec3 specular   = spec * lightColor;

    vec3 finalColor = clamp((1.2*ambient + 0.2*diffuse + 0.1*specular) * objectColor,
                            0.0, 1.0);

    FragColor = vec4(finalColor, 1.0);
}
