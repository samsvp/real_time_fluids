#version 330 core

out vec4 FragColor;

in vec2 coords;

uniform vec2 offset;
uniform float scale;
uniform sampler2D p;
uniform sampler2D x;

void boundary()
{
    FragColor = scale * texture(x, coords + offset);
}