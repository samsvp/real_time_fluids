#version 330 core

out vec4 FragColor;

in vec2 coords;

uniform float sigma_x_inv;
uniform sampler2D p;
uniform sampler2D w;

void main()
{
    float pL = texture(p, coords - vec2(1, 0)).x;
    float pR = texture(p, coords + vec2(1, 0)).x;
    float pT = texture(p, coords - vec2(0, 1)).x;
    float pB = texture(p, coords + vec2(0, 1)).x;
    
    FragColor = texture(w, coords);
    FragColor.xy -= 0.5 * sigma_x_inv * vec2(pR - pL, pT - pB);
}
