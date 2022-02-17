#version 330 core

out float FragColor;

in vec2 coords;

uniform float sigma_x_inv;
uniform sampler2D w;

void divergence()
{
    vec4 wL = texture(w, coords - vec2(1, 0));
    vec4 wR = texture(w, coords + vec2(1, 0));
    vec4 wT = texture(w, coords - vec2(0, 1));
    vec4 wB = texture(w, coords + vec2(0, 1));

    vec2 C = texture(w, coords).xy;
    if (wL.x < 0.0) 
        wL.x = -C.x;
    if (wR.x > 1.0) 
        wR.x = -C.x;
    if (wT.y > 1.0) 
        wT.x = -C.y;
    if (wB.y < 0.0) 
        wB.x = -C.y;
    
    div = 0.5 * sigma_x_inv * ((wR.x - wL.x) + (wT.y - wB.y));
    FragColor = vec4(div, 0, 0, 1.0);
}