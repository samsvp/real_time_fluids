#version 330 core

out vec4 FragColor;

in vec2 coords;

uniform float dt;
uniform float sigma_x_inv;
uniform sampler2D u;
uniform sampler2D x;

// performs bilinear interpolation
vec4 bilinear_interp(sampler2D sam, vec2 coords, vec2 tsize) 
{
    vec2 st = coords / tsize - 0.5;
    vec2 iuv = floor(st);
    vec2 fuv = fract(st);
    vec4 a = texture(sam, (iuv + vec2(0.5, 0.5)) * tsize);
    vec4 b = texture(sam, (iuv + vec2(1.5, 0.5)) * tsize);
    vec4 c = texture(sam, (iuv + vec2(0.5, 1.5)) * tsize);
    vec4 d = texture(sam, (iuv + vec2(1.5, 1.5)) * tsize);
    return mix(mix(a, b, fuv.x), mix(c, d, fuv.x), fuv.y);
}


void main()
{
    coords = aPos.xy;
    // change texture2D to texture
    vec2 pos = coords - dt * sigma_x_inv * texture(u, coords).xy;
    FragColor = texture(x, pos);
}