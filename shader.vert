#version 440

layout(location = 0)in vec3 v_pos;
layout(location = 1)in vec3 v_col;

out vec3 col;

void main()
{
    col = v_col;
    gl_Position = vec4(v_pos, 1.0);
}
