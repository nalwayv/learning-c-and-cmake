#version 440
 
in vec3 col;
out vec4 frag_col;

void main()
{
    frag_col = vec4(col, 1.0);
}
