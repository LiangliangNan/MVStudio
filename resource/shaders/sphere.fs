#version 120

//uniform vec3 lightDir; 

float Ns = 250.0;
vec4 mat_specular=vec4(1); 
vec4 light_specular=vec4(1); 

void main(void)
{
	vec3 lightDir = vec3(0.0, 0.0, 1.0);

    // calculate normal from texture coordinates
    vec3 N;
    N.xy = gl_PointCoord * 2.0 - vec2(1.0);    
    float mag = dot(N.xy, N.xy);
    if (mag > 1.0) discard;   // kill pixels outside circle
    N.z = sqrt(1.0-mag);

    // calculate lighting
    float diffuse = max(0.0, dot(lightDir, N));
 
    vec3 eye = vec3(0.0, 0.0, 1.0);
    vec3 halfVector = normalize(eye + lightDir);
    float spec = max(pow(dot(N, halfVector), Ns), 0.0); 
    vec4 S = light_specular * mat_specular * spec;
    
	gl_FragColor = gl_Color * diffuse + S;
}
