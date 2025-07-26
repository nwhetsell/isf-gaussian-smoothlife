/*{
    "CATEGORIES": [
        "Generator"
    ],
    "CREDIT": "cornusammonis <https://www.shadertoy.com/user/cornusammonis>",
    "DESCRIPTION": "Variant of Stephan Rafler’s SmoothLife that uses a separable Gaussian kernel to compute inner and outer fullness, converted from <https://www.shadertoy.com/view/XtVXzV>",
    "INPUTS": [
        {
            "NAME": "restart",
            "LABEL": "Restart",
            "TYPE": "event"
        },
        {
            "NAME": "addCells",
            "TYPE": "bool",
            "DEFAULT": false
        },
        {
            "NAME": "mouse",
            "TYPE": "point2D",
            "DEFAULT": [0.5, 0.5],
            "MIN": [0, 0],
            "MAX": [1, 1]
        }
    ],
    "ISFVSN": "2",
    "PASSES": [
        {
            "TARGET": "bufferA",
            "PERSISTENT": true,
            "FLOAT": true
        },
        {
            "TARGET": "bufferB",
            "PERSISTENT": true,
            "FLOAT": true
        },
        {
            "TARGET": "bufferC",
            "PERSISTENT": true,
            "FLOAT": true
        },
        {

        }
    ]
}
*/

//
// ShaderToy Buffer A
//

// ---------------------------------------------
const float or = 18.0;         // outer gaussian std dev
const float ir = 6.0;         // inner gaussian std dev
const float b1 = 0.19;         // birth1
const float b2 = 0.212;        // birth2
const float s1 = 0.267;        // survival1
const float s2 = 0.445;        // survival2
const float dt = 0.2;          // timestep
const float alpha_n = 0.017;   // sigmoid width for outer fullness
const float alpha_m = 0.112;   // sigmoid width for inner fullness
// ---------------------------------------------

bool reset() {
    return restart;
}

// the logistic function is used as a smooth step function
float sigma1(float x,float a,float alpha)
{
    return 1.0 / ( 1.0 + exp( -(x-a)*4.0/alpha ) );
}

float sigma2(float x,float a,float b,float alpha)
{
    return sigma1(x,a,alpha)
        * ( 1.0-sigma1(x,b,alpha) );
}

float sigma_m(float x,float y,float m,float alpha)
{
    return x * ( 1.0-sigma1(m,0.5,alpha) )
        + y * sigma1(m,0.5,alpha);
}

// the transition function
// (n = outer fullness, m = inner fullness)
float s(float n,float m)
{
    return sigma2( n, sigma_m(b1,s1,m,alpha_m),
        sigma_m(b2,s2,m,alpha_m), alpha_n );
}

#define T(d) texture(bufferA, fract(uv+d)).x


//
// ShaderToy Buffer B
//

#define SQRT_2_PI 2.50662827463

// ---------------------------------------------
// const float or = 18.0;         // outer gaussian std dev
// const float ir = 6.0;          // inner gaussian std dev
const int   oc = 50;           // sample cutoff
// ---------------------------------------------


vec2 gaussian(float i, vec2 a, vec2 d) {
     return a * exp( -(i*i) / d );
}


#define fragColor gl_FragColor
#define fragCoord gl_FragCoord
#define iFrame FRAMEINDEX
#define iResolution RENDERSIZE


void main()
{
    if (PASSINDEX == 0) // ShaderToy Buffer A
    {
        vec2 tx = 1.0 / iResolution.xy;
        vec2 uv = fragCoord.xy * tx;

        const float _K0 = -20.0/6.0; // center weight
        const float _K1 = 4.0/6.0;   // edge-neighbors
        const float _K2 = 1.0/6.0;   // vertex-neighbors

        // We can optionally add a laplacian to the update rule
    	// to decrease the appearance of aliasing, but we also
    	// introduce subtle anisotropy by doing so.
        // vec4 t = vec4(tx, -tx.y, 0.0);
        // float u =    T( t.ww); float u_n =  T( t.wy); float u_e =  T( t.xw);
        // float u_s =  T( t.wz); float u_w =  T(-t.xw); float u_nw = T(-t.xz);
        // float u_sw = T(-t.xy); float u_ne = T( t.xy); float u_se = T( t.xz);
        // float lapl  = _K0*u + _K1*(u_n + u_e + u_w + u_s) + _K2*(u_nw + u_sw + u_ne + u_se);

        vec4 current = IMG_NORM_PIXEL(bufferA, uv);
        vec2 fullness = IMG_NORM_PIXEL(bufferC, uv).xy;

        float delta =  2.0 * s( fullness.x, fullness.y ) - 1.0;
        float new = clamp( current.x + dt * delta, 0.0, 1.0 );

        if(addCells) {
            // from chronos' SmoothLife shader https://www.shadertoy.com/view/XtdSDn
            float dst = length(fragCoord.xy - mouse.xy * RENDERSIZE);
            if(dst <= or) {
            	new = step((ir+1.5), dst) * (1.0 - step(or, dst));
            }
        }

        float initialCellCount = min(iResolution.x, iResolution.y) / 50.;
        if(iFrame < initialCellCount || reset()) {
            float angle = 8. * atan(1) * // 2 pi
                          iFrame / initialCellCount;
           	vec2 initialCoordinate = 0.25 * vec2(cos(angle), sin(angle)) + 0.5;
            float dst = length(fragCoord.xy - initialCoordinate * RENDERSIZE);
            if(dst <= or) {
            	new = step((ir+1.5), dst) * (1.0 - step(or, dst));
            }
        }
        fragColor = vec4(new, fullness, current.w);
    }
    else if (PASSINDEX == 1) // ShaderToy Buffer B
    {
        vec2 tx = 1.0 / iResolution.xy;
        vec2 uv = fragCoord.xy * tx;
        tx = (mod(float(iFrame),2.0) < 1.0) ? vec2(tx.x,0) : vec2(0,tx.y);

        vec2 sigma = vec2(or, ir);
        vec2 a = vec2(1.0 / (sigma * SQRT_2_PI));
        vec2 d = vec2(2.0 * sigma * sigma);
        vec2 acc = vec2(0.0);
        vec2 sum = vec2(0.0);

        // centermost term
        acc += a * IMG_NORM_PIXEL(bufferA, uv).x;
        sum += a;

        // sum up remaining terms symmetrically
        for (int i = 1; i <= oc; i++) {
            float fi = float(i);
            vec2 g = gaussian(fi, a, d);
            vec2 posL = fract(uv - tx * fi);
            vec2 posR = fract(uv + tx * fi);
            acc += g * (IMG_NORM_PIXEL(bufferA, posL).x + IMG_NORM_PIXEL(bufferA, posR).x);
            sum += 2.0 * g;
        }

        vec2 x_pass = acc / sum;

        fragColor = vec4(x_pass,0,0);
    }
    else if (PASSINDEX == 2) // ShaderToy Buffer C
    {
        vec2 tx = 1.0 / iResolution.xy;
        vec2 uv = fragCoord.xy * tx;
        tx = (mod(float(iFrame),2.0) < 1.0) ? vec2(0,tx.y) : vec2(tx.x,0);

        vec2 sigma = vec2(or, ir);
        vec2 a = vec2(1.0 / (sigma * SQRT_2_PI));
        vec2 d = vec2(2.0 * sigma * sigma);
        vec2 acc = vec2(0.0);
        vec2 sum = vec2(0.0);

        // centermost term
        acc += a * IMG_NORM_PIXEL(bufferB, uv).x;
        sum += a;

        // sum up remaining terms symmetrically
        for (int i = 1; i <= oc; i++) {
            float fi = float(i);
            vec2 g = gaussian(fi, a, d);
            vec2 posL = fract(uv - tx * fi);
            vec2 posR = fract(uv + tx * fi);
            acc += g * (IMG_NORM_PIXEL(bufferB, posL).xy + IMG_NORM_PIXEL(bufferB, posR).xy);
            sum += 2.0 * g;
        }

        vec2 y_pass = acc / sum;

        fragColor = vec4(y_pass,0,0);
    }
    else if (PASSINDEX == 3) // ShaderToy Image
    {
        vec2 uv = fragCoord.xy / iResolution.xy;
    	vec4 col = IMG_NORM_PIXEL(bufferA, uv);

        fragColor = col.x*vec4(1.0) + col.y*vec4(1,0.5,0,0) + col.z*vec4(0,0.5,1,0);
    }
}
