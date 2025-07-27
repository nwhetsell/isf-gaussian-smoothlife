/*{
    "CATEGORIES": [
        "Generator"
    ],
    "CREDIT": "cornusammonis <https://www.shadertoy.com/user/cornusammonis>",
    "DESCRIPTION": "Variant of Stephan Rafler’s SmoothLife that uses a separable Gaussian kernel to compute inner and outer fullness, converted from <https://www.shadertoy.com/view/XtVXzV>",
    "INPUTS": [
        {
            "NAME": "addCells",
            "LABEL": "Add cells",
            "TYPE": "event"
        },
        {
            "NAME": "addCellsWithMouse",
            "LABEL": "Add cells with mouse",
            "TYPE": "bool",
            "DEFAULT": false
        },
        {
            "NAME": "mouse",
            "TYPE": "point2D",
            "DEFAULT": [0.5, 0.5],
            "MIN": [0, 0],
            "MAX": [1, 1]
        },
        {
            "NAME": "or",
            "LABEL": "Outer Gaussian std. dev.",
            "TYPE": "float",
            "DEFAULT": 18,
            "MAX": 100,
            "MIN": 0
        },
        {
            "NAME": "ir",
            "LABEL": "Inner Gaussian std. dev.",
            "TYPE": "float",
            "DEFAULT": 6,
            "MAX": 100,
            "MIN": 0
        },
        {
            "NAME": "b1",
            "LABEL": "Birth 1",
            "TYPE": "float",
            "DEFAULT": 0.19,
            "MAX": 1,
            "MIN": 0
        },
        {
            "NAME": "b2",
            "LABEL": "Birth 2",
            "TYPE": "float",
            "DEFAULT": 0.212,
            "MAX": 1,
            "MIN": 0
        },
        {
            "NAME": "s1",
            "LABEL": "Survival 1",
            "TYPE": "float",
            "DEFAULT": 0.267,
            "MAX": 1,
            "MIN": 0
        },
        {
            "NAME": "s2",
            "LABEL": "Survival 2",
            "TYPE": "float",
            "DEFAULT": 0.445,
            "MAX": 1,
            "MIN": 0
        },
        {
            "NAME": "dt",
            "LABEL": "Time step",
            "TYPE": "float",
            "DEFAULT": 0.2,
            "MAX": 10,
            "MIN": 0
        },
        {
            "NAME": "alpha_n",
            "LABEL": "Sigmoid width for outer fullness",
            "TYPE": "float",
            "DEFAULT": 0.017,
            "MAX": 1,
            "MIN": 0
        },
        {
            "NAME": "alpha_m",
            "LABEL": "Sigmoid width for inner fullness",
            "TYPE": "float",
            "DEFAULT": 0.112,
            "MAX": 10,
            "MIN": 0
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

// the logistic function is used as a smooth step function
float sigma1(float x, float a, float alpha)
{
    return 1. / (1. + exp(-(x-a) * 4. / alpha));
}

float sigma2(float x, float a, float b, float alpha)
{
    return sigma1(x, a, alpha) * (1. - sigma1(x, b, alpha));
}

float sigma_m(float x, float y, float m, float alpha)
{
    return x * (1. - sigma1(m, 0.5, alpha)) + y * sigma1(m, 0.5, alpha);
}

// the transition function
// (n = outer fullness, m = inner fullness)
float s(float n, float m)
{
    return sigma2(
        n,
        sigma_m(b1, s1, m, alpha_m),
        sigma_m(b2, s2, m, alpha_m),
        alpha_n
    );
}


//
// ShaderToy Buffer B
//

// ---------------------------------------------
const int   oc = 50;           // sample cutoff
// ---------------------------------------------

struct GaussianSummation {
    vec2 a;
    vec2 d;
    vec2 acc;
    vec2 sum;
};

GaussianSummation GaussianSummation_create(float or, float ir)
{
    vec2 sigma = vec2(or, ir);

    GaussianSummation gaussianSummation;
    gaussianSummation.a = vec2(1. / (sigma * 2.50662827463)); // sqrt(2 * pi)
    gaussianSummation.d = vec2(2. * sigma * sigma);
    gaussianSummation.acc = vec2(0.);
    gaussianSummation.sum = vec2(0.);
    return gaussianSummation;
}

vec2 GaussianSummation_computeGaussian(GaussianSummation gaussianSummation, float i)
{
    return gaussianSummation.a * exp(-i * i / gaussianSummation.d);
}

#define fragColor gl_FragColor
#define fragCoord gl_FragCoord
#define iFrame FRAMEINDEX
#define iResolution RENDERSIZE


void addCell(inout float new, vec2 nomralizedCoordinate)
{
    // from chronos' SmoothLife shader https://www.shadertoy.com/view/XtdSDn
    float dst = length(fragCoord.xy - nomralizedCoordinate * RENDERSIZE);
    if (dst <= or) {
    	new = step((ir + 1.5), dst) * (1. - step(or, dst));
    }
}

void main()
{
    vec2 tx = 1. / iResolution.xy;
    vec2 uv = fragCoord.xy * tx;

    if (PASSINDEX == 0) // ShaderToy Buffer A
    {
        vec4 current = IMG_NORM_PIXEL(bufferA, uv);
        vec2 fullness = IMG_NORM_PIXEL(bufferC, uv).xy;

        float delta =  2. * s(fullness.x, fullness.y) - 1.;
        float new = clamp(current.x + dt * delta, 0., 1.);

        if (addCellsWithMouse) {
            addCell(new, mouse.xy);
        }

        // For unclear reasons, iFrame must be strictly less than 2, not 1,
        // here.
        if (iFrame < 2 || addCells) {
#ifdef VIDEOSYNC
            float initialCellCount = min(iResolution.x, iResolution.y) / 50.;
#else
            const float initialCellCount = 20.;
#endif
            for (float i = 0.; i < initialCellCount; i++) {
                float angle = 6.2831853072 * // 2 pi
                              i / initialCellCount;
               	vec2 initialCoordinate = 0.25 * vec2(cos(angle), sin(angle)) + 0.5;
                addCell(new, initialCoordinate);
            }
        }

        fragColor = vec4(new, fullness, current.w);
    }
    else if (PASSINDEX == 1) // ShaderToy Buffer B
    {
        if (mod(iFrame, 2) < 1) {
            tx.y = 0.;
        } else {
            tx.x = 0.;
        }

        GaussianSummation gaussianSummation = GaussianSummation_create(or, ir);

        // Incredibly, GLSL and ISF have so many limitations that avoiding code
        // duplication for the Gaussian summation seems to be impossible. The
        // only things that vary between the two summations are the image name
        // (bufferA vs. bufferB) and the component access in the for-loops
        // (x vs. xy), but there doesn’t seem to be a way to encapsulate this in
        // GLSL.

        // centermost term
        gaussianSummation.acc += gaussianSummation.a * IMG_NORM_PIXEL(bufferA, uv).x;
        gaussianSummation.sum += gaussianSummation.a;

        // sum up remaining terms symmetrically
        for (int i = 1; i <= oc; i++) {
            float fi = float(i);
            vec2 g = GaussianSummation_computeGaussian(gaussianSummation, fi);
            vec2 posL = fract(uv - tx * fi);
            vec2 posR = fract(uv + tx * fi);
            gaussianSummation.acc += g * (IMG_NORM_PIXEL(bufferA, posL).x + IMG_NORM_PIXEL(bufferA, posR).x);
            gaussianSummation.sum += 2. * g;
        }

        vec2 x_pass = gaussianSummation.acc / gaussianSummation.sum;

        fragColor = vec4(x_pass, 0., 0.);
    }
    else if (PASSINDEX == 2) // ShaderToy Buffer C
    {
        if (mod(iFrame, 2) < 1) {
            tx.x = 0.;
        } else {
            tx.y = 0.;
        }

        GaussianSummation gaussianSummation = GaussianSummation_create(or, ir);

        // centermost term
        gaussianSummation.acc += gaussianSummation.a * IMG_NORM_PIXEL(bufferB, uv).x;
        gaussianSummation.sum += gaussianSummation.a;

        // sum up remaining terms symmetrically
        for (int i = 1; i <= oc; i++) {
            float fi = float(i);
            vec2 g = GaussianSummation_computeGaussian(gaussianSummation, fi);
            vec2 posL = fract(uv - tx * fi);
            vec2 posR = fract(uv + tx * fi);
            gaussianSummation.acc += g * (IMG_NORM_PIXEL(bufferB, posL).xy + IMG_NORM_PIXEL(bufferB, posR).xy);
            gaussianSummation.sum += 2. * g;
        }

        vec2 y_pass = gaussianSummation.acc / gaussianSummation.sum;

        fragColor = vec4(y_pass, 0., 0.);
    }
    else if (PASSINDEX == 3) // ShaderToy Image
    {
    	vec4 col = IMG_NORM_PIXEL(bufferA, uv);

        fragColor = vec4(col.x * vec3(1.) + col.y * vec3(1., 0.5, 0.) + col.z * vec3(0., 0.5, 1.), 1.);
    }
}
