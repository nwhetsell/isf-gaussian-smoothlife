/*{
    "CATEGORIES": [
        "Filter",
        "Generator"
    ],
    "CREDIT": "cornusammonis <https://www.shadertoy.com/user/cornusammonis>",
    "DESCRIPTION": "Variant of Stephan Rafler’s SmoothLife that uses a separable Gaussian kernel to compute inner and outer fullness, converted from <https://www.shadertoy.com/view/XtVXzV>",
    "INPUTS": [
        {
            "NAME" : "inputImage",
            "TYPE" : "image"
        },
        {
            "NAME": "inputImageAmount",
            "LABEL": "Input image amount",
            "TYPE": "float",
            "DEFAULT": 0,
            "MIN": 0,
            "MAX": 1
        },
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
            "TARGET": "cells",
            "PERSISTENT": true,
            "FLOAT": true
        },
        {
            "TARGET": "summation1",
            "PERSISTENT": true,
            "FLOAT": true
        },
        {
            "TARGET": "summation2",
            "PERSISTENT": true,
            "FLOAT": true
        },
        {

        }
    ]
}
*/

#define SQRT_TWO_PI 2.5066282746310005024157652848110

// Constants and functions from LYGIA <https://github.com/patriciogonzalezvivo/lygia>
#define TWO_PI 6.2831853071795864769252867665590

vec2 polar2cart(in vec2 polar) {
    return vec2(cos(polar.x), sin(polar.x)) * polar.y;
}

//
// ShaderToy Buffer A
//

// the logistic function is used as a smooth step function
float logistic(float x, float midpoint, float quarterInverseSteepness)
{
    float steepness = 4. / quarterInverseSteepness;
    return 1. / (1. + exp(-steepness * (x - midpoint)));
}

float logisticPulse(float x, float midpoint1, float midpoint2, float quarterInverseSteepness)
{
    // In general, midpoint1 should be less than midpoint2.
    return logistic(x, midpoint1, quarterInverseSteepness) * (1. - logistic(x, midpoint2, quarterInverseSteepness));
}

float transformedLogistic(float x, float yShift, float maxValue, float quarterInverseSteepness)
{
    // The original ShaderToy shader effectively uses:
    //    return yShift * (1. - logistic(x, 0.5, quarterInverseSteepness)) + maxValue * logistic(x, 0.5, quarterInverseSteepness);
    // Simplifying this shows that this isn’t really two logistic functions.
    // It’s actually one logistic function shifted vertically by yShift with a
    // maximum value of maxValue:
    //
    //          ╭              1          ╮                       1
    // yShift × │ 1 - ------------------- │ + maxValue × -------------------
    //          ╰     1 + exp(-(x - 0.5)) ╯              1 + exp(-(x - 0.5))
    //
    //
    //                yShift                maxValue
    // yShift - ------------------- + -------------------
    //          1 + exp(-(x - 0.5))   1 + exp(-(x - 0.5))
    //
    //           maxValue - yShift
    // yShift + -------------------
    //          1 + exp(-(x - 0.5))
    return yShift + (maxValue - yShift) * logistic(x, 0.5, quarterInverseSteepness);
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
    gaussianSummation.a = vec2(1. / (sigma * SQRT_TWO_PI));
    gaussianSummation.d = vec2(2. * sigma * sigma);
    gaussianSummation.acc = vec2(0);
    gaussianSummation.sum = vec2(0);

    return gaussianSummation;
}

vec2 GaussianSummation_computeGaussian(GaussianSummation gaussianSummation, float i)
{
    return gaussianSummation.a * exp(-i * i / gaussianSummation.d);
}


void addCell(inout float new, vec2 nomralizedCoordinate)
{
    // from chronos' SmoothLife shader https://www.shadertoy.com/view/XtdSDn
    float dst = length(gl_FragCoord.xy - nomralizedCoordinate * RENDERSIZE);
    if (dst <= or) {
    	new = step((ir + 1.5), dst) * (1. - step(or, dst));
    }
}

void main()
{
    vec2 tx = 1. / RENDERSIZE;
    vec2 uv = gl_FragCoord.xy * tx;

    if (PASSINDEX == 0) // ShaderToy Buffer A
    {
        vec4 current = IMG_NORM_PIXEL(cells, uv);
        current += inputImageAmount * IMG_NORM_PIXEL(inputImage, uv);

        vec2 fullness = IMG_NORM_PIXEL(summation2, uv).xy;
        float outerFullness = fullness.x;
        float innerFullness = fullness.y;
        float midpoint1 = transformedLogistic(innerFullness, b1, s1, alpha_m);
        float midpoint2 = transformedLogistic(innerFullness, b2, s2, alpha_m);
        float delta = 2. * logisticPulse(outerFullness, midpoint1, midpoint2, alpha_n) - 1.;

        float new = clamp(current.x + dt * delta, 0., 1.);

        if (addCellsWithMouse) {
            addCell(new, mouse.xy);
        }

        // For unclear reasons, FRAMEINDEX must be strictly less than 2, not 1,
        // here.
        if (FRAMEINDEX < 2 || addCells) {
#ifdef VIDEOSYNC
            float initialCellCount = min(RENDERSIZE.x, RENDERSIZE.y) / 50.;
#else
            const float initialCellCount = 20.;
#endif
            for (float i = 0.; i < initialCellCount; i++) {
               	vec2 initialCoordinate = polar2cart(vec2(TWO_PI * i / initialCellCount, 0.25)) + 0.5;
                addCell(new, initialCoordinate);
            }
        }

        gl_FragColor = vec4(new, fullness, current.w);
    }
    else if (PASSINDEX == 1) // ShaderToy Buffer B
    {
        if (mod(float(FRAMEINDEX), 2.) < 1.) {
            tx.y = 0.;
        } else {
            tx.x = 0.;
        }

        GaussianSummation gaussianSummation = GaussianSummation_create(or, ir);

        // Incredibly, GLSL and ISF have so many limitations that avoiding code
        // duplication for the Gaussian summation seems to be impossible. The
        // only things that vary between the two summations are the image name
        // (cells vs. summation1) and the component access in the for-loops
        // (x vs. xy), but there doesn’t seem to be a way to encapsulate this in
        // GLSL.

        // centermost term
        gaussianSummation.acc += gaussianSummation.a * IMG_NORM_PIXEL(cells, uv).x;
        gaussianSummation.sum += gaussianSummation.a;

        // sum up remaining terms symmetrically
        for (int i = 1; i <= oc; i++) {
            float fi = float(i);
            vec2 g = GaussianSummation_computeGaussian(gaussianSummation, fi);
            vec2 posL = fract(uv - tx * fi);
            vec2 posR = fract(uv + tx * fi);
            gaussianSummation.acc += g * (IMG_NORM_PIXEL(cells, posL).x + IMG_NORM_PIXEL(cells, posR).x);
            gaussianSummation.sum += 2. * g;
        }

        vec2 x_pass = gaussianSummation.acc / gaussianSummation.sum;

        gl_FragColor = vec4(x_pass, 0, 1);
    }
    else if (PASSINDEX == 2) // ShaderToy Buffer C
    {
        if (mod(float(FRAMEINDEX), 2.) < 1.) {
            tx.x = 0.;
        } else {
            tx.y = 0.;
        }

        GaussianSummation gaussianSummation = GaussianSummation_create(or, ir);

        // centermost term
        gaussianSummation.acc += gaussianSummation.a * IMG_NORM_PIXEL(summation1, uv).x;
        gaussianSummation.sum += gaussianSummation.a;

        // sum up remaining terms symmetrically
        for (int i = 1; i <= oc; i++) {
            float fi = float(i);
            vec2 g = GaussianSummation_computeGaussian(gaussianSummation, fi);
            vec2 posL = fract(uv - tx * fi);
            vec2 posR = fract(uv + tx * fi);
            gaussianSummation.acc += g * (IMG_NORM_PIXEL(summation1, posL).xy + IMG_NORM_PIXEL(summation1, posR).xy);
            gaussianSummation.sum += 2. * g;
        }

        vec2 y_pass = gaussianSummation.acc / gaussianSummation.sum;

        gl_FragColor = vec4(y_pass, 0, 1);
    }
    else // ShaderToy Image
    {
    	vec4 color = IMG_NORM_PIXEL(cells, uv);

        gl_FragColor = vec4(color.r * vec3(1) + color.g * vec3(1, 0.5, 0) + color.b * vec3(0, 0.5, 1), 1);
    }
}
