precision mediump float;
uniform sampler2D tDiffuse;
uniform sampler2D tBackground;
uniform float uBlobThreshold;
varying vec2 vUv;
uniform vec2 uTexelSize;
uniform vec2 uPixelCount;
uniform float uSquarePixels;
uniform float uGradScale;

vec3 softchecker(in vec2 uv, float smoothness) {
    vec2 stripesi = mod(uv - 0.5, 2.);
    vec2 stripes = smoothstep(0.5 - smoothness, 0.5 + smoothness, abs(stripesi - 1.));
    return vec3(mix(vec3(0.2, 0.5, 0.7), vec3(0.2, 0.2, 0.4), float(abs(stripes.x - stripes.y))));
}

vec3 checker(in vec2 uv) {
    ivec2 id = ivec2(floor(uv));
    ivec2 stripes = id % 2;
    return vec3(mix(vec3(0.2, 0.3, 0.5), vec3(0.2, 0.1, 0.1), float(abs(stripes.x - stripes.y))));
}

void main() {
    vec3 c = texture2D(tDiffuse, vUv).rgb;
    vec2 bgUv = vUv / uPixelCount / uSquarePixels;
    vec3 bg = checker(bgUv);
    vec3 finalColor = bg * 0.95;

    if (c.g > 0. || c.r > 0.) {
        vec3 accum = vec3(0.);
        for (int x = 0; x <= 3; x++) {
            for (int y = 0; y <= 3; y++) {
                accum += texture2D(tDiffuse, vUv + uTexelSize * vec2(float(x), float(y))).rgb;
            }
        }
        accum /= 16.;
        c.rg = accum.rg;
    }

    if (c.g > uBlobThreshold || c.r > uBlobThreshold) {
        vec2 rGrad = vec2(dFdx(c.r), dFdy(c.r)) * uGradScale;
        vec2 gGrad = vec2(dFdx(c.g), dFdy(c.g)) * uGradScale;
        finalColor = (c.g > c.r) ? vec3(0.2, 0.7, 0.2) : vec3(0.9, 0.2, 0.);
        // float w = c.r / (c.r + c.g);
        // finalColor = mix(vec3(0, 0.6, 0.2), vec3(0.8, 0.2, 0.0), w);
        vec3 normal = (c.g > c.r) ? normalize(vec3(-gGrad, 1. * c.g * c.g)) : normalize(vec3(-rGrad, 1. * c.r * c.r)); 
        // vec3 normal = mix(normalize(vec3(-gGrad, 0.5 * c.g * c.g)), normalize(vec3(-rGrad, 0.5 * c.r * c.r)), 0.); 

        vec3 halfvec = normalize(vec3(vec2(2.), 1.));
        finalColor += 0.8 * smoothstep(0.8, 0.9, dot(normal, halfvec));
        halfvec = normalize(vec3(-vec2(2.), 1.));
        finalColor += 0.6 * vec3(0.2, 0.2, 0.2) * smoothstep(0.2, 0.8, abs(dot(normal, halfvec)));
        finalColor = mix(finalColor, softchecker(bgUv - normal.xy * 1., 0.05) * 1., 0.5);
        finalColor *= 1.1;
        // finalColor = normal;
        // finalColor = vec3(c.r);
        gl_FragColor = vec4(pow(finalColor, vec3(1.0)), 1.0);
    } else {
        vec2 w = vUv;
        float b = w.x / ( 0.05 + w.x) * ( 1.0 - w.x) / ( 1.05 - w.x) * w.y / ( 0.1 + w.y);
        finalColor *= 0.6 + 1. * b;
        gl_FragColor = vec4(pow(finalColor, vec3(1.0 / 1.2)), 1.0);
    }
}