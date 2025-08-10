precision mediump float;
uniform sampler2D tDiffuse;
uniform sampler2D tBackground;
varying vec2 vUv;
uniform vec2 uPixelCount;
uniform vec2 uTexelSize;
uniform vec3 uFluidColor;
uniform float uFluidThreshold;
uniform float uBlobThreshold;
uniform float uSquarePixels;
uniform float uGradScale;
uniform float t;

vec3 softchecker(in vec2 uv, float smoothness) {
    vec2 stripesi = mod(uv - 0.5, 2.);
    vec2 stripes = smoothstep(0.5 - smoothness, 0.5 + smoothness, abs(stripesi - 1.));
    return vec3(mix(vec3(0.2, 0.5, 0.7), vec3(0.2, 0.2, 0.4), float(abs(stripes.x - stripes.y))));
}

void main() {
    vec3 c = texture2D(tDiffuse, vUv).rgb;
    vec3 back = texture2D(tBackground, vUv).rgb;
    vec3 finalColor = back;
    vec2 bgUv = vUv / uPixelCount / uSquarePixels;


    if (c.b > 0.) {
        vec3 accum = vec3(0.);
        for (int x = 0; x <= 2; x++) {
            for (int y = 0; y <= 2; y++) {
                accum += texture2D(tDiffuse, vUv + uTexelSize * vec2(float(x), float(y))).rgb;
            }
        }
        accum /= 9.;
        c.b = accum.b;
    }

    if (c.r > 0. || c.g > 0.) {
        vec3 accum = vec3(0.);
        for (int x = 0; x <= 3; x++) {
            for (int y = 0; y <= 3; y++) {
                accum += texture2D(tDiffuse, vUv + uTexelSize * vec2(float(x), float(y))).rgb;
            }
        }
        accum /= 16.;
        c.rg = accum.rg;
    }

    vec2 grad = vec2(dFdx(c.b), dFdy(c.b)) * uGradScale;

    if (c.b > uFluidThreshold) {
        finalColor = uFluidColor;

        vec2 v = smoothstep(uFluidThreshold * 1.7, uFluidThreshold, c.b) * grad;
        vec2 v_ = 0.03 * grad;
        // From Maeda Mameo's "aka 6" : https://www.mameson.com/experiment/glsl/aka_6/aka_6.html
        vec2 w = 0.5 + ( 0.5 + length(vUv - 0.5)) * (vUv - 0.5);
        vec2 refrac = 2.0 * sin( t + 30.0 * w.yx) * ( sin( 70.0 * w.yx) + sin( 30.0 * w));
        float edgeFlag = smoothstep(uFluidThreshold * 1.5, uFluidThreshold, c.b);

        if (c.r < uBlobThreshold && c.g < uBlobThreshold) {
            finalColor = mix(finalColor, softchecker(bgUv + (v + v_) * 30. + 0.01 * refrac, 0.0).rgb, (edgeFlag > 0.) ? 0.8 : 0.6);
        } else {
            finalColor = mix(finalColor, back.rgb, (edgeFlag > 0.) ? 0.8 : 0.6);
        }   

        vec3 normal = 1. * normalize(vec3(-grad, 0.4 * c.b * c.b)); 
        finalColor += 0.5 * smoothstep(0.94, 0.96, dot(normal, normalize(vec3( -0.6, 0.55, 0.577))));
        finalColor -= 0.2 * smoothstep(0.6, 1., dot(normal, vec3( 0, -1, 0)));
        finalColor = pow(finalColor, vec3(1.0 / 1.));

        // From Maeda Mameo's "aka 6" : https://www.mameson.com/experiment/glsl/aka_6/aka_6.html
        w = vUv;
        float b = w.x / ( 0.03 + w.x) * ( 1.0 - w.x) / ( 1.03 - w.x) * w.y / ( 0.1 + w.y);
        finalColor *= 0.6 + 0.9 * b;
        finalColor *= 1.1;
    } 
    gl_FragColor = vec4(pow(finalColor, vec3(1.0 / 1.)), 1.0);
    gl_FragColor.b *= 0.97;
}