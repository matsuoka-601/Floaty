precision mediump float;
in float colorFlag;

void main() {
    float dist = length(gl_PointCoord - vec2(0.5));
    float alpha = (0.6) * smoothstep(pow(0.5, 1./4.), 0., pow(dist, 1./4.)); // ここに係数をつける
    if (int(floor(colorFlag)) == 1) {
        gl_FragColor = vec4(vec3(1, 0, 0), alpha);
    } else {
        gl_FragColor = vec4(vec3(0, 1, 0), alpha);
    }
}