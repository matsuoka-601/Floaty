precision mediump float;
uniform float uStrength;
void main() {
    float dist = length(gl_PointCoord - vec2(0.5));
    float alpha = uStrength * smoothstep(pow(0.5, 1./2.), 0., pow(dist, 1./2.)); // ここに係数をつける
    gl_FragColor = vec4(vec3(0, 0, 1), alpha);
}