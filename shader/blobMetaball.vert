uniform float uSize;
uniform float pointsPerBlob;
out float colorFlag;

void main() {
    gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0);
    gl_PointSize = uSize;
    colorFlag = float((gl_VertexID / int(pointsPerBlob)) % 2); 
}