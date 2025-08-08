import * as THREE from 'three';

export const maxParticles = 10000;
export const maxNeighbors = 64;

export interface MouseState {
    screen: THREE.Vector2;  // NDC座標 (-1〜1)
    world: Vec2;            // z=0平面上のワールド座標
    isDragging: boolean;    // マウスボタンが押されているか
}

export class Particle {
    pos: Vec2;
    prevPos: Vec2;
    v: Vec2;

    constructor(v: Vec2) {
        // そのまま x にすると soa2aos でおかしくなる（調べる）
        this.pos = new Vec2(v.x, v.y);
        this.prevPos = new Vec2(v.x, v.y);
        this.v = new Vec2(0, 0);
    }
}

export class SimulationParams {
    distStiffness: number;
    pairStiffness: number;
    volumeStiffness: number;
    volumeScale: number;
    substeps: number;
    eps: number;
    cfl: number;
    restDensity: number; 
    running: boolean;
}
export class BlobParams {
    particlesPerBlob: number;
	blobRadius: number;
	blobMass: number;
}

export default class Vec2 {
    constructor(public x: number, public y: number) {}
    
    add(v: Vec2): Vec2 {
        return new Vec2(this.x + v.x, this.y + v.y);
    }
    
    subtract(v: Vec2): Vec2 {
        return new Vec2(this.x - v.x, this.y - v.y);
    }
    
    multiply(scalar: number): Vec2 {
        return new Vec2(this.x * scalar, this.y * scalar);
    }
    
    divide(scalar: number): Vec2 {
        return new Vec2(this.x / scalar, this.y / scalar);
    }

    rotate90CW(): Vec2 {
        return new Vec2(this.y, -this.x);
    }
    
    length(): number {
        return Math.sqrt(this.x * this.x + this.y * this.y);
    }

    dot(v: Vec2): number {
        return this.x * v.x + this.y * v.y;
    }

    cross(v: Vec2): number {
        return this.x * v.y - this.y * v.x;
    }

    squaredLength(): number {
        return this.x * this.x + this.y * this.y;
    }
    
    normalize(): Vec2 {
        const len = this.length();
        if (len === 0) return new Vec2(0, 0);
        return new Vec2(this.x / len, this.y / len);
    }
}