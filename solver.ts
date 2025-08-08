import { RGBA_ASTC_10x10_Format } from 'three';
import Vec2 from './common';
import { maxParticles, maxNeighbors, MouseState, Particle } from './common';
import { SimulationParams, BlobParams } from './common';
import { solve_fluid_jacobi, alloc_buffer, NeighborLinkedList } from './pkg/out.js';

export type BoundingBox = { minX: number; maxX: number; minY: number, maxY: number };

export class DistanceConstraint {
    restDistance: number;
    id0: number;
    id1: number;

    constructor(restDistance: number, id0: number, id1: number) {
        this.restDistance = restDistance;
        this.id0 = id0;
        this.id1 = id1;
    }

    project(particles: Particle[], stiffness: number, invM: number) {
        let dx = particles[this.id1].pos.subtract(particles[this.id0].pos);
        let dist = dx.length();
        let w0 = invM;
        let w1 = invM;
        let corr = 0;
        if (dist > 0) {
            corr = (this.restDistance - dist) / dist / (w0 + w1) * stiffness;
        }
        particles[this.id0].pos = particles[this.id0].pos.subtract(dx.multiply(w0 * corr));
        particles[this.id1].pos = particles[this.id1].pos.add(dx.multiply(w1 * corr));
    }

    static mouseProject(particles: Particle[], mouseCoord: Vec2, draggedParticleIdx: number, blobToParticleIdx: number[][], blobIdx: number) {
        let target = mouseCoord;
        let dx = target.subtract(particles[draggedParticleIdx].pos);
        let halfDraggedParticleCnt = 5;
        for (var i = 0; i < blobToParticleIdx[blobIdx].length; i++) {
            let currentParticleIdx = blobToParticleIdx[blobIdx][i];
            if (Math.abs(currentParticleIdx - draggedParticleIdx) <= halfDraggedParticleCnt) {
                particles[currentParticleIdx].pos = particles[currentParticleIdx].pos.add(dx.multiply(0.002));
            }
        }
    }
}

export class VolumeConstraint {
    restVolume: number;
    ids: number[];

    constructor(restVolume: number, ids: number[]) {
        this.restVolume = restVolume;
        this.ids = ids;
    }

    project(particles: Particle[], stiffness: number, volumeScale: number, invM: number) {
        let n = this.ids.length;
        let V = 0.;
        for (var i = 0; i < n; i++) {
            let id0 = this.ids[i];
            let id1 = this.ids[(i+1) % n];
            let pos0 = particles[id0].pos;
            let pos1 = particles[id1].pos;
            V += (pos0.x * pos1.y - pos0.y * pos1.x);
        }
        V *= 0.5;

        let gradients: Vec2[] = [];
        let gradSquaredSum = 0;
        for (var i = 0; i < n; i++) {
            let idPrev = this.ids[(i-1+n) % n];
            let idNext = this.ids[(i+1) % n];
            let prev = particles[idPrev].pos;
            let next = particles[idNext].pos;
            let gradX = 0.5 * (next.y - prev.y);
            let gradY = 0.5 * (prev.x - next.x);
            let gradient = new Vec2(gradX, gradY);
            gradSquaredSum += gradient.squaredLength() * invM;
            gradients.push(gradient);
        }

        if (gradSquaredSum > 0.) {
            let C = V - volumeScale * this.restVolume;
            let s = C / gradSquaredSum;
            for (var i = 0; i < n; i++) {
                let id = this.ids[i];
                let dx = gradients[i].multiply(-stiffness * s * invM);
                particles[id].pos = particles[id].pos.add(dx);
            }
        }
    }
}

export default class PBD {
    distanceConstraints: DistanceConstraint[] = [];
    volumeConstraints: VolumeConstraint[] = [];
    blobParticles: Particle[] = [];
    fluidParticles: Particle[] = [];
    initPos: Vec2[] = [];
    prevIsDragging: boolean;
    draggedIdx: number;
    width: number;
    height: number;
    initWidth: number;
    initHeight: number;

    particleCount = 0;
    firstNeighbor: Int32Array;
    firstNeighborPtr: any;
    neighborsArray: Int32Array;
    neighborsArrayPtr: any;
    blobFirstNeighbor: Int32Array; 
    blobFirstNeighborPtr: any;
    blobNeighborsArray: Int32Array; 
    blobNeighborsArrayPtr: any;

    posX: Float32Array;
    posY: Float32Array;
    posXptr: any;
    posYptr: any;
    
    lambdas: Float32Array;
    lambdasPtr: any;
    gradsX: Float32Array;
    gradsXPtr: any; 
    gradsY: Float32Array;
    gradsYPtr: any;
    fluidDeltaX: Float32Array;
    fluidDeltaXPtr: any;
    fluidDeltaY: Float32Array;
    fluidDeltaYPtr: any;

    neighborLinkedList: NeighborLinkedList;

    constructor(width: number, height: number, initSpacing: number, wasm: any) {
        let numX = 60;
        let numY = 50;
        let left = 0.1 * width;
        let bottom = 0.1 * height;

        this.posXptr = alloc_buffer(4 * maxParticles);
        this.posYptr = alloc_buffer(4 * maxParticles);
        this.posX = new Float32Array(wasm.memory.buffer, this.posXptr, maxParticles);
        this.posY = new Float32Array(wasm.memory.buffer, this.posYptr, maxParticles);

        let cnt = 0;
        for (var i = 0; i < numY; i++) {
            for (var j = 0; j < numX; j++) {
                let pos = new Vec2(left + j * 1.0 * initSpacing + 0.0 * Math.random() + 10, bottom + i * 1.0 * initSpacing + 0.0 * Math.random());
                this.fluidParticles.push(new Particle(pos));
                this.posX[cnt] = pos.x;
                this.posY[cnt] = pos.y;
                cnt++;
            }
        }

        this.width = width;
        this.height = height;
        this.initWidth = width;
        this.initHeight = height;

        this.firstNeighborPtr = alloc_buffer(4 * (maxParticles + 1));
        this.firstNeighbor = new Int32Array(wasm.memory.buffer, this.firstNeighborPtr, maxParticles + 1);
        const maxNeighborParticles = maxNeighbors * maxParticles;
        this.neighborsArrayPtr = alloc_buffer(4 * maxNeighborParticles); 
        this.neighborsArray = new Int32Array(wasm.memory.buffer, this.neighborsArrayPtr, maxNeighborParticles);

        this.blobFirstNeighborPtr = alloc_buffer(4 * (maxParticles + 1));
        this.blobFirstNeighbor = new Int32Array(wasm.memory.buffer, this.blobFirstNeighborPtr, maxParticles + 1);
        this.blobNeighborsArrayPtr = alloc_buffer(4 * maxNeighborParticles);
        this.blobNeighborsArray = new Int32Array(wasm.memory.buffer, this.blobNeighborsArrayPtr, maxNeighborParticles);

        this.lambdasPtr = alloc_buffer(4 * maxParticles);
        this.lambdas = new Float32Array(wasm.memory.buffer, this.lambdasPtr, maxParticles);
        this.gradsXPtr = alloc_buffer(4 * maxParticles);
        this.gradsX = new Float32Array(wasm.memory.buffer, this.gradsXPtr, maxParticles);
        this.gradsYPtr = alloc_buffer(4 * maxParticles);
        this.gradsY = new Float32Array(wasm.memory.buffer, this.gradsYPtr, maxParticles);
        this.fluidDeltaXPtr = alloc_buffer(4 * maxParticles);
        this.fluidDeltaX = new Float32Array(wasm.memory.buffer, this.fluidDeltaXPtr, maxParticles);
        this.fluidDeltaYPtr = alloc_buffer(4 * maxParticles);
        this.fluidDeltaY = new Float32Array(wasm.memory.buffer, this.fluidDeltaYPtr, maxParticles);

        this.neighborLinkedList = new NeighborLinkedList(maxParticles);
    }

    addParticle(x: Vec2) {
        this.blobParticles.push(new Particle(x));
        this.initPos.push(x);
        this.particleCount++;
    }

    addDistanceConstraint(restDistance: number, id0: number, id1: number) {
        this.distanceConstraints.push(new DistanceConstraint(restDistance, id0, id1));
    }

    addVolumeConstraint(restVolume: number, ids: number[]) {
        this.volumeConstraints.push(new VolumeConstraint(restVolume, ids));
    }

    distanceToVector(v: Vec2, p: Vec2, start: Vec2) {
        let cross = Math.abs((p.x - start.x) * v.y - (p.y - start.y) * v.x);
        return cross / Math.sqrt(v.x * v.x + v.y * v.y); 
    }

    isPointInsidePolygon(p: Vec2, indices: number[]) {
        let inside = false;
        let n = indices.length;
        // just a simple raycasting
        for (var i = 0; i < n; i++) {
            let j = (i + 1) % n;
            let xi = this.blobParticles[indices[i]].pos.x;
            let yi = this.blobParticles[indices[i]].pos.y;
            let xj = this.blobParticles[indices[j]].pos.x;
            let yj = this.blobParticles[indices[j]].pos.y;

            if ((yi > p.y) != (yj > p.y)) {
                let intersectX = xi + (p.y - yi) * (xj - xi) / (yj - yi);
                if (intersectX > p.x) {
                    inside = !inside;
                }
            }
        }
        return inside;
    }

    solveFluidBoundaries() {
        for (var i = 0; i < this.fluidParticles.length; i++) {
            let px = this.posX[i];
            let py = this.posY[i];
            let epsilon = 1e-2;
            if (px < 0.) { this.posX[i] = epsilon * Math.random(); }
            if (px > this.width) { this.posX[i] = this.width - epsilon * Math.random(); }
            if (py < 0.) { this.posY[i] = epsilon * Math.random() }
            if (py > this.height) { this.posY[i] = this.height - epsilon * Math.random(); }
        }
    }
      
    boundingBoxOverlapCheck(a: BoundingBox, b: BoundingBox) {
        let xOverlap = !((a.maxX < b.minX) || (b.maxX < a.minX));
        let yOverlap = !((a.maxY < b.minY) || (b.maxY < a.minY));
        return xOverlap && yOverlap;
    }

    aos2soa() {
        for (var i = 0; i < this.fluidParticles.length + this.blobParticles.length; i++) {
            if (i < this.fluidParticles.length) {
                this.posX[i] = this.fluidParticles[i].pos.x;
                this.posY[i] = this.fluidParticles[i].pos.y;
            } else {
                this.posX[i] = this.blobParticles[i - this.fluidParticles.length].pos.x;
                this.posY[i] = this.blobParticles[i - this.fluidParticles.length].pos.y;
            }
        }
    }

    soa2aos() {
        for (var i = 0; i < this.fluidParticles.length + this.blobParticles.length; i++) {
            if (i < this.fluidParticles.length) {
                this.fluidParticles[i].pos.x = this.posX[i];
                this.fluidParticles[i].pos.y = this.posY[i];
            } else {
                this.blobParticles[i - this.fluidParticles.length].pos.x = this.posX[i];
                this.blobParticles[i - this.fluidParticles.length].pos.y = this.posY[i];
            }
        }
    }

    makeParticleCollisionConstraintsForBlob(
        distThresh: number, 
        isInsideForBlob: boolean[][], pairConstraints: DistanceConstraint[], 
        particleIdToBlob: { [key: number]: number }
    ) {

        // just assigning [] does not initialize the array
        pairConstraints.length = 0;

        for (var i = 0; i < this.blobParticles.length; i++) {
            this.posX[i] = this.blobParticles[i].pos.x;
            this.posY[i] = this.blobParticles[i].pos.y;
        }

        // neighborhood search (only for blob particles)
        this.neighborLinkedList.find_neighbors(
            this.posXptr, this.posYptr, 
            this.blobFirstNeighborPtr, this.blobNeighborsArrayPtr, 
            this.blobParticles.length, 
            distThresh, maxNeighbors
        );

        for (var i = 0; i < this.blobParticles.length; i++) {
            let first = this.blobFirstNeighbor[i];
            let num = this.blobFirstNeighbor[i + 1] - first;
            for (let j = 0; j < num; j++) {
                let id = this.blobNeighborsArray[first + j];
                // continue when particle i is inside a blob corresponding to id
                if (isInsideForBlob[i][particleIdToBlob[id]]) {   
                    continue;
                }
                if (i < id && particleIdToBlob[i] != particleIdToBlob[id]) {
                    pairConstraints.push(new DistanceConstraint(distThresh, i, id));
                }
            }
        }

        for (var i = 0; i < this.blobParticles.length; i++) {
            this.blobParticles[i].pos.x = this.posX[i];
            this.blobParticles[i].pos.y = this.posY[i];
        }
    }

    blobBoundingBoxes(blobToParticleIdx: number[][]) {
        let boundingBoxes: BoundingBox[] = [];
        for (var i = 0; i < blobToParticleIdx.length; i++) {
            var minX = 1e9;
            var maxX = -1e9;
            var minY = 1e9;
            var maxY = -1e9;
            for (const idx of blobToParticleIdx[i]) {
                minX = Math.min(minX, this.blobParticles[idx].pos.x);
                maxX = Math.max(maxX, this.blobParticles[idx].pos.x);
                minY = Math.min(minY, this.blobParticles[idx].pos.y);
                maxY = Math.max(maxY, this.blobParticles[idx].pos.y);
            }
            boundingBoxes.push({minX, maxX, minY, maxY});
        }
        return boundingBoxes;
    }

    makeIsInsideForBlob(blobToParticleIdx: number[][], isOverlap: boolean[][], isInsideForBlob: boolean[][]) {
        let boundingBoxes = this.blobBoundingBoxes(blobToParticleIdx);
        for (var i = 0; i < blobToParticleIdx.length; i++) {
            for (var j = 0; j < blobToParticleIdx.length; j++) {
                if (i == j) {
                    continue;
                }

                if (this.boundingBoxOverlapCheck(boundingBoxes[i], boundingBoxes[j])) {
                    for (const idx of blobToParticleIdx[i]) {
                        if (this.isPointInsidePolygon(this.blobParticles[idx].pos, blobToParticleIdx[j])) {
                            isInsideForBlob[idx][j] = true;
                            isOverlap[i][j] = true;
                        }
                    }
                }
            }
        }
    }

    // check if each fluid particle is inside blobs
    makeIsInsideForFluid(blobToParticleIdx: number[][], isInsideForFluid: boolean[][]) {
        let boundingBoxes = this.blobBoundingBoxes(blobToParticleIdx);

        for (var i = 0; i < this.fluidParticles.length; i++) {
            for (var j = 0; j < boundingBoxes.length; j++) {
                let boundingBox = boundingBoxes[j];
                if (
                    boundingBox.minX <= this.fluidParticles[i].pos.x && 
                    this.fluidParticles[i].pos.x <= boundingBox.maxX && 
                    boundingBox.minY <= this.fluidParticles[i].pos.y &&
                    this.fluidParticles[i].pos.y <= boundingBox.maxY
                ) {
                    if (this.isPointInsidePolygon(this.fluidParticles[i].pos, blobToParticleIdx[j])) {
                        isInsideForFluid[i][j] = true;
                    }
                }
            }
        }
    }

    enforceRepulsiveForceForBlob(blobToParticleIdx: number[][], isOverlap: boolean[][], blobParams: BlobParams) {
        for (var i = 0; i < blobToParticleIdx.length; i++) { 
            for (var j = 0; j < blobToParticleIdx.length; j++) { 
                if (i == j) {
                    continue;
                }

                if (isOverlap[i][j]) {
                    let ci = new Vec2(0, 0);
                    let cj = new Vec2(0, 0);

                    for (const k of blobToParticleIdx[i]) { ci = ci.add(this.blobParticles[k].pos); }
                    for (const k of blobToParticleIdx[j]) { cj = cj.add(this.blobParticles[k].pos); }
                    ci = ci.multiply(1.0 / blobParams.particlesPerBlob);
                    cj = cj.multiply(1.0 / blobParams.particlesPerBlob);

                    let d = ci.subtract(cj);
                    let C = d.length() - blobParams.blobRadius * 4;  // just an emperical parameter
                    d = d.normalize();
                    let stiffness = 0.005;

                    for (const k of blobToParticleIdx[i]) { 
                        this.blobParticles[k].pos = this.blobParticles[k].pos.add(d.multiply(-C * stiffness));
                    }
                    for (const k of blobToParticleIdx[j]) { 
                        this.blobParticles[k].pos = this.blobParticles[k].pos.add(d.multiply(C * stiffness));
                    }
                }
            }
        }
    }

    solveFluid(isInsideForFluid: boolean[][], kernelRadius: number, 
        restDensity: number, eps: number, blobParams: BlobParams) {
        this.aos2soa();
        // if inside[i] = true, i'th fluid particle is inside a blob
        const inside = isInsideForFluid.map(row => row.some(v => v));
        for (var i = 0; i < this.fluidParticles.length; i++) {
            this.lambdas[i] = inside[i] ? 1. : 0;
        }
        solve_fluid_jacobi(
            this.posXptr, this.posYptr, 
            this.lambdasPtr, this.gradsXPtr, this.gradsYPtr, 
            this.fluidDeltaXPtr, this.fluidDeltaYPtr, 
            this.firstNeighborPtr, this.neighborsArrayPtr, 
            this.fluidParticles.length + this.blobParticles.length, this.fluidParticles.length, 
            kernelRadius, restDensity, eps, 
            1.0 / blobParams.blobMass, maxNeighbors
        );
        this.solveFluidBoundaries();
        for (var i = 0; i < this.fluidParticles.length + this.blobParticles.length; i++) { 
            this.posX[i] += this.fluidDeltaX[i];
            this.posY[i] += this.fluidDeltaY[i];
        }
        this.soa2aos();
    }

    findNeighbors(kernelRadius: number) {
        this.aos2soa();
        let gridSpacing = 1.3 * kernelRadius;
        this.neighborLinkedList.find_neighbors(
            this.posXptr, this.posYptr, 
            this.firstNeighborPtr, this.neighborsArrayPtr, 
            this.fluidParticles.length + this.blobParticles.length, 
            gridSpacing, maxNeighbors
        );
        this.soa2aos();
    }

    enforceCFL(maxVel: number, particle: Particle, dt: number) {
        let dx = particle.pos.x - particle.prevPos.x;
        let dy = particle.pos.y - particle.prevPos.y;
        const v = dx * dx + dy * dy;
        if (v > maxVel * maxVel) { // CFL
            dx *= maxVel / Math.sqrt(v);
            dy *= maxVel / Math.sqrt(v);
            particle.pos.x = particle.prevPos.x + dx;
            particle.pos.y = particle.prevPos.y + dy;
        }
        particle.v.x = dx / dt;
        particle.v.y = dy / dt;
    }

    solveFloorCollision(pos: Vec2, prevPos: Vec2, floorY: number) {
        if (pos.y < floorY) {
            const d = Math.abs(pos.y - floorY); // penetration distance (assuming radius 0)
            pos.y = floorY; 

            const deltaTangent = pos.x - prevPos.x;
            const muS = 0.5;
            const muK = 0.3;
            let deltaX = deltaTangent;
            if (Math.abs(deltaX) >= muS * d) {
                deltaX = deltaTangent * Math.min(muK * d / Math.abs(deltaTangent), 1);
            }
    
            pos.x -= deltaX;
        }
    }

    solveBlobBoundaries(particle: Particle) {
        let pos = particle.pos;
        let prevPos = particle.prevPos;

        const minX = 0;
        const maxX = this.width;
        const minY = 0;
        const maxY = this.height;
            
        if (pos.x < minX) { pos.x = minX; }
        if (pos.x > maxX) { pos.x = maxX; }
        this.solveFloorCollision(pos, prevPos, minY);
        if (pos.y > maxY) { pos.y = maxY; }
    }

    setNearestBlobIdToMouse(mouseState: MouseState) {
        let curMinDist = 1e32;
        for (var i = 0; i < this.blobParticles.length; i++) { 
            let dist = this.blobParticles[i].pos.subtract(mouseState.world).length();
            if (dist < curMinDist) {
                curMinDist = dist;
                this.draggedIdx = i;
            }
        }
    }

    simulate(
        dtFrame: number, gravity: Vec2, mouseState: MouseState, 
        params: SimulationParams, blobParams: BlobParams, 
        blobToParticleIdx: number[][], kernelRadius: number, 
        width: number, height: number
    ) {

        let dt = dtFrame / params.substeps;

        let pairConstraints: DistanceConstraint[] = [];
        // isInsideForBlob[i][j] = whether blob particle i is inside blob j 
        let isInsideForBlob: boolean[][] = Array.from({ length: this.blobParticles.length }, () => Array(blobToParticleIdx.length).fill(false));
        // isInsideForFluid[i][j] = whether fluid particle i is inside blob j 
        let isInsideForFluid: boolean[][] = Array.from({ length: this.fluidParticles.length }, () => Array(blobToParticleIdx.length).fill(false));
        let isOverlap: boolean[][] = Array.from({ length: blobToParticleIdx.length }, () => Array(blobToParticleIdx.length).fill(false));
        let particleIdToBlob = {};
        for (var i = 0; i < blobToParticleIdx.length; i++) {
            for (const j of blobToParticleIdx[i]) {
                particleIdToBlob[j] = i;
            }
        }

        this.makeIsInsideForBlob(blobToParticleIdx, isOverlap, isInsideForBlob);
        this.makeIsInsideForFluid(blobToParticleIdx, isInsideForFluid);
        this.findNeighbors(kernelRadius);

        this.width = width;
        this.height = height;

        for (var substep = 0; substep < params.substeps; substep++) {
            for (var i = 0; i < this.fluidParticles.length; i++) { 
                this.fluidParticles[i].v.x += gravity.x * dt;
                this.fluidParticles[i].v.y += gravity.y * dt;
                // ここはめちゃくちゃ注意する必要がある
                // pos をそのままコピーすると、参照がコピーされることになるため、意図しない挙動になってしまう
                // あとで詳しく分析
                this.fluidParticles[i].prevPos.x = this.fluidParticles[i].pos.x;
                this.fluidParticles[i].prevPos.y = this.fluidParticles[i].pos.y;
                this.fluidParticles[i].pos.x += this.fluidParticles[i].v.x * dt;
                this.fluidParticles[i].pos.y += this.fluidParticles[i].v.y * dt;
            }
            for (var i = 0; i < this.blobParticles.length; i++) {
                this.blobParticles[i].v.x += gravity.x * dt;
                this.blobParticles[i].v.y += gravity.y * dt;
                this.blobParticles[i].prevPos.x = this.blobParticles[i].pos.x;
                this.blobParticles[i].prevPos.y = this.blobParticles[i].pos.y;
                this.blobParticles[i].pos.x += this.blobParticles[i].v.x * dt;
                this.blobParticles[i].pos.y += this.blobParticles[i].v.y * dt;
            }

    
            this.solveFluid(isInsideForFluid, kernelRadius, params.restDensity, params.eps, blobParams);

            if (substep % 5 == 0) {
                const collisionDistThresh = 0.4;
                this.makeParticleCollisionConstraintsForBlob(collisionDistThresh, isInsideForBlob, pairConstraints, particleIdToBlob);
                this.enforceRepulsiveForceForBlob(blobToParticleIdx, isOverlap, blobParams);
            }

            for (var i = 0; i < this.distanceConstraints.length; i++) { 
                this.distanceConstraints[i].project(this.blobParticles, params.distStiffness, 1. / blobParams.blobMass);
            }
            for (var i = 0; i < this.volumeConstraints.length; i++) { 
                this.volumeConstraints[i].project(this.blobParticles, params.volumeStiffness, params.volumeScale, 1. / blobParams.blobMass);
            }
            for (var i = 0; i < pairConstraints.length; i++) {
                pairConstraints[i].project(this.blobParticles, params.pairStiffness, 1. / blobParams.blobMass);
            }

            if (!this.prevIsDragging && mouseState.isDragging) {
                this.setNearestBlobIdToMouse(mouseState);
            }

            if (mouseState.isDragging) {
                DistanceConstraint.mouseProject(this.blobParticles, mouseState.world, this.draggedIdx, blobToParticleIdx, particleIdToBlob[this.draggedIdx]);
            }

            let maxVel = params.cfl * kernelRadius;

            for (var i = 0; i < this.fluidParticles.length; i++) { 
                this.enforceCFL(0.9 * maxVel, this.fluidParticles[i], dt);
            }
            for (var i = 0; i < this.blobParticles.length; i++) { 
                this.solveBlobBoundaries(this.blobParticles[i]);
                this.enforceCFL(maxVel, this.blobParticles[i], dt);
            }
        }
        this.prevIsDragging = mouseState.isDragging;
    }
}