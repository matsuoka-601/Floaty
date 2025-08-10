import * as THREE from 'three';
import PBD from './solver';
import GUI from 'lil-gui';
import Vec2 from './common';
import { BlobParams, SimulationParams } from './common';
import { MouseState, Particle } from './common';
import init, { initThreadPool } from './pkg/out.js';
import blobVert from './shader/blobMetaball.vert'
import blobFrag from './shader/blobMetaball.frag'
import fluidVert from './shader/fluidMetaball.vert'
import fluidFrag from './shader/fluidMetaball.frag'
import blobThresholdVert from './shader/blobThreshold.vert';
import blobThresholdFrag from './shader/blobThreshold.frag';
import fluidThresholdVert from './shader/fluidThreshold.vert';
import fluidThresholdFrag from './shader/fluidThreshold.frag';

const worldSize = 16;
let fluidThresholdMaterial: THREE.ShaderMaterial;

function pixelToLength(pixelCount: number, pixelHeight: number, worldHeight: number) {
	return pixelCount * (worldHeight / pixelHeight);
}
function lengthToPixel(worldLength: number, pixelHeight: number, worldHeight: number) {
	return Math.floor(worldLength * (pixelHeight / worldHeight));
}
function calcSquarePixels() {
	const squareCount = 6;
	if (window.innerWidth > window.innerHeight) {
		return window.innerHeight / squareCount;
	} else {
		return window.innerWidth / squareCount;
	}
}
function calcGradScale() {
	return Math.min(window.innerWidth, window.innerHeight) / 671;
}
function makeMouseState(camera: THREE.OrthographicCamera) {
	const mouseState: MouseState = {
		screen: new THREE.Vector2(),
		world: new Vec2(0, 0),
		isDragging: false,
	};

	const raycaster = new THREE.Raycaster();
	const planeZ = new THREE.Plane(new THREE.Vector3(0, 0, 1), 0);
	const intersection = new THREE.Vector3();

	function updateMouseState(event: MouseEvent) {
		mouseState.screen.x = (event.clientX / window.innerWidth) * 2 - 1;
		mouseState.screen.y = - (event.clientY / window.innerHeight) * 2 + 1;

		raycaster.setFromCamera(mouseState.screen, camera);
		if (raycaster.ray.intersectPlane(planeZ, intersection)) {
			mouseState.world.x = intersection.x;
			mouseState.world.y = intersection.y;
		}
	}

	// const label = document.getElementById('mouseLabel')!;

	window.addEventListener('pointerdown', (event) => {
		updateMouseState(event);
		mouseState.isDragging = true;
		// label.style.display = 'block';
	});

	// prevent swipe back (only on iOS?)
	window.addEventListener('touchstart', (e) => {
		e.preventDefault();
	});

	window.addEventListener('pointermove', (event) => {
		updateMouseState(event);
	});

	window.addEventListener('pointerup', () => {
		mouseState.isDragging = false;
		// label.style.display = 'none';
	});

	return mouseState;
}

function initScene(
	fluidParticlesCount: number, blobCnt: number, particlesPerBlob: number
) {
	const scene = new THREE.Scene();
	// const camera = new THREE.PerspectiveCamera(75, window.innerWidth / window.innerHeight, 0.1, 1000);

	const camera = new THREE.OrthographicCamera(
		0,      // left
		0,      // right
		0,      // top
		0,      // bottom
		-1000,     // near
		+1000      // far
	);

	const renderer = new THREE.WebGLRenderer();
	document.body.appendChild(renderer.domElement);

	// camera.position.set(width / 2, height / 2, 10);
	// camera.lookAt(new THREE.Vector3(width / 2, height / 2, 0));

	let mouseState = makeMouseState(camera);

	const gui = new GUI();
	let folder = gui.addFolder('render params');

	// 三角形
	const triGeometry = new THREE.BufferGeometry();
	const triInitPosition = new Float32Array(blobCnt * particlesPerBlob * 3);
	triGeometry.setAttribute('position', new THREE.BufferAttribute(triInitPosition, 3));


	const blobMaterial = new THREE.ShaderMaterial({
		transparent: true,
		blending: THREE.AdditiveBlending,
		depthTest: false,
		uniforms: {
			uSize: { value: {} },
			pointsPerBlob: { value: particlesPerBlob } 
		},
		vertexShader: blobVert,
		fragmentShader: blobFrag
	});

	const triMesh = new THREE.Points(triGeometry, blobMaterial);
	scene.add(triMesh);

	// 点
	const pointGeometry = new THREE.BufferGeometry();
	const pointInitPosition = new Float32Array(fluidParticlesCount * 3);
	pointGeometry.setAttribute('position', new THREE.BufferAttribute(pointInitPosition, 3));

	// 3) PointsMaterial にセット
	const fluidMaterial = new THREE.ShaderMaterial({
		transparent: true,
		blending: THREE.AdditiveBlending,
		depthTest: false,
		uniforms: {
			uSize: { value: 16. },
			uStrength: { value: 0.4 }
		},
		vertexShader: fluidVert,
		fragmentShader: fluidFrag
	});
	const points = new THREE.Points(pointGeometry, fluidMaterial);
	scene.add(points);

	folder.add(fluidMaterial.uniforms.uSize, 'value', 0., 100).name('fluid radius');
	folder.add(fluidMaterial.uniforms.uStrength, 'value', 0., 2.0).name('strength');

	const downsampleFactor = 0.5;
	const renderTarget = new THREE.WebGLRenderTarget(window.innerWidth * downsampleFactor, window.innerHeight * downsampleFactor, {
		minFilter: THREE.LinearFilter,
		magFilter: THREE.LinearFilter,
		format: THREE.RGBAFormat
	});
	const temporaryTarget = new THREE.WebGLRenderTarget(window.innerWidth, window.innerHeight, {
		minFilter: THREE.LinearFilter,
		magFilter: THREE.LinearFilter,
		format: THREE.RGBAFormat
	});

	// ポストプロセス
	const quadScene = new THREE.Scene();
	const quadFinalScene = new THREE.Scene();
	const quadCamera = new THREE.OrthographicCamera(-1, 1, 1, -1, 0, 1);

	let thresholdParams = {
		blobThreshold: 0.3 * 0.25, 
		fluidThreshold: 0.15, 
	}

	const blobThresholdMaterial = new THREE.ShaderMaterial({
		uniforms: {
			tDiffuse: { value: renderTarget.texture },
			uBlobThreshold: { value: thresholdParams.blobThreshold },
			uTexelSize: { value: {} },
			uPixelCount: { value: {} },
			uSquarePixels: { value: {} },
			uGradScale: { value: {} }
		},
		vertexShader: blobThresholdVert,
		fragmentShader: blobThresholdFrag
	});

	const quad = new THREE.Mesh(new THREE.PlaneGeometry(2, 2), blobThresholdMaterial);
	quadScene.add(quad);

	fluidThresholdMaterial = new THREE.ShaderMaterial({
		uniforms: {
			tDiffuse: { value: renderTarget.texture },
			tBackground: { value: temporaryTarget.texture },
			uFluidThreshold: { value: thresholdParams.fluidThreshold },
			uBlobThreshold: { value: thresholdParams.blobThreshold },
			uTexelSize: { value: {} },
			uPixelCount: { value: {} },
			uSquarePixels: { value: {} },
			uGradScale: { value: {} }, 
			uFluidColor: { value: {} }, 
			t: { value: 0 }
		},
		vertexShader: fluidThresholdVert,
		fragmentShader: fluidThresholdFrag
	});

	let colorParams = {
		// colorObj: [25, 120, 210],
		// colorObj: [50, 100, 180]
		colorObj: [50, 100, 180],
	};
	fluidThresholdMaterial.uniforms.uFluidColor.value = 
		new THREE.Vector3(colorParams.colorObj[0] / 255, colorParams.colorObj[1] / 255, colorParams.colorObj[2] / 255);

	folder.addColor(colorParams, 'colorObj', 255).onChange(value => { fluidThresholdMaterial.uniforms.uFluidColor.value = new THREE.Vector3(value[0] / 255, value[1] / 255, value[2] / 255); });

	const quadFinal = new THREE.Mesh(new THREE.PlaneGeometry(2, 2), fluidThresholdMaterial);
	quadFinalScene.add(quadFinal);

	function onResize() {
		renderer.setSize(window.innerWidth, window.innerHeight);
		let [w, h] = initWorldSize(worldSize);
		camera.left = 0;
		camera.right = w;
		camera.top = h;
		camera.bottom = 0;
		camera.updateProjectionMatrix();
		renderTarget.setSize(
			Math.floor(window.innerWidth * downsampleFactor),
			Math.floor(window.innerHeight * downsampleFactor)
		);
		temporaryTarget.setSize(
			Math.floor(window.innerWidth),
			Math.floor(window.innerHeight)
		)

		if (window.innerWidth > window.innerHeight) {
			fluidMaterial.uniforms.uSize.value = lengthToPixel(pixelToLength(25, 671, 15), window.innerHeight, h);
			blobMaterial.uniforms.uSize.value = lengthToPixel(pixelToLength(30, 671, 15), window.innerHeight, h);
			console.log(blobMaterial.uniforms.uSize.value);
		} else { 
			fluidMaterial.uniforms.uSize.value = lengthToPixel(pixelToLength(25, 671, 15), window.innerWidth, w);
			blobMaterial.uniforms.uSize.value = lengthToPixel(pixelToLength(30, 671, 15), window.innerWidth, w);
		}

		fluidThresholdMaterial.uniforms.uTexelSize.value = new THREE.Vector2(1.0 / renderTarget.width, 1.0 / renderTarget.height);
		fluidThresholdMaterial.uniforms.uPixelCount.value = new THREE.Vector2(1.0 / window.innerWidth, 1.0 / window.innerHeight);
		fluidThresholdMaterial.uniforms.uSquarePixels.value = calcSquarePixels();
		fluidThresholdMaterial.uniforms.uGradScale.value = calcGradScale();
		blobThresholdMaterial.uniforms.uTexelSize.value = new THREE.Vector2(1.0 / renderTarget.width, 1.0 / renderTarget.height);
		blobThresholdMaterial.uniforms.uPixelCount.value = new THREE.Vector2(1.0 / window.innerWidth, 1.0 / window.innerHeight);
		blobThresholdMaterial.uniforms.uSquarePixels.value = calcSquarePixels();
		blobThresholdMaterial.uniforms.uGradScale.value = calcGradScale();

		renderer.render(scene, camera);
	}

	window.addEventListener('resize', onResize);
	onResize();

	folder.close();

	return { scene, quadScene, quadFinalScene, camera, quadCamera, renderer, mouseState, triGeometry, pointGeometry, renderTarget, temporaryTarget, gui };
}

function setPositionsBlobs(
	particles: Particle[], geometry: THREE.BufferGeometry, 
	particlesPerBlob: number, blobCnt: number
) {
	const positions = geometry.getAttribute('position') as THREE.BufferAttribute;

	let offset = 0;
	const scaleFromCenter = 0.9;
	for (let blob = 0; blob < blobCnt; blob++) {
		let center = new Vec2(0, 0);
		for (let i = 0; i < particlesPerBlob; i++) {
			center = center.add(particles[offset + i].pos);
		}
		center = center.multiply(1.0 / particlesPerBlob);
		for (let i = 0; i < particlesPerBlob; i++) {
			let diff = particles[offset + i].pos.subtract(center);
			let afterPos = center.add(diff.multiply(scaleFromCenter));
			positions.array[3 * (offset + i)] = afterPos.x;
			positions.array[3 * (offset + i) + 1] = afterPos.y;
		}
		offset += particlesPerBlob;
	}
	positions.needsUpdate = true;
}

function setPositionsFluid(particles: Particle[], geometry: THREE.BufferGeometry) {
	const positions = geometry.getAttribute('position') as THREE.BufferAttribute;
	for (var i = 0; i < particles.length; i++) {
		positions.array[3 * i] = particles[i].pos.x;
		positions.array[3 * i + 1] = particles[i].pos.y;
	}
	positions.needsUpdate = true;
}

function initWorldSize(basis: number) {
	if (window.innerHeight > window.innerWidth) { 
		return [basis, basis * (window.innerHeight / window.innerWidth)];
	} else { 
		return [basis * (window.innerWidth / window.innerHeight), basis];
	}
}

function initBlobParticles(
	pbdSimulator: PBD, blobCnt: number, blobParams: BlobParams
) {
	let vertexIndices: number[][] = [];

	for (var blob = 0; blob < blobCnt; blob++) {
		let indices: number[] = [];
		let dist = 0;
		let offset = [3.5 * (blob % 6) + 5, 3.5 * Math.floor(blob / 6) + 15];
		for (var i = 0; i < blobParams.particlesPerBlob; i++) {
			let pos0 = new Vec2(
				blobParams.blobRadius * Math.cos(2 * Math.PI / blobParams.particlesPerBlob * i) + offset[0], 
				blobParams.blobRadius * Math.sin(2 * Math.PI / blobParams.particlesPerBlob * i) + offset[1]
			);
			let pos1 = new Vec2(
				blobParams.blobRadius * Math.cos(2 * Math.PI / blobParams.particlesPerBlob * (i + 1)) + offset[0], 
				blobParams.blobRadius * Math.sin(2 * Math.PI / blobParams.particlesPerBlob * (i + 1)) + offset[1]
			);
			dist = pos0.subtract(pos1).length();
			pbdSimulator.addParticle(pos0);
			indices.push(i + blob * blobParams.particlesPerBlob);
		}

		for (var i = 0; i < blobParams.particlesPerBlob; i++) {
			pbdSimulator.addDistanceConstraint(dist, indices[i], indices[(i + 1) % blobParams.particlesPerBlob]);
		}

		vertexIndices.push(indices);
		pbdSimulator.addVolumeConstraint(
			blobParams.particlesPerBlob / 2 * Math.sin(2 * Math.PI / blobParams.particlesPerBlob) * Math.pow(blobParams.blobRadius, 2.), 
			vertexIndices[blob]
		);
	}

	return vertexIndices;
}

async function main() {
	let particleRadius = 0.15;
	let particleDiameter = 2.0 * particleRadius; 
	let kernelRadius = 3.0 * particleRadius;

	console.log("waiting");
	const wasm = await init();
	const numThreads = Math.min(8, navigator.hardwareConcurrency);
	await initThreadPool(numThreads)
		.then(() => console.log('Thread pool initialized'))
		.catch(e => console.error('Thread pool init failed', e));
	console.log("end waiting");
	document.getElementById("numThreads").textContent = `Using ${numThreads} threads`;

	let [initWidth, initHeight] = initWorldSize(worldSize);
	let pbdSimulator = new PBD(initWidth, initHeight, 0.8 * particleDiameter, wasm);
	let blobCnt = 10;
	let blobParams: BlobParams = {
		particlesPerBlob: 48, 
		blobRadius: 1.5, 
		blobMass: 2., 
	};
	let vertexIndices = initBlobParticles(pbdSimulator, blobCnt, blobParams);

	const { scene, quadScene, quadFinalScene, camera, quadCamera, renderer, mouseState, triGeometry, pointGeometry, renderTarget, temporaryTarget, gui } =
		initScene(pbdSimulator.fluidParticles.length, blobCnt, blobParams.particlesPerBlob);

	const params: SimulationParams = {
		distStiffness: 0.6,
		pairStiffness: 1.,
		volumeStiffness: 1.,
		volumeScale: 1.,
		substeps: 10,
		eps: 50,
		restDensity: 1.0 / (particleDiameter * particleDiameter) * 1.5,
		cfl: 0.12,
		running: true,
	};
	const folder = gui.addFolder('sim params');
	folder.add(params, 'distStiffness', 0., 1.0).name('distStiffness');
	folder.add(params, 'pairStiffness', 0., 1.0).name('pairStiffness');
	folder.add(params, 'volumeStiffness', 0., 1.0).name('volumeStiffness');
	folder.add(params, 'volumeScale', 0., 10.0).name('volumeScale');
	folder.add(params, 'substeps', 1, 50).name('substeps').step(1);
	folder.add(params, 'eps', 1, 10000).name('eps').step(1);
	folder.add(params, 'cfl', 0.1, 1.).name('cfl');
	folder.add(params, 'running');
	folder.close();

	console.log("restDensity: ", params.restDensity);

	gui.destroy();

	let [width, height] = [initWidth, initHeight];
	let frameCount = 0;
	let timeSum = 0;
	function animate() {
		requestAnimationFrame(animate);

		[width, height] = initWorldSize(worldSize);

		let gravity = new Vec2(0, -30.);
		let start = performance.now();
		if (params.running) {
			pbdSimulator.simulate(1.0 / 60, gravity, mouseState, params, blobParams, vertexIndices, kernelRadius, width, height);
		}

		setPositionsFluid(pbdSimulator.fluidParticles, pointGeometry);
		setPositionsBlobs(pbdSimulator.blobParticles, triGeometry, blobParams.particlesPerBlob, blobCnt);

		renderer.setRenderTarget(renderTarget);
		renderer.render(scene, camera);

		renderer.setRenderTarget(temporaryTarget);
		renderer.render(quadScene, quadCamera);

		renderer.setRenderTarget(null);
		renderer.render(quadFinalScene, quadCamera);

		let end = performance.now();
		timeSum += end - start;
		if (frameCount % 100 == 0) {
		    console.log("all: ", timeSum / 100);
		    timeSum = 0;
		}

		frameCount++;
		fluidThresholdMaterial.uniforms.t.value += 0.2;
	}

	document.getElementById('loading-wrapper').style.display = 'none';
	animate();
}
main();
