/*
Check the encoder WebAssembly modules with the native decoder

Copyright 2026 Ahmet Inan <inan@aicodix.de>
*/

import fs from "fs";
import crypto from "crypto";
import { spawnSync } from "child_process";

const modulations = ["BPSK", "QPSK", "8PSK", "QAM16", "QAM64", "QAM256", "QAM1024", "QAM4096"];
const code_rates = ["1/2", "2/3", "3/4", "5/6"];
const frame_sizes = ["short", "normal"];

function base40(str) {
	let acc = 0n;
	for (const c of str.toUpperCase()) {
		acc *= 40n;
		if (c == "/")
			acc += 3n;
		else if (c >= "0" && c <= "9")
			acc += BigInt(c.charCodeAt(0) - 48 + 4);
		else if (c >= "A" && c <= "Z")
			acc += BigInt(c.charCodeAt(0) - 65 + 14);
		else if (c != " ")
			return -1n;
	}
	return acc;
}

function wav_file(samples, rate) {
	const header = Buffer.alloc(44);
	header.write("RIFF", 0);
	header.writeUInt32LE(36 + 4 * samples.length, 4);
	header.write("WAVE", 8);
	header.write("fmt ", 12);
	header.writeUInt32LE(16, 16);
	header.writeUInt16LE(3, 20);
	header.writeUInt16LE(1, 22);
	header.writeUInt32LE(rate, 24);
	header.writeUInt32LE(4 * rate, 28);
	header.writeUInt16LE(4, 32);
	header.writeUInt16LE(32, 34);
	header.write("data", 36);
	header.writeUInt32LE(4 * samples.length, 40);
	return Buffer.concat([header, Buffer.from(samples.buffer, samples.byteOffset, 4 * samples.length)]);
}

function encode(wasm, rate, mode, payload, call_sign) {
	const chunks = [new Float32Array(rate)];
	const take = len => chunks.push(new Float32Array(wasm.memory.buffer, wasm.audio.value, len).slice());
	let len = wasm.start(mode);
	if (len < 0)
		return null;
	take(len);
	new Uint8Array(wasm.memory.buffer, wasm.payload.value, payload.length).set(payload);
	if (wasm.begin(payload.length, call_sign) < 0)
		return null;
	while ((len = wasm.next()) > 0)
		take(len);
	take(wasm.stop());
	chunks.push(new Float32Array(rate));
	const samples = new Float32Array(chunks.reduce((sum, chunk) => sum + chunk.length, 0));
	let pos = 0;
	for (const chunk of chunks) {
		samples.set(chunk, pos);
		pos += chunk.length;
	}
	return samples;
}

let passed = 0, failed = 0;
const call_sign = base40("ANONYMOUS");
for (const rate of [44100, 48000]) {
	const module = new WebAssembly.Module(fs.readFileSync(`temp/encoder_${rate}.wasm`));
	const wasm = new WebAssembly.Instance(module).exports;
	wasm.__wasm_call_ctors();
	for (let mode = 0; mode < 128; ++mode) {
		const name = `${rate} Hz ${modulations[(mode >> 4) & 7]} ${code_rates[(mode >> 1) & 7]} ${frame_sizes[mode & 1]}`;
		if (((mode >> 1) & 7) > 3) {
			if (wasm.start(mode) != -1) {
				console.log(`FAILED ${mode}: unsupported mode accepted`);
				++failed;
			}
			continue;
		}
		if (wasm.start(mode) < 0) {
			console.log(`FAILED ${name}: start`);
			++failed;
			continue;
		}
		const length = wasm.payload_length();
		const payload = crypto.randomBytes(length - (mode & 3));
		const samples = encode(wasm, rate, mode, payload, call_sign);
		if (!samples) {
			console.log(`FAILED ${name}: encode`);
			++failed;
			continue;
		}
		fs.writeFileSync("temp/check.wav", wav_file(samples, rate));
		fs.rmSync("temp/check.dat", { force: true });
		spawnSync("./decode", ["temp/check.wav", "temp/check.dat"], { stdio: "ignore" });
		const expected = Buffer.concat([payload, Buffer.alloc(length - payload.length)]);
		const decoded = fs.existsSync("temp/check.dat") ? fs.readFileSync("temp/check.dat") : null;
		if (decoded && expected.equals(decoded)) {
			++passed;
		} else {
			console.log(`FAILED ${name}: decoded data differs`);
			++failed;
		}
	}
}
fs.rmSync("temp/check.wav", { force: true });
fs.rmSync("temp/check.dat", { force: true });
console.log(`${passed} passed, ${failed} failed`);
process.exit(failed ? 1 : 0);
