/*
Check the encoder WebAssembly modules with the native decoder

Copyright 2026 Ahmet Inan <inan@aicodix.de>
*/

import fs from "fs";
import crypto from "crypto";
import { spawnSync } from "child_process";

/* the seconds whose frame lasts closest to 4 << code payload symbols */
const seconds = [1, 2, 3, 5, 9, 18, 35];

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
	/* the payload of every valid mode, the symbol count in the top three bits, the efficiency index in the low five */
	const payload_of = new Map();
	for (let mode = 0; mode < 256; ++mode) {
		const ok = wasm.start(mode) >= 0;
		if (mode >= 224 && ok) {
			console.log(`FAILED ${rate} Hz mode ${mode}: unsupported mode accepted`);
			++failed;
		}
		if (ok)
			payload_of.set(mode, wasm.payload_length());
	}
	/* every valid mode up to 32 symbols, the lowest and the highest valid index of the longer ones */
	const modes = [];
	for (let code = 0; code < 7; ++code) {
		const valid = [...payload_of.keys()].filter(mode => mode >> 5 == code);
		modes.push(...(code < 4 ? valid : [valid[0], valid[valid.length - 1]]));
	}
	for (const mode of modes) {
		const name = `${rate} Hz mode ${mode} (${4 << (mode >> 5)} symbols, index ${mode & 31})`;
		const length = payload_of.get(mode);
		/* the most robust mode for exactly its payload is the mode itself */
		if (wasm.pick_mode(seconds[mode >> 5], length) != mode) {
			console.log(`FAILED ${name}: pick_mode`);
			++failed;
			continue;
		}
		if (wasm.start(mode) < 0) {
			console.log(`FAILED ${name}: start`);
			++failed;
			continue;
		}
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
	/* no mode carries more than the largest payload */
	const largest = Math.max(...payload_of.values());
	for (const s of seconds) {
		if (wasm.pick_mode(s, largest + 1) != -1) {
			console.log(`FAILED ${rate} Hz ${s} seconds: payload too long accepted`);
			++failed;
		}
	}
}
fs.rmSync("temp/check.wav", { force: true });
fs.rmSync("temp/check.dat", { force: true });
console.log(`${passed} passed, ${failed} failed`);
process.exit(failed ? 1 : 0);
