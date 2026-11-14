/*
Check the decoder WebAssembly modules with the native encoder

Copyright 2026 Ahmet Inan <inan@aicodix.de>
*/

import fs from "fs";
import crypto from "crypto";
import { spawnSync } from "child_process";

/* the seconds whose frame lasts closest to 2 << code payload symbols, zero for the pings of code zero */
const seconds = [0, 1, 2, 3, 5, 9, 18, 35];
const decoded_frame = 1 << 11, ping = 1 << 13;

function wav_samples(buffer) {
	const pos = buffer.indexOf("data") + 8;
	const data = buffer.subarray(pos);
	return new Float32Array(data.buffer.slice(data.byteOffset, data.byteOffset + data.length - data.length % 4));
}

function decode(wasm, samples) {
	const block = new Float32Array(wasm.memory.buffer, wasm.samples.value, 1024);
	const frames = [];
	for (let pos = 0; pos < samples.length || wasm.busy(); pos += block.length) {
		block.fill(0);
		if (pos < samples.length)
			block.set(samples.subarray(pos, pos + block.length));
		for (let offset = 0; offset < block.length;) {
			offset += wasm.feed(offset, block.length - offset);
			if (wasm.events() & ping)
				frames.push({ ping: wasm.oper_mode() & 31 });
			if (!(wasm.events() & decoded_frame))
				continue;
			const length = wasm.payload_length();
			frames.push({
				mode: wasm.oper_mode(),
				call_sign: new TextDecoder().decode(new Uint8Array(wasm.memory.buffer, wasm.call_sign.value, 9)).trim(),
				payload: Buffer.from(new Uint8Array(wasm.memory.buffer, wasm.payload.value, length)),
			});
		}
	}
	return frames;
}

let passed = 0, failed = 0;
for (const rate of [44100, 48000]) {
	const module = new WebAssembly.Module(fs.readFileSync(`temp/decoder_${rate}.wasm`));
	const wasm = new WebAssembly.Instance(module).exports;
	wasm.__wasm_call_ctors();
	/* the encoder module knows the payload of every mode */
	const encoder = new WebAssembly.Instance(new WebAssembly.Module(fs.readFileSync(`temp/encoder_${rate}.wasm`))).exports;
	encoder.__wasm_call_ctors();
	const payload_of = new Map();
	for (let mode = 32; mode < 256; ++mode)
		if (encoder.start(mode) >= 0)
			payload_of.set(mode, encoder.payload_length());
	/* pings carry their value in the mode and nothing else */
	for (const value of [0, 7, 31]) {
		spawnSync("./encode", ["temp/check.wav", String(rate), "32", "ANONYMOUS", "0", String(value)], { stdio: "ignore" });
		const frames = decode(wasm, wav_samples(fs.readFileSync("temp/check.wav")));
		if (frames.length == 1 && frames[0].ping === value) {
			++passed;
		} else {
			console.log(`FAILED ${rate} Hz ping ${value}`);
			++failed;
		}
	}
	/* every valid mode up to 32 symbols, the lowest and the highest valid index of the longer ones */
	const modes = [];
	for (let code = 1; code < 8; ++code) {
		const valid = [...payload_of.keys()].filter(mode => mode >> 5 == code);
		modes.push(...(code < 5 ? valid : [valid[0], valid[valid.length - 1]]));
	}
	for (const mode of modes) {
		const name = `${rate} Hz mode ${mode} (${2 << (mode >> 5)} symbols, index ${mode & 31})`;
		/* one byte more than the next lower index carries, so the native encoder picks this mode */
		const lower = [...payload_of.keys()].filter(m => m >> 5 == mode >> 5 && m < mode);
		const payload = crypto.randomBytes(lower.length ? payload_of.get(lower[lower.length - 1]) + 1 : 1);
		fs.writeFileSync("temp/check.dat", payload);
		spawnSync("./encode", ["temp/check.wav", String(rate), "32", "ANONYMOUS", String(seconds[mode >> 5]), "temp/check.dat"], { stdio: "ignore" });
		const frames = decode(wasm, wav_samples(fs.readFileSync("temp/check.wav")));
		const ok = frames.length == 1 && frames[0].mode == mode && frames[0].call_sign == "ANONYMOUS" &&
			frames[0].payload.length == payload_of.get(mode) && frames[0].payload.subarray(0, payload.length).equals(payload);
		if (ok) {
			++passed;
		} else {
			console.log(`FAILED ${name}: ${frames.length} frames decoded`);
			++failed;
		}
	}
}
fs.rmSync("temp/check.wav", { force: true });
fs.rmSync("temp/check.dat", { force: true });
console.log(`${passed} passed, ${failed} failed`);
process.exit(failed ? 1 : 0);
