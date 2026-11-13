/*
Check the decoder WebAssembly modules with the native encoder

Copyright 2026 Ahmet Inan <inan@aicodix.de>
*/

import fs from "fs";
import crypto from "crypto";
import { spawnSync } from "child_process";

const modulations = ["BPSK", "QPSK", "8PSK", "QAM16", "QAM64", "QAM256", "QAM1024", "QAM4096"];
const code_rates = ["1/2", "2/3", "3/4", "5/6"];
const frame_sizes = ["short", "normal"];
const decoded_frame = 1 << 11;

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
	for (let mode = 0; mode < 128; ++mode) {
		if (((mode >> 1) & 7) > 3)
			continue;
		const modulation = modulations[(mode >> 4) & 7], code_rate = code_rates[(mode >> 1) & 7], frame_size = frame_sizes[mode & 1];
		const name = `${rate} Hz ${modulation} ${code_rate} ${frame_size}`;
		const payload = crypto.randomBytes(7000);
		fs.writeFileSync("temp/check.dat", payload);
		spawnSync("./encode", ["temp/check.wav", String(rate), "32", "ANONYMOUS", modulation, code_rate, frame_size, "temp/check.dat"], { stdio: "ignore" });
		const frames = decode(wasm, wav_samples(fs.readFileSync("temp/check.wav")));
		const ok = frames.length == 1 && frames[0].mode == mode && frames[0].call_sign == "ANONYMOUS" &&
			frames[0].payload.equals(payload.subarray(0, frames[0].payload.length));
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
