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

function wav_file(samples, rate, channels) {
	const header = Buffer.alloc(44);
	header.write("RIFF", 0);
	header.writeUInt32LE(36 + 4 * samples.length, 4);
	header.write("WAVE", 8);
	header.write("fmt ", 12);
	header.writeUInt32LE(16, 16);
	header.writeUInt16LE(3, 20);
	header.writeUInt16LE(channels, 22);
	header.writeUInt32LE(rate, 24);
	header.writeUInt32LE(4 * rate * channels, 28);
	header.writeUInt16LE(4 * channels, 32);
	header.writeUInt16LE(32, 34);
	header.write("data", 36);
	header.writeUInt32LE(4 * samples.length, 40);
	return Buffer.concat([header, Buffer.from(samples.buffer, samples.byteOffset, 4 * samples.length)]);
}

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
			const count = wasm.values_length();
			frames.push({
				mode: wasm.oper_mode(),
				call_sign: new TextDecoder().decode(new Uint8Array(wasm.memory.buffer, wasm.call_sign.value, 9)).trim(),
				payload: length < 0 ? null : Buffer.from(new Uint8Array(wasm.memory.buffer, wasm.payload.value, length)),
				values: count < 0 ? null : new Float32Array(wasm.memory.buffer, wasm.values.value, 2 * count).slice(),
			});
		}
	}
	return frames;
}

let passed = 0, failed = 0;
/* with sync detection in the time and in the frequency domain */
for (const [rate, sync] of [[44100, 0], [44100, 1], [48000, 0], [48000, 1]]) {
	const module = new WebAssembly.Module(fs.readFileSync(`temp/decoder_${rate}.wasm`));
	const wasm = new WebAssembly.Instance(module).exports;
	wasm.__wasm_call_ctors();
	wasm.frequency_sync(sync);
	const sync_name = sync ? "frequency sync" : "time sync";
	for (let mode = 0; mode < 128; ++mode) {
		if (((mode >> 1) & 7) > 3)
			continue;
		const modulation = modulations[(mode >> 4) & 7], code_rate = code_rates[(mode >> 1) & 7], frame_size = frame_sizes[mode & 1];
		const name = `${rate} Hz ${sync_name} ${modulation} ${code_rate} ${frame_size}`;
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
	for (const count of [1, 256, 1000, 32768]) {
		const name = `${rate} Hz ${sync_name} analog with ${count} complex values`;
		const values = new Float32Array(2 * count).map(() => 2 * Math.random() - 1);
		fs.writeFileSync("temp/check.dat", wav_file(values, wasm.analog_rate(), 2));
		spawnSync("./encode", ["temp/check.wav", String(rate), "32", "ANONYMOUS", "ANALOG", "temp/check.dat"], { stdio: "ignore" });
		const frames = decode(wasm, wav_samples(fs.readFileSync("temp/check.wav")));
		/* analog frames are lossy, the encoder clips and filters, so the worst error grows with their length */
		const expected = new Float32Array(2 * 256 * Math.ceil(count / 256));
		expected.set(values);
		const decoded = frames.length == 1 && frames[0].mode == 128 + (count - 1 >> 8) ? frames[0].values : null;
		const error = decoded && decoded.length == expected.length ? expected.map((v, i) => Math.abs(v - decoded[i])) : null;
		const ok = error && error.every(e => e < 0.1) && error.reduce((a, e) => a + e) / error.length < 0.002;
		if (ok) {
			++passed;
		} else {
			console.log(`FAILED ${name}: decoded values differ`);
			++failed;
		}
	}
}
fs.rmSync("temp/check.wav", { force: true });
fs.rmSync("temp/check.dat", { force: true });
console.log(`${passed} passed, ${failed} failed`);
process.exit(failed ? 1 : 0);
