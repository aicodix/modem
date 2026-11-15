/*
Check that the native decoder and the decoder WebAssembly modules still decode the recordings of the frozen air interface

Copyright 2026 Ahmet Inan <inan@aicodix.de>
*/

import fs from "fs";
import { spawnSync } from "child_process";

const decoded_frame = 1 << 11, ping = 1 << 13;

/* the samples of a mono WAV file with 16 bit integer or 32 bit float samples */
function wav_samples(buffer) {
	let format = 0, bits = 0, pos = 12;
	while (pos + 8 <= buffer.length) {
		const id = buffer.toString("ascii", pos, pos + 4), size = buffer.readUInt32LE(pos + 4);
		if (id == "fmt ") {
			format = buffer.readUInt16LE(pos + 8);
			bits = buffer.readUInt16LE(pos + 22);
		} else if (id == "data") {
			const data = buffer.subarray(pos + 8, pos + 8 + size);
			if (format == 1 && bits == 16) {
				const samples = new Float32Array(data.length / 2);
				for (let i = 0; i < samples.length; ++i)
					samples[i] = data.readInt16LE(2 * i) / 32768;
				return samples;
			}
			if (format == 3 && bits == 32)
				return new Float32Array(data.buffer.slice(data.byteOffset, data.byteOffset + data.length));
			break;
		}
		pos += 8 + size + (size & 1);
	}
	throw new Error("unsupported WAV file");
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
				call_sign: new TextDecoder().decode(new Uint8Array(wasm.memory.buffer, wasm.call_sign.value, 9)).trim(),
				payload: Buffer.from(new Uint8Array(wasm.memory.buffer, wasm.payload.value, length)),
			});
		}
	}
	return frames;
}

let passed = 0, failed = 0;
const modules = {};
for (const rate of [44100, 48000]) {
	modules[rate] = new WebAssembly.Instance(new WebAssembly.Module(fs.readFileSync(`temp/decoder_${rate}.wasm`))).exports;
	modules[rate].__wasm_call_ctors();
}
/* the payload of a frame is the text padded with zeros to the payload of its mode */
function padded(text, decoded) {
	return decoded.length >= text.length && text.equals(decoded.subarray(0, text.length)) && decoded.subarray(text.length).every(b => b == 0);
}
/* each recording with its sample rate, call sign and text, or ping value */
for (const frame of JSON.parse(fs.readFileSync("frozen/frames.json"))) {
	const name = `frozen/${frame.wav}`;
	const text = frame.payload ? fs.readFileSync(`frozen/${frame.payload}`) : null;
	const frames = decode(modules[frame.rate], wav_samples(fs.readFileSync(name)));
	const ok = frames.length == 1 && (text ? frames[0].call_sign == frame.call_sign && padded(text, frames[0].payload) : frames[0].ping === frame.ping);
	if (ok) {
		++passed;
	} else {
		console.log(`FAILED ${name}: decoder module, ${frames.length} frames decoded`);
		++failed;
	}
	fs.rmSync("temp/check.dat", { force: true });
	const result = spawnSync("./decode", [name, "temp/check.dat"], { encoding: "utf8" });
	/* the native decoder right aligns the call sign to nine characters */
	const call_sign = result.stderr.split("\n").some(line => /^call sign: +(\S+)$/.exec(line)?.[1] == frame.call_sign);
	const native = text ? fs.existsSync("temp/check.dat") && padded(text, fs.readFileSync("temp/check.dat")) && call_sign : result.stderr.includes(`ping ${frame.ping} received`);
	if (native) {
		++passed;
	} else {
		console.log(`FAILED ${name}: native decoder`);
		++failed;
	}
}
fs.rmSync("temp/check.dat", { force: true });
console.log(`${passed} passed, ${failed} failed`);
process.exit(failed ? 1 : 0);
