/*
Service worker that keeps the page and its modules for offline use

Copyright 2026 Ahmet Inan <inan@aicodix.de>
*/

const cache_name = "modem";
const files = [
	"./", "index.html", "manifest.webmanifest", "icon.svg", "icon-192.png", "icon-512.png",
	"temp/encoder_44100.wasm", "temp/encoder_48000.wasm", "temp/decoder_44100.wasm", "temp/decoder_48000.wasm"
];
self.addEventListener("install", event => {
	event.waitUntil(caches.open(cache_name).then(cache => cache.addAll(files)).then(() => self.skipWaiting()));
});
self.addEventListener("activate", event => {
	event.waitUntil(self.clients.claim());
});
/* network first, so the page and its modules stay up to date while online, and the cache when offline */
self.addEventListener("fetch", event => {
	const request = event.request;
	if (request.method != "GET" || new URL(request.url).origin != self.location.origin)
		return;
	event.respondWith(fetch(request).then(response => {
		if (response.ok) {
			const copy = response.clone();
			event.waitUntil(caches.open(cache_name).then(cache => cache.put(request, copy)));
		}
		return response;
	}).catch(() => caches.match(request, { ignoreSearch: true }).then(response => response || Response.error())));
});
