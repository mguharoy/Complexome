const addResourcesToCache = async (resources) => {
  const cache = await caches.open("v1");
  await cache.addAll(resources);
};

async function putInCache(request, response) {
  if (response.ok) {
    const cache = await caches.open("v1");
    await cache.put(request, response);
  } else {
    console.error(`Cannot cache response: ${response}.`);
  }
}

async function cacheFirst({ request, preloadResponse }) {
  const responseFromCache = await caches.match(request);
  if (responseFromCache) {
    console.log("Returning cached response");
    return responseFromCache;
  }
  // not in the cache

  const preloaded = await preloadResponse;
  if (preloaded) {
    console.log("Returning pre-loaded response");
    putInCache(request, preloaded.clone());
    return preloaded;
  }
  // not already preloaded

  try {
    const responseFromNetwork = await fetch(request.clone());
    putInCache(request, responseFromNetwork.clone());
    return responseFromNetwork;
  } catch (error) {
    console.error(`Network request failed with: ${error}`);
    return new Response("Network error", {
      status: 408,
      headers: { "Content-Type": "text/plain" },
    });
  }
}

self.addEventListener("activate", (event) => {
  event.waitUntil(self.registration?.navigationPreload.enable());
  console.log("Service worker activated");
});

self.addEventListener("install", (event) => {
  console.log("Caching complexome species...");
  event.waitUntil(
    addResourcesToCache([
      "https://ftp.ebi.ac.uk/pub/databases/intact/complex/current/complextab/10090.tsv",
      "https://ftp.ebi.ac.uk/pub/databases/intact/complex/current/complextab/559292.tsv",
      "https://ftp.ebi.ac.uk/pub/databases/intact/complex/current/complextab/83333.tsv",
    ]),
  );
});

self.addEventListener("fetch", (event) => {
  event.respondWith(
    cacheFirst({
      request: event.request,
      preloadResponse: event.preloadResponse,
    }),
  );
});
