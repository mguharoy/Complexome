/// <declare lib="esnext" />

export {};
declare global {
  interface WorkerGlobalScope {
    cache: Map<string, [Map<string, Set<string>>, Map<string, string>, Map<string, string[]>, Map<string, string[]>]>;
  }
}

export {};
