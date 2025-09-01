declare global {
  interface Window {
    complexome?: [Map<string, Set<string>>, Map<string, string>, Map<string, string[]>, Map<string, string[]>];
    userdata?: Map<string, [number, number]>;
  }
}

export {};
