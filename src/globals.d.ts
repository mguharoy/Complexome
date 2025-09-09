declare global {
  interface CustomEventMap {
    "custom:table-sort": CustomEvent
  }

  interface Window extends EventTarget{
    complexome?: [Map<string, Set<string>>, Map<string, string>, Map<string, string[]>, Map<string, string[]>];
    userdata?: Map<string, [number, number]>;

    addEventListener<K extends keyof CustomEventMap>(event: K, listener: ((this: Window, event: CustomEventMap[K]) => any) | null, options?: AddEventListenerOptions | boolean): void;
    addEventListener(type: string, callback: EventListenerOrEventListenerObject | null, options?: AddEventListenerOptions | boolean): void;
  }
}

export {};
