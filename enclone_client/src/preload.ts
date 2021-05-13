/* eslint-disable  @typescript-eslint/no-explicit-any */

import { contextBridge, ipcRenderer } from 'electron';

export type EncloneApi = typeof encloneApi

// Because of updated Electron security policies, it's better to expose api methods
// individually rather than expose the full main process api.

const encloneApi = {
  // FIXME: Maybe should be integer.
  getClonotypesById: (clonotypeId: string) => {
    ipcRenderer.send('get-clonotypes-by-id', clonotypeId);
  },
  onClonotypesByIdResponse: (cb: (data: any) => void) => {
    ipcRenderer.on('get-clonotypes-by-id-response', (event, ...args) => cb(args[0]));
  },
  enclone: (encloneArgs: string) => {
    ipcRenderer.send('enclone', encloneArgs);
  },
  onEncloneResponse: (cb: (data: any) => void) => {
    ipcRenderer.on('enclone-response', (event, ...args) => cb(args[0]));
  },
};

// Create the api that the render process can use to tell the main process what to do.
export function createRenderProcessApi(): void {
  contextBridge.exposeInMainWorld('encloneApi', encloneApi);
}

createRenderProcessApi();

