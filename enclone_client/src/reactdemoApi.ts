/* eslint-disable @typescript-eslint/no-explicit-any, @typescript-eslint/explicit-module-boundary-types */
import path from 'path';

import * as grpc from '@grpc/grpc-js';
import { ipcMain } from 'electron';
import * as protoLoader from '@grpc/proto-loader';

const PROTO_NAME = 'server.proto';

const PROTO_PATH =
  process.env.NODE_ENV === 'development'
    ? `../enclone_server_proto/${PROTO_NAME}`
    : path.join(process.resourcesPath, PROTO_NAME);

const packageDefinition = protoLoader.loadSync(PROTO_PATH, {
  keepCase: true,
  longs: String,
  enums: String,
  defaults: true,
  oneofs: true,
});

export const server = grpc.loadPackageDefinition(packageDefinition).enclone_server;


const getClient = () => new server.Analyzer('127.0.0.1:7000', grpc.credentials.createInsecure());

export const clonotypesById = (msg: string, callback: any): void => {
  getClient().GetClonotype({ clonotypeNumber: msg }, (err: any, response: any) => {
    callback(response);
  });
};

export const enclone = (msg: string, callback: any): void => {
  getClient().Enclone({ args: msg }, (err: any, response: any) => {
    callback(response);
  });
};

// Until we get errors back from server, catch them here
function encloneError() {
  return new Error('Enclone error');
}

// Create the event listeners that the render process can dispatch to the main process
// to do the grpc work.
export function initGrpcProcessListeners(): void {
  ipcMain.on('get-clonotypes-by-id', (event, clonotypeId) => {
    clonotypesById(clonotypeId, function (response: any, err: any) {
      if (err) {
        console.error(err);
      } else {
        event.sender.send('get-clonotypes-by-id-response', response);
      }
    });
  });

  ipcMain.on('enclone', (event, encloneArgs) => {
    enclone(encloneArgs, function (response: any, err: any) {
      if (err || !response) {
        event.sender.send('enclone-response', err || encloneError());
      } else {
        event.sender.send('enclone-response', response);
      }
    });
  });
}
